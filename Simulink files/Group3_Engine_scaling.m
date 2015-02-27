
%================================================================================================
% Willans Line Model & Scaling
%======================================================================

% ICE displacement
ICE.scale.Vd=1.3e-3;
%======================================================================

ICE.data.w_RPM=xlsread('engine_map_fuel_consumption','C11:O11'); % [rad/s]
ICE.data.T=xlsread('engine_map_fuel_consumption','C12:O23')'; % [Nm]
ICE.data.m_dot=xlsread('engine_map_fuel_consumption','C31:O42')'*32/44.4/3600; % [kg/s] 
ICE.data.w=ICE.data.w_RPM*2*pi/60;    % [rad/s]
ICE.data.Vd=1.498e-3;               % m^3
ICE.data.B=78e-3;                   % m
ICE.data.S=78.4e-3;                 % m

ICE.scale.B=(ICE.scale.Vd/1.1/pi)^(1/3);
ICE.scale.S=1.1*ICE.scale.B;
% ICE.scale.B=60e-3;                % m
% ICE.scale.S=70e-3;                 % m
% ICE.scale.Vd=pi*ICE.scale.B^2*ICE.scale.S;               % m^3
ICE.data.QLHV=44e6;                 % lower heating value [J/kg]
nr=2;
% % aa=zeros(13,12)
% % for i=1:13
% %     aa(i,:)=ICE.data.w_RPM(i)*ones(1,12)/60*2*pi;
% % end
% % maxp=max(max(aa.*ICE.data.T))
% Maximum Torque [Nm] - Speed [rad/s] 

ICE.data.maxT_w = ICE.data.w;
ICE.data.maxT_T= ICE.data.T(:,12)';

%===============
A=zeros(7,7);   
B=zeros(7,1);
for i=1:length(ICE.data.w)                
    for j=1:length(ICE.data.T(1,:))
        ICE.data.cm(i)=ICE.data.S/pi*ICE.data.w(i);
        ICE.data.pma(i,j)=2*pi*nr/ICE.data.Vd*ICE.data.QLHV*ICE.data.m_dot(i,j)/ICE.data.w(i); 
        ICE.data.pme(i,j)=2*pi*nr/ICE.data.Vd*ICE.data.T(i,j);         
        if ICE.data.T(i,j)~=0
            cm=ICE.data.cm(i);
            pma=ICE.data.pma(i,j);
            pme=ICE.data.pme(i,j);            
            A(1,:)=A(1,:)+[pma  cm*pma  cm*cm*pma  -pma*pma  -cm*pma*pma -1  -cm*cm];           
            A(2,:)=A(2,:)+[pma*cm  cm*pma*cm  cm*cm*pma*cm  -pma*pma*cm  -cm*pma*pma*cm  -1*cm  -cm*cm*cm];
            A(3,:)=A(3,:)+[pma*cm^2  cm*pma*cm^2  cm*cm*pma*cm^2  -pma*pma*cm^2  -cm*pma*pma*cm^2  -1*cm^2  -cm*cm*cm^2];
            A(4,:)=A(4,:)+[pma*cm^3  cm*pma*cm^3  cm*cm*pma*cm^3  -pma*pma*cm^3  -cm*pma*pma*cm^3  -1*cm^3  -cm*cm*pme^3];
            A(5,:)=A(5,:)+[pma*pma  cm*pma*pma  cm*cm*pma*pma  -pma*pma*pma  -cm*pma*pma*pma  -1*pma  -cm*cm*pma];
            A(6,:)=A(6,:)+[pma*pma^2  cm*pma*pma^2  cm*cm*pma*pma^2  -pma*pma*pma^2  -cm*pma*pma*pma^2  -1*pma^2  -cm*cm*pma^2];
            A(7,:)=A(7,:)+[pma*pma^3  cm*pma*pma^3  cm*cm*pma*pma^3  -pma*pma*pma^3  -cm*pma*pma*pma^3  -1*pma^3  -cm*cm*pma^3];
            B(:,1)=B(:,1)+[pme  pme*cm  pme*cm^2  pme*cm^3  pme*pma  pme*pma^2  pme*pma^3]';
        end
    end
end

C=inv(A)*B;
w=ICE.data.w;
wRPM=ICE.data.w_RPM;
cm=ICE.scale.S/pi.*w;
T=0:200;
pme=4*pi*T/ICE.scale.Vd;
for i=1:length(w)
      for j=1:length(T)
            b=(C(1)+C(2)*cm(i)+C(3)*cm(i)^2);
            a=(C(4)+C(5)*cm(i));
            pml(i)=C(6)+C(7)*cm(i)^2;
            if (b^2-4*a*(pme(j)+pml(i)))<0
                pma(i,j)=NaN;
            else
                pma(i,j)=(b-sqrt(b^2-4*a*(pme(j)+pml(i))))/(2*a);
            end
          ICE.scale.m_dot(i,j)=pma(i,j)*w(i)*ICE.scale.Vd/(2*pi*nr*ICE.data.QLHV)*1000;
          ICE.scale.eff(i,j)=pme(j)/pma(i,j);
      end 
end

% Eliminate NaN's by extrapolation
ICE_m_f_dot_map_check1 = isnan(ICE.scale.m_dot);
ICE_m_f_dot_map_check2 = isnan(ICE.scale.eff);
for i = 1:length(w)
    for j = 1:length(T)
        if ICE_m_f_dot_map_check1(i,j) == 1
            ICE.scale.m_dot(i,j) = ICE.scale.m_dot(i,j-1) + (ICE.scale.m_dot(i,j-1)-ICE.scale.m_dot(i,j-2));
        end
        if ICE_m_f_dot_map_check2(i,j) == 1
            ICE.scale.eff(i,j) = ICE.scale.eff(i,j-1) + (ICE.scale.eff(i,j-1)-ICE.scale.eff(i,j-2));
        end
    end
end

for i = 1:length(w)
    for j = 1:length(T)
        if ICE.scale.eff(i,j) < 0.01
            ICE.scale.eff(i,j) = 0.01;
        end
    end
end

% Fuel consumption
figure(1)
clf
cont = [0.5 1 2 3 4 5 6 7 8 9 10 11 12 13];
[cc,h] = contour(wRPM,T,ICE.scale.m_dot',cont);
clabel(cc,h);
axis([0 7000 0 200])
grid on;
hold on;

% Max Torque Line fit
% ICE.data.cm_max=ICE.data.L*w/pi;
% p=polyfit(ICE.data.cm_max,ICE.data.pme_max,2);
% pme_max=p(1)*cm.^2+p(2)*cm+p(3);        
% Max torque Line scale
% ICE.maxT_T= pme_max*1e5*ICE.scale.Vd/(2*pi*ICE.nr);

ICE.data.pme_max=4*pi*ICE.data.maxT_T/ICE.data.Vd/1e5; %bar
ICE.maxT_T=ICE.data.pme_max*1e5*ICE.scale.Vd/(2*pi*nr);
maxw_RPM=wRPM;
maxw=max(maxw_RPM);
minw=min(maxw_RPM);
maxT=200;
minT=200;
xp = [minw maxw_RPM maxw ];
yp = [minT ICE.maxT_T maxT ];
patch(xp,yp,'white')
plot(wRPM,ICE.maxT_T,'LineWidth',2);
% plot(4400,174,'o','LineWidth',2,'color','r');

% Engine Idle/redline
idleRPM=750;
redlineRPM=6300;
maxT=200;
h = line([idleRPM idleRPM],[0 maxT]);
set(h,'LineWidth',2,'Color','red')
h = line([redlineRPM redlineRPM],[0 maxT]);
set(h,'LineWidth',2,'Color','red')
text(redlineRPM+200,maxT/2.5,'Red Line','Rotation',90,'HorizontalAlignment','center', ...
    'Color','red')
text(idleRPM-200,maxT/2.5,'Idle','Rotation',90,'HorizontalAlignment','center', ...
    'Color','red')
 xlabel('Engine Speed (RPM)')
 ylabel('Torque (Nm)')
 title('Willans line model with curve fits - Fuel consumption map (g/sec)')
 
 
% Efficiency
figure(2)
clf
cont = [0.05 0.07 0.1 0.14 0.16 0.2 0.24 0.26 0.28 2.9 0.3 0.31 0.32 0.33 0.35];
[cc,h] = contour(wRPM,T,ICE.scale.eff',cont);
clabel(cc,h);
axis([0 7000 0 200])
grid on;
hold on;
      
% Max torque Line scale
ICE.data.pme_max=4*pi*ICE.data.maxT_T/ICE.data.Vd/1e5; %bar
ICE.maxT_T=ICE.data.pme_max*1e5*ICE.scale.Vd/(2*pi*nr);
maxw_RPM=wRPM;
maxw=max(maxw_RPM);
minw=min(maxw_RPM);
maxT=200;
minT=200;
xp = [minw maxw_RPM maxw ];
yp = [minT ICE.maxT_T maxT ];
patch(xp,yp,'white')
plot(wRPM,ICE.maxT_T,'LineWidth',2);
% plot(4400,174,'o','LineWidth',2,'color','r');

% Engine Idle/redline
idleRPM=750;
redlineRPM=6300;
maxT=200;
h = line([idleRPM idleRPM],[0 maxT]);
set(h,'LineWidth',2,'Color','red')
h = line([redlineRPM redlineRPM],[0 maxT]);
set(h,'LineWidth',2,'Color','red')
text(redlineRPM+200,maxT/2.5,'Red Line','Rotation',90,'HorizontalAlignment','center', ...
    'Color','red')
text(idleRPM-200,maxT/2.5,'Idle','Rotation',90,'HorizontalAlignment','center', ...
    'Color','red')
 xlabel('Engine Speed (RPM)')
 ylabel('Torque (Nm)')
 title('Willans line model with curve fits - Efficiency map')

 % maximum power vs engine speed
 figure(3)
 plot(wRPM,ICE.maxT_T.*w/1000,'LineWidth',2);
 hold on
%  plot(6000,98,'o','LineWidth',2,'color','r');
 xlabel('Engine Speed (RPM)')
 ylabel('Power (kW)')
 title('Willans line model with curve fits - Maximum Power vs Engine speed')
 grid on
% % ICE.scale.w=w;
% % ICE.scale.T=T;
% % ICE.scale.maxw_T=ICE.maxT_T;
% % ICE.scale.maxw_w=maxw_RPM;
% % 
% % load Copy_of_Engine_map.mat
% % ICEe.map_T=T;
% % ICEe.maxT_T=ICE.scale.maxw_T;
% % ICEe.map_m_dot=ICE.scale.m_dot/1000;
% % ICEe.map_eff=ICE.scale.eff;
% % ICEe.idle_m_dot=ICEe.map_m_dot(1,1);
% % ICE=ICEe;
% % save('Engine2_map.mat','ICE')
% % 
% % aa=zeros(13,1);
% % for i=1:13
% %     aa(i)=ICE.maxT_w(i)*ICE.maxT_T(i);
% % end
% % maxp_Eng=max(aa)
% % ICE.m=(maxp_Eng/1000/74.5)^(1.5)*160;
% % save('Engine2_map.mat','ICE')
%============================================================================================



