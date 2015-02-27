% AuE 817
% Project #1: Conventional Vehicle Analysis
%
% By: Pierluigi Pisu
%==========================================================================
% This m-file will generate figures based on the aquired during the last
% simulation. Code can be added to this file to generate additinoal
% figures. 
%==========================================================================

clc
% Drive Cyce
if DRIV_cycle == 1
    name = 'Federal Urban Driving Schedule (FUDS) Velocity Profile';
elseif DRIV_cycle == 1
    name = 'Federal Highway Driving Schedule (FHDS) Velocity Profile';
elseif DRIV_cycle == 3
    name = 'US06 Driving Cycle Velocity Profile';
elseif DRIV_cycle == 4
    name = 'Acceleration Test Velocity Profile';
elseif DRIV_cycle == 5
    name = 'NEDC Velocity Profile';
end


% figure(11)
% eff=zeros(1,length(t));
% for i=1:(length(t))
%    eff(i)=interp2(ICE.map_w,ICE.map_T,ICE.map_eff',ICE_w(i),ICE_T(i)); 
% end
% [f,x]=hist(eff,50);%# create histogram from a normal distribution.
% bar(x,f/sum(f));
% 
% figure(12)
% eff=zeros(1,length(t));
% for i=1:(length(t))
%    eff(i)=interp2(EM.map_w,EM.map_T,EM.map_eff',EM_w(i),EM_T(i)); 
% end
% [f,x]=hist(eff,50);%# create histogram from a normal distribution.
% bar(x,f/sum(f));



% 3-a
% figure(10)
%     plot(t_cyc, v_cyc)  % Disired velocity
%     xlabel('Time [s]');
%     ylabel('Vehicle Velocity [m/s]');
%     grid on
%     hold on;
%     plot(t, VEH_v, 'r');    % Actual velocity
% xlabel('Time [s]')
% ylabel('Velocity [m/s]')
% title('Actual velocity vs Desired velocity (UDDS)')
%         hold on;
% plot(t_cyc+597, v_cyc)
% plot(t_cyc+2*597, v_cyc)
%     plot(t, VEH_v, 'r')
   
% 3-b
%     figure(2)
%     plot(t, ICE_w/2/pi*60)
%     title('Engine speed vs. Time')
%     xlabel('Time [s]')
%     ylabel('Engine speed [RPM]')
%     grid on
%     ylim([700 2500])
    
% 3-c
%     figure(3)
%     plot(t, ICE_T)
%     title('Engine torque vs. Time')
%     xlabel('Time [s]')
%     ylabel('Engine torque [Nm]')
%     grid on

% 3-d
%     figure(4)
%     plot(t, TRAN_gear)
%     title('Gear selection vs. Time')
%     xlabel('Time [s]')
%     ylabel('Gear selection')



% 3-e (please refer the last line of this code)
    
% 3-f
%     figure(6)
%     plot(t, EM_T)
%     title('Electric motor torque vs. Time')
%     xlabel('Time [s]')
%     ylabel('Electric motor torque [Nm]')
%     grid on

% 3-g
%     figure(7)
%     plot(t, EM_w/2/pi*60)
%     title('Electric motor speed vs. Time')
%     xlabel('Time [s]')
%     ylabel('Electric motor speed [RPM]')
%     grid on

% 3-h
%     figure(8)
%     plot(t, BATT_V_oc)
%     title('Battery voltage vs. Time')
%     xlabel('Time [s]')
%     ylabel('Battery voltage [V]')
%     grid on
%     ylim([785 805])
    
%     figure(9)
%     plot(t, BATT_i)
%     title('Battery current vs. Time')
%     xlabel('Time [s]')
%     ylabel('Battery current [A]')
%     grid on

% 3-i    
%     figure(10)
%     plot(t, BATT_SOC)
%     title('Battery state of charge vs. Time')
%     xlabel('Time [s]')
%     ylabel('Battery SOC')
%     grid on
%     
% 3-j
EM=EM2;
figure(8)
cont2 = [0.39 0.45 0.50 0.55 0.60 0.65 0.7 0.75 0.80 0.83 0.85 0.9 0.95];
axis1 = ([50 625 0 350]);
[c h] = contour(EM.map_w, EM.map_T, EM.map_eff', cont2);
hold on;
grid on;
clabel(c,h);
plot(EM.maxT_w, EM.maxT_T, 'r-', 'LineWidth', 2)
    min_w = min(EM.maxT_w);
    max_w = max(EM.maxT_w);
    xp=[min_w EM.maxT_w' max_w];
    min_T = min(EM.maxT_T);
    max_T = max(EM.maxT_T);
    yp=[800 EM.maxT_T' 800]; 
patch(xp,yp,'white','Edgecolor','r')
    min_w = min(EM.maxT_w);
    max_w = max(EM.maxT_w);
    xp=[min_w EM.maxT_w' max_w];
    min_T = min(EM.minT_T);
    max_T = max(EM.minT_T);
    yp=[-800 EM.minT_T' -800]; 
patch(xp,yp,'white','Edgecolor','r')
%axis([ICE.idle_w ICE.red_line_w 0 (max(ICE.maxT_T))])
hold off
% axis = axis1;
xlabel('EM Speed [rad/s]')
ylabel('EM Torque [Nm]')
title('EM Efficiency Map')
legend('Iso-Efficiency Contours', 'Maximum Torque Curve', 4)
hold on;
plot(EM_w,EM_T,'o');
%==========================================================================
% 4 Run the simulator as it is for the gasoline engine over one
% % acceleration test

%  4-1 Top speed of the vehicle
% v_top=max(VEH_v)

%  4-2 0-60 MPH time
% qqq=VEH_v*3.6/1.609;
%    for i=1:101 
%         if VEH_v(i)>=60*1.609/3.6
%         t_60=i;
%         break;
%         end
%    end
%    t_60=t_60-3+(60-qqq(t_60-1))/(qqq(t_60)-qqq(t_60-1))

% %  4-3 1/4mile time, exit speed
%    for i=1:101
%       if VEH_x(i)>=0.25*1609
%           t_qm=i;
%            v_qm=VEH_v(i);
%           break; 
%       end
%    end
% t_qm=t_qm-2+(0.25*1609-VEH_x(t_qm-1))/(VEH_x(t_qm)-VEH_x(t_qm-1))
% v_qm=VEH_v(i-1)+(0.25*1609-VEH_x(i-1))/(VEH_x(i)-VEH_x(i-1))*(VEH_v(i)-VEH_v(i-1));
% v_qm=v_qm*3.6/1.609
% 
% % 4-4 0-20 20-40 30-50 50-70 mph time intervals
% t1=0;
% t2=0;
% t3=0;
% t4=0;
% t5=0;
% for i=1:101
%    if t1==0
%        if VEH_v(i)>=20*1.609/3.6
%            t1=i;
%        end 
%    end
%    if (t1>0&&t2==0)
%        if VEH_v(i)>=30*1.609/3.6
%            t2=i;
%        end
%    end
%    if (t1>0&&t2>0&&t3==0)
%        if VEH_v(i)>=40*1.609/3.6
%            t3=i;
%        end
%    end
%    if (t1>0&&t2>0&&t3>0&&t4==0)
%        if VEH_v(i)>=50*1.609/3.6
%            t4=i;
%        end
%    end
%    if (t1>0&&t2>0&&t3>0&&t4>0&&t5==0)
%        if VEH_v(i)>=70*1.609/3.6
%            t5=i;
%            break;
%        end
%    end
% end
% 
% t1=t1-2+(20-qqq(t1-1))/(qqq(t1)-qqq(t1-1));
% t2=t2-2+(30-qqq(t2-1))/(qqq(t2)-qqq(t2-1));
% t3=t3-2+(40-qqq(t3-1))/(qqq(t3)-qqq(t3-1));
% t4=t4-2+(50-qqq(t4-1))/(qqq(t4)-qqq(t4-1));
% t5=t5-2+(70-qqq(t5-1))/(qqq(t5)-qqq(t5-1));
% t1
% t3-t1
% t4-t2
% t5-t4

%==========================================================================
% (5)

%==========================================================================
% Engine Maps
if ICE_choice == 1
    label1 = '3.4L V6 Gasoliine Engine Efficiency Map';
    label2 = '3.4L V6 Gasoline Engine Fuel Consumption Map [g/s]';
    cont1 = [0:0.5:11];
    cont2 = [0.0 0.05 0.10 0.15 0.20 0.25 0.27 0.28 0.30 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39];
    axis1 = ([50 625 0 350]);
    min_w = min(ICE.maxT_w);
    max_w = max(ICE.maxT_w);
    xp=[min_w ICE.maxT_w max_w];
    min_T = min(ICE.maxT_T);
    max_T = max(ICE.maxT_T);
    yp=[800 ICE.maxT_T 800]; 
elseif ICE_choice == 2
    label1 = '1.9L I4 Diesel Engine Efficiency Map';
    label2 = '1.9L I4 Diesel Engine Fuel Consumption Map [g/s]';
    cont1 = [0:0.5:7];
    cont2 = [0.0 0.05 0.10 0.15 0.20 0.25 0.30 0.33 0.35 0.37 0.39 0.40 0.41 0.415];
    axis1 = ([100 475 0 400]);
    min_w = min(ICE.maxT_w);
    max_w = max(ICE.maxT_w);
    xp=[min_w ICE.maxT_w max_w];
    min_T = min(ICE.maxT_T);
    max_T = max(ICE.maxT_T);
    yp=[800 ICE.maxT_T 800];   
elseif ICE_choice == 3
    label1 = 'Engine Efficiency Map';
    label2 = 'Engine Fuel Consumption Map [g/s]';
    cont1 = [0:0.5:11];
    cont2 = [0.0 0.05 0.10 0.15 0.20 0.22 0.25 0.27 0.29 0.30 0.31 0.315 0.32 0.33 ];
    axis1 = ([100 475 0 400]);
    min_w = min(ICE.maxT_w);
    max_w = max(ICE.maxT_w);
    xp=[min_w ICE.maxT_w max_w];
    min_T = min(ICE.maxT_T);
    max_T = max(ICE.maxT_T);
    yp=[800 ICE.maxT_T 800];   
end
     
figure(2)
[c,h] = contour(ICE.map_w, ICE.map_T, ICE.map_m_dot'.*1000, cont1);
%scatter(ICE_w,ICE_T,ICE_m_f_dot);
grid on;
clabel(c,h);
hold on;
plot(ICE.maxT_w, ICE.maxT_T, 'r-', 'LineWidth', 2)
patch(xp,yp,'white','Edgecolor','r')
%axis([ICE.idle_w ICE.red_line_w 0 (max(ICE.maxT_T))])
hold off
axis = axis1;
xlabel('Engine Speed [rad/s]')
ylabel('Engine Torque [Nm]')
title(label2)
legend('Iso-Efficiency Contours', 'Maximum Torque Curve', 4)
hold on
plot(ICE_w,ICE_T,'o');
ylim([0 140])

figure(3)
[c,h] = contour(ICE.map_w, ICE.map_T, ICE.map_eff', cont2);
hold on;
grid on;
clabel(c,h);
plot(ICE.maxT_w, ICE.maxT_T, 'r-', 'LineWidth', 2)
patch(xp,yp,'white','Edgecolor','r')
%axis([ICE.idle_w ICE.red_line_w 0 (max(ICE.maxT_T))])
hold off
axis = axis1;
xlabel('Engine Speed [rad/s]')
ylabel('Engine Torque [Nm]')
title(label1)
title('Engine Efficiency Map ');
legend('Iso-Efficiency Contours', 'Maximum Torque Curve', 4)
% 3-e
hold on;
% plot(ICE_w(1:25,1),ICE_T(1:25,1),'o');
hold on
plot(ICE_w,ICE_T,'o');
ylim([0 140])

%==========================================================================



    