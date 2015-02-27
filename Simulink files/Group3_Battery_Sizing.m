% Battery modelling
%=============================================================
clear all
clc


%% Part 1 
% Perkaut equation fitting
I=[1 5 20 30];
C=[7450 7420 7320 6900];
figure(1)
plot(I,C,'*','color','r');
hold on
% General model:
%      f(x) = a*x^b
% Coefficients (with 95% confidence bounds):
       a =        7521;  %(6718, 8324)
       b =    -0.01692;  %(-0.0621, 0.02827)
c=1:0.1:30;
plot(c,a*c.^b,'LineWidth',2);
% xlim([0 2.2])
xlabel('Current (A)')
ylabel('Capacity (mAh)')
title('Curve fitting for Peukert equation')

 





%% Part 2 
% Discharge fitting
load('Original_data.mat')

SOC=1-[520/7450 2006/7450 3061/7450  3889/7450 4931/7450 6115/7450 5361/7420 1825/7420 2594/7420 3521/7320 4245/7320  1154/7320 5534/7320 5648/6900 3024/6900 4064/6900  4999/6900 2126/6900];
I=[1 1 1 1 1 1 5 5 5 20 20  20 20 30 30  30 30 30];
V=[ 3.2907 3.2618 3.2462 3.2351 3.2083 3.1683 3.07 3.1683 3.1572 2.8501 2.8145 2.9013 2.73 2.5096 2.665 2.6186 2.5675 2.6943];
T=273+[ 25.2 25.2 25.2 25.2 25.2 25.2 28.6 28.6 28.6 40.6 40.6 40.6 40.6 49.1 49.1  49.1  49.1 49.1];

C=zeros(length(T),12);
A=zeros(length(T),12);

for i=1:length(T)
    C(i,:)=[1 SOC(i) SOC(i)^2 T(i) SOC(i)*T(i) SOC(i)^2*T(i) -I(i) -I(i)*SOC(i) -I(i)*SOC(i)^2 -I(i)*T(i) -I(i)*SOC(i)*T(i) -I(i)*SOC(i)^2*T(i)]; 
    A(i,:)=[0 0 0 0 0 0 -1 -SOC(i) -SOC(i)^2 -T(i) -SOC(i)*T(i) -SOC(i)^2*T(i)]; % the Coefficient matrix of resistance
end 
d=V';
b=-0.001*ones(length(T),1);
X = lsqlin(C,d,A,b);


figure(2)
I= 10;
T=32.5+273;
soc=0:0.01:1;
VV=X(1)+X(2)*soc+X(3)*soc.^2+X(4)*T+X(5)*T*soc+X(6)*T*soc.^2-I*(X(7)+X(8)*soc+X(9)*soc.^2+X(10)*T+X(11)*T*soc+X(12)*T*soc.^2);
plot(1-cap_10/7370,vol_10,'b')
hold on
plot(soc,VV,'g')
xlabel('SOC')
ylabel('Voltage')
legend('Data','fitted')
title('V10 vs SOC discharging')
xlim([0.2 0.8])

% % % % % % % % % % 
% % E0&R


SOC=[0:1:100]/100;

T252=25.2+273;
T286=28.6+273;
T325=32.5+273;
T406=40.6+273;
T491=49.1+273;

figure(3)
E252=X(1)+X(2).*SOC+X(3).*SOC.^2+...
    X(4).*T252+X(5).*SOC.*T252+X(6).*SOC.^2.*T252;
E286=X(1)+X(2).*SOC+X(3).*SOC.^2+...
    X(4).*T286+X(5).*SOC.*T286+X(6).*SOC.^2.*T286;
E325=X(1)+X(2).*SOC+X(3).*SOC.^2+...
    X(4).*T325+X(5).*SOC.*T325+X(6).*SOC.^2.*T325;
E406=X(1)+X(2).*SOC+X(3).*SOC.^2+...
    X(4).*T406+X(5).*SOC.*T406+X(6).*SOC.^2.*T406;
E491=X(1)+X(2).*SOC+X(3).*SOC.^2+...
    X(4).*T491+X(5).*SOC.*T491+X(6).*SOC.^2.*T491;



plot(SOC,E252,'g',SOC,E286,'m',SOC,E325,'c',SOC,E406,'y',SOC,E491,'b')
xlabel('Normalized SOC [%]')
ylabel('E_0 [V]')
title('Open Circuit Voltage - Normalized SOC')
legend('25.2deg_C','28.6deg_C','32.5deg_C','40.6deg_C','49.1deg_C')
xlim([0.2 0.8])


figure(4)
R252=X(7)+X(8).*SOC+X(9).*SOC.^2+...
    X(10).*T252+X(11).*SOC.*T252+X(12).*SOC.^2.*T252;
R286=X(7)+X(8).*SOC+X(9).*SOC.^2+...
    X(10).*T286+X(11).*SOC.*T286+X(12).*SOC.^2.*T286;
R325=X(7)+X(8).*SOC+X(9).*SOC.^2+...
    X(10).*T325+X(11).*SOC.*T325+X(12).*SOC.^2.*T325;
R406=X(7)+X(8).*SOC+X(9).*SOC.^2+...
    X(10).*T406+X(11).*SOC.*T406+X(12).*SOC.^2.*T406;
R491=X(7)+X(8).*SOC+X(9).*SOC.^2+...
    X(10).*T491+X(11).*SOC.*T491+X(12).*SOC.^2.*T491;

plot(SOC,R252,'g',SOC,R286,'m',SOC,R325,'c',SOC,R406,'y',SOC,R491,'b')
xlabel('Normalized SOC [%]')
ylabel('R [ohm]')
title('Resistance - Normalized SOC')
legend('25.2deg_C','28.6deg_C','32.5deg_C','40.6deg_C','49.1deg_C')
xlim([0.2 0.8])





%% Part 3
% Charge fitting
load('Original_data.mat')
% with 03 05 2 3
SOC=1-[ 939/7300 2061/7300 2951/7300 4015/7300 5058/7300 6050/7300 939/7300 2061/7300 2951/7300 4015/7300 5058/7300 6050/7300 939/7300 2061/7300 2951/7300 4015/7300 5058/7300 6050/7300 1083/7300 2104/7300 3067/7300 4080/7300  5174/7300 6086/7300];
I=[0.3 0.3 0.3 0.3 0.3 0.3 0.5 0.5 0.5 0.5 0.5 0.5 2 2 2 2 2 2 3 3 3 3 3 3];
V=[ 3.26 3.3289 3.342 3.3596 3.3684 3.3859  3.26 3.3289 3.342 3.3596 3.3684 3.3859 3.26 3.3289 3.342 3.3596 3.3684 3.3859 3.3049 3.3655 3.3787 3.3854 3.4137 3.4355];
T=[ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

% % wirh 05 2 3
% SOC=1-[ 939/7300 2061/7300 2951/7300 4015/7300 5058/7300 6050/7300 939/7300 2061/7300 2951/7300 4015/7300 5058/7300 6050/7300 939/7300 2061/7300 2951/7300 4015/7300 5058/7300 6050/7300 1083/7300 2104/7300 3067/7300 4080/7300  5174/7300 6086/7300];
% I=[0.3 0.3 0.3 0.3 0.3 0.3 0.5 0.5 0.5 0.5 0.5 0.5 2 2 2 2 2 2 3 3 3 3 3 3];
% V=[  3.26 3.3289 3.342 3.3596 3.3684 3.3859 3.26 3.3289 3.342 3.3596 3.3684 3.3859 3.3049 3.3655 3.3787 3.3854 3.4137 3.4355];
% T=[  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

% % % with 05 3
% SOC=1-[ 939/7300 2061/7300 2951/7300 4015/7300 5058/7300 6050/7300 939/7300 2061/7300 2951/7300 4015/7300 5058/7300 6050/7300 939/7300 2061/7300 2951/7300 4015/7300 5058/7300 6050/7300 1083/7300 2104/7300 3067/7300 4080/7300  5174/7300 6086/7300];
% I=[0.3 0.3 0.3 0.3 0.3 0.3 0.5 0.5 0.5 0.5 0.5 0.5 2 2 2 2 2 2 3 3 3 3 3 3];
% V=[ 3.26 3.3289 3.342 3.3596 3.3684 3.3859  3.26 3.3289 3.342 3.3596 3.3684 3.3859 3.26 3.3289 3.342 3.3596 3.3684 3.3859 3.3049 3.3655 3.3787 3.3854 3.4137 3.4355];
% T=[ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];


C=zeros(length(T),12);
for i=1:length(T)
    C(i,:)=[1 SOC(i) SOC(i)^2 T(i) SOC(i)*T(i) SOC(i)^2*T(i) -I(i) -I(i)*SOC(i) -I(i)*SOC(i)^2 -I(i)*T(i) -I(i)*SOC(i)*T(i) -I(i)*SOC(i)^2*T(i)]; 
    A(i,:)=[0 0 0 0 0 0 -1 -SOC(i) -SOC(i)^2 -T(i) -SOC(i)*T(i) -SOC(i)^2*T(i)]; % the Coefficient matrix of resistance
end 
% d=V'
% X=C\d;
d=V'
b=-0.001*ones(length(T),1);
X = lsqlin(C,d,A,b);


% C=zeros(length(T),12);
% for i=1:length(T)
%     C(i,:)=[1 SOC(i) SOC(i)^2 T(i) SOC(i)*T(i) SOC(i)^2*T(i) -I(i) -I(i)*SOC(i) -I(i)*SOC(i)^2 -I(i)*T(i) -I(i)*SOC(i)*T(i) -I(i)*SOC(i)^2*T(i)]; 
% %     A(i,:)=[0 0 0 0 0 0 -1 -SOC(i) -SOC(i)^2 -T(i) -SOC(i)*T(i) -SOC(i)^2*T(i)]; % the Coefficient matrix of resistance
% end 
% d=V';
% % b=-1*ones(length(T),1);
% % X = lsqlin(C,d,A,b);
% X = C\d;
% 
% SOC_1 = 1-cap_10/7450;
T = 0;

figure(5)
% I= 0.3;
% soc=0:0.01:1;
% VV03=X(1)+X(2)*soc+X(3)*soc.^2+X(4)*T+X(5)*T*soc+X(6)*T*soc.^2-I*(X(7)+X(8)*soc+X(9)*soc.^2+X(10)*T+X(11)*T*soc+X(12)*T*soc.^2);
% plot(soc,VV03,'g')
% hold on
% I= 0.5;
% soc=0:0.01:1;
% VV05=X(1)+X(2)*soc+X(3)*soc.^2+X(4)*T+X(5)*T*soc+X(6)*T*soc.^2-I*(X(7)+X(8)*soc+X(9)*soc.^2+X(10)*T+X(11)*T*soc+X(12)*T*soc.^2);
% plot(soc,VV05,'b')
% hold on
I= 1;
soc=0.2:0.01:0.85;
VV1=X(1)+X(2)*soc+X(3)*soc.^2+X(4)*T+X(5)*T*soc+X(6)*T*soc.^2-I*(X(7)+X(8)*soc+X(9)*soc.^2+X(10)*T+X(11)*T*soc+X(12)*T*soc.^2);
plot(1-soc,VV1,'r')
hold on
% I= 2;
% soc=0:0.01:1;
% VV2=X(1)+X(2)*soc+X(3)*soc.^2+X(4)*T+X(5)*T*soc+X(6)*T*soc.^2-I*(X(7)+X(8)*soc+X(9)*soc.^2+X(10)*T+X(11)*T*soc+X(12)*T*soc.^2);
% plot(soc,VV2,'g')
% hold on
% I= 3;
% soc=0:0.01:1;
% VV3=X(1)+X(2)*soc+X(3)*soc.^2+X(4)*T+X(5)*T*soc+X(6)*T*soc.^2-I*(X(7)+X(8)*soc+X(9)*soc.^2+X(10)*T+X(11)*T*soc+X(12)*T*soc.^2);
% plot(soc,VV3,'k')
% hold on

% plot(1-ch_cap_2/7300,ch_vol_2,'b')
% plot(1-ch_cap_03/7300,ch_vol_03,'b')
% plot(1-ch_cap_05/7300,ch_vol_05,'b')
% plot(1-ch_cap_3/7300,ch_vol_3,'k')
plot(ch_cap_1/7300,ch_vol_1,'b')

xlim([0.2 0.8])
legend('fitted','data')
title('V2 vs SOC charging')
xlabel('SOC')
ylabel('Voltage')
hold off

%%%
%%E&R
SOC=[0:1:100]/100;


figure(6)
E=X(1)+X(2).*SOC+X(3).*SOC.^2;
plot(1-SOC,E)
xlabel('Normalized SOC')
ylabel('E_0 [V]')
title('Open Circuit Voltage - Normalized SOC')
xlim([0.2 0.8])

figure(7)
R=X(7)+X(8).*SOC+X(9).*SOC.^2
plot(1-SOC,R)
xlabel('Normalized SOC')
ylabel('R [ohm]')
title('Resistance - Normalized SOC')
xlim([0.2 0.8])
V_oc_SOC=E(20:1:88);




