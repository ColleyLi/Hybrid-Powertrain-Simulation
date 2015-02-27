% AuE 817
% Project #2: Hybrid-Electric Vehicle Analysis
%Modified by P. Pisu On 02/02/2010


%==========================================================================
% Locate the parameter(s) of interest within this m-file. Once you have
% made the necessary changes, run this m-file prior to running the
% Conventional Vehicle Simulator. Everytime you make a change within this
% m-file, you must run it to load all of the necessary vehicle parameters
% to the Workspace for use in the Simulator. 
%==========================================================================

clc
close all
clear all
global ENV VEH TC TRAN ICE EM  


Ah = 7.3;        % Battery capacity, A.h
Thermal=1;
load EM2_map
load Mat_files/BATT_Data
BATT.Ah=7.3;
BATT.V0=3.2;
N=250;
BATT.m=BATT.Ah*BATT.V0*N/60*1.8;
BATT.num_series=N;
save('Mat_files/BATT_Data.mat','BATT')
T_amb = 20+273;     % Ambient temperature, deg. C
hA_0 = 0.07;    % Effective convective heat transfer per cell, no fan, W/deg. C
fs_0 = 1;       % Initial fan setting
mc = 50;        % Effective heat capacity per cell, J/deg. C
fs_act = 1;         % No fan control
fs(1) = fs_0; 


%===========================
load('ecms_table_n5.mat')
%==========================================================================
%  Environmentla Data

ENV.grade = 0;            
% Road Grade [%]
ENV.g = 9.81;               % Acceleration Due to Gravity [m/s^2]
ENV.den_air = 1.225;         % Air Density [kg/m^3]
ENV.C_fric_tire_gnd = 0.2;  % Dynamic Coefficient of Friction Between Tire and Ground (unitless)

save('Mat_files/ENV_Data.mat', 'ENV');

%==========================================================================
% Transmission Data
% TRAN.ratio_gear=[ 1st   2nd   3rd   4th   ]; % Gear
TRAN.ratio_gear = [2.847 1.552 1.000 0.700  ]; % Transmission Gear Ratios

TRAN.ratio_diff = 4.130;                                  % Differential Ratio
TRAN.eff = 0.92;                                         % Transmission Efficiency

save('Mat_files/TRAN_Data.mat', 'TRAN');

%==========================================================================
% Vehicle Parameters
VEH.M = 1300+75;                   % Mass [kg]
VEH.A_frontal = 2.48;           % Frontal Area [m^2]
VEH.C_drag = 0.32;             % Coefficient of Drag (unitless)
VEH.C_roll_res = 0.008;         % Coefficient of Rolling Resistance (unitless)
VEH.r_tire = 0.357;            % Tire Radius [m]
VEH.brake_F_propor = 0.60;                      % Front Brake Proportioning
VEH.brake_R_propor = 1 - VEH.brake_F_propor;    % Rear Brake Proportioning

load Mat_files/TRAN_Data.mat

for i = 1:length(TRAN.ratio_gear)   % Mass Factor
    VEH.mass_factor_index(i) = i;
    VEH.mass_factor(i) = 1 + 0.04*TRAN.ratio_gear(i)*TRAN.ratio_diff + 0.0025*(TRAN.ratio_gear(i)*TRAN.ratio_diff)^2;
end

save('Mat_files/VEH_Data.mat', 'VEH');

%==========================================================================
% Torque Converter Data
% 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % DO NOT CHANGE ANYTHING BELOW THIS LINE %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load Mat_files/ICE_Data.mat;
load Mat_files/EM_Data.mat;

TC.c1 = 5.7656E-03;
TC.c2 = 0.3107E-03;
TC.c3 = -5.4323E-03;
TC.c4 = 3.4325E-03;
TC.c5 = 2.2210E-03;
TC.c6 = -4.6041E-03;
TC.c7 = -6.7644E-03;
TC.c8 = 32.0084E-03;
TC.c9 = -25.2441E-03;

% Threshold Speed Ratio (at which mode transfers from Torque Multiplication to Torque Coupling
TC.thresh_w_ratio = ((TC.c2 - TC.c5) + sqrt((TC.c5 - TC.c2)^2 - 4*(TC.c6 - TC.c3)*(TC.c4 - TC.c1)))/(2*(TC.c6 - TC.c3));

% Threshold Pump Torque Factor (at which mode transfers from Torque Multiplication to Torque Coupling
TC.thresh_T_p = (TC.c7/TC.thresh_w_ratio^2 + TC.c8/TC.thresh_w_ratio + TC.c9);

% Speed Ratio vs. Torque Ratio Curve
% For POSITIVE Torus Flow: w_pump > w_turbine
% TC.w_ratio = w_turbine / w_pump
% TC.T_ratio = T_turbine / T_pump

i_index = 0;
i_step = 0.01;
for i = i_step:i_step:floor(TC.thresh_w_ratio/i_step)*i_step      % Torque Multiplication Stage for Positive Torus Flow
    i_index = i_index + 1;
    TC.pos_w_ratio(i_index) = i;
    TC.pos_T_ratio(i_index) = TC.c1/(TC.c4 + TC.c5*TC.pos_w_ratio(i_index) + TC.c6*TC.pos_w_ratio(i_index)^2)...
                    + TC.c2/(TC.c4/TC.pos_w_ratio(i_index) + TC.c5 + TC.c6*TC.pos_w_ratio(i_index))...
                    + TC.c3/(TC.c4/TC.pos_w_ratio(i_index)^2 + TC.c5/TC.pos_w_ratio(i_index) + TC.c6);
end
for i = ceil(TC.thresh_w_ratio/i_step)*i_step:i_step:1      % Torque Coupling Stage for Positive Torus Flow
    i_index = i_index + 1;
    TC.pos_w_ratio(i_index) = i;
    TC.pos_T_ratio(i_index) = 1;
end

% Speed Ratio vs. Torque Ratio Curve
% For NEGATIVE Torus Flow: w_turbine > w_pump
% TC.w_ratio = w_pump / w_turbine
% TC.T_ratio = T_pump / T_turbine

i_index = 0;
i_step = 0.01;
for i = 1:i_step:2                                          % Torque Coupling Stage for Negative Torus Flow (No Torque Multiplication Stage for Negative Torus Flow)
    i_index = i_index + 1;
    TC.neg_w_ratio(i_index) = i;
    TC.neg_T_ratio(i_index) = 1;
end

% Torque Converter Efficiency
Pump_w_step = 5;
Turbine_w_step = 25;
i_index = 0; j_index = 0;
for i = floor(ICE.idle_w/Pump_w_step)*Pump_w_step:Pump_w_step:ceil(ICE.red_line_w/Pump_w_step)*Pump_w_step
    i_index = i_index + 1;
    TC.map_w_p(i_index) = i;
    for j = 0:Turbine_w_step:ceil(max(EM.map_w)/EM.gearbox_ratio*TRAN.ratio_diff*max(TRAN.ratio_gear)/Turbine_w_step)*Turbine_w_step
        j_index = j_index + 1;
        TC.map_w_t(j_index) = j;
        if TC.map_w_p(i_index) >= TC.map_w_t(j_index)
            TC.map_eff(i_index, j_index) = TC.map_w_t(j_index)*interp1(TC.pos_w_ratio, TC.pos_T_ratio, (TC.map_w_t(j_index)/TC.map_w_p(i_index)), 'linear', 'extrap')/TC.map_w_p(i_index);
        else TC.map_eff(i_index, j_index) = TC.map_w_p(i_index)/(interp1(TC.neg_w_ratio, TC.neg_T_ratio, (TC.map_w_t(j_index)/TC.map_w_p(i_index)), 'linear', 'extrap')*TC.map_w_t(j_index));
        end
    end
    j_index = 0;
end

save('Mat_files/TC_Data.mat', 'TC')

%==========================================================================
% Electric Motor Data

load Mat_files/BATT_Data.mat
VEH.M = VEH.M +BATT.num_series*BATT.V0*BATT.Ah/60*1.8;                   % Mass [kg]
save('Mat_files/BATT_Data.mat', 'BATT')
% Gearbox Ratio and Efficiency
EM.gearbox_ratio = 10.946;
EM.gearbox_eff = 0.90;
    
% Maximum Torque (Nm) - Speed (rad/s) Curve
EM.maxT_w = [500 996 2999 5002 7005 8998 11001]*2*pi/60;
EM.maxT_T = [132.1 132.8 128.9 121.8 109.0 80.9 62.3]*1.355818;
EM.minT_T = -EM.maxT_T;

% Efficiency vs. Speed (rad/s) and Torque (Nm)
EM_eff_w_temp = [500*ones(1,11) 996*ones(1,11) 2999*ones(1,11) 5002*ones(1,11)...
    7005*ones(1,11) 8998*ones(1,11) 11001*ones(1,11)]*2*pi/60;

EM_eff_T_temp = [1.2 7.5 15.2 24.2 33.7 43.7 57.5 62.7 73.2 82.7 132.1...
    1.2 7.1 15.4 24.3 33.6 43.4 57.3 62.5 72.3 82.0 132.8...
    0.9 5.4 14.1 22.7 31.5 41.1 54.5 59.2 68.9 78.1 128.9...
    0.7 4.3 11.7 19.8 28.4 37.8 50.7 55.1 64.5 73.0 121.8...
    0.5 3.5 9.7 17.0 25.1 34.0 46.3 50.3 58.7 64.9 109.0...
    0.5 2.9 8.1 14.6 22.1 29.4 39.2 42.6 49.3 54.0 80.9...
    0.5 2.5 7.0 12.7 17.8 23.3 32.1 35.0 39.8 44.1 62.3]*1.355818;

EM_invert_eff = [35.1 63.0 72.9 73.4 71.4 73.3 74.7 75.0 74.6 76.5 74.7...
    47.2 76.6 71.3 81.5 81.4 91.7 73.8 77.5 83.7 80.7 80.5...
    71.7 86.9 90.1 89.1 89.9 89.5 89.1 89.9 89.4 88.3 87.6...
    76.9 88.1 91.8 91.2 91.8 91.7 91.9 91.8 91.1 91.3 89.7...
    77.5 89.8 93.2 92.6 92.5 92.7 92.9 92.9 92.4 91.8 91.2...
    79.9 90.3 93.6 93.5 93.7 94.3 93.8 93.4 92.9 92.2 93.7...
    81.6 91.1 94.2 94.4 94.5 94.5 94.3 93.9 93.1 93.0 93.3]/100;

EM_motor_eff = [89.8 82.3 71.6 73.9 74.7 71.5 67.9 66.6 64.7 60.9 52.7...
    100.6 84.7 92.4 81.6 81.8 71.4 86.8 82.0 74.0 75.3 67.1...
    86.3 86.8 88.2 90.0 89.7 89.7 89.6 88.6 88.0 88.3 83.9...
    83.8 87.8 89.3 91.6 91.6 91.8 91.4 91.2 91.3 90.4 88.3...
    80.8 86.2 89.0 91.5 92.5 92.8 92.6 92.3 92.3 91.9 89.8...
    80.9 85.7 88.2 91.3 92.6 92.6 92.7 93.1 92.6 92.0 88.2...
    96.5 89.0 89.6 91.4 92.5 92.6 92.3 92.3 91.8 90.9 88.0]/100;

EM_overall_eff = EM_invert_eff.*EM_motor_eff;

EM.map_w = [500:500:11000]*2*pi/60;
EM.map_T = [5:5:180]';

EM_eff_map_temp = griddata(EM_eff_w_temp, EM_eff_T_temp, EM_overall_eff, EM.map_w, EM.map_T)';

EM.map_T = [-fliplr(EM.map_T') EM.map_T'];
EM.map_eff = horzcat(fliplr(EM_eff_map_temp), EM_eff_map_temp);

% Maximum Torque (Nm) - Speed (rad/s) Curve while Battery Limited
j_index = 0;
for i = 1:length(EM.map_w)
    for j = (length(EM.map_T)/2 + 1):length(EM.map_T)
        j_index = j_index + 1;
        if EM.map_T(j) <= interp1(EM.maxT_w, EM.maxT_T, EM.map_w(i), 'linear', 'extrap') & ...
            EM.map_w(i)*EM.map_T(j)/EM.map_eff(i,j) <= BATT.P_max_disch
            EM_maxTcurve_batt_lim(i,j_index) = EM.map_T(j);
        else EM_maxTcurve_batt_lim(i,j_index) = NaN;
        end
    end
    j_index = 0;
end

% Minimum Torque (Nm) - Speed (rad/s) Curve while Battery Limited
for i = 1:length(EM.map_w)
    for j = 1:length(EM.map_T)/2
        if EM.map_T(j) >= interp1(EM.maxT_w, EM.minT_T, EM.map_w(i), 'linear', 'extrap') & ...
            EM.map_w(i)*EM.map_T(j)*EM.map_eff(i,j) >= BATT.P_max_ch
            EM_minTcurve_batt_lim(i,j) = EM.map_T(j);
        else EM_minTcurve_batt_lim(i,j) = NaN;
        end
    end
end

for i = 1:length(EM.map_w)
    EM.batt_lim_maxT_w(i) = EM.map_w(i);
    EM.batt_lim_maxT_T(i) = max(EM_maxTcurve_batt_lim(i,:));
    EM.batt_lim_minT_T(i) = min(EM_minTcurve_batt_lim(i,:));
end

% Eliminate NaN's by extrapolation
EM_eff_map_check = isnan(EM_eff_map_temp);
for i = 1:length(EM.map_w)
    for j = 1:length(EM.map_T)/2
        if EM_eff_map_check(i,j) == 1
            EM_eff_map_temp(i,j) = EM_eff_map_temp(i,j-1) + (EM_eff_map_temp(i,j-1)-EM_eff_map_temp(i,j-2));
        end
    end
end
EM.map_eff = horzcat(fliplr(EM_eff_map_temp), EM_eff_map_temp);

save('Mat_files/EM_Data.mat','EM')

%==========================================================================
% Driving Cycle Information
Sport=0;
disp('Available Drivng Cycles for Simulation:')
disp('     1 - FUDS: Federal Urban Driving Schedule')
disp('     2 - FHDS: Federal Highway Driving Schedue')
disp('     3 - US06: More Agressive Urban Driving Schedule')
disp('     4 - Acceleration Test')
disp('     5 - UDDS: Urban Dynamometer Driving Schedule')
disp('     6 - HWFET: Highway Fuel Economy Test')
disp('     7 - Cycle4.87 miles')
disp('     8 - Cycle10.6 miles')
DRIV_cycle = input('Enter the number corresponding to the driving cycle you wish to analyze and Press ENTER: ');

if DRIV_cycle == 1
    load Velocity_Profiles/FUDS;
elseif DRIV_cycle == 2
   	load Velocity_Profiles/FHDS;
elseif DRIV_cycle == 3
    load Velocity_Profiles/US06;
    v_cyc = v_cyc.*0.27778;
elseif DRIV_cycle ==4
    Sport=1
    disp('Sport is on now !!!')
    t_cyc = [1:1:1000];
    for i=1:length(t_cyc)
        v_cyc(i) = 200;
    end
    t_cyc = t_cyc';
    v_cyc = v_cyc';
elseif DRIV_cycle == 5
    load Velocity_Profiles/UDDS;
    v_cyc = v_cyc.*1.609/3.6;
elseif DRIV_cycle == 6
    load Velocity_Profiles/HWFET;
    v_cyc = v_cyc.*1.609/3.6;
elseif DRIV_cycle == 7
    load Velocity_Profiles/cycle487;
    v_cyc = v_cyc.*1.609/3.6;
elseif DRIV_cycle == 8
    load Velocity_Profiles/cycle1060;
    v_cyc = v_cyc.*1.609/3.6;
end


cycle_num = input('Enter the number of times you would like this Driving Schedule to be executed and Press ENTER: ');

if cycle_num==0
    disp('Invalid Entry. Re-execute this m-file and enter a value greater than zero.')
else
    Stoptime = length(t_cyc)*cycle_num;
    
    save('Mat_files/Drive_cycle.mat','DRIV_cycle','Stoptime','v_cyc','t_cyc')

end
fprintf('\n')

%==========================================================================
% Engine Selection

disp('Available Engines for Simulation')
disp('     1 - 3.4L V6 Gasoline')
disp('     2 - 1.9L I4 Diesel')
disp('     3 - 1.3L I4 Gasoline 817 FINAL Team 3')
ICE_choice = input('Enter the number corresponding to the engine you wish to use and press ENTER: ');

if ICE_choice == 1
    load Mat_files/ICE_Gas_Data.mat
elseif ICE_choice == 2
    load Mat_files/ICE_Diesel_Data.mat
elseif ICE_choice == 3
    load Engine2_map.mat
    VEH.M = VEH.M+(-160+ICE.m);                   % Mass [kg]
    save('Mat_files/VEH_Data.mat', 'VEH');
else
    disp('Invalid Selection. Re-execute this m-file and enter either (1), (2) or (3)')
end

save('Mat_files/ICE_Data.mat','ICE')

disp('Initialization data loaded successfully.')
disp('You are now ready to run the project simulator.')

%**********************************************************************
soc_h=0.85;
soc_l=0.55;

% load EMSC map

%load('new_em.mat')
%load('gear_map.mat')
%load('ecms_table2.mat')
%load('ecms_table_pem.mat')


% load new EM map
load('EM2_map')
VEH.M = VEH.M+EM2.m;                   % Mass [kg]
save('Mat_files/VEH_Data.mat', 'VEH');




