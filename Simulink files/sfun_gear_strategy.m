function [sys,x0,str,ts] = sfun_gear_shifting(t,x,u,flag)
%**************************************************************************
% This S-function implements the shifting strategy for a 6-gear transmission
%
% Last review: 10/24/2005
% By K. Koprubasi & M. Arnett 
%**************************************************************************


           %%%%%%%%%%%%%%%%%%%%%%%%%%
           % DO NOT CHANGE ANY OF   % 
           % THE CODE ABOVE LINE 106%
           %%%%%%%%%%%%%%%%%%%%%%%%%% 
           
switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes;

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%  
  case 1,  
    sys=[];
    
  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2,    
    sys=[];
    
  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,x,u);
    
  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
  case 4,   
    sys=[];
    
  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9,        
    sys=[];
    
  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    error(['Unhandled flag = ',num2str(flag)]);

end

%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts]=mdlInitializeSizes()

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%

sizes = simsizes;

sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 1;
sizes.NumInputs      = 4;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  = [];

%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [0 0];

%end mdlInitializeSizes


%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)

ICE_w = u(1);
up = u(2);
down = u(3);
TRAN_gear = u(4);

if (TRAN_gear==1) && (ICE_w>=up)
    gear_cmd = 1;    
elseif (TRAN_gear~=1) && (ICE_w > down) && (ICE_w < up)
    gear_cmd = 0;
elseif (TRAN_gear~=1) && (ICE_w <= down)
    gear_cmd = -1;
elseif (TRAN_gear~=1) && (ICE_w >= up)
    gear_cmd = 1;
else 
    gear_cmd = 0;
end

% TRAN_gear = TRAN_gear + gear_cmd;
% 
% if (TRAN_gear> 6)||(TRAN_gear < 1)
%     gear_cmd = 0; 
% end
% sys=TRAN_gear;


if (TRAN_gear+gear_cmd> 4)||(TRAN_gear+gear_cmd < 1)
    gear_cmd = 0; 
end
TRAN_gear = TRAN_gear + gear_cmd;
sys=TRAN_gear;
% end mdlOutputs

