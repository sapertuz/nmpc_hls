%NonLinear model predivtive control for Inverted Pendulum
% 

% clear 
% clc

calc_state = 1;

%% Inverted Pendulum constants and data definition
global invPend
invPend.m = 7.3;    % [Kg] Uniformily distributed mass of the pendulum
invPend.M = 14.6;   % [Kg] Cart mass
invPend.g = 9.81;   % [m/s^2] Gravity
invPend.l = 1.2;    % [m] half the lenght of the pendulum
invPend.b = 14.6;   % [Kg/s] surface friction
invPend.h = 0.0136; % [Kg.m^2/s] rotation friction damping coeficient

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  T    tx   ty   tz
%Input limits
invPend.u_max  = [50.0];
invPend.u_min  = [-50.0];
invPend.state_upper_limits = [1, 1e3, 1e3, 1e3] ;
invPend.state_lower_limits = [-1, -1e3, -1e3, -1e3] ;

Ts = 0.1;
xmeasure      = [0.0, 0.0, 3.1415926536, 0.0];
%% Either perform new NMPC drone implementation
[time, xHistory, uHistory, xRef_out] = InvPendulumNMPC();
     
%% Or load previously generated data
% load data_line.mat
% load data_spiral.mat
% load data_ring.mat
% load data_ring_noInputWeights.mat

%% Or load from c data
% run ../nmpc_pso_xhls/matlab/nmpc_sim.m
% Ts = 0.1;
% uHistory = NMPC_SIM.control_history;
% xHistory = NMPC_SIM.state_history;
% xRef_out = NMPC_SIM.xref;

%% Recalculate state
clear xHistory_new
xHistory_new(1,:) = xmeasure;
for i=1:length(time)
    if i==1
        x_ant = xmeasure;
    else
        x_ant = xHistory_new(i-1,:);
    end
    xHistory_new(i,:) = InvPendulumSystem(0, x_ant, uHistory(i,:), Ts);
end