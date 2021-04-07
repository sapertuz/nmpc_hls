% 
clear 
clc

path_sim_data = '../nmpc_pso_hls/matlab/';

%% Or load previously generated data
run([path_sim_data,'nmpc_sim.m']);

Ts = Model.Ts;
time = 0:Ts:Model.SimulationTime;
% xHistory = [NMPC_SIM.state_history(:,1:2:end), NMPC_SIM.state_history(:,2:2:end)];
% xRef_out = [NMPC_SIM.xref(:,1:2:end), NMPC_SIM.xref(:,2:2:end)];
xHistory = NMPC_SIM.state_history;
xRef_out = NMPC_SIM.xref(:,1:3);

uHistory = NMPC_SIM.control_history;

% axis_lim = [-5 40 -15 15 -1 45];
axis_lim = [4 10 9 17 -1 3];

global drone
drone.l  = .25;

%% Vizualization
% close all
figure(1), clf
for i=1:length(time)-1
    animateQuadrotor(time(i), xHistory(i,:), uHistory(i,:), xRef_out, axis_lim);
    pause(Ts);    
end