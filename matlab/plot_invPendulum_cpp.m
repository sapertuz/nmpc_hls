% 
close all
clear 
clc

path_sim_data = '../nmpc_pso_hls/matlab/';

%% Load previously generated data
run([ path_sim_data, 'inverted_pendulum_kpso.m'])

Ts = Model.Ts;
time = 0:Ts:(Model.SimulationTime-Ts);
n_step = length(time);
xHistory = NMPC_SIM.state_history;
uHistory = NMPC_SIM.control_history;
xRef_out = NMPC_SIM.xref;
uRef_out = NMPC_SIM.uref;

x_max = Model.x_max;
x_min = Model.x_min;
u_max = Model.u_max;
u_min = Model.u_min;

x = xHistory(:,1);
x_ref = xRef_out(:,1);
theta = xHistory(:,3);
theta_ref = xRef_out(:,3);

u_ref = uRef_out(:,1);
%% 
subplot(3,1,1),
    plot(time, theta,'-k')
    hold on
    plot(time, theta_ref, '--b')
    plot(time, theta_ref + 2*pi, '--b')
subplot(3,1,2),
    plot(time, x,'-k')
    hold on
    plot(time, x_ref, '--b')
    plot(time, repmat(x_max(1),n_step,1), '--r')
    plot(time, repmat(x_min(1),n_step,1), '--r')
subplot(3,1,3),
    plot(time, uHistory,'-k')
    hold on
    plot(time, u_ref, '--b')
    plot(time, repmat(u_max(1),n_step,1), '--r')
    plot(time, repmat(u_min(1),n_step,1), '--r')

