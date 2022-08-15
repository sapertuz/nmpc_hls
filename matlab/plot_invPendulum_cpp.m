% 
close all
clear 
clc

path_sim_data = '../nmpc_pso_hls/matlab/';
xmeasure      = [0.0, 0.0, 3.1415926536, 0.0];

%% fmin matlab plot
load invPend_.mat

x_matlab = xHistory(:,1);
theta_matlab = xHistory(:,3);
u_matlab = uHistory;


%% PSO Plot
run([ path_sim_data, 'inverted_pendulum_pso.m'])

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

col = '#696969';	 %[1 1 1]*.15;

figure
subplot(3,1,1),
    plot(time, theta_matlab,'Color',col,'LineWidth',1)
    hold on
    plot(time, theta,'--g')
subplot(3,1,2),
    plot(time, x_matlab,'Color',col,'LineWidth',1)
    hold on
    plot(time, x,'--g')
subplot(3,1,3),
    plot(time, u_matlab,'Color',col,'LineWidth',1)
    hold on
    plot(time, uHistory,'--g')

%% MPSO plot
run([ path_sim_data, 'test.m'])

xHistory = NMPC_SIM.state_history;
uHistory = NMPC_SIM.control_history;

x = xHistory(:,1);
theta = xHistory(:,3);

%%
subplot(3,1,1),
    hold on
    plot(time, theta,'-k','LineWidth',1)
    plot(time, theta_ref, '--b')
    plot(time, theta_ref + 2*pi, '--b')
    axis([0, 50, -.5, 2*pi+.5])
    ylabel('State $\theta$ [rad]','Interpreter','latex','FontSize',12)
subplot(3,1,2),
    hold on
    plot(time, x,'-k','LineWidth',1)
    plot(time, x_ref, '--b')
    plot(time, repmat(x_max(1),n_step,1), '--r')
    plot(time, repmat(x_min(1),n_step,1), '--r')
    axis([0, 75, -1.2, 1.2])
    ylabel('State $x$ [m]','Interpreter','latex','FontSize',12)
subplot(3,1,3),
    hold on
    plot(time, uHistory,'-k','LineWidth',1)
    plot(time, u_ref, '--b')
    plot(time, repmat(u_max(1),n_step,1), '--r')
    plot(time, repmat(u_min(1),n_step,1), '--r')
    axis([0, 75, -55, 55])
    ylabel('Input [N]','Interpreter','latex','FontSize',12)

xlabel('time [s]','Interpreter','latex','FontSize',12)
lgd = legend({'Active-set QP', 'PSO', 'mPSO', 'Reference', 'Limits'},'Interpreter','latex', ...
        'FontSize',12);
lgd.NumColumns = 4;
% lgd.Location = 'southoutside';
