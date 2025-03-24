% 
% close all
clear 
clc

path_sim_data = '../nmpc_pso_hls/matlab/';
xmeasure      = [0.0, 0.0, 3.1415926536, 0.0];

n_step = 700;

figure(2), clf

%% fmin matlab plot
load invPend_.mat

x_matlab = xHistory(1:n_step,1);
theta_matlab = xHistory(1:n_step,3);
u_matlab = uHistory(1:n_step);


%% PSO Plot
run([ path_sim_data, 'inverted_pendulum_pso.m'])

Ts = Model.Ts;
time = 0:Ts:(Model.SimulationTime-Ts);
time = time(1:n_step);

xHistory = NMPC_SIM.state_history;
uHistory = NMPC_SIM.control_history;

xRef_out = NMPC_SIM.xref;
uRef_out = NMPC_SIM.uref;

x_max = Model.x_max;
x_min = Model.x_min;
u_max = Model.u_max;
u_min = Model.u_min;

x_pso = xHistory(1:n_step,1);
theta_pso = xHistory(1:n_step,3);
u_pso = uHistory(1:n_step,1);

theta_ref = xRef_out(1:n_step,3);
x_ref = xRef_out(1:n_step,1);
u_ref = uRef_out(1:n_step,1);

%%
% run([ path_sim_data, 'inverted_pendulum_kpso.m'])
run([ path_sim_data, 'inverted_pendulum_kpso.m'])

xHistory = NMPC_SIM.state_history;
uHistory = NMPC_SIM.control_history;

x_mpso = xHistory(1:n_step,1);
theta_mpso = xHistory(1:n_step,3);
u_mpso = uHistory(1:n_step,1);


%%
col = '#696969';	 %[1 1 1]*.15;
col2 = '#708090';	 %[1 1 1]*.15;


%%
subplot(3,1,1),
    hold on
    plot(time, theta_matlab,'Color',col,'LineWidth',1)
    plot(time, theta_pso,'--','Color',col2)
    plot(time, theta_mpso,'-k','LineWidth',1)
    plot(time, theta_ref, '--b')
    plot(time, theta_ref + 2*pi, '--b')
    axis([0, 60, -.5, 2*pi+.5])
    ylabel('State $\theta$ [rad]','Interpreter','latex','FontSize',14)
    grid on
subplot(3,1,2),
    hold on
    plot(time, x_matlab,'Color',col,'LineWidth',1)
    plot(time, x_pso,'--','Color',col2)
    plot(time, x_mpso,'-k','LineWidth',1)
    plot(time, x_ref, '--b')
    plot(time, repmat(x_max(1),n_step,1), '-.r')
    plot(time, repmat(x_min(1),n_step,1), '-.r')
    axis([0, 60, -1.2, 1.2])
    ylabel('State $x$ [m]','Interpreter','latex','FontSize',14)
    grid on
subplot(3,1,3),
    hold on
    plot(time, u_matlab,'Color',col,'LineWidth',1)
    plot(time, u_pso,'--','Color',col2)
    plot(time, u_mpso,'-k','LineWidth',1)
    plot(time, u_ref, '--b')
    plot(time, repmat(u_max(1),n_step,1), '-.r')
    plot(time, repmat(u_min(1),n_step,1), '-.r')
    axis([0, 60, -55, 55])
    ylabel('Input [N]','Interpreter','latex','FontSize',14)
    grid on

xlabel('Time [s]','Interpreter','latex','FontSize',14)
lgd = legend({'Active-set QP', 'PSO', 'mPSO', 'Reference', 'Limits'},'Interpreter','latex', ...
        'FontSize',14);
lgd.NumColumns = 5;
% lgd.Location = 'southoutside';

%% MSE
index = theta_matlab > pi;
theta_matlab_tmp(index) = (pi*2 - theta_matlab(index));
theta_matlab_tmp(~index) = theta_matlab(~index);
theta_matlab_tmp = theta_matlab_tmp(:);

index = theta_pso > pi;
theta_pso_tmp(index) = (pi*2 - theta_pso(index));
theta_pso_tmp(~index) = theta_pso(~index);
theta_pso_tmp = theta_pso_tmp(:);

index = theta_mpso > pi;
theta_mpso_tmp(index) = (pi*2 - theta_mpso(index));
theta_mpso_tmp(~index) = theta_mpso(~index);
theta_mpso_tmp = theta_mpso_tmp(:);

mse = [
    sum((x_matlab - x_ref).^2)/n_step;
    sum((x_pso - x_ref).^2)/n_step;
    sum((x_mpso - x_ref).^2)/n_step;

    sum((theta_matlab_tmp - theta_ref).^2)/n_step;
    sum((theta_pso_tmp - theta_ref).^2)/n_step;
    sum((theta_mpso_tmp - theta_ref).^2)/n_step;

    sum((u_matlab - u_ref).^2)/n_step;
    sum((u_pso - u_ref).^2)/n_step;
    sum((u_mpso - u_ref).^2)/n_step
];
