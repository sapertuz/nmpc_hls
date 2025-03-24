% 
clear 
clc

path_sim_data = '../nmpc_pso_hls/matlab/';

%% Load data
run([path_sim_data, 'drone_kpso.m']);

n_step = 375;

x_max = Model.x_max;
x_min = Model.x_min;
u_max = Model.u_max;
u_min = Model.u_min;

Ts = Model.Ts;
x_mpso = NMPC_SIM.state_history(1:n_step,:);
u_mpso = NMPC_SIM.control_history(1:n_step,:);
xref = NMPC_SIM.xref(1:n_step,:);
uref = NMPC_SIM.uref(1:n_step,:);

run([path_sim_data, 'drone_pso.m']);
x_pso = NMPC_SIM.state_history(1:n_step,:);
u_pso = NMPC_SIM.control_history(1:n_step,:);

load drone_spiral_dir.mat
x_matlab = xHistory(1:n_step,:);
u_matlab = uHistory(1:n_step,:);

time = 0:Ts:Model.SimulationTime;
time = time(1:n_step);

% axis_lim = [-5 40 -15 15 -1 45];
axis_lim = [3.5 10 9 16 -1 3];

col = '#696969';	 %[1 1 1]*.15;
col2 = '#708090';	 %[1 1 1]*.15;

%% Vizualization
% close all
figure(1), clf

subplot(4,2,[1:4])
    axis(axis_lim), view([30 30])
    title('Drone Position','interpreter','latex','FontSize',14)
    grid on, hold on
    plot3(x_matlab(:,1), x_matlab(:,2), x_matlab(:,3),'-','Color',col,'linewidth',1)
    plot3(x_pso(:,1), x_pso(:,2), x_pso(:,3),'--','Color',col2,'linewidth',1)
    plot3(x_mpso(:,1), x_mpso(:,2), x_mpso(:,3),'-k','linewidth',1.5)
    plot3(xref(:,1), xref(:,2), xref(:,3),'--b','linewidth',2)
    xlabel('x [m]','interpreter','latex','FontSize',14)
    ylabel('y [m]','interpreter','latex','FontSize',14)
    zlabel('z [m]','interpreter','latex','FontSize',14)
% subplot(4,2,5)
%     grid on, hold on
%     plot(time, x_matlab(:,6),'-','Color',col,'linewidth',1)
%     plot(time, x_pso(:,6),'--','Color',col2,'linewidth',1)
%     plot(time, x_kpso(:,6),'-k','linewidth',1)
%     plot(time, xref(:,6),'--b','linewidth',1)
axis_2d = [-inf inf -110 110];
subplot(4,2,5)
    grid on, hold on
    plot(time, u_matlab(:,2),'-','Color',col,'linewidth',1)
    plot(time, u_pso(:,2),'--','Color',col2)
    plot(time, u_mpso(:,2),'-k','linewidth',1)
    plot(time, uref(:,2),'--b','linewidth',1)
    plot(time, repmat(u_max(1),n_step,1), '-.r')
    plot(time, repmat(u_min(1),n_step,1), '-.r')
    axis(axis_2d)
    ylabel('Pitch [\%]','interpreter','latex','FontSize',14)
subplot(4,2,6)
    grid on, hold on
    plot(time, u_matlab(:,3),'-','Color',col,'linewidth',1)
    plot(time, u_pso(:,3),'--','Color',col2)
    plot(time, u_mpso(:,3),'-k','linewidth',1)
    plot(time, uref(:,3),'--b','linewidth',1)
    plot(time, repmat(u_max(1),n_step,1), '-.r')
    plot(time, repmat(u_min(1),n_step,1), '-.r')
    axis(axis_2d)
    ylabel('Roll [\%]','interpreter','latex','FontSize',14)
subplot(4,2,7)
    grid on, hold on
    plot(time, u_matlab(:,1),'-','Color',col,'linewidth',1)
    plot(time, u_pso(:,1),'--','Color',col2)
    plot(time, u_mpso(:,1),'-k','linewidth',1)
    plot(time, uref(:,1),'--b','linewidth',1)
    plot(time, repmat(u_max(1),n_step,1), '-.r')
    plot(time, repmat(u_min(1),n_step,1), '-.r')
    axis(axis_2d)
    ylabel('Throttle [\%]','interpreter','latex','FontSize',14)
    xlabel('Time [s]','interpreter','latex','FontSize',14)
subplot(4,2,8)
    grid on, hold on
    plot(time, u_matlab(:,4),'-','Color',col,'linewidth',1)
    plot(time, u_pso(:,4),'--','Color',col2)
    plot(time, u_mpso(:,4),'-k','linewidth',1)
    plot(time, uref(:,4),'--b','linewidth',1)
    plot(time, repmat(u_max(1),n_step,1), '-.r')
    plot(time, repmat(u_min(1),n_step,1), '-.r')
    axis(axis_2d)
    ylabel('Yaw [\%]','interpreter','latex','FontSize',14)
    xlabel('Time [s]','interpreter','latex','FontSize',14)

lgd = legend({'Active-set QP', 'PSO', 'mPSO', 'Reference', 'Limits'},'Interpreter','latex', ...
        'FontSize',14);
lgd.NumColumns = 5;

%% MSE
uref = uref(:,1:4);

mse = [
    sum(sum((x_matlab(:,1:3) - xref(:,1:3)).^2))/n_step,...
    sum(sum((x_pso(:,1:3) - xref(:,1:3)).^2))/n_step,...
    sum(sum((x_mpso(:,1:3) - xref(:,1:3)).^2))/n_step,...
...
    sum(sum((x_matlab(:,4:6) - xref(:,4:6)).^2))/n_step,...
    sum(sum((x_pso(:,4:6) - xref(:,4:6)).^2))/n_step,...
    sum(sum((x_mpso(:,4:6) - xref(:,4:6)).^2))/n_step,...
...
    sum(sum((u_matlab - uref).^2))/n_step,...
    sum(sum((u_pso - uref).^2))/n_step,...
    sum(sum((u_mpso - uref).^2))/n_step
];

