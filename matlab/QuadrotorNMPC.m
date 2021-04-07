function [time, xHistory, uHistory, xRef_out] = QuadrotorNMPC()
nx = 12; % 12 states [x,y,z,phi,theta,psi,... diff ant]
ny = 12; % 12 outputs [x,y,z,phi,theta,psi,... diff ant]
nu = 4;  % 4 inputs [w_1, w_2, w_3, w_4].^2 (velocities are squared)
tot_time      = 20;

%% NMPC configuration data
Ts            = 0.1;        % Sampling interval
mpciterations = tot_time/Ts;% Number of MPC iterations to be performed
N_h           = round(3/Ts);% Length of optimization horizon
tmeasure      = 0.0;        % Time measurement of initial value
                            % State measurement of initial value
xmeasure      = [7,10,0,0,0,0,0,0,0,0,0,0];
u0            = zeros(nu,N_h); % Initial guess of open loop control

 %% NMPC implementation
[time, xHistory, uHistory] = nmpc(@runningcosts, @terminalcosts, @constraints, ...
         @terminalconstraints, @linearconstraints, @system_drone, ...
         mpciterations, N_h, Ts, tmeasure, xmeasure, u0);
xRef_out = x_ref(time);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the NMPC functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cost = runningcosts(t, x, u)
    % Input weights
    u_w = [0.02; 0.00; 0.00; 0.00];
    % Output weights
    x_w = [1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0];
    
    x_ref_out = x_ref(t, u);
    
    % Calculates from reference trayectory and control input ideally is 0
    cost = sum((x_w.*(x - x_ref_out)).^2) + sum(u_w.*u.^2);

function cost = terminalcosts(t, x)
    cost = 0.0;

function [c,ceq] = constraints(t, x, u)
    c   = [];
    ceq = [];

function [c,ceq] = terminalconstraints(t, x)
    c   = [];
    ceq = [];

function [A, b, Aeq, beq, lb, ub] = linearconstraints(t, x, u)
    global drone;
    
    A   = [];
    b   = [];
    Aeq = [];
    beq = [];
    lb = drone.norm_min;
    ub = drone.norm_max;

