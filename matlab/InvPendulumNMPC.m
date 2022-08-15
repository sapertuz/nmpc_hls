function [time, xHistory, uHistory, xRef_out] = InvPendulumNMPC()
nx = 4; % 4 states [x,x_dot,th,th_dot]
ny = 4; % 4 states [x,x_dot,th,th_dot]
nu = 1;  % 1 input
tot_time      = 80;

%% NMPC configuration data
Ts            = 0.1;        % Sampling interval
mpciterations = tot_time/Ts;% Number of MPC iterations to be performed
N_h           = round(3/Ts);% Length of optimization horizon
tmeasure      = 0.0;        % Time measurement of initial value
                            % State measurement of initial value
xmeasure      = [0,0,pi,0];
u0            = zeros(nu,N_h); % Initial guess of open loop control

 %% NMPC implementation
[time, xHistory, uHistory] = nmpc(@runningcosts, @terminalcosts, @constraints, ...
         @terminalconstraints, @linearconstraints, @InvPendulumSystem, ...
         mpciterations, N_h, Ts, tmeasure, xmeasure, u0);
xRef_out = InvPendulum_xref(time);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the NMPC functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cost = runningcosts(t, x, u)    
    global invPend;
    
    % Input weights
    u_w = [0.001];
    % Output weights
    x_w = [1, 0.015, 1, 0];
    
    x_ref_out = InvPendulum_xref(t, u);
    
    penality = ones(1,4);
    index = x > invPend.state_upper_limits ;
    penality(index) = 1e4;
    index = x < invPend.state_lower_limits;
    penality(index) = 1e4;
    
    % Calculates from reference trayectory and control input ideally is 0
    J_x = x_w.*(x - x_ref_out);
    if ((x(3) > pi) && (x_ref_out(3) < pi))
        J_x(3) = ((pi*2 - x(3)) - x_ref_out(3));
    end
    J_x = penality.*J_x;
    
    cost = sum((J_x).^2) + sum(u_w.*u.^2);

function cost = terminalcosts(t, x)
    cost = 0.0;

function [c,ceq] = constraints(t, x, u)
    c   = [];
    ceq = [];

function [c,ceq] = terminalconstraints(t, x)
    c   = [];
    ceq = [];

function [A, b, Aeq, beq, lb, ub] = linearconstraints(t, x, u)
    global invPend;
    
    A   = [];
    b   = [];
    Aeq = [];
    beq = [];
    lb = invPend.u_min;
    ub = invPend.u_max;

