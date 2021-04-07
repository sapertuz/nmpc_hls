function y = system_drone(~, x, u, Ts)
    global  counter
    if isempty( counter )
        counter=0; %Initializing counter
    end
    counter = counter + 1;
    
    u = u(:);
    global drone
    u_min = drone.u_min; u_max = drone.u_max;
    norm_min = drone.norm_min; norm_max = drone.norm_max;
    u_real = map(u, norm_min, norm_max, u_min, u_max);
    w = mma(u_real);
    y = QuadrotorStateFcn(x',w);
    if (nargin>3)
        y = x + (y*Ts)';
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = mma(u)
%MMA Motor Mixing algotithm. This is a simple motor mixing algorithm that
%can convert between the intuitive roll, pitch, yaw, and thrust, and the
%less intuitive motor speeds.
%   W = MMA(CMD) returns the squared motor speeds in U with CMD commands
% 
%INPUT ARGUMENTS
%   U     array of commands in this format: [thrust, roll, pitch, yaw]
% 
%OUTPUT ARGUMENTS
%   W       motor speeds array of quadrotor from -1 to 1 in this format: 
%           [m_1, m_2, m_3, m_4]

    global drone

    kf = drone.kf ;  % Thrust (lift) coefficient [Ns^2]
    kM = drone.kM;   % Moment (drag) coefficient [Nms^2]
    l = drone.l;    % Arm length [m]
    b = drone.b; % Speed offset [rad^2/sec^2]
    M_mma = drone.M_mma;

    u = u(:);

    param = [1/kf, 1/(l*kf), 1/(l*kf), 1/kM];
    w = inv(M_mma)*(u.*param') + b;

function u = Inv_mma(w)
%MMA Motor Mixing algotithm. This is a simple motor mixing algorithm that
%can convert between motor speeds to roll, pitch, yaw, and thrust.
%   U = MMA(CMD) returns roll, pitch, yaw, and thrust
% 
%INPUT ARGUMENTS
%   W     motor speeds array of quadrotor from -1 to 1 in this format: 
%           [m_1, m_2, m_3, m_4]
% 
%OUTPUT ARGUMENTS
%   U       array of commands in this format: [thrust, roll, pitch, yaw]

    global drone

    kf = drone.kf ;  % Thrust (lift) coefficient [Ns^2]
    kM = drone.kM;   % Moment (drag) coefficient [Nms^2]
    l = drone.l;    % Arm length [m]
    b = drone.b; % Speed offset [rad^2/sec^2]
    M_mma = drone.M_mma;

    w = w(:);
    w = w - b;

    param = [kf, l*kf, l*kf, kM];
    u = (M_mma*w).*param';
    
function output = map(input, in_min, in_max, out_min, out_max)
    output = (input - in_min) .* (out_max - out_min) ./ (in_max - in_min) + out_min;