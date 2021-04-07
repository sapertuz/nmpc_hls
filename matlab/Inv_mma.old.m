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
end

