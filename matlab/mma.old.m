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
% 
%            (y)
%           
%            m_1
%             |
%             |
%    m_4 ---- o ---- m_2 (x)
%             |
%             |
%            m_3
% 
% 
global drone

kf = drone.kf ;  % Thrust (lift) coefficient [Ns^2]
kM = drone.kM;   % Moment (drag) coefficient [Nms^2]
l = drone.l;    % Arm length [m]
b = drone.b; % Speed offset [rad^2/sec^2]
M_mma = drone.M_mma;

u = u(:);
      
param = [1/kf, 1/(l*kf), 1/(l*kf), 1/kM];
w = inv(M_mma)*(u.*param') + b;
end

