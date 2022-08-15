%NonLinear model predivtive control for Quadrotor
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
clear 
clc

calc_state = 1;

%% Drone model constants and data definition
global drone
drone.I  = [1.2, 1.2, 2.3];	% Moment of Inertia (Ixx Iyy Izz)[kg.m^2] 
drone.l  = .25;             % Arm length [m]
drone.kf = 1;               % Thrust (lift) coefficient [Ns^2]
drone.kM = 0.2;             % Moment (drag) coefficient [Nms^2]
drone.m  = 2;               % Mass of quadcopter [kg]
drone.g  = 9.81;            % Gravity [m/s^2]
drone.b  = 4.9050;          % Speed offset [rad^2/sec^2]4.9050
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Motor Mixin Algorithm matrix
%                  w1   w2   w3   w4
drone.M_mma = [     1    1    1    1; % thrust (T)
                    0   -1    0    1; % roll  (tx)
                   -1    0    1    0; % pitch (ty)
                   -1    1   -1    1];% yaw   (tz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  T    tx   ty   tz
%Input limits
drone.u_min =    [-15;  -3;  -3;  -3];
drone.u_max =    [ 15;   3;   3;   3];
%Normalization of inputs
drone.norm_min= [-100;-100;-100;-100];
drone.norm_max= [ 100; 100; 100; 100];

Ts = 0.05;
xmeasure      = [7,10,0,0,0,0,0,0,0,0,0,0];

%% Either perform new NMPC drone implementation
[time, xHistory, uHistory, xRef_out] = QuadrotorNMPC();
     
%% Or load previously generated data
% load data_line.mat
% load data_spiral.mat
% load data_ring.mat
% load data_ring_noInputWeights.mat

%% Or load from c data
% run ../nmpc_pso_xhls/matlab/nmpc_sim.m
% Ts = 0.1;
% uHistory = NMPC_SIM.control_history;
% xHistory = NMPC_SIM.state_history;
% xRef_out = NMPC_SIM.xref;

%% 
% if calc_state
% clear xHistory
% xHistory(1,:) = xmeasure;
% for i=1:size(time)
%     xHistory(i+1,:) = QuadrotorSystem(0, xHistory(i,:), uHistory(i,:), Ts);
% end
% end

%% Vizualization
% close all
vidfile = VideoWriter('testmovie.mp4','MPEG-4');
vidfile.FrameRate = 10;
open(vidfile);
figure(1), clf
for i=1:size(time)
    animateQuadrotor(time(i), xHistory(i,:), uHistory(i,:), xRef_out);
    set(gcf, 'Position', [10 10 1350 900]);
    F = getframe(figure(1));
    writeVideo(vidfile, F);
    pause(Ts);    
end
close(vidfile)