clear
clc

addpath ../

%% Parameters
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

%%
n_U = 4;
Ts = 0.1;
xmeasure      = [7,10,0,0,0,0,0,0,0,0,0,0];

M_mma_inv = eye(4)/(drone.M_mma);

axis_lim = [4 10 9 17 -1 3];

Ixx = drone.I(1);
Iyy = drone.I(2);
Izz = drone.I(3);
I = drone.I;
l  = drone.l;
kf = drone.kf;
kM = drone.kM;
m  = drone.m;
g  = drone.g;
b  = drone.b;

u_min =    drone.u_min;
u_max =    drone.u_max;
norm_min = drone.norm_min;
norm_max = drone.norm_max;

%%

load ../data_line.mat
init = xHistory(1,:);
clear xHistory
%Input limits


xHistory(1,:) = init;

for i=1:size(uHistory,1)
%     control = (uHistory(i,:)' - norm_min) .* (u_max - u_min) ./ (norm_max - norm_min) + u_min;
%     state_dot(i,:) = model_drone(xHistory(i,:), control,Ixx,Iyy,Izz,kM,kf,l,m,g,b,M_mma_inv,n_U)

    state_dot(i,:) = drone_model(xHistory(i,:), uHistory(i,:));

    xHistory(i+1,:) = xHistory(i,:) + state_dot(i,:)*Ts;
end

%% Vizualization
% close all
figure(1), clf
for i=1:length(time)-1
    animateQuadrotor(time(i), xHistory(i,:), uHistory(i,:), xRef_out, axis_lim);
    pause(Ts);    
end

rmpath ../