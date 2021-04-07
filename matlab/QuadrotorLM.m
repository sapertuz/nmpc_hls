function x_out = QuadrotorLM(x_input, w)
%QUADROTORLM Linear model of quadrotor according to conference paper
%"Quadrotor Control" Master Project by Franchesco Sabatino.
%   X_OUT = QUADROTORLM(X,U) returns the 12 states of the system 
%                    [x,y,z,diff(x,y,z),phi,theta,psi,diff(phi,theta,psi)]
% 
%INPUT ARGUMENTS
%   X_INPUT array of input states
%   U       array of input commands [][thrusth, roll, pitch, yaw]
% 
%OUTPUT ARGUMENTS
%   X_OUT   array of output states 
%           
%       w_2       w_1
%          \     /
%           \   / 
%             o
%           /   \
%          /     \
%       w_3       w_4

% Parameters and Initial Conditions (Blade Element Theory)
I = [7.5e-3;7.5e-3;7.5e-3];         % Moment of Inertia [kg.m^2] 
l = .23;                            % Arm length [m]
Ir = 6e-5;                          % Inertia of Motor [kg.m^2]
kf = 3.13e-5;                       % Thrust (lift) coefficient [Ns^2]
kM = 7.5e-7;                        % Moment (drag) coefficient [Nms^2]
m = 0.65;                           % Mass of quadcopter [kg]
g = 9.81;                           % Gravity [m/s^2]
kt = [0.1 0 0; 0 0.1 0; 0 0 0.1];   % Aerodynamic thrust drag coef [Ns/m]
kr = [0.1 0 0; 0 0.1 0; 0 0 0.1];   % Aerodynamic moment drag coef [Nm.s]

% Input states
x =     x_input(1);
y =     x_input(2);
z =     x_input(3);
phi =   x_input(4);
theta = x_input(5);
psi =   x_input(6);
x_p =   x_input(7);
y_p =   x_input(8);
z_p =   x_input(9);
phi_p =   x_input(10);
theta_p = x_input(11);
psi_p =   x_input(12);

ft = u(1);
tx = u(2);
ty = u(3);
tz = u(4);

% Motor Mixing 
M_mma = [  1  1  1  1;
          -1 -1  1  1;
           1 -1 -1  1;
           1 -1  1 -1];

% Matrices
%   [   x   y   z   phi     theta   psi     xp      yp      zp      phip    thetap  psip]
A = [   0   0   0   0       0       0       1       0       0       0       0       0;  % xp
        0   0   0   0       0       0       0       1       0       0       0       0;  % yp
        0   0   0   0       0       0       0       0       1       0       0       0;  % zp
        0   0   0   0       0       0       0       0       0       1       0       0;  % phip
        0   0   0   0       0       0       0       0       0       0       1       0;  % thetap
        0   0   0   0       0       0       0       0       0       0       0       1;  % psip
        0   0   0	-g      0       0       0       0       0       0       0       0;  % xpp
        0	0	0	0       g       0       0       0       0       0       0       0;  % ypp
        0	0	0   0       0       0       0       0       0       0       0       0;  % zpp
        0   0   0   0       0       0       0       0       0       0       0       0;  % phipp
        0   0   0   0       0       0       0       0       0       0       0       0;  % thetapp
        0   0   0   0       0       0       0       0       0       0       0       0]  % psipp
%   [   ft  tx  ty  tz ]
B = [   0   0   0   0;
        0   0   0   0;
        0   0   0   0;
        0   0   0   0;
        0   0   0   0;
        0   0   0   0;
        0   0   0   0;
        0   0   0   0;
        0   0   0   0;
        0   0   0   0;
        0   0   0   0;
        0   0   0   0];