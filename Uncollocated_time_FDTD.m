% Uncollocated_FDTD.m
%
% Author : Derek Swanson
% Date   : 5 July 2011
%
% This program uses the Uncollocated FDTD technique to solve for the
% voltage and current along a 1-D transmission line

clear; close all;
scrnsz = get(0,'ScreenSize');
figure('Position',[scrnsz(1) scrnsz(2) scrnsz(3) scrnsz(4)])

ustep=@(t) 0.5*(sign(t)+1);
pulse=@(t) ustep(t+.5) - ustep(t-.5);

% Transmission Line Parameters
R = 0;
G = 0;
L = 250e-9;
C = 1e-9;

% Simulation Parameters
M = 100;  % Number of Nodes
N = 1000; % Number of Time Steps
Len = 1;
A = 1;

% Boundary Conditions
RS = 10;
RL = sqrt(L/C)/100;

% Wave Propegation Speed
up = sqrt(1/L/C);

% Distance between adjacent nodes
dz = Len/(M-1);


% "Magic" time step (Courant-Fredrichs-Lewy stability requirement)
dt = dz/up;

% z vector for plotting
z = linspace(-Len,0,2*M);
t = 0:dt:(N-1)*dt;

% % Sinusoidal Source
% f = 400e6;
% periods = 2;
% vg = A*sin(f*2*pi*t).*(1-ustep(t - 1/f*periods));

% Pulse Source
vg = 1-ustep(t - 2.5e-9);

% % DC Source
% vg = ones(length(t));
% vg = vg*A;

% % Sawtooth Source
% f = 400e6;
% periods = 2;
% vg = (A+A*sawtooth(f*2*pi*t)).*(1-ustep(t-1/f*periods));

% % Triangle Wave Source
% f = 400e6;
% periods = 2;
% vg = (A+A*sawtooth(f*2*pi*t,.5)).*(1-ustep(t-1/f*periods));

% % Square Wave Source
% f = 400e6;
% periods = 2;
% vg = (A*square(f*2*pi*t)).*(1-ustep(t-1/f*periods));


% Initial Conditions
v  = zeros(1,2*M);
i  = zeros(1,2*M);
vn = zeros(1,2*M);
in = zeros(1,2*M);

for n = 1:N

    % generate plots
    subplot(2,1,1);
    plot(z,v)
    axis([-Len,0,-2,2]);
    xlabel('distance (m)');
    ylabel('voltage (V)');
    title(sprintf('t = %.3f ns',(n-1)*dt*1e9));
    subplot(2,1,2);
    plot(z,i)
    axis([-Len,0,-.3,.3]);
    xlabel('distance (m)');
    ylabel('current (A)');
    shg;
    
    % Update Voltages
    for m = 2:2*M-1
        vn(m) = (v(m)*(C/dt - G/2) + (i(m-1) - i(m+1))/dz)/(C/dt + G/2);
    end
    
    % Boundary conditions
    vn(1) = vg(n);
    vn(2*M) = (-2*vn(2*M-1)*RL/dz + i(2*M)*(R*RL - 2*RL*L/dt))/(-2*RL/dz - 2*L/dt);
    
    % Update Currents
    for m = 2:2*M-1
        in(m) = (i(m)*(L/dt - R/2) + (vn(m-1) - vn(m+1))/dz)/(R/2 + L/dt);
    end
    
    % Boundary Conditions
    in(1) = in(2);
    in(2*M) = (-2*in(2*M-1)/dz + vn(2*M)*(G - 2*C/dt))/(-2/dz - 2*C*RL/dt);
    
    % Update arrays
    v = vn;
    i = in;

end


