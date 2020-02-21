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
RS = sqrt(L/C)/10;
RL = sqrt(L/C);

% Wave Propegation Speed
up = sqrt(1/L/C);

% Distance between adjacent nodes
dz = Len/(M-1);

% "Magic" time step (Courant-Fredrichs-Lewy stability requirement)
dt =  dz/up;

% z vector for plotting
z = -Len:dz:0;
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
v = zeros(1,M);
i = zeros(1,M-1);
va = zeros(1,M);
ia = zeros(1,M-1);
vn = zeros(1,M);
in = zeros(1,M-1);

vg(1) = 0;

for n = 1:N

    % generate plots
    subplot(2,1,1);
    plot(z,v)%,z,vv);
    axis([-Len,0,-2,2]);
    xlabel('distance (m)');
    ylabel('voltage (V)');
    title(sprintf('t = %.3f ns',(n-1)*dt*1e9));
    subplot(2,1,2);
    plot(z(1:M-1)+dz/2,i)%,z,ii);
    axis([-Len,0,-.3,.3]);
    xlabel('distance (m)');
    ylabel('current (A)');
    shg;

    for k = 1:2
    
    for m = 1:M-2
        
        vn(m+1) = dt/C*((ia(m) - ia(m+1))/dz - G*va(m+1)) + v(m+1);
        in(m)   = dt/L*((va(m) - va(m+1))/dz - R*ia(m)) + i(m);
        
    end
    
    in(M-1)   = dt/L*((va(M-1) - va(M))/dz - R*ia(M-1)) + i(M-1);
    
    % Boundary Condition at the termination
    vn(M) = (v(M)*(RL*G/2 - RL*C/dt + 1/dz) - 2*RL*ia(M-1)/dz)/(-1/dz - RL*G/2 - RL*C/dt);
    %vn(1) = (-C*v(1)*RS/dt - 2*vg(n)/dz + 2*in(1)*RS/dz)/(-2/dz - G*RS - RS*C/dt);
    %vn(1) = (-2*C*RS*va(1)/dz - 2*vg(n)/dz + 2*RS*in(1)/dz)/(-2/dz - 2*RS*C/dz - G*RS);
    vn(1) = ((vg(n) - va(1))/RS - ia(1) + C*dz*v(1)/2/dt)/(C*dz/2);
    
    % Update arrays
    v = va;
    i = ia;
    va = vn;
    ia = in;
    
    end
    
end


