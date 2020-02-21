% Uncollocated_FDTD.m
%
% Author : Derek Swanson
% Date   : 5 July 2011
%
% This program uses the Uncollocated FDTD technique to solve for the
% voltage and current along a 1-D transmission line

clear;

ustep=@(t) 0.5*(sign(t)+1);
pulse=@(t) ustep(t+.5) - ustep(t-.5);

% Transmission Line Parameters
R = 0;
G = 0;
C = 250e-9;
L = 1e-9;

% Simulation Parameters
M = 100;  % Number of Nodes
N = 1000; % Number of Time Steps
Len = 1;
A = 1;

% Boundary Conditions
RS = 0;
RL = 1e9;

% Wave Propegation Speed
up = sqrt(1/L/C);

% Distance between adjacent nodes
dz = Len/(M-1);

% "Magic" time step (Courant-Fredrichs-Lewy stability requirement)
dt = dz/up;

% z vector for plotting
z = -Len:dz:0;
t = 0:dt:(N-1)*dt;

% Sinusoidal Source
f = 400e6;
periods = 2;
vg = A*sin(f*2*pi*t).*(1-ustep(t - 1/f*periods));

% % Pulse Source
% vg = 1-ustep(t - 2.5e-9);

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
vn = zeros(1,M);
in = zeros(1,M-1);

for n = 1:N

    % generate plots
    subplot(2,1,1);
    plot(z,v)%,z,vv);
    axis([-Len,0,-2,2]);
    xlabel('distance (m)');
    ylabel('voltage (V)');
    title(sprintf('t = %.3f ns',(n-1)*dt*1e9));
    subplot(2,1,2);
    plot(z(1:M-1)+dt/2,i)%,z,ii);
    axis([-Len,0,-.3,.3]);
    xlabel('distance (m)');
    ylabel('current (A)');
    shg;

    % Update Currents
    for m = 1:M-1
        in(m) = dt/L*((v(m) - v(m+1))/dz - i(m)*R) + i(m);
    end
    
    % Boundary Condition at the Source
    
    %vn(1) = (v(1)*(-G*RS + C/dt*RS - 1/2/dz) + in(1)*RS/dz + vg(n)/dz)/(1/2/dz + C*RS/dt);
    vn(1) =(v(1)*(G*RS - C*RS/dt + 1/2/dz) - vg(n)/dz + in(1)*RS/dz)/(-1/2/dz - C/dt*RS);q 
    
    % Update Voltages
    for m = 1:M-2
        vn(m+1) = dt/C*((in(m) - in(m+1))/dz - v(m+1)*G) + v(m+1);
    end
    
    % Boundary Condition at the termination
    vn(M) = v(M)*(-G*RL + C*RL/dt - 1/2/dz)/(1/2/dz + C*RL/dt);
    
    % Update arrays
    v = vn;
    i = in;
    
    pause
    
end


