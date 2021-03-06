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
L = @(m) 100e-9*(1 - ustep(m - 500)) + 250e-9*ustep(m - 500);


% Transmission Line Parameters
R = 0;
G = 0;
%L = 250e-9;
C = 1e-9;

% Simulation Parameters
M = 1000;  % Number of Nodes
N = 2000; % Number of Time Steps
Len = 10;
A = 1;

% Boundary Conditions
RS = 0;
RL = 0;

% Wave Propegation Speed
up = sqrt(1/L(1)/C);

% Distance between adjacent nodes
dz = Len/(M-1);

% "Magic" time step (Courant-Fredrichs-Lewy stability requirement)
dt = dz/up;

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
% f = 400e7;
% periods = 2;
% vg = (A/2+A/2*sawtooth(f*2*pi*t)).*(1-ustep(t-1/f*periods));

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

vg(1) = 0;

for n = 2:N

    % generate plots
    subplot(2,1,1);
    plot(z,v)%,z,vv);
    axis([-Len,0,-1.5,1.5]);
    xlabel('distance (m)');
    ylabel('voltage (V)');
    title(sprintf('t = %.3f ns',(n-1)*dt*1e9));
    subplot(2,1,2);
    plot(z(1:M-1)+dz/2,i)%,z,ii);
    axis([-Len,0,-.3,.3]);
    xlabel('distance (m)');
    ylabel('current (A)');
    shg;

    % Update Currents
    for m = 1:M-1
        %in(m) = dt/L*((v(m) - v(m+1))/dz - i(m)*R) + i(m);
        in(m) = ((v(m) - v(m+1))/dz + i(m)*(-R/2 + L(m)/dt))/(R/2 + L(m)/dt);
    end
    

    %vn(1) = vg(n);
    
    
    % Update Voltages
    for m = 1:M-2
        %vn(m+1) = dt/C*((in(m) - in(m+1))/dz - v(m+1)*G) + v(m+1);
        vn(m+1) = ((in(m) - in(m+1))/dz + v(m+1)*(-G/2 + C/dt))/(G/2 + C/dt);
    end
    
    vn(1) = (v(1)*(1/dz + G*RS/2 - C*RS/dt) + 2*RS*in(1)/dz - vg(n-1)/dz - vg(n)/dz)/(-1/dz - G*RS/2 - C*RS/dt);
        
    % Boundary Condition at the termination
    vn(M) = (v(M)*(RL*G/2 - RL*C/dt + 1/dz) - 2*RL*in(M-1)/dz)/(-1/dz - RL*G/2 - RL*C/dt);
    
    % Update arrays
    v = vn;
    i = in;
    
end


