% Name   : fdtline.m
% Author : Derek Swanson
% Date   : 24 MAY 2011
% 
% Description:
%   Simulates the transient response of a transmission line using the
%   finite difference time domain (FDTD) technique.

clear;

% Line inputs
Len=1;%input('Length = ');
M=100;%input('Nodes = ');
dz=Len/(M-1);
z=linspace(-Len,0,M);

% Line parameters 
C=1e-9;
L=250e-9;
R = 0;
G = 0;
A=1;%A = input('Amplitude = ');

% Velocity calc
up=1/sqrt(L*C);

% Time increment and no. of steps
dt=dz/up; %Courant condition: up*dt/dz <= 1

% Number of time steps and time vector
N=1000;%N=input('Time Steps = ');
tmax=(N-1)*dt;
t=linspace(0,tmax,N);

% Prompt the user to enter the boundary condition
RL = 1e90;%input('RL = ');
Rs = 0;%input('Rs = ');

% Prompt the user to select the voltage source
% display('Select voltage source:');
% display('  1 - 400MHz Sine Wave');
% display('  2 - 2.5ns pulse');
% display('  3 - DC');
source = 2;%input('Enter a number 1-3: ');

% Set up the voltage source
if(source == 1)
    % Sinusoidal
    vg = A*cos(400e6*2*pi*t);
elseif(source == 2)
    % Pulse
    vg = zeros(length(t));
    tt = 0;
    i = 1;
    while(tt <= 2.5e-9)
       vg(i) = A;
       tt = tt + dt;
       i = i + 1;
    end
else
    % DC
    vg = ones(length(t));
    vg = vg*A;
end

% Initial conditions
v=zeros(1,M); i=zeros(1,M);
vn=zeros(1,M); in=zeros(1,M);
vp=zeros(1,M); ip=zeros(1,M);

% Time loop t=0:dt:tmax
for n=1:N-1
    
    % generate plots
    subplot(2,1,1);
    plot(z,v);
    axis([-Len,0,-2,2]);
    xlabel('distance (m)');
    ylabel('voltage (V)');
    title(sprintf('t = %.3f ns',(n-2)*dt*1e9));
    subplot(2,1,2);
    plot(z,i);
    axis([-Len,0,-.3,.3]);
    xlabel('distance (m)');
    ylabel('current (A)');
%     subplot(3,1,3);
%     plot(z,v./i);
%     axis([-Len,0,-2*Rs,2*Rs]);
%     xlabel('distance (m)');
%     ylabel('input impedance (ohm)');
    shg;
    
    % "Next" values @ z=-Len+dz:dz:0
    for m=2:M-1
%       vn(m)=(v(m+1)+v(m-1))/2-dt/C*(i(m+1)-i(m-1))/2/dz;
% 	  in(m)=(i(m+1)+i(m-1))/2-dt/L*(v(m+1)-v(m-1))/2/dz;
        
%          in(m) = ip(m) + 2*dt/L*((v(m-1) - v(m+1))/2/dz - i(m)*R);
%          vn(m) = vp(m) + 2*dt/C*((i(m-1) - i(m+1))/2/dz - v(m)*G);

%       vn(m) = v(m) + dt/C*((i(m-1) - i(m+1))/2/dz - G*v(m));
%       in(m) = i(m) + dt/L*((v(m-1) - v(m+1))/2/dz - R*i(m));

%      vn(m) = dt/C*((i(m) - i(m+1))/dz - G*v(m)) + v(m);
%      in(m) = dt/L*((v(m) - v(m+1))/dz - R*i(m)) + i(m);




    end
    
    % Boundary condition at the source
    vn(1) = ((vg(n+1) - Rs*in(2))/dz + C*Rs/dt*v(1))/(1/dz + G*Rs + C*Rs/dt);
    in(1) = ((vg(n+1) - vn(2))/dz + L*i(1)/dt)/(Rs/dz + R + L/dt);
    
    % Boundary condition at the termination
    vn(M) = (RL*in(M-1)/dz + RL*C/dt*v(M))/(1/dz + G*RL + C*RL/dt);
    in(M) = (vn(M-1)/dz + L*i(M)/dt)/(RL/dz + R + L/dt);

%     in(M) =    ((-4*vn(M-1) + vn(M-2))/2/dz + L/2/dt*(-4*i(M) + ip(M)))/(-3*RL/2/dz - R - 3*L/2/dt);
%     vn(M) = RL*((-4*in(M-1) + in(M-2))/2/dz + C/2/dt*(-4*v(M) + vp(M)))/(-3/2/dz - G*RL -3*C*RL/2/dt);
    
%pause;    

    % Updating arrays
    vp = v;
    ip = i;
    v=vn;
    i=in;
    
end