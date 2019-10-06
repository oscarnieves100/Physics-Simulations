% This program simulates projectiles in 2D by using random wind forces
% (white noise)
clear all; clc; close all;

% Inputs
g = 9.8; % gravitational constant
theta0 = 15; % initial launch angle
v0 = 1; % initial launch speed
y0 = 0; % initial elevation
T_flight = 1.5*2*v0*sind(theta0)/g; % time of simulation
t = linspace(0,T_flight,1000); % time vector
m = 1; % mass of projectile
sigma = 0.75; % noise strength for wind forces
dt = t(2)-t(1); % step
Nt = length(t); 
Nsims = 100; % number of independent simulations (projectiles)
X = zeros(Nt,Nsims);
Y = X;

% Main loop
for nn = 1:Nsims
    % initialize random vectors, initial conditions
    Rx = sqrt(dt).*randn(Nt,1);
    Ry = sqrt(dt).*randn(Nt,1);
    vx(1) = v0*cosd(theta0);
    vy(1) = v0*sind(theta0);
    x(1) = 0;
    y(1) = 0;
    
    % Implement Euler-Mayurama scheme
    for kk = 1:Nt-1
        vx(kk+1) = vx(kk) - 1/m*sqrt(sigma)*Rx(kk);
        vy(kk+1) = vy(kk) - g*dt - 1/m*sqrt(sigma)*Ry(kk);
        x(kk+1) = x(kk) + dt*vx(kk);
        y(kk+1) = y(kk) + dt*vy(kk);
    end
    X(:,nn) = x;
    Y(:,nn) = y;
    disp(nn);
end

% Export to video
vid = VideoWriter([num2str(Nsims) ' projectiles random.mp4'],'MPEG-4');
open(vid);

figure(1);
set(gcf,'color','w');
axis tight manual
for kk = 1:10:Nt
    plot(X(1:kk,:),Y(1:kk,:),'LineWidth',2);
    axis([0,max(max(X)),0,max(max(Y))]);
    xlabel('x (m)');
    ylabel('y (m)');
    drawnow;
    frame = getframe(gcf); 
    writeVideo(vid,frame); % export frame to video
end
close(vid);