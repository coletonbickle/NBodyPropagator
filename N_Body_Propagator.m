%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title: N-Body Propagator
% Author: Coleton C. Bickle
%
% Description:
% This Script Numerically solves the N-Body problem using the Runge-Kutta
% Numerical Integration method. This script has no limit on the number of
% celestial body interactions. The performace of the scipt is a direct
% result on the number of inputted bodies. The bodies all have the same
% mass, but can be manually changed to investigate the interaction between
% many bodies of different masses.
%
% User Input:
%              N = Positive Integer Number of Bodies
%         ThreeD = 1 or 0
%             t0 = Initial Time of Simulation
%             tf = Final Time of Simulation
%
% Outputs:
%         Figure 1 = Projected Trajectory of Bodies
%         Figure 2 = Animated Plot of Trajectories.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preamble
% Clears and Resets Workspace
close all; clear; clc;

% Sets Interpreter to Latex for Plots
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% User Input
N = 3; % Number of Bodies
ThreeD = 0; % Use (1) To Plot in Three Dimensions
t0 = 0;
tf = 8000;

%% Constants
global G m
G = 6.67430e-20; % [km^3/kgs^2]
mass = 0.01*5.9722e22*ones(1,N);
m = mass;

%% Initial Conditions
rng(42);

initial = GetICs(N,ThreeD);
% initial = [-200 0 0 0 -0.1 0 200 0 0 0 0.1 0];

%% Integrator

[t,x] = OrbPert(@(t,y) NBodyPerts(t,y),[t0 tf],initial);

%% Plots N-Body Propagations

plotter(t, x, mass, ThreeD)


%-------------------------------Functions---------------------------------%
%% Initial Condition Function
function [y] = GetICs(N,ThreeD)
% This Function creates the initial conditions for the problem

y = [];

for j = 1:N
    if ThreeD == 1
        pos = 0.05*4e3*randn(1,3); vel = [0 0 0];%0.1*rand(1,3);
        y = [y pos vel];
    else
%         s = rng;
        pos = [0.05*4e3*randn(1,2) 0]; vel = [0 0 0];%[0.1*rand(1,2) 0]; %
        y = [y pos vel];
    end
end

end

%% Runge-Kutta Function
function [t,y] = OrbPert(f,tspan,y_0) 
% This Function is a numerical integrator that utilizes Runge-Kutta method
% to integrate the inputted Ordinary Differential Equation

if length(tspan) == 2
    h = round((tspan(end)-tspan(1))/250,2);
else
    h = tspan(2) - tspan(1);
end
t = tspan(1):h:tspan(end);

y = nan(numel(y_0),numel(t));

y(:,1) = y_0;

for i = 1:length(t)-1
    k1 = f(t(i),y(:,i));
    k2 = f(t(i)+.5*h, y(:,i)+h*.5*k1);
    k3 = f(t(i)+.5*h, y(:,i)+h*.5*k2);
    k4 = f(t(i)+h, y(:,i)+h*k3);

    y(:,i+1) = y(:,i) + 1/6*(k1 + 2*k2 + 2*k3 + k4)*h;
end
y = y';

end

%% N-Body Perturbation ODE Function
function [up] = NBodyPerts(t,u)
global G m

up = [];

% This function is Intended to be used in conjuntion with ODE113 to
% numerically integrate the acceleration caused by N-Body Perturbations
N = length(u)/6;
pos = zeros(N,3);
vel = pos;
for j = 1:N
    temp = u(6*j-5:6*j,:);

    pos(j,:) = temp(1:3,:);
    vel(j,:) = temp(4:6,:);
end

for i = 1:N
    a = zeros(1,3);
    for j = 1:N
        if i ~= j
            r = pos(j,:) - pos(i,:);
            r_mag = norm([r 100]);
            dvdt = G*m(j)/r_mag^3*r;
            a = a + dvdt;
        end
    end
    up = [up vel(i,:) a];
end
up = up';

end

%% Animated Plotting Function
function [] = plotter(t, x, mass, ThreeD)
% This Function Animates a Plot for the trajectory of the N-Body
% interactions
N = length(x(1,:))/6;
% Plot of Projected Trajectory
figure(1)
subplot(1,2,1)
Body = [];
hold on
for i = 1:N
    plot(x(:,i*6-5),x(:,i*6-4),'linewidth',1.1)
    Body = [Body "Body " + i];
end
title("N-Body Projected Trajectory")
legend(Body,'location','best')
axis equal
xlabel("X-Axis"); ylabel("Y-Axis")


dt = t(2) - t(1);
if ceil(t(end)/dt/250) > 1
    split = ceil(t(end)/dt/250);
else
    split = 1;
end

subplot(1,2,2)
for j = 1:split:length(x(:,1))
    x1 = []; y1 = x1; z1 = x1;
    if ThreeD == 1
        for k = 1:length(mass)

            if k == 1
                plot3(x(max(j-500,1):j,k*6-5),x(max(j-500,1):j,k*6-4),x(max(j-500,1):j,k*6-3),'color','#A9A9A9','LineWidth',1.25)
                hold on
            else
                hold on
                plot3(x(max(j-500,1):j,k*6-5),x(max(j-500,1):j,k*6-4),x(max(j-500,1):j,k*6-3),'color','#A9A9A9','LineWidth',1.25)
            end
            x1 = [x1 x(j,k*6-5)*mass(1)/sum(mass)]; y1 = [y1 x(j,k*6-4)*mass(1)/sum(mass)]; z1 = [z1 x(j,k*6-3)*mass(1)/sum(mass)];
            plot3(x(j,k*6-5),x(j,k*6-4),x(j,k*6-3),'.k','MarkerSize',10)
            hold off
        end
        com = [sum(x1) sum(y1) sum(z1)];
        axis square
        grid on
        axis([-5e2+com(1) 5e2+com(1) -5e2+com(2) 5e2+com(2) -5e2+com(3) 5e2+com(3)])
        title("3D N-Body Animated Trajectory")
        xlabel("X-Axis"); ylabel("Y-Axis"); zlabel("Z-Axis")
        drawnow
    else
        for k = 1:length(mass)

            if k == 1
                plot(x(max(j-500,1):j,k*6-5),x(max(j-500,1):j,k*6-4),'color','#A9A9A9','LineWidth',1.25)
                hold on
            else
                hold on
                plot(x(max(j-500,1):j,k*6-5),x(max(j-500,1):j,k*6-4),'color','#A9A9A9','LineWidth',1.25)
            end
            x1 = [x1 x(j,k*6-5)*mass(1)/sum(mass)]; y1 = [y1 x(j,k*6-4)*mass(1)/sum(mass)];
            plot(x(j,k*6-5),x(j,k*6-4),'.k','MarkerSize',10)
            hold off
        end
        com = [sum(x1) sum(y1)];
        axis square
        
%         axis([-5e2+com(1) 5e2+com(1) -5e2+com(2) 5e2+com(2)])
        axis equal
        title("2D N-Body Animated Trajectory")
        xlabel("X-Axis"); ylabel("Y-Axis")
        drawnow
    end
    
end

end
