%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Mars Re-entry Simulation                            %
%                                                                           %
% This script simulates the re-entry trajectory of a spacecraft on Mars,    %
% considering the effects of gravity, drag, and lift forces. The dynamics   %
% are calculated using numerical integration, and the re-entry path is      %
% plotted in 3D.                                                            %
%                                                                           %
% Input Parameters:                                                         %
%   - Initial velocity (v0) in m/s                                          %
%   - Mass (m) of the spacecraft in kg                                      %
%   - Initial altitude (h) in meters                                        %
%   - Mars gravitational acceleration (gm) in m/s^2                         %
%   - Mars radius (Rm) in meters                                            %
%   - Mars gravitational constant (Mum) in m^3/s^2                          %
%   - Cross-sectional area (Scross) in m^2                                  %
%   - Drag coefficient (Cd)                                                 %
%   - Lift coefficient (Cl)                                                 %
%                                                                           %
% The code performs numerical integration using a time step (Dt) to update  %
% the spacecraft's position and velocity, calculating the forces at each    %
% time step. The trajectory is terminated when the spacecraft reaches the   %
% surface of Mars (impact condition).                                       %
%                                                                           %
% Author: [Florian Calatayud]                                               %
% Date: [05/2018]                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;  % Remove all variables from the workspace
clc;        % Clear the command window


clear all  % Remove all variables from the workspace
clc        % Clear the command window

% Initial parameters
i = 0;
v0 = 3.35e3;   % Initial velocity in m/s
m = 3;         % Mass in kg
h = 300e3;     % Initial altitude in m
gm = 3.711;    % Mars gravitational acceleration in m/s^2
Rm = 3.3895e6; % Mars radius in m
Mum = 42828e9; % Mars gravitational constant in m^3/s^2
Scross = 0.01; % Cross-sectional area in m^2
fpa = 0;       % Flight path angle in degrees
Cd = 2;        % Drag coefficient
Cl = 0;        % Lift coefficient (no lift in this case)

% Initial position in spherical coordinates
pos_theta(1) = 0;
pos_phi(1) = 0;
pos_mag(1) = Rm + h;
[pos_x(1), pos_y(1), pos_z(1)] = sph2cart(pos_theta(1), pos_phi(1), pos_mag(1));

% Initial velocity in spherical coordinates
vit_theta(1) = (pi/2 + (pi*fpa/180));
vit_phi(1) = i;
vit_mag(1) = v0;
[vit_x(1), vit_y(1), vit_z(1)] = sph2cart(vit_theta(1), vit_phi(1), vit_mag(1));

% Initial gravity calculation
grav_theta(1) = pos_theta(1);
grav_phi(1) = pos_phi(1);
grav_mag(1) = -(Mum*m) / (pos_mag(1)^2);

Dt = 1; % Time step in seconds

% Function for Mars atmospheric density based on altitude
rho = @(x) 0.02045 * exp(-x * 0.1038e-3);

% Simulation loop
i_max = 60000 / Dt;
for i = 2:i_max
    % Gravity force calculation
    grav_theta(i) = pos_theta(i-1);
    grav_phi(i) = pos_phi(i-1);
    grav_mag(i) = -(Mum*m) / (pos_mag(i-1)^2);
    [grav_x(i), grav_y(i), grav_z(i)] = sph2cart(grav_theta(i), grav_phi(i), grav_mag(i));
    
    % Drag force calculation
    drag_theta(i) = vit_theta(i-1);
    drag_phi(i) = vit_phi(i-1);
    drag_mag(i) = -0.5 * rho(pos_mag(i-1) - Rm) * Scross * (vit_mag(i-1)^2) * Cd;
    [drag_x(i), drag_y(i), drag_z(i)] = sph2cart(drag_theta(i), drag_phi(i), drag_mag(i));
    
    % Lift force calculation (here Cl=0, so lift_mag = 0)
    lift_theta(i) = vit_theta(i-1);
    lift_phi(i) = vit_phi(i-1);
    lift_mag(i) = 0.5 * rho(pos_mag(i-1) - Rm) * Scross * (vit_mag(i-1)^2) * Cl;
    [lift_x(i), lift_y(i), lift_z(i)] = sph2cart(lift_theta(i), lift_phi(i), lift_mag(i));
    
    % Total forces
    Fdx(i) = lift_x(i) + drag_x(i);
    Fdy(i) = lift_y(i) + drag_y(i);
    Fdz(i) = lift_z(i) + drag_z(i);
    
    % Net accelerations
    net_x(i) = grav_x(i) + Fdx(i);
    net_y(i) = grav_y(i) + Fdy(i);
    net_z(i) = grav_z(i) + Fdz(i);
    
    % Update velocity
    vit_x(i) = vit_x(i-1) + (net_x(i)/m) * Dt;
    vit_y(i) = vit_y(i-1) + (net_y(i)/m) * Dt;
    vit_z(i) = vit_z(i-1) + (net_z(i)/m) * Dt;
    [vit_theta(i), vit_phi(i), vit_mag(i)] = cart2sph(vit_x(i), vit_y(i), vit_z(i));
    
    % Update position
    pos_x(i) = pos_x(i-1) + vit_x(i);
    pos_y(i) = pos_y(i-1) + vit_y(i);
    pos_z(i) = pos_z(i-1) + vit_z(i);
    [pos_theta(i), pos_phi(i), pos_mag(i)] = cart2sph(pos_x(i), pos_y(i), pos_z(i));
    
    % Impact condition
    if pos_mag(i) < Rm
        disp('Impact')
        break;
    end
end

% Plot trajectory in 3D
figure;
hold on;
% Plot Mars as a sphere
[X, Y, Z] = sphere(50);
surf(Rm*X/1000, Rm*Y/1000, Rm*Z/1000, 'FaceColor', [1 0.3 0], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
colormap([1 0.3 0]);

% Plot re-entry trajectory
plot3(pos_x/1000, pos_y/1000, pos_z/1000, 'b', 'LineWidth', 1.5);
xlabel('x (km)');
ylabel('y (km)');
zlabel('z (km)');
title('Mars Re-entry Trajectory');
grid on;
hold off;
