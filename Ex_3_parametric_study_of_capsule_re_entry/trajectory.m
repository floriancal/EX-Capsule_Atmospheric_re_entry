function [Qdot, Heatload, gload, isimpact, timebeforeimpact, Dt] = trajectory(fpa, Cd, v0, m, Cl, i0, h, Scross, Rayoncourburecaps)
%% REENTRY TRAJECTORY SIMULATION
% This function simulates the atmospheric reentry of a vehicle, computing 
% key aerodynamic and thermal parameters along the trajectory.
%
% Inputs:
%   - fpa: Flight path angle at entry (degrees)
%   - Cd: Drag coefficient
%   - v0: Initial velocity (m/s)
%   - m: Vehicle mass (kg)
%   - Cl: Lift coefficient (typically zero for ballistic entry)
%   - i0: Initial angle of attack (degrees)
%   - h: Initial altitude (meters)
%   - Scross: Cross-sectional area of the vehicle (m²)
%   - Rayoncourburecaps: Radius of curvature for the capsule (meters)
%
% Outputs:
%   - Qdot: Heat flux during reentry (W/m²)
%   - Heatload: Total heat load accumulated during reentry (J/m²)
%   - gload: Acceleration experienced (in g)
%   - isimpact: Flag indicating ground impact (1 if impact occurs)
%   - timebeforeimpact: Time before ground impact (s)
%   - Dt: Simulation time step (s)
%
% Simulation steps:
% 1. **Initialization**: Defines initial conditions (altitude, velocity, etc.).
% 2. **Force Computation**: Computes gravitational, aerodynamic, and lift forces.
% 3. **Integration**: Updates velocity and position at each time step.
% 4. **Thermal Analysis**: Estimates heat flux and accumulated heat load.
% 5. **Impact Detection**: Stops the simulation upon surface contact.
%
% Assumptions:
% - The atmosphere is modeled using the NRLMSISE-00 density model.
% - The Earth is a perfect sphere with constant gravity.
% - No control forces are applied (purely ballistic trajectory).

isimpact = 0; % Impact flag
Rm = 6378e3; % Earth radius (meters)
Mum = 3.98600e14; % Earth's gravitational parameter (m^3/s^2)

Dt = 1;   % Time step (seconds)
sizecalcul = 60000 / Dt;  % Number of iterations based on the time step

%% Preallocate memory for speed and force calculations
timevector = zeros(sizecalcul, 1);
grav_mag = zeros(sizecalcul, 1);
grav_x = zeros(sizecalcul, 1);
grav_y = zeros(sizecalcul, 1);
grav_z = zeros(sizecalcul, 1);
drag_mag = zeros(sizecalcul, 1);
drag_x = zeros(sizecalcul, 1);
drag_y = zeros(sizecalcul, 1);
drag_z = zeros(sizecalcul, 1);
lift_mag = zeros(sizecalcul, 1);
lift_x = zeros(sizecalcul, 1);
lift_y = zeros(sizecalcul, 1);
lift_z = zeros(sizecalcul, 1);
Fdx = zeros(sizecalcul, 1);
Fdy = zeros(sizecalcul, 1);
Fdz = zeros(sizecalcul, 1);
vit_x = zeros(sizecalcul, 1);
vit_y = zeros(sizecalcul, 1);
vit_z = zeros(sizecalcul, 1);
vit_mag = zeros(sizecalcul, 1);
pos_x = zeros(sizecalcul, 1);
pos_y = zeros(sizecalcul, 1);
pos_z = zeros(sizecalcul, 1);
Qdot = zeros(sizecalcul, 1);
net_x = zeros(sizecalcul, 1);
net_y = zeros(sizecalcul, 1);
net_z = zeros(sizecalcul, 1);
net_mag = zeros(sizecalcul, 1);
pos_mag = zeros(sizecalcul, 1);
Heatload = zeros(sizecalcul, 1);
Altitudevector = zeros(sizecalcul, 1);
gload = zeros(sizecalcul, 1);
drag_theta = zeros(sizecalcul, 1);
drag_phi = zeros(sizecalcul, 1);
vit_theta = zeros(sizecalcul, 1);
vit_phi = zeros(sizecalcul, 1);
pos_theta = zeros(sizecalcul, 1);
pos_phi = zeros(sizecalcul, 1);
grav_theta = zeros(sizecalcul, 1);
grav_phi = zeros(sizecalcul, 1);

%% Initial conditions
pos_mag(1) = Rm + h;

pos_theta(1) = 0;
pos_phi(1) = 0;
[pos_x(1), pos_y(1), pos_z(1)] = sph2cart(pos_theta(1), pos_phi(1), pos_mag(1));

vit_theta(1)=(pi/2+(pi*fpa/180));
vit_phi(1)=pi*i0/180;
vit_mag(1) = v0;
[vit_x(1), vit_y(1), vit_z(1)] = sph2cart(vit_theta(1), vit_phi(1), vit_mag(1));

grav_mag(1) = -(Mum * m) / (pos_mag(1)^2);
Qdot(1) = 0;
Altitudevector(1) = h;

%% Main loop
for i = 2:sizecalcul
    timevector(i) = i * Dt;
    
    % Atmospheric conditions (density)
    if (pos_mag(i - 1) - Rm) < 100000 && (pos_mag(i - 1) - Rm) > -100000
        [~, rho] = atmosnrlmsise00((pos_mag(i - 1) - Rm), 0, 0, 2007, 1, 0);
        rho = rho(6); % Extract density
    else
        rho = 0;
    end
   
    % Gravitational force calculation
    grav_theta(i) = pos_theta(i-1);
    grav_phi(i) = pos_phi(i-1);
    grav_mag(i) = -(Mum * m) / (pos_mag(i - 1)^2);
    [grav_x(i), grav_y(i), grav_z(i)] = sph2cart(grav_theta(i), grav_phi(i), grav_mag(i));
    
    % Drag force calculation
    drag_theta(i)=vit_theta(i-1);
    drag_phi(i)=vit_phi(i-1);
    drag_mag(i)=-1/2*rho*Scross*(vit_mag(i-1)^2)*Cd;
    [drag_x(i),drag_y(i),drag_z(i)]=sph2cart(drag_theta(i),drag_phi(i),drag_mag(i));
    
    % Lift force calculation
    lift_mag(i) = 0.5 * rho * Scross * (vit_mag(i - 1)^2) * Cl;  % Lift coefficient is 0 (not used)
    [lift_x(i), lift_y(i), lift_z(i)] = sph2cart(0, 0, lift_mag(i)); % 0 lift 
    
    % Net forces (gravity + drag + lift)
    Fdx(i) = lift_x(i) + drag_x(i);
    Fdy(i) = lift_y(i) + drag_y(i);
    Fdz(i) = lift_z(i) + drag_z(i);
    
    net_x(i) = grav_x(i) + Fdx(i);
    net_y(i) = grav_y(i) + Fdy(i);
    net_z(i) = grav_z(i) + Fdz(i);
    [~, ~, net_mag(i)] = cart2sph(net_x(i), net_y(i), net_z(i));
    
    % Velocity update
    vit_x(i) = vit_x(i - 1) + (net_x(i) / m) * Dt;
    vit_y(i) = vit_y(i - 1) + (net_y(i) / m) * Dt;
    vit_z(i) = vit_z(i - 1) + (net_z(i) / m) * Dt;
    [vit_theta(i),vit_phi(i),vit_mag(i)]=cart2sph(vit_x(i),vit_y(i),vit_z(i));
    
    % Position update
    pos_x(i) = pos_x(i - 1) + vit_x(i);
    pos_y(i) = pos_y(i - 1) + vit_y(i);
    pos_z(i) = pos_z(i - 1) + vit_z(i);
    [pos_theta(i), pos_phi(i), pos_mag(i)] = cart2sph(pos_x(i), pos_y(i), pos_z(i));
    
    % Heat flux and load
    Qdot(i) = 1.9027e-4 * sqrt(rho / Rayoncourburecaps) * vit_mag(i - 1)^3;
    Heatload(i) = Qdot(i) * (i * Dt);
    Altitudevector(i) = pos_mag(i - 1) - Rm;
    gload(i) = net_mag(i) / (9.80665 * m);
    
    % Impact condition
    if pos_mag(i) < Rm
        timebeforeimpact = i * Dt;
        isimpact = 1;
        break
    end
end

%% Post-processing: Trim excess data after impact
if isimpact == 1
    idx = find(pos_mag < Rm+1100, 1);
    Heatload(idx:end) = [];
    gload(idx:end) = [];
    Qdot(idx:end) = [];
else
    timebeforeimpact = 0;
end
end
