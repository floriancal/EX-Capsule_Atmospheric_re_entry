function [epaisseur] = shield_thickness(Qdot, timebeforeimpact, Dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Heat Shield Thickness Calculation                    %
%                                                                           %
% This function calculates the thickness of a heat shield based on the      %
% instantaneous heat flux and temperature profile. The code simulates the    %
% thermal recession of the heat shield material (PICA3) as it experiences   %
% heat flux during re-entry, using the finite difference method to solve    %
% the temperature distribution and compute the resulting recession.         %
%                                                                           %
% Input Parameters:                                                         %
%   - Qdot: Vector of heat flux values (W/m^2)                              %
%   - timebeforeimpact: Time before impact in seconds                       %
%   - Dt: Time step for the simulation (seconds)                            %
%                                                                           %
% Output:                                                                   %
%   - epaisseur: Calculated thickness of the heat shield (in meters)        %
%                                                                           %
% Material properties for PICA3 are considered, including the maximum      %
% temperature (Tabla), specific heat capacity (cp), thermal conductivity    %
% (k), and density (rho). The function calculates the recession over time   %
% as the heat shield absorbs heat and begins to erode.                      %
%                                                                           %
% Author: [Florian Calatayud]                                                       %
% Date: [05/2018]                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Mission parameters and material properties for PICA3
Tabla = 2500; % K (Maximum temperature of PICA3)
cp = 1884;    % J/kg.K (Specific heat capacity of PICA3)
k = 0.72;     % W/m.K (Thermal conductivity of PICA3)
rho = 270;    % kg/m^3 (Density of PICA3)

% Initialization parameters
dt = 0.01;       % Time step in seconds
dx = 0.1 / 1000; % Spatial step size (m)

% Initial temperature setup
T0 = 300;        % Initial temperature (K)
n = timebeforeimpact; % Time before impact in seconds
K = k / (rho * cp);  % Thermal diffusivity
S = K * dt / dx^2;   % Stability criterion
spacestep = 1000;     % Number of spatial steps
timestep = n / dt;    % Number of time steps

% Preallocation of variables
Qdotreshape = zeros(timestep, 1);  % Reshaped Qdot
Recessinstantanee = zeros(timestep, 1);  % Instantaneous recession

% Reshape Qdot to match the time steps
coef = Dt / dt;  % Coefficient for Qdot reshaping
j = 1;
k1 = 0;
for i = 1:timestep
    k1 = k1 + 1;
    Qdotreshape(i) = Qdot(j);
    if k1 == coef
        j = j + 1;
        k1 = 0;
    end
    if j == length(Qdot)
        j = j - 1; % Prevent exceeding the length of Qdot
    end
end

% Initialize temperature matrix (Timestep x Space step)
a = ones(timestep + 2, spacestep) * T0;

% Main loop for solving the temperature profile
for i = 1:timestep
    % Check if temperature exceeds the limit and set to Tabla
    if a(i, 1) >= 2000
        a(i, 1) = Tabla;
        
    end
    
    % Calculate the temperature distribution using the finite difference method
    for j = 2:spacestep - 1
        a(i + 1, j) = S * (a(i, j + 1) + a(i, j - 1)) + (1 - 2 * S) * a(i, j);
    end
    
    % Apply heat flux at the boundary
    a(i + 2, 1) = (Qdotreshape(i)) * dx + a(i + 1, 2);
    a(i + 2, spacestep) = a(i + 1, 1);
    
    % Calculate recession when temperature reaches 1500K
    if a(i, 1) >= 1500
        Recessinstantanee(i) = 0.08;  % Recession rate in mm/s
    end
end

% Calculate total recession (in meters)
Recesstotale = sum(Recessinstantanee) * dt; % Total recession in mm
epaisseur = Recesstotale * 10^-3; % Convert to meters

end
