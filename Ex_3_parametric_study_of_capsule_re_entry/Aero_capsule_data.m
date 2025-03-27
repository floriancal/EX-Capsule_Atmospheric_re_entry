function [Resultante, Scross] = Aero_capsule_data(gamma, Rayoncourburecaps, Cutting_shape)
%% AERODYNAMIC COEFFICIENT COMPUTATION FOR A HEMISPHERICAL SHIELD
% This function calculates the aerodynamic forces on a hemispherical reentry shield
% under hypersonic flow conditions.
%
% Inputs:
%   - gamma: Specific heat ratio (e.g., 1.4 for air)
%   - Rayoncourburecaps: Radius of the hemispherical shield (in meters)
%   - Cutting_shape: Height at which the shield is cut off (in meters)
%
% Outputs:
%   - Resultante: Aerodynamic force components [Fx, Fy, Fz] (in Newtons)
%   - Scross: Cross-sectional area of the shield (in square meters)
%
% Methodology:
% 1. Computes the pressure coefficient (Cp) using hypersonic theory.
% 2. Constructs a mesh grid representation of the hemispherical shield.
% 3. Calculates surface normals and pressure coefficients for each triangular element.
% 4. Integrates aerodynamic forces over the shield surface.
%
% Assumptions:
% - The flow is hypersonic, allowing the use of a simplified Cp formulation.
% - The shield is perfectly hemispherical with a radius of 5.5 meters.
% - The aerodynamic forces are computed based on a normal pressure distribution.
%
% Dependencies:
% - None

% Calculate Cpmax (Pressure coefficient) for hypersonic flow (Mach -->
% +inf)
Cpmax = ((gamma+1)^2/(4*gamma))^(gamma/(gamma-1))*(4/(gamma+1));
% Define grid for the mesh
debut = -5; % Starting value for mesh grid
fin = 5; % Ending value for mesh grid
pas = 0.05; % Step size for mesh grid

% Create mesh grid for surface definition
[x1, y1] = meshgrid(debut:pas:fin, debut:pas:fin);
tri = delaunay(x1, y1);  % Delaunay triangulation to form the surface




% Initialize variables
j = 0;
 

vent = [0 0 -1]; % Wind direction vector (negative z-direction)

% Pre-allocate arrays for Cp and normal vector components
Cp1 = zeros(1, length(tri));
Cp2 = zeros(1, length(tri));
Cp3 = zeros(1, length(tri));
surface = zeros(1, length(tri));

% Preallocate z array
z = NaN(size(x1));

%% Code for calculating the aerodynamic forces
% Generate the shape of the hemispherical shield
for x = debut:pas:fin
    j = j + 1;
    i = 0;
    for y = debut:pas:fin
        i = i + 1;
        z(j, i) = real(sqrt(Rayoncourburecaps^2 - x^2 - y^2));  % Define the hemisphere shape
        if z(j, i) <= Cutting_shape 
            z(j, i) = 0;  % Cut off below the height of 4.5 meters
        end
    end
end

%% Compute cross-sectional area (Scross)
R_eff = sqrt(Rayoncourburecaps^2 - Cutting_shape^2);  % Efecie Scross radius
Scross = pi * R_eff^2; % Cross-sectional area

% Plot the mesh surface
figure;
hold on;
trimesh(tri, x1, y1, z);


%% Calculating the normal vectors and pressure coefficients
for i = 1:length(tri)
    % Extract the coordinates of the three vertices of each triangle
    P1 = [x1(tri(i, 1)), y1(tri(i, 1)), z(tri(i, 1))];
    P2 = [x1(tri(i, 2)), y1(tri(i, 2)), z(tri(i, 2))];
    P3 = [x1(tri(i, 3)), y1(tri(i, 3)), z(tri(i, 3))];
    
    % Calculate the vectors of the triangle edges
    %V1 = P2 - P3;
    V2 = P3 - P1;
    V3 = P2 - P1;
    
    % Compute the normal vector to the triangle
    n = cross(V3, V2);
    
    % Store the components of the normal vector
    n1 = n(1);
    n2 = n(2);
    n3 = n(3);

    % Compute the surface area of the triangle
    surface(i) = norm(n) / 2;
    
    % Calculate the angle of attack between the triangle and the wind direction
    alpha = (pi/2) - acos(dot(n, vent) / (norm(n) * norm(vent)));
    
    % If the triangle lies on the equatorial plane (z=0), there is no contribution to Cp
    if P1(3) == 0 && P2(3) == 0 && P3(3) == 0
        Cp1(i) = 0;
        Cp2(i) = 0;
        Cp3(i) = 0;
    else
        % Calculate the pressure coefficient for each vertex
        Cp1(i) = Cpmax * sin(alpha)^2 * (n1 / norm(n));
        Cp2(i) = Cpmax * sin(alpha)^2 * (n2 / norm(n));
        Cp3(i) = Cpmax * sin(alpha)^2 * (n3 / norm(n));
    end
end

% Calculate the total aerodynamic forces in each direction (x, y, z)
Resultante1 = sum(Cp1 .* surface) / Scross;
Resultante2 = sum(Cp2 .* surface) / Scross;
Resultante3 = sum(Cp3 .* surface) / Scross;




% Output the total aerodynamic forces as a vector
Resultante = [Resultante1, Resultante2, Resultante3];

end
