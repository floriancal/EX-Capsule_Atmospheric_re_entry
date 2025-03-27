%% AERODYNAMIC ANALYSIS OF A SPHERICAL CAPSULE
% This script computes the aerodynamic force coefficient on a hemispherical 
% heat shield exposed to hypersonic flow. The main steps include:
%
% 1. **Triangle Base Setup**: Defines a reference triangular surface 
%    and computes its normal vector and area.
% 2. **Mesh Generation**: Creates a grid to represent the capsule surface.
% 3. **Aerodynamic Computation**:
%    - Determines the normal vectors for each surface element.
%    - Computes the local pressure coefficient (Cp) based on the angle of attack.
%    - Integrates the pressure forces to obtain the resultant aerodynamic force.
% 4. **Visualization**: Plots the triangulated surface mesh.
%
% Output:
% - **Resultant Force Coefficient**: Displays the non-dimensional force coefficient.
%
% Assumptions:
% - Hypersonic limit with constant maximum Cp.
% - Wind direction aligned with the -Z axis.
% - Perfectly hemispherical shape.

clear all; % Remove all variables from the workspace
clc; % Clear the command window

% Define wind direction
wind = [0 0 -1];

% Define maximum pressure coefficient
Cp_max = 1.83;

% Define grid for the 3D surface
start_point = -5;
end_point = 5;
step = 0.05;
[x1, y1] = meshgrid(start_point:step:end_point, start_point:step:end_point);
tri = delaunay(x1, y1);

% Define sphere radius
R = 5;

% Compute Z-coordinates for the spherical surface
z = NaN(size(x1));
j = 0;
for x = start_point:step:end_point
    j = j + 1;
    i = 0;
    for y = start_point:step:end_point
        i = i + 1;
        z(j, i) = real(sqrt(R^2 - x^2 - y^2));
    end 
end

% Plot the mesh surface
figure;
hold on;
trimesh(tri, x1, y1, z);

Cp = zeros(1,length(tri(:,1)));
for i = 1:length(tri)
    % Extract triangle vertices
    P1 = [x1(tri(i,1)), y1(tri(i,1)), z(tri(i,1))];
    P2 = [x1(tri(i,2)), y1(tri(i,2)), z(tri(i,2))];
    P3 = [x1(tri(i,3)), y1(tri(i,3)), z(tri(i,3))];
    
    % Compute triangle vectors
    V1 = P2 - P3;
    V2 = P3 - P1;
    V3 = P2 - P1;
    
    % Compute normal vector
    n = cross(V3, V2);
%     n1(i) = n(1);
%     n2(i) = n(2);
%     n3(i) = n(3);
    
    % Compute surface area
    surface(i) = norm(n) / 2;
    
    % Compute angle of attack
    alpha = (pi / 2) - acos(dot(n, wind) / (norm(n) * norm(wind)));
    alpha_deg = rad2deg(alpha);
    
    % Compute pressure coefficient
    if P1(3) == 0 && P2(3) == 0 && P3(3) == 0
        Cp(i) = 0;
    else
        Cp(i) = Cp_max * sin(alpha)^2;
    end 
end

% Compute resultant force
Resultant = sum(Cp .* surface) / (pi * R^2);

disp(['Resultant Force Coefficient: ', num2str(Resultant)]);

hold off;