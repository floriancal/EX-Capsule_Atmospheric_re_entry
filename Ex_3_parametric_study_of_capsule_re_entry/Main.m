%% MERIT DIAGRAM COMPUTATION FOR REENTRY TRAJECTORY
% This script evaluates key parameters of a reentry profile as a function 
% of the Flight Path Angle (FPA) for a simple capsule.
%
% Parameters analyzed:
% - Maximum heat flux (Qdot) during reentry
% - Total heat load experienced by the capsule
% - Maximum G-load encountered
% - Required thermal protection shield thickness (assuming PICA 3 material)
%
% Methodology:
% 1. Computes the aerodynamic drag coefficient of the capsule.
% 2. Simulates the reentry trajectory for different FPA values.
% 3. Extracts relevant performance parameters for each case.
% 4. Plots the merit diagrams to visualize trade-offs.
%
% Assumptions:
% - The aerodynamic coefficients are computed using a predefined function.
% - The shield thickness is estimated based on the heat flux and exposure time.
% - The analysis considers only no lift entry (Cl = 0).
%
% Dependencies:
% - `trajectory.m`: Simulates the reentry trajectory.
% - `Aerocoeff_capsule.m`: Computes aerodynamic coefficients and key datas.
% - `shield_thickness.m`: Estimates the required thermal protection thickness (ablative type).

% Author: [Florian Calatayud]                                                      
% Date: [05/2018] 


clear all 
clc
hold on 
warning('off', 'all');

%% Mission Parameters
Cl = 0;
v0 = 7.557e3; % Initial velocity in m/s
m = 6100;   % Mass in kg
% Other relevant Parameters are located in trajectorymerite.m and
% Aerocoeff_capsule.m

%% Mission Parameters (Constants)
i0 = 0.;  % Orbital angle (degrees)
h = 408e3;  % Altitude (meters)
Rayoncourburecaps = 5.5;  % Curvature radius of the capsule (m)
Cutting_shape = 4.5; % parameter to define capsule shape (m)


%% First we compute the drag coefficient of the capsule like in Example 2.
gamma = 1.4;% Typical value for air 
[Resultante, Scross] = Aero_capsule_data(gamma, Rayoncourburecaps, Cutting_shape);
Cd = Resultante(3);  % Aerodynamic drag coefficient 


%% Merit Diagram: FPA (in degrees)
% Define FPA range from 15 to 30 degrees
fpa_vec = linspace(0, 10, 20);

%Prealocation of key data 
Qdotmerite = zeros(length(fpa_vec),1);
Heatloadmerite = zeros(length(fpa_vec),1);
isimpactmerite = zeros(length(fpa_vec),1);
gloadmerite = zeros(length(fpa_vec),1);
epaisseurmerite  = zeros(length(fpa_vec),1);


var = 0;
% Loop over each FPA value
for fpa = fpa_vec
    var = var + 1;
    % Call the trajectory function to compute various parameters
    [Qdot, Heatload, gload, isimpact, timebeforeimpact, Dt] = trajectory(fpa, Cd, v0, m, Cl, i0, h, Scross, Rayoncourburecaps);
    
    % Store the results for merit diagram plotting
    Qdotmerite(var) = max(Qdot);
    Heatloadmerite(var) = Heatload(end);
    isimpactmerite(var) = isimpact;
    gloadmerite(var) = max(gload);
    
   %  Calculate shield thickness 
    epaisseurmerite(var) = shield_thickness(Qdot, timebeforeimpact, Dt);
    
end

figure
set(gcf, 'Color', 'k')  % Fond de la figure en noir

%% Subplot 1: Max Qdot and Heat Load
set(gcf, 'Color', 'k')  % Fond de la figure en noir

%% Subplot 1: Max Qdot and Heat Load
subplot(2,2,1)
set(gca, 'Color', 'k', 'XColor', 'w')  % Fond noir pour chaque subplot, axes X en blanc
yyaxis left
plot(fpa_vec, Qdotmerite, 'cyan', 'LineWidth', 2)  % Couleur cyan pour le premier graphe
ylabel('Max Qdot [W/m^2]', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'w')
set(gca, 'YColor', 'cyan')  % Couleur de l'axe gauche (cyan)

yyaxis right
plot(fpa_vec, Heatloadmerite, 'magenta', 'LineWidth', 2)  % Couleur magenta pour le deuxième graphe
ylabel('Heat Load [W/m^2]', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'w')
set(gca, 'YColor', 'magenta')  % Couleur de l'axe droit (magenta)

xlabel('FPA [deg]', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'w')
title('Qdot & Heat Load', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'w')
grid on
grid minor

%% Subplot 2: G-load
subplot(2,2,2)
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w')  % Fond noir pour chaque subplot, axes en blanc
plot(fpa_vec, gloadmerite, 'yellow', 'LineWidth', 2)  % Couleur jaune pour le graphe
ylabel('G-load', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'w')
xlabel('FPA [deg]', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'w')
title('G-load', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'w')
grid on
grid minor
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w')  

%% Subplot 3: Shield Thickness
subplot(2,2,3)
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w')  % Fond noir pour chaque subplot, axes en blanc
plot(fpa_vec, epaisseurmerite, 'yellow', 'LineWidth', 2)  % Couleur jaune pour le graphe
ylabel('Shield Thickness [m]', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'w')
xlabel('FPA [deg]', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'w')
title('Shield Thickness', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'w')
grid on
grid minor
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w')  

%% Global title
sgtitle('Merit Diagrams for FPA', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial', 'Color', 'w')
