function [aoa_var, cl_confMap] = Theoretical
%% Parameter settings.

% Import the wind tunnel experiment setup.
%load('Data/setup.mat')

% General parameters.
c       = 0.45;  % Chordwide length of the wing [m].
nPanel  = 200;          % Number of panels used.
v_inf   = 10.83;        % Free stream velocity [m/s].
aoa     = 0;            % Angle of attack [°].

% Airfoil parameters.
NACA_id = '0018';
eps = str2double(NACA_id(1))  /100;  % Maximal camber ratio.
p   = str2double(NACA_id(2))  /10;   % Location of maximal camber from LE.
tau = str2double(NACA_id(3:4))/100;  % Thickness ratio.

% NACA definition of tickness.
T = @(x) 10 * tau * c * ( ...
   0.2969 * sqrt(x/c)    ...
 - 0.1260 *     (x/c)    ...
 - 0.3537 *     (x/c).^2 ...
 + 0.2843 *     (x/c).^3 ...
 - 0.1015 *     (x/c).^4 ...
);

%% Conformal mapping.

% x-coord sample along the chord.
x = linspace(0, c, 200);
% Mean of the thickness.
T_bar = mean(T(x));
% aoa sample.
aoa_var = -10:1:25;

% Implementation of Conformal mapping cl formula for thick airfoil (see report
% and slide 76, lesson 3).
cl_confMap = 2*pi*(1+(4*T_bar/(c*3*sqrt(3))))*sind(aoa_var);

% Plot the results.
% plot(aoa_var, cl_confMap); grid;
% xlabel("Angle of attack (°)");
% ylabel("cl computed from conformal mapping");
end