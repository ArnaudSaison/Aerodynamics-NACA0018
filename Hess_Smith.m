%     __  __                    _____           _ __  __
%    / / / /__  __________     / ___/____ ___  (_) /_/ /_
%   / /_/ / _ \/ ___/ ___/_____\__ \/ __ `__ \/ / __/ __ \
%  / __  /  __(__  |__  )_____/__/ / / / / / / / /_/ / / /
% /_/ /_/\___/____/____/     /____/_/ /_/ /_/_/\__/_/ /_/.m
%
% This matlab script implements the Hess & Smith straight panel method. The code
% below was written as part of the aerodynamics course project (AERO0001-1),
% academic year 2021-2022.
%
% 1. Parameters settings.
% 2. Airfoil geometry.
% 3. Panel definition.
% 4. Construction of the linear system.
% 5. System resolution.
% 6. Results: aerodynamic forces and coefficients.

clear; close all;

%% 1. Parameter settings. {{{1

% Import the wind tunnel experiment setup.
load('Lab_data/setup.mat')

% General parameters.
c       = setup.chord;  % Chordwide length of the wing [m].
rho     = setup.rho;    % Air density [kg/m³].
nPanel  = 200;          % Number of panels used.
v_inf   = 10.83;        % Free stream velocity [m/s].
aoa     = 5;            % Angle of attack [°].
aoa_rad = deg2rad(aoa);

% Airfoil parameters.
NACA_id = '0018';
eps = str2double(NACA_id(1))  /100;  % Maximal camber ratio.
p   = str2double(NACA_id(2))  /10;   % Location of maximal camber from LE.
tau = str2double(NACA_id(3:4))/100;  % Thickness ratio.

% Assertions on parameters value.
assert(~mod(nPanel, 2), ...
	"Number of panels must be even.");
assert(all([c, nPanel, v_inf] >= 0), ...
	"Aberrant input value detected (negative quantity).");
assert(length(NACA_id) == 4, ...
	"This script only supports NACA 4-digits airfoils.");

%% 2. Airfoil geometry. {{{1

% panel definition: x-component boundary points.
xi = linspace(0, 2*pi, nPanel+1);
x  = c/2 * (cos(xi) + 1);

% x-half : x coordinate along the x axis.
xh = x(nPanel/2 +1 : end);

% NACA definition of thickness.
T = 10 * tau * c * ( ...
	  0.2969 * sqrt(xh/c)    ...
	- 0.1260 *     (xh/c)    ...
	- 0.3537 *     (xh/c).^2 ...
	+ 0.2843 *     (xh/c).^3 ...
	- 0.1015 *     (xh/c).^4 ...
);

% NACA definition of camber line.
inf_pc = xh < p*c;
Y_bar = [ ...
	eps*     xh( inf_pc)  /      p ^2 .* (  - xh( inf_pc)/c + 2*p), ...
	eps*(c - xh(~inf_pc)) / (1 - p)^2 .* (1 + xh(~inf_pc)/c - 2*p)  ...
];

% First-order approximation of the camber line derivative.
dY_bar = diff(Y_bar)./diff(xh);
% Add the backward difference for the last point (@ TE).
dY_bar = [dY_bar, dY_bar(end)];
theta_camber = atan(dY_bar);

% Coordinates of the lower and upper surface.
X_u = xh    - T/2 .* sin(theta_camber);
X_l = xh    + T/2 .* sin(theta_camber);
Y_u = Y_bar + T/2 .* cos(theta_camber);
Y_l = Y_bar - T/2 .* cos(theta_camber);
% Gather upper and lower surface coordinate points. Unlike the x-component
% boundary points, the (X, Y) coordinate points run clockwise around the
% airfoil, starting and finishing at the TE. We trim the first value of X_u
% (i.e. the coordinate point of the LE) to not repeat it, as this point is
% already taken in account by the last element of X_l.
X = [flip(X_l), X_u(2:end)];
Y = [flip(Y_l), Y_u(2:end)];

% Airfoil preview.
figure('WindowStyle', 'docked'); hold on;
plot(X, Y, 'marker', '*');
plot(xh, Y_bar);
title(['NACA-', NACA_id, ' airfoil preview']);
xlabel('X/m'); ylabel('Y/m');
axis equal; grid;

%% 3. Panel definition. {{{1

% Control points (x_bar(i), y_bar(i)), located @ mid-point of each panel.
x_bar = movmean(X, [0, 1]);
y_bar = movmean(Y, [0, 1]);
% Remove irrelevant last value.
x_bar = x_bar(1:end-1);
y_bar = y_bar(1:end-1);

% Visualization of the control point placements.
plot(x_bar, y_bar, 'linestyle', 'none', 'marker', 'o');
legend( ...
	'Panel boundary points.', ...
	'Camber line.',           ...
	'Panel control points.'   ...
);

% Diffs between successive boundary points.
dx = diff(X);
dy = diff(Y);

% Panel length.
l = sqrt(dx.^2 + dy.^2);

% Panel inclination to the x-axis.
theta = atan2(dy, dx);

%% 4. Construction of the linear system. {{{1

% Matrices preallocation.
A = zeros(nPanel, nPanel);
B = zeros(nPanel, nPanel);
C = zeros(nPanel, nPanel);
D = zeros(nPanel, nPanel);
E = zeros(nPanel, nPanel);
F = zeros(nPanel, nPanel);
G = zeros(nPanel, nPanel);

% Matrices construction.
for i = 1:nPanel
	for j = 1:nPanel
		A(i,j) = -(x_bar(i)-X(j)) * cos(theta(j)) ...
		         -(y_bar(i)-Y(j)) * sin(theta(j));
		E(i,j) =  (x_bar(i)-X(j)) * sin(theta(j)) ...
		         -(y_bar(i)-Y(j)) * cos(theta(j));
		C(i,j) = sin(theta(i)-theta(j));
		D(i,j) = cos(theta(i)-theta(j));
		B(i,j) = (x_bar(i)-X(j))^2 + (y_bar(i)-Y(j))^2;
		F(i,j) = log(1 + (l(j)^2 + 2*A(i,j)*l(j)) / B(i,j));
		G(i,j) = atan2((E(i,j)*l(j)), A(i,j)*l(j)+B(i,j));
	end
end

% Contruction of matrix M.

M_11 = C.*F/2 - D.*G;
M_12 = D.*F/2 + C.*G;

% Force the values on the M_11 and M_12 diagonals.
M_11 = M_11 - diag(diag(M_11)) + pi*eye(size(M_11));  % pi on M_11 diagonal
M_12 = M_12 - diag(diag(M_12));                       % 0  on M_12 diagonal

M_12 = sum(M_12,2);

M_21 = ...
	  D(end,:).*F(end,:)/2 + C(end,:).*G(end,:) ...
	+ D(1,  :).*F(1,  :)/2 + C(1,  :).*G(1,  :);
M_22 = - sum( ...
	  C(end,:).*F(end,:)/2 - D(end,:).*G(end,:) ...
	+ C(1  ,:).*F(1,  :)/2 - D(1,  :).*G(1,  :) ...
);

M = 1/(2*pi) * [M_11, M_12; M_21, M_22];

% Construction of matrix N.

N = v_inf * [ ...
	sin(theta-aoa_rad),                               ...
	(cos(theta(1)-aoa_rad) + cos(theta(end)-aoa_rad)) ...
]';

%% 5. System resolution. {{{1

res = M\N;

% Sources.
q = res(1:end-1);
% Vorticity.
gamma = res(end);

%% 6. Results: aerodynamic forces and coefficients. {{{1

% TODO: Try to vectorize v_tan.
% Tangential velocity on each panel.
v_tan = zeros(nPanel,1);
% Second value of the tangential velocity.
val_2 = (D.*F)/2 + C.*G;
term2 = zeros(nPanel, 1);
for i = 1 : nPanel
	val_2(i, i) = 0;
	for j = 1 : nPanel
		term2(i) = term2(i) + q(j) * val_2(i, j) / (2*pi);
	end
end
% third value of the tangential velocity.
val_3 = sum(M_11, 2);
for i = 1 : nPanel
	v_tan(i, 1) = v_inf*cos(theta(i)-aoa_rad) - term2(i) + gamma/(2*pi)*val_3(i);
end

% Pressure coefficient for each panel.
cp = 1 - (v_tan/v_inf).^2;
diff_cp = cp(1:end/2) - flip(cp(end/2+1:end));
diff_pressure = diff_cp * (0.5*rho*v_inf^2);

% Axial and normal forces.
Nprime = trapz(x_bar(end/2+1:end), flip(diff_pressure));
Aprime = trapz(x_bar(end/2+1:end),flip(-diff_pressure'.*(dy(1:end/2)./dx(1:end/2))));

% Lift and drag forces.
Lprime = Nprime*cosd(aoa) - Aprime*sind(aoa);
Dprime = Nprime*sind(aoa) + Aprime*cosd(aoa);

% Lift and drag coefficients.
cl = Lprime / (0.5*rho*v_inf^2*c);
cd = Dprime / (0.5*rho*v_inf^2*c);

% Lift coefficient from Kutta-Joukowski theorem.
% To be compared with previously computed cl.
cl_KJ = 2*gamma*sum(l)/(v_inf*c);

% Pressure coefficient graph.
figure('WindowStyle', 'docked'); hold on; grid;
plot( ...
	x_bar(1:nPanel/2)./c, ...
	cp(1:nPanel/2),       ...
	'color', 'red',       ...
	'Linewidth', 1.5      ...
);
plot( ...
	x_bar(nPanel/2:nPanel)./c, ...
	cp(nPanel/2:end),          ...
	'color', 'blue',           ...
	'Linewidth', 1.5           ...
);
xlabel('x/c');
ylabel('Cp');
legend('Lower surface','Upper surface');
title('Distribution of the pressure coefficient Cp along the chord');
set(gca, 'YDir', 'reverse')