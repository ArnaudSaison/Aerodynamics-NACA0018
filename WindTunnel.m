close all; clear all; clc;

% parameters
fig.print_figures = 1;
fig.print_folder = './Figures/';

%%
% =========================================================================
% Data
% =========================================================================
% -------------------------------------------------------------------------
% Importing the data
% -------------------------------------------------------------------------
load('./Data/DataEP.mat');
load('./Data/DataNoEP.mat');
load('./Data/setup.mat');
load('./Data/num_cpx.mat');
load('./Data/xfoil_polar.mat');
xfoilpolar = table2array(xfoilpolar);
xfoilpolar_shift = 2;
xfoilpolar_0 = 16+xfoilpolar_shift;
xfoilpolar = xfoilpolar(11-xfoilpolar_shift:42+xfoilpolar_shift, :);
xfoilpolar(6+xfoilpolar_shift,:) = [];

rho = setup.rho;
taps = setup.coord_taps;        % taps coordinates
sp = setup.span;
ch = setup.chord;

tap_nb = size(taps, 2);         % nb of taps
meas_EP_nb = length(DataEP);    % nb of measurements with endplate
meas_NoEP_nb = length(DataNoEP);% nb of measurements without endplate

ref.angles = [0, 5, 10, 15];
ref.Uinf = 10.83; % m/s


% -------------------------------------------------------------------------
% Checking position of pressure tabs
% -------------------------------------------------------------------------
% indices for position of taps
taps_cr_up = 1:19;  % crosswise taps (up)
taps_cr_do = 29:49; % crosswise taps (down)
taps_cr = [taps_cr_up, taps_cr_do];

taps_sp_up = 20:28; % spanwise taps (up)
taps_sp_do = 50:58; % spanwise taps (down)
taps_sp = [taps_sp_up, taps_sp_do];

fig.press_tabs = figure('Name', 'Position of pressure tabs', 'WindowStyle', 'docked');
ax.press_tabs = axes('Parent', fig.press_tabs);

plot3(taps(1,taps_cr_up), taps(2,taps_cr_up), taps(3,taps_cr_up), 'Color', 'red', 'Marker', '.', 'MarkerSize', 12); hold on;
plot3(taps(1,taps_cr_do), taps(2,taps_cr_do), taps(3,taps_cr_do), 'Color', 'green', 'Marker', '.', 'MarkerSize', 12); hold on;
plot3(taps(1,taps_sp_up), taps(2,taps_sp_up), taps(3,taps_sp_up), 'Color', 'blue', 'Marker', '.', 'MarkerSize', 12); hold on;
plot3(taps(1,taps_sp_do), taps(2,taps_sp_do), taps(3,taps_sp_do), 'Color', 'magenta', 'Marker', '.', 'MarkerSize', 12); hold on;
set(gca,'DataAspectRatio',[1 1 1])
xlabel('X');
ylabel('Y');
zlabel('Z');


% -------------------------------------------------------------------------
% Correction for wind tunnel (solid blockage)
% -------------------------------------------------------------------------
K_1 = 0.52;                     % for vertical model
S = 2.5 * 1.8;                  % surface of the wind tunnel section where the model is located
V_b = 0.0247118638062730 * sp;  % body volume (surface of airfoil taken from panel code)

eps_b = K_1 * V_b / S^(3/2);    % correction factor for velocity

% Correcting the effective wind tunnel velocity
for i = 1:meas_EP_nb
    DataEP(i).Uinf = DataEP(i).Uinf * (1 + eps_b);
end

for i = 1:meas_NoEP_nb
	DataNoEP(i).Uinf = DataNoEP(i).Uinf * (1 + eps_b);
end



%%
% =========================================================================
% Coefficients
% =========================================================================
% -------------------------------------------------------------------------
% Pressure coefficients
% -------------------------------------------------------------------------
C_p.EP = zeros(meas_EP_nb, tap_nb);
C_p.NoEP = zeros(meas_NoEP_nb, tap_nb);

for i = 1:meas_EP_nb
    C_p.EP(i,:) = DataEP(i).p ./ (1/2 * rho * DataEP(i).Uinf^2);
end

for i = 1:meas_NoEP_nb
	C_p.NoEP(i,:) = DataNoEP(i).p ./ (1/2 * rho * DataNoEP(i).Uinf^2);
end

% -------------------------------------------------------------------------
% C_L and C_D
% -------------------------------------------------------------------------
% with endplate
c_l.EP = zeros(meas_EP_nb, 1);
c_d.EP = zeros(meas_EP_nb, 1);

for i = 1:meas_EP_nb
    % Normal force' by integrating pressure along chord (fricton neglected)
    % crosswise pressures up and down
    % sign changed because upper side must have positive contribution
    N_p = (abs(trapz(taps(1,taps_cr_up), DataEP(i).p(taps_cr_up))) - abs(trapz(taps(1,taps_cr_do), DataEP(i).p(taps_cr_do))));
    
    % Lift' and Drag' by using AoA
    c_l.EP(i) = N_p * cos(DataEP(i).AoA / 180 * pi) / (1/2 * rho * DataEP(i).Uinf^2 * ch);
    c_d.EP(i) = N_p * sin(DataEP(i).AoA / 180 * pi) / (1/2 * rho * DataEP(i).Uinf^2 * ch);
end


% without endplate
c_l.NoEP = zeros(meas_NoEP_nb, 1);
c_d.NoEP = zeros(meas_NoEP_nb, 1);

for i = 1:meas_NoEP_nb
    % Normal force' by integrating pressure along chord (fricton neglected)
    % crosswise pressures up and down
    % sign changed because upper side must have positive contribution
    N_p = (abs(trapz(taps(1,taps_cr_up), DataNoEP(i).p(taps_cr_up))) - abs(trapz(taps(1,taps_cr_do), DataNoEP(i).p(taps_cr_do))));
    
    % Lift' and Drag' by using AoA
    c_l.NoEP(i) = N_p * cos(DataNoEP(i).AoA / 180 * pi) / (1/2 * rho * DataNoEP(i).Uinf^2 * ch);
    c_d.NoEP(i) = N_p * sin(DataNoEP(i).AoA / 180 * pi) / (1/2 * rho * DataNoEP(i).Uinf^2 * ch);
end


% H-S
meas_HS_nb = size(num_cpx.hessSmith, 2) - 1;
c_l.HS = zeros(meas_HS_nb, 1);
c_d.HS = zeros(meas_HS_nb, 1);

for i = (1:meas_HS_nb)
    % Normal force' by integrating pressure along chord (fricton neglected)
    % crosswise pressures up and down
    % sign changed because upper side must have positive contribution
    N_p = - trapz(num_cpx.hessSmith(:,1), num_cpx.hessSmith(:,i+1));
    
    % Lift' and Drag' by using AoA
    c_l.HS(i) = N_p * cos(ref.angles(i) / 180 * pi) / (1/2 * rho * ref.Uinf^2 * ch);
    c_d.HS(i) = N_p * sin(ref.angles(i) / 180 * pi) / (1/2 * rho * ref.Uinf^2 * ch);
end



%%
% =========================================================================
% Plots
% =========================================================================
% -------------------------------------------------------------------------
% creating figures
% -------------------------------------------------------------------------
fig.C_p.EP_cr = figure('Name', 'Pressure coefficient EP, crosswise', 'WindowStyle', 'docked');
ax.C_p.EP_cr = axes('Parent', fig.C_p.EP_cr);

fig.C_p.EP_sp = figure('Name', 'Pressure coefficient EP, spanwise', 'WindowStyle', 'docked');
ax.C_p.EP_sp = axes( 'Parent', fig.C_p.EP_sp);

fig.C_p.NoEP_cr = figure('Name', 'Pressure coefficient no EP, crosswise', 'WindowStyle', 'docked');
ax.C_p.NoEP_cr = axes('Parent', fig.C_p.NoEP_cr);

fig.C_p.NoEP_sp = figure('Name', 'Pressure coefficient no EP, spanwise', 'WindowStyle', 'docked');
ax.C_p.NoEP_sp = axes('Parent', fig.C_p.NoEP_sp);

fig.C_p.spanwise = figure('Name', 'Spanwise pressure coefficient', 'WindowStyle', 'docked');
ax.C_p.spanwise = axes('Parent', fig.C_p.spanwise);

% -------------------------------------------------------------------------

fig.stall = figure('Name', 'Stall curve', 'WindowStyle', 'docked');
ax.stall = axes('Parent', fig.stall);

fig.drag_polar = figure('Name', 'Drag polar', 'WindowStyle', 'docked');
ax.drag_polar = axes('Parent', fig.drag_polar);

fig.cobra_span = figure('Name', 'Cobra span v', 'WindowStyle', 'docked');
ax.cobra_span = axes('Parent', fig.cobra_span);

fig.C_p.separ = figure('Name', 'Effect of separation on pressure coefficient', 'WindowStyle', 'docked');
ax.C_p.separ = axes('Parent', fig.C_p.separ);

fig.cobra_sp_u = figure('Name', 'Cobra span u', 'WindowStyle', 'docked');
ax.cobra_sp_u = axes('Parent', fig.cobra_sp_u);

% -------------------------------------------------------------------------
% Plotting data
% -------------------------------------------------------------------------
% C_p with EP crosswise # GRAPH 1: Cp(x/c)
plot_dataset_EP = 8;
plot_angle = 3;

plot(ax.C_p.EP_cr,     num_cpx.hessSmith(:,1), num_cpx.hessSmith(:,plot_angle+1), ...
                       'LineStyle', '--', ...
                       'LineWidth', 2, ...
                       'Color', 'blue', ...
                       'Marker', 'none', ...
                       'DisplayName', 'Hess-Smith');
hold(ax.C_p.EP_cr,     'on');

plot(ax.C_p.EP_cr,     num_cpx.xfoil(:,1), num_cpx.xfoil(:,plot_angle+1), ...
                       'LineStyle', '-', ...
                       'Color', 'red', ...
                       'Marker', 'none', ...
                       'DisplayName', 'XFOIL');
hold(ax.C_p.EP_cr,     'on');

plot(ax.C_p.EP_cr,     [(taps(1,taps_cr_up) + 0.19) / ch, (taps(1,taps_cr_do) + 0.19) / ch], ...
                       [C_p.EP(plot_dataset_EP,taps_cr_up),  C_p.EP(plot_dataset_EP,taps_cr_do)], ...
                       'LineStyle', 'none', ...
                       'Marker', '*', ...
                       'MarkerEdgeColor', 'black', ...
                       'DisplayName', 'Experimental');
hold(ax.C_p.EP_cr,     'on');

% C_p with EP spanwise
plot_dataset_EP = 12;

plot(ax.C_p.EP_sp,     taps(3,taps_sp_up) / sp, C_p.EP(plot_dataset_EP,taps_sp_up)); hold(ax.C_p.EP_sp, 'on');
plot(ax.C_p.EP_sp,     taps(3,taps_sp_do) / sp, C_p.EP(plot_dataset_EP,taps_sp_do)); hold(ax.C_p.EP_sp, 'on');

% C_p with no EP crosswise
plot_dataset_NoEP = 1:4;

plot(ax.C_p.NoEP_cr,   taps(1,taps_cr_up) / sp, C_p.NoEP(plot_dataset_NoEP,taps_cr_up)); hold(ax.C_p.NoEP_cr, 'on');
plot(ax.C_p.NoEP_cr,   taps(1,taps_cr_do) / sp, C_p.NoEP(plot_dataset_NoEP,taps_cr_do)); hold(ax.C_p.NoEP_cr, 'on');

% C_p with no EP spanwise
plot_dataset_NoEP = 1:4;

plot(ax.C_p.NoEP_sp,   taps(3,taps_sp_up) / sp, C_p.NoEP(plot_dataset_NoEP,taps_sp_up)); hold(ax.C_p.NoEP_sp, 'on');
plot(ax.C_p.NoEP_sp,   taps(3,taps_sp_do) / sp, C_p.NoEP(plot_dataset_NoEP,taps_sp_do)); hold(ax.C_p.NoEP_sp, 'on');

% Stall
plot_dataset_Stall = 2;
plot_dataset_Stall = (1:3:12) + plot_dataset_Stall - 1;

[conf_map.aoa_var, conf_map.cl_confMap] = Theoretical;

plot(ax.stall,         conf_map.aoa_var, conf_map.cl_confMap, ...
                       'LineStyle', '-', ...
                       'LineWidth', 1, ...
                       'Color', 'cyan', ...
                       'Marker', 'none', ...
                       'DisplayName', 'Conformal mapping'); hold(ax.stall, 'on');

plot(ax.stall,         xfoilpolar(xfoilpolar_0:end,1), xfoilpolar(xfoilpolar_0:end,2), ...
                       'LineStyle', '-', ...
                       'Color', 'red', ...
                       'Marker', 'none', ...
                       'DisplayName', 'XFOIL'); hold(ax.stall, 'on');

plot(ax.stall,         ref.angles, c_l.EP(plot_dataset_Stall), ...
                       'LineStyle', 'none', ...
                       'Marker', '*', ...
                       'MarkerEdgeColor', 'black', ...
                       'DisplayName', 'Experimental (with EP)'); hold(ax.stall, 'on');
                   
plot(ax.stall,         ref.angles, c_l.NoEP(:), ...
                       'LineStyle', 'none', ...
                       'Marker', 'o', ...
                       'MarkerEdgeColor', 'green', ...
                       'DisplayName', 'Experimental (no EP)'); hold(ax.stall, 'on');

% Drag polar
plot(ax.drag_polar,    xfoilpolar(:,3), xfoilpolar(:,2), ...
                       'LineStyle', '-', ...
                       'Color', 'red', ...
                       'Marker', 'none', ...
                       'DisplayName', 'XFOIL'); hold(ax.drag_polar, 'on');

plot(ax.drag_polar,    c_d.EP(plot_dataset_Stall), c_l.EP(plot_dataset_Stall), ...
                       'LineStyle', 'none', ...
                       'Marker', '*', ...
                       'MarkerEdgeColor', 'black', ...
                       'DisplayName', 'Experimental'); hold(ax.drag_polar, 'on');

drag_polar_intervals = 5;
for j = (1:drag_polar_intervals:size(xfoilpolar, 1)) + xfoilpolar_shift
    plot(ax.drag_polar,xfoilpolar(j,3), xfoilpolar(j,2), ...
                       'LineStyle', 'none', ...
                       'Marker', '.', ...
                       'MarkerEdgeColor', 'black', ...
                       'HandleVisibility','off'); 
    hold(ax.drag_polar,'on');
    
    text(ax.drag_polar,xfoilpolar(j,3), xfoilpolar(j,2), ...
                       [' ' num2str(xfoilpolar(j,1)) '°']); 
    hold(ax.drag_polar,'on');
end

% spanwise
plot_dataset_EP = 3;
plot_dataset_EP = (1:3:12) + plot_dataset_EP - 1;
plot_dataset_NoEP = 1:4;

plot(ax.C_p.spanwise,  taps(3,taps_sp_up) / sp, C_p.EP(plot_dataset_EP,taps_sp_up), ...
                       'LineStyle', '-', ...
                       'LineWidth', 1.5, ...
                       'Color', 'black', ...
                       'DisplayName', 'Experimental (with EP)'); hold(ax.C_p.spanwise, 'on');

plot(ax.C_p.spanwise,  taps(3,taps_sp_do) / sp, C_p.EP(plot_dataset_EP,taps_sp_do), ...
                       'LineStyle', '-', ...
                       'Color', 'black', ...
                       'HandleVisibility','off'); hold(ax.C_p.spanwise, 'on');

plot(ax.C_p.spanwise,  taps(3,taps_sp_up) / sp, C_p.NoEP(plot_dataset_NoEP,taps_sp_up), ...
                       'LineStyle', '-', ...
                       'LineWidth', 1.5, ...
                       'Color', 'green', ...
                       'DisplayName', 'Experimental (no EP)'); hold(ax.C_p.spanwise, 'on');

plot(ax.C_p.spanwise,  taps(3,taps_sp_do) / sp, C_p.NoEP(plot_dataset_NoEP,taps_sp_do), ...
                       'LineStyle', '-', ...
                       'Color', 'green', ...
                       'HandleVisibility','off'); hold(ax.C_p.spanwise, 'on');
                   
% cobra span v
cobra_comp = 2; % v component
cobra_NoEP_datasets = 3:3:12;
cobra_textAoA_pos = 3;

plot(ax.cobra_span,    [0, sp], [0,0], ...
                       'LineStyle', '-', ...
                       'Linewidth', 1.5, ...
                       'Color', 'magenta', ...
                       'Marker', '|', ...
                       'DisplayName','Wing span'); 
hold(ax.cobra_span, 'on');

for i = 1:meas_NoEP_nb
    
    text(ax.cobra_span,DataNoEP(i).Pos_Span(3,cobra_textAoA_pos), DataNoEP(i).Cobra_Span(cobra_comp,cobra_textAoA_pos), ...
                       [' ' num2str(DataNoEP(i).AoA) '°'], ...
                       'VerticalAlignment', 'bottom'); hold(ax.cobra_span, 'on');
    
    if i ~= meas_NoEP_nb
    
    plot(ax.cobra_span,DataEP(cobra_NoEP_datasets(i)).Pos_Span(3,:), DataEP(cobra_NoEP_datasets(i)).Cobra_Span(cobra_comp,:), ...
                       'LineStyle', '-', ...
                       'Color', 'black', ...
                       'HandleVisibility','off'); hold(ax.cobra_span, 'on');
    
    plot(ax.cobra_span,DataNoEP(i).Pos_Span(3,:), DataNoEP(i).Cobra_Span(cobra_comp,:), ...
                       'LineStyle', '-', ...
                       'Color', 'green', ...
                       'HandleVisibility','off'); hold(ax.cobra_span, 'on');
    
    else
    
    plot(ax.cobra_span,DataEP(cobra_NoEP_datasets(i)).Pos_Span(3,:), DataEP(cobra_NoEP_datasets(i)).Cobra_Span(cobra_comp,:), ...
                       'LineStyle', '-', ...
                       'Color', 'black', ...
                       'DisplayName','with EP'); hold(ax.cobra_span, 'on');
    
    plot(ax.cobra_span,DataNoEP(i).Pos_Span(3,:), DataNoEP(i).Cobra_Span(cobra_comp,:), ...
                       'LineStyle', '-', ...
                       'Color', 'green', ...
                       'DisplayName','no EP'); hold(ax.cobra_span, 'on');
    
    end
end

% Cp(x) for stall (with EP)
plot_dataset_EP_list = 10:12;
C_p_separ_colors = copper(3);
C_p_separ_colors([2,3],:) = C_p_separ_colors([3,2],:);

for plot_dataset_EP = flip(plot_dataset_EP_list)
    plot(ax.C_p.separ, [(taps(1,taps_cr_up) + 0.19) / ch, flip((taps(1,taps_cr_do) + 0.19) / ch)], ...
                       [C_p.EP(plot_dataset_EP,taps_cr_up),  flip(C_p.EP(plot_dataset_EP,taps_cr_do))], ...
                       'LineStyle', '-', ...
                       'Color', C_p_separ_colors(plot_dataset_EP-9,:), ...
                       'LineWidth', (plot_dataset_EP-9), ...
                       'DisplayName', [' ' num2str(DataEP(plot_dataset_EP).Uinf, '%.2f') ' m/s']);
    hold(ax.C_p.separ, 'on');
end

% cobra span u
cobra_comp = 1; % u component
cobra_sp_u_angle = 0:3:9;
plot_crobra_span_u_list = 1:3;
k = 1;

for j = cobra_sp_u_angle
    
    for i = plot_crobra_span_u_list
        if i ~= plot_crobra_span_u_list(end)
            plot(ax.cobra_sp_u, DataEP(i+j).Cobra_Cross(cobra_comp,:), DataEP(i+j).Pos_Cross(2,:), ...
                       'LineStyle', '-', ...
                       'Color', 'black', ...
                       'LineWidth', k, ...
                       'HandleVisibility','off'); 
            hold(ax.cobra_sp_u, 'on');
        else
            plot(ax.cobra_sp_u, DataEP(i+j).Cobra_Cross(cobra_comp,:), DataEP(i+j).Pos_Cross(2,:), ...
                       'LineStyle', '-', ...
                       'Color', 'black', ...
                       'LineWidth', k, ...
                       'DisplayName', [' ' num2str(DataEP(j+1).AoA) '°']); 
            hold(ax.cobra_sp_u, 'on');
        end
    end
    
    k = k + 1;
end


% -------------------------------------------------------------------------
% Plot style
% -------------------------------------------------------------------------
set(ax.C_p.EP_cr,      'Ydir', 'reverse')
grid(ax.C_p.EP_cr,     'on');
xlabel(ax.C_p.EP_cr,   '$x/c$','interpreter','latex');
ylabel(ax.C_p.EP_cr,   '$C_p$','interpreter','latex');
box(ax.C_p.EP_cr,      'on');
legend(ax.C_p.EP_cr);

set(ax.C_p.EP_sp,      'Ydir', 'reverse')
grid(ax.C_p.EP_sp,     'on');
xlabel(ax.C_p.EP_sp,   '$\displaystyle z$','interpreter','latex');
ylabel(ax.C_p.EP_sp,   '$C_p$','interpreter','latex');
box(ax.C_p.EP_sp,      'on');

set(ax.C_p.NoEP_cr,    'Ydir', 'reverse')
grid(ax.C_p.NoEP_cr,   'on');
xlabel(ax.C_p.NoEP_cr, '$\displaystyle z$','interpreter','latex');
ylabel(ax.C_p.NoEP_cr, '$C_p$','interpreter','latex');
box(ax.C_p.NoEP_cr,    'on');


set(ax.C_p.NoEP_sp,    'Ydir', 'reverse')
grid(ax.C_p.NoEP_sp,   'on');
xlabel(ax.C_p.NoEP_sp, '$\displaystyle z$','interpreter','latex');
ylabel(ax.C_p.NoEP_sp, '$C_p$','interpreter','latex');
box(ax.C_p.NoEP_sp,    'on');

grid(ax.stall,         'on');
xlabel(ax.stall,       'Angle of attack $(\alpha)$ [$\circ$]','interpreter','latex');
ylabel(ax.stall,       '$c_\ell$','interpreter','latex');
box(ax.stall,          'on');
legend(ax.stall,       'Location', 'northwest');
xlim(ax.stall,         [0, 18])

grid(ax.drag_polar,    'on');
xlabel(ax.drag_polar,  '$c_d$','interpreter','latex');
ylabel(ax.drag_polar,  '$c_\ell$','interpreter','latex');
box(ax.drag_polar,     'on');
legend(ax.drag_polar,  'Location', 'east');

grid(ax.C_p.spanwise,  'on');
xlabel(ax.C_p.spanwise,'span ($z$) [m]','interpreter','latex');
ylabel(ax.C_p.spanwise,'$C_p$','interpreter','latex');
box(ax.C_p.spanwise,   'on');
legend(ax.C_p.spanwise,'Location', 'west');

grid(ax.cobra_span,    'on');
xlabel(ax.cobra_span,  'span ($z$) [m]','interpreter','latex');
ylabel(ax.cobra_span,  '$v$ velocity component [m/s]','interpreter','latex');
box(ax.cobra_span,     'on');
legend(ax.cobra_span,  'Location', 'northwest');
xlim(ax.cobra_span,    [min(DataNoEP(1).Pos_Span(3,:)) - 0.05, max(DataNoEP(1).Pos_Span(3,:)) + 0.05])

set(ax.C_p.separ,      'Ydir', 'reverse')
grid(ax.C_p.separ,     'on');
xlabel(ax.C_p.separ,   '$x/c$','interpreter','latex');
ylabel(ax.C_p.separ,   '$C_p$','interpreter','latex');
box(ax.C_p.separ,      'on');
legend(ax.C_p.separ);

grid(ax.cobra_sp_u,    'on');
xlabel(ax.cobra_sp_u,  'crosswise velocity component $u$ [m/s]','interpreter','latex');
ylabel(ax.cobra_sp_u,  'normal position along $y$ [m]','interpreter','latex');
box(ax.cobra_sp_u,     'on');
legend(ax.cobra_sp_u);

% -------------------------------------------------------------------------
% Prining figures
% -------------------------------------------------------------------------
fig.height = 1;
fig.width = fig.height * 1.5;

if fig.print_figures
    make_fig(fig.press_tabs, fig.print_folder, 'press_tabs', fig.height, fig.width, 0);
    make_fig(fig.C_p.EP_cr, fig.print_folder, 'C_p_EP_cr', fig.height, fig.width, 0);
    make_fig(fig.C_p.EP_sp, fig.print_folder, 'C_p_EP_sp', fig.height, fig.width, 0);
    make_fig(fig.C_p.NoEP_cr, fig.print_folder, 'C_p_NoEP_cr', fig.height, fig.width, 0);
    make_fig(fig.C_p.NoEP_sp, fig.print_folder, 'C_p_NoEP_sp', fig.height, fig.width, 0);
    make_fig(fig.stall, fig.print_folder, 'stall', fig.height, fig.width, 0);
    make_fig(fig.drag_polar, fig.print_folder, 'drag_polar', fig.height, fig.width, 0);
    make_fig(fig.C_p.spanwise, fig.print_folder, 'Cp_spanwise', fig.height, fig.width, 0);
    make_fig(fig.cobra_span, fig.print_folder, 'cobra_span', fig.height, fig.width, 0);
    make_fig(fig.C_p.separ, fig.print_folder, 'C_p_separ', fig.height, fig.width, 0);
    make_fig(fig.cobra_sp_u, fig.print_folder, 'cobra_sp_u', fig.height, fig.width, 0);
end






mu = 1.8e-5;
Re = rho * ref.Uinf * ch / mu









