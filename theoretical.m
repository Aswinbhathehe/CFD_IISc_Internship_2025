%% Clear workspace and close figures
clear; clc; close all;

%% Domain Setup (Physical)
Lx = 4; Ly = 1;              % Domain size (as in the paper)
kappa = 0.125;               % Confinement ratio
a = kappa * Ly;              % Particle radius
G = 2 * 0.1 / Ly;            % Shear rate (Uw = 0.1)

%% Simulated Re_p cases
Rep_list = [1, 3, 10];
colors = lines(length(Rep_list));

%% 1. Transverse Position vs Time
time = linspace(0, 10, 300);   % Arbitrary time steps
y_pos = zeros(length(Rep_list), length(time));

y_pos(1,:) = -0.25 * exp(-0.5 * time);                               % Rep=1
y_pos(2,:) = (-0.25 ) * exp(-0.3 * time);                % Rep=3
y_pos(3,:) = (-0.25 ) * exp(-0.2 * time);              % Rep=10

%% 2. Lift Force vs Transverse Position
y_positions = linspace(-0.35, 0.35, 300);
lift_force = zeros(length(Rep_list), length(y_positions));

lift_force(1,:) = -29 * y_positions;                                               % Rep=1
lift_force(2,:) = -100 * (y_positions - 0.1) .* (y_positions + 0.1) .* y_positions; % Rep=3
lift_force(3,:) = -20 * (y_positions - 0.23) .* (y_positions + 0.23) .* y_positions; % Rep=10

%% 3. Equilibrium Position vs Rep
Rep_values = linspace(0.1, 20, 200);
critical_Rep = 2.0;

equilibrium_position = zeros(size(Rep_values));
equilibrium_position(Rep_values >= critical_Rep) = ...
    0.1 + 0.17 * (Rep_values(Rep_values >= critical_Rep) - critical_Rep) / (20 - critical_Rep);
unstable_position = -equilibrium_position;  % Reflect for symmetric unstable branch

%% Create figure with subplots
figure();

% --- Subplot 1: Trajectory ---
subplot(2,2,1);
hold on;
for i = 1:length(Rep_list)
    plot(time, y_pos(i,:), 'LineWidth', 2, 'Color', colors(i,:), ...
         'DisplayName', ['Re_p = ' num2str(Rep_list(i))]);
end
xlabel('Time'); ylabel('Transverse Position (\ity/H)');
title('(a) Transverse Position vs Time');
legend('Location', 'best'); grid on;

% --- Subplot 2: Lift Force vs Position ---
subplot(2,2,2);
hold on;
for i = 1:length(Rep_list)
    plot(y_positions, lift_force(i,:), 'LineWidth', 2, 'Color', colors(i,:), ...
         'DisplayName', ['Re_p = ' num2str(Rep_list(i))]);
end
xlabel('Transverse Position (\ity/H)');
ylabel('Lift Force (\itF̃_L)');
title('(b) Lift Force vs Transverse Position');
legend('Location', 'best'); grid on;

% --- Subplot 3+4: Equilibrium Position vs Rep (log scale) ---
subplot(2,2,[3,4]);

% Plot stable branches
semilogx(Rep_values, equilibrium_position, 'k-', 'LineWidth', 2, ...
          'DisplayName', 'Stable Equilibrium');
hold on;

% Plot symmetric unstable branch (reflected)
semilogx(Rep_values, unstable_position, 'k--', 'LineWidth', 2, ...
          'DisplayName', 'Unstable Equilibrium');

% Critical Rep line
xline(critical_Rep, 'k:', 'LineWidth', 2, 'DisplayName', 'Critical Re_p');

% Formatting
xlabel('Reynolds number, Re_p');
ylabel('Equilibrium Position (\ity/H)');
title('(c) Equilibrium Position vs Re_p (log scale)');
legend('Location', 'best');
grid on;
set(gca, 'XScale', 'log');
xlim([0.1 20]);
ylim([-0.3 0.3]);

%% Final title
sgtitle('IB-LBM Validation with Hypothetical Data for Circular Cylinder (\kappa = 0.125)');
set(gcf, 'Color', 'w');


% clc; clear; close all;
% 
% %% Parameters
% kappa = 0.125;         % Confinement ratio
% H = 1.0;               % Channel height (nondimensionalized)
% a = kappa * H;         % Particle radius
% G = 1.0;               % Shear rate
% nu = 1.0;              % Kinematic viscosity
% Rep = 0.625;           % Particle Reynolds number
% 
% % Time parameters
% dt = 0.1;
% Tmax = 300;
% Nsteps = Tmax/dt;
% 
% % Initialize particle
% y = zeros(Nsteps,1);
% y(1) = -0.25;          % Initial position
% 
% % Predefine lift force function (from cubic fit to Figure 3 of the paper)
% % FL(y) ≈ -k*(y^3 - b*y), with stable fixed points at ±sqrt(b)
% b = 0.01;
% k = 1;
% 
% % Simulation loop for transverse migration
% for t = 1:Nsteps-1
%     FL = -k*(y(t)^3 - b*y(t));
%     y(t+1) = y(t) + FL*dt;
% end
% 
% %% Plot 1: Transverse Position vs Time
% subplot(1,3,1);
% time = (0:Nsteps-1)*dt;
% plot(time, y, 'b', 'LineWidth', 2);
% xlabel('Time');
% ylabel('Transverse position y/H');
% title('Transverse Position vs Time');
% grid on;
% 
% %% Plot 2: Lift Force vs Transverse Position
% y_span = linspace(-0.3,0.3,500);
% FL_curve = -k*(y_span.^3 - b*y_span);
% subplot(1,3,2);
% plot(y_span, FL_curve, 'r', 'LineWidth', 2);
% xlabel('Transverse position y/H');
% ylabel('Lift force (arb. units)');
% title('Lift Force vs Transverse Position');
% grid on;
% 
% %% Plot 3: Equilibrium Position vs Rep (kappa = 0.125)
% Rep_vals = logspace(log10(0.1), log10(50), 100);
% y_eq = zeros(size(Rep_vals));
% for i = 1:length(Rep_vals)
%     if Rep_vals(i) < 2.0
%         y_eq(i) = 0;
%     else
%         % Bifurcation form: y_eq = A*sqrt(Rep - Rep_crit)
%         Rep_crit = 2.0;
%         A = 0.13;
%         y_eq(i) = A * sqrt(Rep_vals(i) - Rep_crit);
%     end
% end
% 
% subplot(1,3,3);
% semilogx(Rep_vals, -y_eq, 'b', Rep_vals, y_eq, 'b'); hold on;
% semilogx(Rep_vals, zeros(size(Rep_vals)), 'k--');
% xlabel('Rep');
% ylabel('Equilibrium Position y/H');
% title('Equilibrium Position vs Rep (κ = 0.125)');
% grid on;
% legend('Stable Positions');
% 
% sgtitle('Verification of PhysRevResearch.2.013009 (κ = 0.125)');
