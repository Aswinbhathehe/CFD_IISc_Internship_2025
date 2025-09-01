clc; clear; close all

%% Domain and LBM Parameters
Nx = 256; Ny = 64;                 % Domain size
AR = Nx / Ny;
dx = 1; dt = .01;

tau = 0.8;                         % Relaxation time
nu = (tau - 0.5) / 3;              % Kinematic viscosity
Uw = 0.1;                          % Wall velocity
G = 2 * Uw / Ny;                   % Shear rate

%% Particle Parameters
kappa = 0.125;                     % Confinement ratio
a = kappa * Ny;                    % Radius in lattice units
Rep = 0.625;
G = Rep * nu / a^2;                % Enforce Rep condition
Uw = G * Ny / 2;

rho_f = 1.0;                       % Fluid density
mass = pi * a^2 * rho_f;          % Mass per unit length
I = 0.5 * mass * a^2;             % Moment of inertia

x0 = Nx / 2;                       % Initial x
y0 = Ny / 2 - 0.25 * Ny;           % Initial y (−0.25 H)
theta = linspace(0, 2*pi, 100);   % Lagrangian boundary
Xb = x0 + a * cos(theta);
Yb = y0 + a * sin(theta);
N_b = length(Xb);

%% LBM Setup
w = [4 1 1 1 1 1 1 1 1]/9;
cx = [ 0 1 0 -1  0 1 -1 -1  1];
cy = [ 0 0 1  0 -1 1  1 -1 -1];
opp = [1 4 5 2 3 8 9 6 7];

ux = zeros(Ny,Nx);
uy = zeros(Ny,Nx);
rho = ones(Ny,Nx);
f = zeros(Ny,Nx,9);
feq = @(rho,ux,uy) arrayfun(@(i) ...
    w(i)*rho .* (1 + 3*(cx(i)*ux + cy(i)*uy) ...
    + 4.5*(cx(i)*ux + cy(i)*uy).^2 - 1.5*(ux.^2 + uy.^2)), ...
    1:9, 'UniformOutput', false);

feq_vals = feq(rho,ux,uy);
for i = 1:9
    f(:,:,i) = feq_vals{i};
end

%% Delta function
delta = @(r) (abs(r)<1).*(1 - abs(r)) + ...
             ((abs(r)>=1)&(abs(r)<2)).*(2 - abs(r))/2;

%% Time-stepping
tmax = 1000; tplot = 1;
traj = zeros(tmax/10, 1);
Up = [0, 0]; Omega = 0;

for t = 1:tmax
    %% Macroscopic Variables
    rho = sum(f,3);
    ux = sum(f .* reshape(cx,1,1,9), 3) ./ rho;
    uy = sum(f .* reshape(cy,1,1,9), 3) ./ rho;

    %% Apply top/bottom wall velocities (shear)
    ux(1,:) = -Uw; uy(1,:) = 0;
    ux(end,:) = Uw; uy(end,:) = 0;

    %% Interpolate fluid velocity to Lagrangian points
    Ub = zeros(N_b,2);
    for k = 1:N_b
        % Clamp inside domain
        Xb(k) = min(max(Xb(k), 2), Nx-2);
        Yb(k) = min(max(Yb(k), 2), Ny-2);
        i = floor(Xb(k)); j = floor(Yb(k));
        dxr = Xb(k) - i; dyr = Yb(k) - j;

        u_interp = [0, 0];
        for ii = -1:2
            for jj = -1:2
                phi = delta(dxr - ii) * delta(dyr - jj);
                iidx = mod(i+ii-1,Nx)+1;
                jidx = min(max(j+jj,1),Ny);
                u_interp = u_interp + phi * [ux(jidx,iidx), uy(jidx,iidx)];
            end
        end
        Ub(k,:) = u_interp;
    end

    %% Desired boundary velocity = rigid-body velocity
    Ub_des = zeros(N_b,2);
    Fb = zeros(N_b,2);
    torque_sum = 0;
    force_sum = [0, 0];
    for k = 1:N_b
        rp = [Xb(k) - x0, Yb(k) - y0];
        v_rigid = Up + Omega * [-rp(2), rp(1)];
        % Fb(k,:) = (v_rigid - Ub(k,:)) / dt;
        Fb(k,:) = (Ub(k,:) - v_rigid) / dt;
        force_sum = force_sum + Fb(k,:);
        torque_sum = torque_sum + (rp(1)*Fb(k,2) - rp(2)*Fb(k,1));
    end

    %% Update translational and angular velocity
    Up = Up + force_sum / mass * dt;
    Omega = Omega + torque_sum / I * dt;

    %% Spread force to Eulerian grid
    Fx = zeros(Ny,Nx); Fy = zeros(Ny,Nx);
    for k = 1:N_b
        i = floor(Xb(k)); j = floor(Yb(k));
        dxr = Xb(k) - i; dyr = Yb(k) - j;
        for ii = -1:2
            for jj = -1:2
                phi = delta(dxr - ii) * delta(dyr - jj);
                iidx = mod(i+ii-1,Nx)+1;
                jidx = min(max(j+jj,1),Ny);
                Fx(jidx,iidx) = Fx(jidx,iidx) + phi * Fb(k,1);
                Fy(jidx,iidx) = Fy(jidx,iidx) + phi * Fb(k,2);
            end
        end
    end

    %% Guo forcing term
    for i = 1:9
        eu = (cx(i)*Fx + cy(i)*Fy);
        f(:,:,i) = f(:,:,i) + (1 - 0.5/tau)*3*w(i)*eu;
    end

    %% Collision step
    feq_vals = feq(rho, ux, uy);
    feq_now = cat(3, feq_vals{:});
    f = f - (1/tau)*(f - feq_now);

    %% Streaming
    for i = 1:9
        f(:,:,i) = circshift(f(:,:,i), [cy(i), cx(i)]);
    end

    %% Update Lagrangian markers
    for k = 1:N_b
        rp = [Xb(k)-x0, Yb(k)-y0];
        vel = Up + Omega * [-rp(2), rp(1)];
        Xb(k) = Xb(k) + vel(1)*dt;
        Yb(k) = Yb(k) + vel(2)*dt;
    end

    %% Update particle center
    x0 = mean(Xb);
    y0 = mean(Yb);

    %% Record position every 10 steps
    if mod(t,10) == 0
        traj(tplot) = (y0 - Ny/2) / Ny;  % dimensionless y~
        tplot = tplot + 1;
    end
end

%% Plotting
time = (0:10:tmax-1)*dt;
plot(time, traj, 'b-', 'LineWidth', 2)
xlabel('Time'); ylabel('Transverse position y~')
title('Migration of Neutrally Buoyant Cylinder (Rep=0.625, \kappa=0.125)')
grid on
