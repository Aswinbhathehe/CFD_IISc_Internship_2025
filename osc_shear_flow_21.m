% Oscillating Shear Flow Simulation with MG, Rhie-Chow
clc; clear; close all;

%% Parameters
Nx =64; 
Ny=64;
Lx = 1; 
Ly =.2;
dx = Lx/Nx; 
dy = Ly/Ny;
x =linspace(0,Lx, Nx);
y =linspace(0,Ly, Ny);
[X, Y] =meshgrid(x, y);

nu = 0.05;             % Kinematic viscosity
U0 = 5;               % Wall velocity amplitude
f = 10;
omega = 2*pi*f;
rho =1;

% Time settings
CFL = 0.35;
dt = CFL*min(dx, dy)^2/nu;
Tf= 2/f;
Nt = ceil(Tf/dt);
dt = Tf/Nt;
time = linspace(0,Tf,Nt);


u = zeros(Ny,Nx);
v = zeros(Ny,Nx);
p = zeros(Ny,Nx);
relL2 = zeros(1, Nt);

alpha = sqrt(omega/(2*nu));
yc = linspace(dy/2,Ly -dy/2,Ny)';
% Figure
figure('Name', 'Oscillating Shear Flow', 'NumberTitle', 'off');

subplot(2,2,1);
hQuiver = quiver(X, Y, zeros(size(X)), zeros(size(Y)), 'AutoScaleFactor', 3);
title('Velocity Field'); xlabel('x'); ylabel('y'); axis equal tight;


subplot(2,2,2);
hAnalytical = plot(zeros(Ny,1), yc, 'r-', 'LineWidth', 2); hold on;
hNumerical = plot(zeros(Ny,1), yc, 'bo--');
title('Velocity Profile'); xlabel('u velocity'); ylabel('y'); legend('Analytical', 'Numerical','Location','south'); grid on;

subplot(2,2,3);
hError = plot(zeros(Ny,1), yc, 'k');
title('Error Distribution'); xlabel('Error'); ylabel('y'); grid on;

subplot(2,2,4);
hConv = plot(0, 0, 'b');
title('Convergence of Relative L2 Error'); xlabel('Time Step'); ylabel('Relative L2 Error'); grid on;



% %% Setup video writer
% k = VideoWriter('CouetteFlowSimulation.mp4', 'MPEG-4'); % Name of video file
% k.FrameRate = 10; % Frames per second
% open(k);

%% Time loop
for n = 1:Nt
    t = time(n);


    u(1,:) = 0;
    u(end,:) = U0 * sin(omega * t);
    v([1 end], :) = 0;
    v(:, [1 end]) = 0;

    [u_star, v_star] = explicit_predictor(u, v, p, rho, nu, dt, dx, dy);

    rhs = divergence(u_star, v_star, dx, dy) / dt;
    

    p_corr = poisson_solver(rhs, dx, dy);

    [u, v] = rhie_chow_projection(u_star, v_star, p_corr, rho, dt, dx, dy);
    p = p + p_corr;


    u_analytical = U0 * exp(-alpha * flip(yc'))' .* sin(omega*t - alpha * flip(yc'))';
    u_center = u(:, round(Nx/2));
    err = abs(u_center - u_analytical);
    relL2(n) = norm(err) / norm(u_analytical);

    if mod(n, 10) == 0 || n == 1 || n == Nt
        set(hQuiver, 'UData', u, 'VData', v);
        set(hAnalytical, 'YData', yc, 'XData', u_analytical);
        set(hNumerical, 'YData', yc, 'XData', u_center);
        set(hError, 'XData', err, 'YData', yc);
        set(hConv, 'XData', 1:n, 'YData', relL2(1:n));
        drawnow;

        % % - Capture frame for video -
        % frame = getframe(gcf);
        % writeVideo(k, frame);
    end
end



% %% Close video file
% close(k);
% disp('Video saved successfully as CouetteFlowSimulation.mp4');

% Helper Functions 
function [u_corr, v_corr] = rhie_chow_projection(u_star, v_star, p_corr, rho, dt, dx, dy)
    [Ny, Nx] = size(u_star);
    u_corr = u_star;
    v_corr = v_star;

    dpdx = zeros(Ny, Nx);
    dpdy = zeros(Ny, Nx);

    dpdx(:,2:Nx-1) = (p_corr(:,3:Nx) - p_corr(:,1:Nx-2)) / (2*dx);
    dpdy(2:Ny-1,:) = (p_corr(3:Ny,:) - p_corr(1:Ny-2,:)) / (2*dy);

    u_corr = u_star - dt / rho * dpdx;
    v_corr = v_star - dt / rho * dpdy;
end

function p = poisson_solver(rhs, dx, dy)
    % V-cycle solver with gs smoothing
    maxLevel = floor(log2(min(size(rhs)))) - 1;
    p = fmg_vcycle(rhs, dx, dy, maxLevel);
end

function p = fmg_vcycle(rhs, dx, dy, level)
    if level == 0
        p = zeros(size(rhs));
        p = gauss_seidel(p, rhs, dx, dy, 150);
    else
        coarse_rhs = restrict(rhs);
        coarse_p = fmg_vcycle(coarse_rhs, 2*dx, 2*dy, level - 1);
        fine_p = prolong(coarse_p);
        fine_p = gauss_seidel(fine_p, rhs, dx, dy, 150);
        p = fine_p;
    end
end

function out = gauss_seidel(p, rhs, dx, dy, iterations)
    [Ny, Nx] = size(p);
    dx2 = dx^2; dy2 = dy^2;
    denom = 2*(dx2 + dy2);
    for iter = 1:iterations
        for j = 2:Ny-1
            for i = 2:Nx-1
                p(j,i) = ((p(j,i+1) + p(j,i-1))*dy2 + (p(j+1,i) + p(j-1,i))*dx2 - rhs(j,i)*dx2*dy2) / denom;
            end
        end
    end
    out = p;
end

function coarse = restrict(fine)
    coarse = fine(1:2:end, 1:2:end);
end

function fine = prolong(coarse)
    [Ny, Nx] = size(coarse);
    fine = zeros(2*Ny, 2*Nx);
    fine(1:2:end, 1:2:end) = coarse;
    fine(2:2:end, 1:2:end) = coarse;
    fine(1:2:end, 2:2:end) = coarse;
    fine(2:2:end, 2:2:end) = coarse;
end
function div = divergence(u, v, dx, dy)
    div = (u(:,[2:end end]) - u(:,[1 1:end-1])) / (2*dx) + ...
          (v([2:end end],:) - v([1 1:end-1],:)) / (2*dy);
end

function [u_star, v_star] = explicit_predictor(u, v, p, rho, nu, dt, dx, dy)
    [Ny, Nx] = size(u);
    u_star = u;
    v_star = v;

    %Laplacians
    uxx = zeros(Ny, Nx); uyy = zeros(Ny, Nx);
    vxx = zeros(Ny, Nx); vyy = zeros(Ny, Nx);

    uxx(:,2:Nx-1) = (u(:,3:Nx) - 2*u(:,2:Nx-1) + u(:,1:Nx-2)) / dx^2;
    uyy(2:Ny-1,:) = (u(3:Ny,:) - 2*u(2:Ny-1,:) + u(1:Ny-2,:)) / dy^2;

    vxx(:,2:Nx-1) = (v(:,3:Nx) - 2*v(:,2:Nx-1) + v(:,1:Nx-2)) / dx^2;
    vyy(2:Ny-1,:) = (v(3:Ny,:) - 2*v(2:Ny-1,:) + v(1:Ny-2,:)) / dy^2;

    %Pressure grads
    px = zeros(Ny, Nx); py = zeros(Ny, Nx);
    px(:,2:Nx-1) = (p(:,3:Nx) - p(:,1:Nx-2)) / (2*dx);
    py(2:Ny-1,:) = (p(3:Ny,:) - p(1:Ny-2,:)) / (2*dy);

    u_star = u + dt * (-px/rho + nu * (uxx + uyy));
    v_star = v + dt * (-py/rho + nu * (vxx + vyy));
end
