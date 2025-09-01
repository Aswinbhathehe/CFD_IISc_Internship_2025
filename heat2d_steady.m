clc; clear; close all;
% Note:Code works only for temperatures applied to top or bottom wall. we have to change analytical solution if different BCs is applied

% Domain and Grid parameters
Lx = 1;
Ly = 1;
Nx = 400;
Ny = 400;
dx = Lx/Nx;
dy = Ly/Ny;

x = linspace(dx/2, Lx-dx/2, Nx);
y = linspace(dy/2, Ly-dy/2, Ny);
[X, Y] = meshgrid(x, y);

% Boundary Conditions (change according to user input)
BC.T_top = 100;
BC.T_bottom = 0;
BC.T_left = 0;
BC.T_right = 0;

% Source Term (Poisson Equation)
% Set f = 0 for Laplace
f = zeros(Ny, Nx);

% Solver Options
solver_type = 'Sparse';   % SOR or Sparse_solver

% Solver Parameters
omega = 2 / (1 + sin(pi/max(Nx, Ny))); % Optimal omega
max_iter = 10000;
tol = 1e-12;

% Solve
if strcmpi(solver_type, 'SOR')
    [numerical_temp, iteration_count, error_history] = SOR_solver(f, BC, dx, dy, omega, max_iter, tol);
elseif strcmpi(solver_type, 'Sparse')
    numerical_temp = Sparse_solver(f, BC, dx, dy);
    error_history = [];
    iteration_count = 1;
else
    error('Unknown solver type. Choose "SOR" or "Sparse".');
end


% Analytical Solution
analytical_temp = zeros(Ny, Nx);
N_terms = 100; % Number of sine terms

for n = 1:2:(2*N_terms-1) % Only odd terms
    lambda= n*pi;
    term = (4*BC.T_top)/(n*pi) *sinh(n*pi*Y) ./ sinh(n*pi*Ly) .*sin(n*pi*X);
    analytical_temp =analytical_temp  + term;
end

%Error Calc

error_abs= abs(numerical_temp- analytical_temp);
L2_error =sqrt(sum(error_abs(:).^2)) /sqrt(sum(analytical_temp(:).^2));

fprintf('L2 Relative Error = %.6e\n', L2_error);

% Plotting (Includes error plot)

figure('Position',[100 100 1200 400]);

subplot(1,3,1);
contourf(X, Y, numerical_temp, 50, 'LineColor','none');
colorbar;
title(['Numerical Solution (' solver_type ')']);
xlabel('X');
ylabel('Y');
axis equal tight;

subplot(1,3,2);
contourf(X, Y, analytical_temp, 50, 'LineColor','none');
colorbar;
title('Analytical Solution');
xlabel('X');
ylabel('Y');
axis equal tight;

subplot(1,3,3);
contourf(X, Y, error_abs ./ abs(analytical_temp + eps), 50, 'LineColor','none'); 
%+eps to avoid divide by zero
colorbar;
title('L2 Relative Error Distribution');
xlabel('X');
ylabel('Y');
axis equal tight;


% Plot Convergence History
if strcmpi(solver_type, 'SOR')
    figure;
    semilogy(error_history,'-o');
    grid on;
    xlabel('Iteration');
    ylabel('Max Error');
    title('SOR Convergence History');
end



% Helper Functions

% Dirichlet BCs (for neumann different analytical solution)
function T = apply_BC(T, BC)
    T(1,:)   = BC.T_top;
    T(end,:) = BC.T_bottom;
    T(:,1)   = BC.T_left;
    T(:,end) = BC.T_right;
end

% SOR Solver 
function [T, iter, error_history] = SOR_solver(f, BC, dx, dy, omega, max_iter, tol)
    [Ny, Nx] = size(f);
    T = zeros(Ny, Nx);
    T = apply_BC(T, BC);

    dx2 = dx^2;
    dy2 = dy^2;
    coeff = 1/(2*(dx2 + dy2));

    error_history = [];

    for iter = 1:max_iter
        T_old = T;
        T(2:end-1,2:end-1) = (1-omega)*T(2:end-1,2:end-1) +omega *coeff *((T(2:end-1,3:end)+ T(2:end-1,1:end-2))*dy2 +(T(3:end,2:end-1)+T(1:end-2,2:end-1))*dx2 -f(2:end-1,2:end-1)*dx2*dy2 );
        T = apply_BC(T, BC);

        % Error Check
        err = max(max(abs(T - T_old)));
        error_history = [error_history; err];

        if err < tol
            fprintf('SOR converged in %d iterations with error %.3e\n', iter, err);
            break
        end
    end

    if iter == max_iter
        fprintf('SOR reached max iterations (%d) with error %.3e\n', iter, err);
    end
end

%Sparse Matrix Solver
function T = Sparse_solver(f, BC, dx, dy)
    [Ny, Nx] = size(f);
    N = Ny * Nx;
    dx2 = dx^2;
    dy2 = dy^2;

    % Sparse Matrix Assembly
    main_diag = -2*(1/dx2 + 1/dy2) * ones(N,1);
    off_diag_x = 1/dx2 * ones(N,1);
    off_diag_y = 1/dy2 * ones(N,1);

    A = spdiags([off_diag_y, off_diag_x, main_diag, off_diag_x, off_diag_y], [-Nx,-1,0,1,Nx], N, N);
    b = -reshape(f,[],1);

    % applying Dirichlet BCs
    top_idx    = (Ny-1)*Nx + (1:Nx);
    bottom_idx = (0)*Nx + (1:Nx);
    left_idx   = ((Ny-1):-1:0)*Nx + 1;
    right_idx  = ((Ny-1):-1:0)*Nx + Nx;

    A(top_idx,:)    = 0; A(sub2ind(size(A), top_idx, top_idx)) = 1; b(top_idx) = BC.T_top;
    A(bottom_idx,:) = 0; A(sub2ind(size(A), bottom_idx, bottom_idx)) = 1; b(bottom_idx) = BC.T_bottom;
    A(left_idx,:)   = 0; A(sub2ind(size(A), left_idx, left_idx)) = 1; b(left_idx) = BC.T_left;
    A(right_idx,:)  = 0; A(sub2ind(size(A), right_idx, right_idx)) = 1; b(right_idx) = BC.T_right;

    T_vec = A\b;

    % Reshape to 2D
    T = reshape(T_vec, [Nx, Ny])';
end
