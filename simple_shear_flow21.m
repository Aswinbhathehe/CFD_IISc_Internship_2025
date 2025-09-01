clc; clear;close all;

% Domain and grid setup
x= 32;
y = 32;
Lx =1;
Ly = 0.2;
dx = Lx/x;
dy =Ly / y;
visc = 0.1;

% Boundary conditions 
Ut = 10; % Top wall velocity
Ub = 0;  % Bottom wall velocity

% Time step based on stability
CFL =0.3;
u_max = max(abs([Ut, Ub]));
dt1 =1e6;
dt2 = CFL * min(dx, dy) / u_max;
dt = min(dt1, dt2)
% dt = 0.000025; % chosen small for stability

% Preallocation
p  = zeros(y+2, x+2);
u = zeros(y+2, x+2);
v  =zeros(y+2, x+2);
ut =zeros(y+2, x+2);
vt= zeros(y+2, x+2);
divut =zeros(y+2, x+2);

% Analytical solution for simple shear flow (linear profile)
ya = linspace(dy/2,Ly -dy/2,y);
ua = (Ut -Ub)/Ly* ya +Ub;


%% Mesh for plotting
[X, Y] = meshgrid(dx/2:dx:Lx - dx/2, dy/2:dy:Ly - dy/2);

%% Setup figure with subplots
figure('Name', 'Simple shear flow between parallel plates', 'NumberTitle', 'off');

subplot(2,3,1);
hQuiver = quiver(X, Y, zeros(size(X)), zeros(size(Y)), 'AutoScaleFactor', 3);
title('Velocity Field');
xlabel('x'); ylabel('y'); axis equal tight;

subplot(2,3,2);
hIm = imagesc(linspace(0,Lx,x), linspace(0,Ly,y), zeros(y,x));
colorbar;
title('Divergence');
xlabel('x'); ylabel('y');
axis equal tight; grid on;

subplot(2,3,3);
hProfile = plot(ua, ya, 'r-', 'LineWidth', 2); hold on;
hNumerical = plot(ua*0, ya, 'bo--');
title('Velocity Profile');
xlabel('u velocity'); ylabel('y');
legend('Analytical','Numerical');
grid on;

subplot(2,3,5);
hError = plot(zeros(y,1), ya, 'k');
title('Error Distribution');
xlabel('Error'); ylabel('y');
grid on;

subplot(2,3,6);
hConv = plot(0, 0, 'b');
title('Convergence of Relative L2 Error');
xlabel('Time Step'); ylabel('Relative L2 Error');
grid on;




% %% Setup video writer
% k = VideoWriter('CouetteFlowSimulation.mp4', 'MPEG-4'); % Name of video file
% k.FrameRate = 10; % Frames per second
% open(k);



%% Time-stepping loop
tsteps =2000;
L2_Err_hist = zeros(tsteps, 1);

for n = 1:tsteps

    [u, v,~] =apply_bcs(u, v, Ut, Ub,p);

    %  X momentum 
    for i = 3:x+1
        for j = 2:y+1
            ue = 0.5 * (u(j, i+1) + u(j, i));
            uw = 0.5 * (u(j, i) + u(j, i-1));
            un = 0.5 * (u(j+1, i) + u(j, i));
            us = 0.5 * (u(j, i) + u(j-1, i));
            vn = 0.5 * (v(j+1, i-1) + v(j+1, i));
            vs = 0.5 * (v(j, i-1) + v(j, i));

            if n > 50
    convection = 0;
else
    convection = -(ue^2 - uw^2)/dx - (un*vn - us*vs)/dy;
end

            diffusion = visc * ((u(j, i-1) - 2*u(j, i) + u(j, i+1))/dx^2 + ...
                                (u(j-1, i) - 2*u(j, i) + u(j+1, i))/dy^2);

            ut(j, i) = u(j, i) + dt * (convection + diffusion);
        end
    end

    %  Y momentum 
    for i = 2:x+1
        for j = 3:y+1
            ve = 0.5 * (v(j, i+1) + v(j, i));
            vw = 0.5 * (v(j, i) + v(j, i-1));
            ue = 0.5 * (u(j, i+1) + u(j-1, i+1));
            uw = 0.5 * (u(j, i) + u(j-1, i));
            vn = 0.5 * (v(j+1, i) + v(j, i));
            vs = 0.5 * (v(j, i) + v(j-1, i));

            if n > 50
    convection = 0;
else
    convection = -(ue*ve - uw*vw)/dx - (vn^2 - vs^2)/dy;
            end
            diffusion = visc * ((v(j, i+1) - 2*v(j, i) + v(j, i-1))/dx^2 + ...
                                (v(j+1, i) - 2*v(j, i) + v(j-1, i))/dy^2);
            vt(j, i) =v(j, i) + dt *(convection +diffusion);
        end
    end

    % pressre correction
    rho = 1;
    divut(2:end-1, 2:end-1) = (ut(2:end-1,3:end) - ut(2:end-1,2:end-1))/dx + ...
                              (vt(3:end,2:end-1) - vt(2:end-1,2:end-1))/dy;

    rhs = rho * divut / dt;
    % [p, ~] = sor_solver(p, rhs, Lx, Ly, x, y);

    % Multigrid_solving=true;
    % 
    % if Multigrid_solving==true
    %     [p, ~] = multigrid_solver(p, rhs, Lx, Ly, x, y);
    % else
    %     [p, ~] = sor_solver(p, rhs, Lx, Ly, x, y);
    % end
    
    [p, ~] = fmg_solver(rhs, Lx, Ly, x, y);


    [~,~,p]=apply_bcs(u, v, Ut, Ub,p);

    %corner smoothining
    u(end,end) = mean([u(end-1,end), u(end,end-1)]);
    v(end,end) = mean([v(end-1,end), v(end,end-1)]);
    
    % Bottom-right corner
    u(1,end) = mean([u(2,end), u(1,end-1)]);
    v(1,end) = mean([v(2,end), v(1,end-1)]);
    
    % Top-left corner
    u(end,1) = mean([u(end-1,1), u(end,2)]);
    v(end,1) = mean([v(end-1,1), v(end,2)]);
    
    % Bottom-left corner
    u(1,1) = mean([u(2,1), u(1,2)]);
    v(1,1) = mean([v(2,1), v(1,2)]);
        
    p(2,2)=0; % Set reference pressure point to zero to prevent accidental existence of pressure gradients by only enforcing neumann equations
  %velocity correction
    [u, v] = rhie_chow_correction(ut, vt, p, dx, dy, dt);
    
    % --- Under-relaxation ---
    alpha = 0.4;

    
    u= alpha* u  +(1- alpha) *ut;
    v = alpha*v + (1 -alpha)* vt;


    % to display the velocity at the geometric center of the cell for aesthetic and accurate plotting.
    uc = 0.5*(u(2:end-1,2:end-1)+ u(2:end-1,  3:end));
    vc = 0.5 *(v(2:end-1,2:end-1) +v(3:end,2:end-1));

    u_profile =mean(uc, 2); % average across x-direction

    % Error
    error_profile =u_profile -ua';
    L2_error= sqrt(sum(error_profile.^2) /length(error_profile)) / ...
        sqrt(sum(ua.^2) /length(ua));
    L2_Err_hist(n) = L2_error;


    %  Update plots every 10 steps 
    if mod(n,10) == 0
        set(hQuiver, 'UData', uc, 'VData', vc);
        set(hIm, 'CData', flipud(divut(2:end-1,2:end-1)));

        set(hNumerical, 'XData', u_profile, 'YData', ya);
        set(hError, 'XData', error_profile, 'YData', ya);
        set(hConv, 'XData', 1:n, 'YData', L2_Err_hist(1:n));

        subplot(2,3,1); title(['Velocity Field Step ', num2str(n)]);
        subplot(2,3,2); title(['Divergence Step ', num2str(n)]);
        subplot(2,3,3); title('Velocity Profile');
        subplot(2,3,5); title('Error Distribution');
        subplot(2,3,6); title('Convergence');

        drawnow;

        % %  Capture frame for video 
        % frame = getframe(gcf);
        % writeVideo(k, frame);
    end

end


% %% Close video file
% close(k);
% disp('Video saved successfully as CouetteFlowSimulation.mp4');

fprintf('Max divergence: %.2e\n',max(abs(divut(:))));


function [u, v,p] = apply_bcs(u, v,Ut,Ub,p)
    % Left wall
    u(:,1) = u(:,2);
    v(:,1) = 0;

    % Right wall
    u(:,end) = u(:,end-1);
    v(:,end) = 0;

    % Top wall (Moving wall)
    u(end,:) = Ut;
    v(end,:) = 0;

    % Bottom wall (Fixed wall)
    u(1,:) = Ub;
    v(1,:) = 0;


  u(:,end) = u(:,end-1);
    v(:,end) = v(:,end-1); % convective (zero-gradient) outflow BCs on velocity
        
    % Apply Neumann BCs on pressure
    p(:,1)   = p(:,2);
    p(:,end) = p(:,end-1);
    p(1,:)   = p(2,:);
    p(end,:) = p(end-1,:);
    
end


function [p,err]=sor_solver(p, S,Lx,Ly,x, y)
    dx=Lx/x ;
    dy=Ly/y ;
    Ae = ones(y+2, x+2) / dx^2;
    Aw = ones(y+2, x+2) / dx^2;
    An = ones(y+2, x+2) / dy^2;
    As = ones(y+2, x+2) / dy^2;
    Ap=-(Ae+Aw+An+As);
    
    it = 0;
    err = 1e10;
    tol = 1e-8;
    maxit=1000;
    B = 1.9; % between 1 and 2

    while err > tol && it < maxit
        pk = p;
        for i =2:x+1
            for j =2:y+1
               ap = Ap(j,i); ae = Ae(j,i); aw = Aw(j,i); an = An(j,i); as = As(j,i);

               pe = p(j,i+1); pw = p(j,i-1); pn = p(j+1,i); ps = p(j-1,i);

               res = S(j,i) - (ae*pe + aw*pw + an*pn + as*ps);
               p(j,i) = B * res / ap + (1-B) * pk(j,i);
            end
        end
        u  = zeros(y+2, x+2);
        v  = zeros(y+2, x+2);
        Ut=10;
        Ub=0;
        [~,~,p]=apply_bcs(u, v, Ut, Ub,p);  % applying pressure BCs to prevent pressure drift in Neumann BCs problems
        err = norm(p(:) - pk(:), 2);
        it = it+1;
    end
end
 


function [p, err] = multigrid_solver(p, rhs, Lx, Ly, Nx, Ny)
    max_iter = 100;
    tol = 1e-8;
    err = 1e10;

    for iter = 1:max_iter
        p_old = p;
        p = V_cycle(p, rhs, Lx, Ly, Nx, Ny);
        res = res_computing(p, rhs, Lx, Ly, Nx, Ny);
        err = norm(res(:), 2);

        if err < tol
            break;
        end
    end
end

% V-Cycle Function
function p = V_cycle(p, rhs, Lx, Ly, Nx, Ny)
    if Nx <= 6 || Ny <= 6
        p = sor_solver_local(p, rhs, Lx, Ly, Nx, Ny, 100);
        return;
    end
    p = sor_solver_local(p, rhs, Lx, Ly, Nx, Ny, 10);
    res = res_computing(p, rhs, Lx, Ly, Nx, Ny);
    res_coarse = restrict(res);

    Nc_x = size(res_coarse,2) - 2;
    Nc_y = size(res_coarse,1) - 2;
    e_coarse = zeros(Nc_y + 2, Nc_x + 2);

    e_coarse = V_cycle(e_coarse, res_coarse, Lx, Ly, Nc_x, Nc_y);

    e_fine = prolong(e_coarse, Nx, Ny);
    p = p + e_fine;

    % Post-smoothing
    p = sor_solver_local(p, rhs, Lx, Ly, Nx, Ny, 7);
end

function p = sor_solver_local(p, rhs, Lx, Ly, Nx, Ny, Niter)
    if Nx < 2 || Ny < 2
        warning('SOR skipped due to small grid');
        return;
    end

    dx = Lx / Nx;
    dy = Ly / Ny;
    B = 1.9;

    for iter = 1:Niter
        p_old = p;
        for i = 2:Nx+1
            for j = 2:Ny+1
                p(j,i) = (1-B)*p(j,i) +B*0.5*((dy^2*(p(j,i+1) +p(j,i-1)) + dx^2*(p(j+1,i)+ p(j-1,i)) -dx^2*dy^2 * rhs(j,i))/ (2*(dx^2 + dy^2)));
            end
        end

        % Neumann BCs
        p(:,1)   = p(:,2);
        p(:,end) = p(:,end-1);
        p(1,:)   = p(2,:);
        p(end,:) = p(end-1,:);

        if norm(p(:) - p_old(:), 2) < 1e-8
            break;
        end
    end
end




function res = res_computing(p, rhs, Lx, Ly, Nx, Ny)
    dx = Lx / Nx;
    dy = Ly / Ny;

    res = zeros(Ny+2, Nx+2);
    for i = 2:Nx+1
        for j = 2:Ny+1
            laplace = (p(j,i+1) - 2*p(j,i) + p(j,i-1)) / dx^2 + ...
                      (p(j+1,i) - 2*p(j,i) + p(j-1,i)) / dy^2;
            res(j,i) = rhs(j,i) - laplace;
        end
    end
end


function coarse = restrict(fine)
    [Nyf, Nxf] = size(fine);
    Nxc = ceil((Nxf - 2)/2);
    Nyc = ceil((Nyf - 2)/2);

    coarse = zeros(Nyc+2, Nxc+2);
    for i = 2:Nxc+1
        for j = 2:Nyc+1
            i_f = 2*(i-1);
            j_f = 2*(j-1);

            neighbors = fine(j_f-1:j_f+1, i_f-1:i_f+1);
            weights = [1 2 1; 2 4 2; 1 2 1];

            % Handle edges
            valid = ~isnan(neighbors);
            w_sum = sum(weights(valid));

            coarse(j,i) = sum(neighbors(valid) .* weights(valid), 'all') / w_sum;
        end
    end
end

% Prolongation 
function fine = prolong(coarse, Nxf, Nyf)
    Nxc = size(coarse, 2) - 2;
    Nyc = size(coarse, 1) - 2;

    fine = zeros(Nyf+2, Nxf+2);

    for i = 2:Nxc+1
        for j = 2:Nyc+1
            i_f = 2 * (i - 1);
            j_f = 2 * (j - 1);

            % Safely assign values to fine grid 
            if j_f <= Nyf && i_f <= Nxf
                fine(j_f,   i_f)   = fine(j_f,   i_f)   + coarse(j, i);
            end
            if j_f + 1 <= Nyf && i_f <= Nxf
                fine(j_f+1, i_f)   = fine(j_f+1, i_f)   + coarse(j, i);
            end
            if j_f <= Nyf && i_f + 1 <= Nxf
                fine(j_f,   i_f+1) = fine(j_f,   i_f+1) + coarse(j, i);
            end
            if j_f + 1 <= Nyf && i_f + 1 <= Nxf
                fine(j_f+1, i_f+1) = fine(j_f+1, i_f+1) + coarse(j, i);
            end
        end
    end

    % average overlapping contributions
    fine(2:end-1, 2:end-1) = fine(2:end-1, 2:end-1) / 4;
end


function [u_corr, v_corr] = rhie_chow_correction(ut, vt, p, dx, dy, dt)
    % appliyng Rhie-Chow interpolation based velocity correction   
    [Ny, Nx] = size(p);
    u_corr = ut;
    v_corr = vt;
    
    %  u  correction 
    u_corr(2:end-1,3:end-1) = ut(2:end-1,3:end-1) - ...
        dt * (p(2:end-1,3:end-1) - p(2:end-1,2:end-2)) / dx;
    
    %  vcorrection 
    v_corr(3:end-1,2:end-1) = vt(3:end-1,2:end-1) - ...
        dt * (p(3:end-1,2:end-1) - p(2:end-2,2:end-1)) / dy;
end

% FMG Solver
function [p, err] = fmg_solver(rhs, Lx, Ly, Nx, Ny)
    levels = floor(log2(min(Nx, Ny))) - 1;
    levels = min(levels, 5);  % Optional hard cap to avoid too deep levels
    if levels < 1
        warning('FMG: Too few grid levels, using regular multigrid.');
        [p, err] = multigrid_solver(zeros(size(rhs)), rhs, Lx, Ly, Nx, Ny);
        return;
    end

    % Coarsest grid size
    Nc_x = floor(Nx / 2^(levels - 1));
    Nc_y = floor(Ny / 2^(levels - 1));

    if Nc_x < 3 || Nc_y < 3
        [p, err] = multigrid_solver(zeros(size(rhs)), rhs, Lx, Ly, Nx, Ny);
        return;
    end

    % Construct RHS at coarsest level
    rhs_c = rhs;
    for l = 1:(levels - 1)
        rhs_c = restrict(rhs_c);
    end
    p_c = zeros(size(rhs_c));
    [Nc_y_full, Nc_x_full] = size(rhs_c);
    p_c = sor_solver_local(p_c, rhs_c, Lx, Ly, Nc_x_full - 2, Nc_y_full - 2, 100);

    for l = (levels - 1):-1:0

        Nxf = floor(Nx / 2^l);
        Nyf = floor(Ny / 2^l);
        p_f = prolong(p_c, Nxf, Nyf);


        rhs_f = rhs;
        for li = 1:l
            rhs_f = restrict(rhs_f);
        end

        p_c =V_cycle(p_f, rhs_f, Lx, Ly, Nxf, Nyf);
    end

    p = p_c;

    res =res_computing(p, rhs, Lx, Ly, Nx, Ny);
    err =norm(res(:), 2);
end
