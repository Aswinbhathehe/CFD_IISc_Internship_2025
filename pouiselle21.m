clc; clear; close all

% grid and domain properties while using staggered grid 
x=32; y=32;
Lx=1; Ly=0.2;
dx=Lx/x; dy=Ly/y;
visc=0.05;

%top and Bottom wall fixed
Ut=0; Ub=0;

% Time step selection based on stability criteria
dt1 = 0.5/(visc*(1/(dx^2) + 1/(dy^2)));
CFL = 0.5;
u_max = max(abs([Ut, Ub]));
if u_max == 0
    dt2 = Inf;
else
    dt2 = CFL * min(dx/u_max, dy/u_max);
end
dt = min(dt1, dt2);
dt = 0.25 * dt; % safety margin
dpx = 0.1; % pressure gradient

% Preallocate fields
p = zeros(y+2, x+2);
u = zeros(y+2, x+2);
v = zeros(y+2, x+2);
ut = zeros(y+2, x+2);
vt = zeros(y+2, x+2);
divut = zeros(y+2, x+2);

% Cell center velocities and meshgrid
uc=  0.5*(u(2:end-1,  2:end-1) +u(2:end-1,3:end));
vc =0.5*(v(2:end-1,2:end-1)+ v(3:end,2:end-1));
[X,Y] = meshgrid(dx/2:dx:Lx-dx/2, dy/2:dy:Ly-dy/2);

fig = figure('Name','Flow Simulation','NumberTitle','off');

subplot(3,2,1);
hQuiver = quiver(X,Y,uc,vc,'AutoScaleFactor',3);
title('Velocity Field'); xlabel('x'); ylabel('y'); axis equal tight;

% Initial divergence
for i=2:x+1
    for j=2:y+1
        divu(j,i) =(u(j,i)-u(j,i-1))/dx   +(v(j,i)-v(j-1,i))/dy;
    end
end
subplot(3,2,2);
hIm = imagesc(linspace(0,Lx,x), linspace(0,Ly,y), flipud(divu(2:end-1,2:end-1)));
colorbar; title('Divergence of Velocity'); xlabel('x'); ylabel('y'); axis equal tight; grid on;

t = 0; tsteps = 8000;
err_hist = []; t_hist = [];
% k = VideoWriter('poiselle_flow_vid.mp4', 'MPEG-4'); k.FrameRate = 20; open(k);
err_initial = NaN;

for n = 1:tsteps
    [u,v,~] = apply_boundary_conditions(u,v,p,dpx,dx);

    % X-momentum
    for i = 3:x+1
        for j = 2:y+1
            ue = 0.5*(u(j,i+1)+u(j,i)); 
            uw = 0.5*(u(j,i)  +u(j,i-1));
            un = 0.5*(u(j+1,i)+u(j,i)); 
            us = 0.5*(u(j,i)  +u(j-1,i));
            vn = 0.5*(v(j+1,i-1)+v(j+1,i));
            vs = 0.5*(v(j,i-1)  +v(j,i));
            convection = -(ue^2- uw^2)/dx-(un*vn- us*vs)/dy;
            diffusion = visc* ((u(j,i-1)-2*u(j,i)+u(j,i+1))/dx^2+ (u(j-1,i)-2*u(j,i)  +u(j+1,i))/dy^2);
            pressure_source = dpx; 
            rho=0.01;
            ut(j,i) = u(j,i)+ dt*(convection+ diffusion+  (pressure_source/ rho));
        end
    end

    % Y-momentum
    for i = 2:x+1
        for j = 3:y+1
            ve = 0.5*(v(j,i+1) +v(j,i)); 
            vw = 0.5*(v(j,i)   +v(j,i-1));
            ue = 0.5*(u(j,i+1) +u(j-1,i+1)); 
            uw = 0.5*(u(j,i)   +u(j-1,i));
            vn = 0.5*(v(j+1,i) +v(j,i)); 
            vs = 0.5*(v(j,i)   +v(j-1,i));
            convection = - (ue*ve-uw*vw)/dx -(vn^2-vs^2)/dy;
            diffusion = visc* ((v(j,i+1)-2*v(j,i)+v(j,i-1))/dx^2+   (v(j+1,i)-2*v(j,i)+ v(j-1,i))/dy^2);
            vt(j,i) = v(j,i) + dt*(convection+diffusion);
        end
    end


    divut(2:end-1,2:end-1) = (ut(2:end-1,3:end)-ut(2:end-1,2:end-1))/dx + (vt(3:end,2:end-1)-vt(2:end-1,2:end-1))/dy;
    rhs = rho * divut / dt;

        [p,~] = sor_solver(p,rhs,Lx,Ly,x,y);


    [~,~,p] = apply_boundary_conditions(u,v,p,dpx,dx);

    % Velocity correction
    u(2:end-1,3:end-1) =ut(2:end-1,3:end-1) -dt*(p(2:end-1,3:end-1)-p(2:end-1,2:end-2))/dx;
    v(3:end-1,2:end-1)=vt(3:end-1,2:end-1) - dt*(p(3:end-1,2:end-1) -p(2:end-2,2:end-1))/dy;

    for i=2:x+1
        for j=2:y+1
            divu(j,i) =(u(j,i)-u(j,i-1))/dx +(v(j,i)-v(j-1,i))/dy;
        end
    end

    uc = 0.5*  (u(2:end-1,2:end-1) +u(2:end-1,3:end));
    vc= 0.5*(v(2:end-1,2:end-1) +v(3:end,2:end-1));

    if mod(n,100)==0 || n==0
        subplot(3,2,1);
        set(hQuiver,'UData',uc,'VData',vc);
        title(['Velocity Field at step ',num2str(n)]);

        subplot(3,2,2);
        set(hIm,'CData',flipud(divu(2:end-1,2:end-1)));
        title(['Divergence at step ',num2str(n)]);

        ya = dy/2 : dy: Ly-dy/2;
        ua = (dpx/(2*visc)) *ya .*(Ly - ya)/rho;
        x_phys = linspace(dx/2,Lx-dx/2, x);
        [~, mid_idx] = min(abs(x_phys - Lx/2));
        u_num = uc(:,mid_idx);

        subplot(3,2,3);
        plot(u_num', ya, 'b-', ua, ya, 'r--', 'LineWidth',2);
        ylabel('y'); xlabel('u(y)'); legend('Numerical','Analytical');
        title('Velocity Profile'); grid on;

        err_profile = u_num - ua';
        subplot(3,2,4);
        plot(err_profile, ya, 'k-', 'LineWidth', 2);
        ylabel('y'); xlabel('Error'); title('Error Profile'); grid on;

        L2_err = norm(err_profile, 2);
        if isnan(err_initial)
            err_initial =L2_err;
        end
        err_rel = L2_err /err_initial;
        err_hist(end+1) = err_rel;
        t_hist(end+1) = t;

        subplot(3,2,[5 6]);
        semilogy(t_hist, err_hist, 'r-o', 'LineWidth', 1.5);
        xlabel('Time'); 
        ylabel('Relative L2 Error');
        title('Convergence History'); grid on;

        drawnow;
        % frame = getframe(gcf);
        % writeVideo(k, frame);
    end
    t = t + dt;
end

% close(k);
% disp('Video saved successfully');

figure('Name','Final Convergence','NumberTitle','off');
semilogy(t_hist, err_hist, 'r-o', 'LineWidth', 1.5);
xlabel('Time'); ylabel('Relative L2 Norm of Error');
title('Final Convergence of Numerical Solution'); grid on;

fprintf('Max divergence: %.2e\n', max(abs(divu(:))));

function [u, v, p] = apply_boundary_conditions(u, v, p,dpx,dx)
    % Left wall
    u(:,1) =u(:,2);
    v(:,1) =0;

    % Right wall
    u(:,end) =u(:,end-1);
    v(:,end) =0;

    % Top wall
    u(end,:) =0;
    v(end,:) =0;

    % Bottom wall
    u(1,:) =0;
    v(1,:) =0;

    % non zero-gradient boundary condition for pressure at boundaries
    p(:,1)     = p(:,2)-dpx*dx;        % left
    p(:,end)   = p(:,end-1)+dpx*dx;    % right

end



% poission solver SOR
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
    tol = 1e-12;
    maxit=20000;
    B = 1.9; % between 1 and 2

    while err > tol && it < maxit
        pk = p;
        for i =2:x+1
            for j =2:y+1
               ap = Ap(j,i); ae = Ae(j,i); aw = Aw(j,i); an = An(j,i); as = As(j,i);
               pe = p(j,i+1); pw = p(j,i-1); pn = p(j+1,i); ps = p(j-1,i);
               res = S(j,i) - (ae*pe + aw*pw + an*pn + as*ps);
               p(j,i) = B * res/ap + (1-B)*pk(j,i);
            end
        end
        err = norm(p(:)-pk(:),2);
        it = it+1;
    end
end
 
