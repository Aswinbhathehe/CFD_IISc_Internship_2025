clc; clear; close all

% Domain and grid setup
x=32; y=32;
Lx=1; Ly=0.2;
dx=Lx/x; dy=Ly/y;
visc=0.01;
Ut=0; Ub=0;
dt=min(1, 0.5*dx^2/visc);

rho=1;               % density
dpx=-1;              % pressure gradient(dp/dx -ve then flow in +ve)
[X, Y]=meshgrid(dx/2:dx:Lx-dx/2, dy/2:dy:Ly-dy/2);
u=zeros(y+2, x+2); v=zeros(y+2, x+2);
ut=u; vt=v; p=zeros(y+2, x+2);


ya=linspace(0, Ly, y);

ua=( -dpx/(2*visc) )*ya .* (Ly - ya); 
Nb=100; s=linspace(0, 1, Nb);

% Time loop
tsteps=300; L2_error_history=zeros(tsteps, 1);
FxL_old=zeros(1, Nb); FyL_old=zeros(1, Nb);
alpha=0.5;

for n=1:tsteps
    [u, v, ~]=applybcs(u,v, Ut, Ub, p);

    fx=ones(y+2, x+2)*(-dpx)/rho;  % uniform body force in x
    fy=zeros(y+2, x+2);            % no vertical force

Nx_u=x+1;   % for u: faces in x
Ny_u=y;

u_staggered=u(2:end-1, 2:end);     % 32×33
fx_staggered=fx(2:end-1, 2:end);   % 32×33

u_new=implicit_diffusion_u(u_staggered+dt*fx_staggered, Lx, Ly, Nx_u, Ny_u, dt, visc);
u(2:end-1, 2:end)=u_new;
ut=u;
    v_center=v(2:end-1, 2:end-1);
    fy_center=fy(2:end-1, 2:end-1);
    v_center_new=implicit_diffusion_v(v_center+dt*fy_center, Lx, Ly, x, y, dt, visc);
    v(2:end-1, 2:end-1)=v_center_new;
    vt=v;

    divut=(ut(2:end-1,3:end) - ut(2:end-1,2:end-1))/dx+...
            (vt(3:end,2:end-1) - vt(2:end-1,2:end-1))/dy;

    rhs_full=zeros(y+2, x+2);
    rhs_full(2:end-1, 2:end-1)=divut/dt;

    for cycle=1:30
        p_old=p;
        p=V_cycle(p, rhs_full, Lx, Ly, x, y);
        if norm(p(:) - p_old(:), 2)/norm(p_old(:), 2) < 1e-8
            break;
        end
    end

    [u, v]=rhie_chow_correction(ut, vt, p, dx, dy, dt);

    uc=0.5*(u(2:end-1,2:end-1)+u(2:end-1,3:end));
    u_profile=mean(uc, 2);
    error_profile=u_profile - ua';
    L2_error=sqrt(sum(error_profile.^2)/length(error_profile))/sqrt(sum(ua.^2)/length(ua));
    L2_error_history(n)=L2_error;

    if mod(n,10) == 0
        fprintf("Step %d, L2 Error=%.2e\n", n, L2_error);

        x_index=round(x/2)+1;
        uc_current=0.5*(u(2:end-1,2:end-1)+u(2:end-1,3:end));

        figure(1); clf;
        subplot(1,2,1);
        plot(uc_current(:, x_index), ya, 'b-', 'LineWidth', 2); hold on;
        plot(ua, ya, 'r--', 'LineWidth', 2);
        xlabel('u'); ylabel('y'); title(['Velocity Profile at x=Lx/2, Step=', num2str(n)]);
        legend('Numerical', 'Analytical','Location','south'); grid on;

        subplot(1,2,2);
        semilogy(1:n, L2_error_history(1:n), 'k-', 'LineWidth', 2);
        xlabel('Time step'); ylabel('Relative L2 Error');
        title('L2 Error Convergence'); grid on;
        drawnow;
    end
end

function [u, v, p]=applybcs(u, v, Ut, Ub, p)
    u(:,1)=u(:,2); 
    u(:,end)=u(:,end-1);

    v(:,1)=0; 
    v(:,end)=0;

    u(end,:)=Ut; 
    u(1,:)=Ub;
    
    v([1 end],:)=0;

    p(:,1)=p(:,2); 
    p(:,end)=p(:,end-1);
    p(1,:)=p(2,:); 
    p(end,:)=p(end-1,:);
end

function [uL, vL]=vel_interpolate(u, v, Xb, Yb, dx, dy)
    Nb=length(Xb); 
    uL=zeros(1,Nb); 
    vL=zeros(1,Nb);
    for k=1:Nb
        xk=Xb(k); yk=Yb(k);
        i0=floor(xk/dx)+1; j0=floor(yk/dy)+1;
        for i=i0-1:i0+2
            for j=j0-1:j0+2
                if i >= 1 && i <= size(u,2) && j >= 1 && j <= size(u,1)
                    xi=(i-1)*dx; yj=(j-1)*dy;
                    phi=delta_kernel((xk - xi)/dx)*delta_kernel((yk - yj)/dy);
                    uL(k)=uL(k)+u(j,i)*phi;
                    vL(k)=vL(k)+v(j,i)*phi;
                end
            end
        end
    end
end

function f=spread_force(FxL, FyL, Xb, Yb, dx, dy, Nx, Ny)
    f=zeros(Ny+2, Nx+2, 2); Nb=length(Xb);
    for k=1:Nb
        xk=Xb(k); yk=Yb(k);
        i0=floor(xk/dx)+1; j0=floor(yk/dy)+1;
        for i=i0-1:i0+2
            for j=j0-1:j0+2
                xi=(i-1)*dx; yj=(j-1)*dy;
                phi=delta_kernel((xk - xi)/dx)*delta_kernel((yk - yj)/dy);
                if i >= 1 && i <= Nx+2 && j >= 1 && j <= Ny+2
                    f(j,i,1)=f(j,i,1)+FxL(k)*phi*dx*dy;
                    f(j,i,2)=f(j,i,2)+FyL(k)*phi*dx*dy;
                end
            end
        end
    end
end

% 4-point kernel (Peskin’s standard) 
function val=delta_kernel(r) 
    r=abs(r);
    if r < 1
        val=0.125*(3 - 2*r+sqrt(1+4*r - 4*r^2));
    elseif r < 2
        val=0.125*(5 - 2*r - sqrt(-7+12*r - 4*r^2));
    else
        val=0;
    end
end



function u_new=implicit_diffusion_u(u_old, Lx, Ly, Nx, Ny, dt, nu)
    dx=Lx/Nx;
    dy=Ly/Ny;
    N=Nx*Ny;
    A=sparse(N, N);
    b=zeros(N, 1);
    coeff_center=1+2*dt*nu*(1/dx^2+1/dy^2);
    coeff_x=-dt*nu/dx^2;
    coeff_y=-dt*nu/dy^2;
    index=@(i,j) (j-1)*Nx+i;
    
    for j=1:Ny
        for i=1:Nx
            n=index(i,j);
            % top and bottom walls  Dirichilet BC
            if j == 1
                A(n,:)=0;
                A(n,n)=1;
                b(n)=0;       % Bottom wall velocity Ub=0
                continue;
            elseif j == Ny
                A(n,:)=0;
                A(n,n)=1;
                b(n)=0;      % Top wall velocity Ut=0
                continue;
            end
            A(n,n)=coeff_center;
            if i > 1
                A(n, index(i-1,j))=coeff_x;
            else
                A(n,n)=A(n,n) - coeff_x;  % Neumann
            end
            if i < Nx
                A(n, index(i+1,j))=coeff_x;
            else
                A(n,n)=A(n,n) - coeff_x;  % Neumann
            end
            A(n, index(i,j-1))=coeff_y;
            A(n, index(i,j+1))=coeff_y;
            b(n)=u_old(j,i);
        end
    end
    u_vec=A \b;
    u_new=reshape(u_vec, [Nx, Ny])';
end

% V-Cycle Function
function p=V_cycle(p, rhs, Lx, Ly, Nx, Ny)
    if Nx <= 4 || Ny <= 4
        p=sor_solve(p, rhs, Lx, Ly, Nx, Ny, 100);
        return;
    end

    p=sor_solve(p, rhs, Lx, Ly, Nx, Ny, 15);
    res=residualcalc(p, rhs, Lx, Ly, Nx, Ny);
    res_coarse=restrict(res);


    Nc_x=size(res_coarse,2) - 2;
    Nc_y=size(res_coarse,1) - 2;
    e_coarse=zeros(Nc_y+2, Nc_x+2);
    e_coarse=V_cycle(e_coarse, res_coarse, Lx, Ly, Nc_x, Nc_y);
    e_fine=prolong(e_coarse, Nx, Ny);
    p=p+e_fine;


    p=sor_solve(p, rhs, Lx, Ly, Nx, Ny, 10);
end


function p=sor_solve(p, rhs, Lx, Ly, Nx, Ny, Niter)
    dx=Lx/Nx;
    dy=Ly/Ny;
    B=1.9;

    for iter=1:Niter
        p_old=p;
        for i=2:Nx+1
            for j=2:Ny+1
                p(j,i)=(1-B)*p(j,i)+B*0.5*((dy^2*(p(j,i+1)+p(j,i-1))+ dx^2*(p(j+1,i)+p(j-1,i)) -dx^2*dy^2*rhs(j,i))/(2*(dx^2+dy^2)));
            end
        end

        % Neumann boundary conditions
        p(:,1)  =p(:,2);
        p(:,end)=p(:,end-1);
        p(1,:)  =p(2,:);
        p(end,:)=p(end-1);

        if norm(p(:) -p_old(:), 2) < 1e-8
            break;
        end
    end
end

function res=residualcalc(p,rhs, Lx, Ly, Nx, Ny)
    dx=Lx/Nx;
    dy=Ly/Ny;
    res=zeros(Ny+2, Nx+2);
    for i=2:Nx+1
        for j=2:Ny+1
            res(j,i)=rhs(j,i) -(p(j,i+1) - 2*p(j,i)+p(j,i-1))/dx^2+ (p(j+1,i) - 2*p(j,i)+p(j-1,i))/dy^2;
        end
    end
end
function coarse=restrict(fine)
    [Nyf, Nxf]=size(fine);
    Nxc=ceil((Nxf - 2)/2);
    Nyc=ceil((Nyf - 2)/2);
    coarse=zeros(Nyc+2, Nxc+2);
    for i=2:Nxc+1
        for j=2:Nyc+1
            i_f=2*(i-1);
            j_f=2*(j-1);
            neighbors=fine(j_f-1:j_f+1, i_f-1:i_f+1);
            weights=[1 2 1; 2 4 2; 1 2 1];
            valid=~isnan(neighbors);
            w_sum=sum(weights(valid));
            coarse(j,i)=sum(neighbors(valid) .* weights(valid), 'all')/w_sum;
        end
    end
end

% Prolongation bilinear
function fine=prolong(coarse, Nxf, Nyf)
    Nxc=size(coarse,2)- 2;
    Nyc=size(coarse,1)- 2;

    fine=zeros(Nyf+2, Nxf+2);

    for i=2:Nxc+1
        for j=2:Nyc+1
            i_f=2*(i-1);
            j_f=2*(j-1);

            fine(j_f, i_f)  =fine(j_f, i_f)  +coarse(j,i);
            fine(j_f+1, i_f)  =fine(j_f+1, i_f)  +coarse(j,i);
            fine(j_f,i_f+1)=fine(j_f, i_f+1)+coarse(j,i);
            fine(j_f+1, i_f+1)=fine(j_f+1, i_f+1)+coarse(j,i);
        end
    end
    fine(2:end-1,2:end-1)=fine(2:end-1,2:end-1)/4;
end


function [u_corr, v_corr]=rhie_chow_correction(ut, vt, p, dx, dy, dt)
    [Ny, Nx]=size(p);
    u_corr=ut;
    v_corr=vt; 
    u_corr(2:end-1,3:end-1)=ut(2:end-1,3:end-1) -dt*(p(2:end-1,3:end-1) - p(2:end-1,2:end-2))/dx;
    v_corr(3:end-1,2:end-1)=vt(3:end-1,2:end-1) -dt*(p(3:end-1,2:end-1) - p(2:end-2,2:end-1))/dy;
end


function v_new=implicit_diffusion_v(v_old, Lx, Ly, Nx, Ny, dt, nu)
    dx=Lx/Nx;
    dy=Ly/Ny;
    N=Nx*Ny;
    A=sparse(N, N);
    b=zeros(N, 1);
    coeff_center=1+2*dt*nu*(1/dx^2+1/dy^2);
    coeff_x=-dt*nu/dx^2;
    coeff_y=-dt*nu/dy^2;
    index=@(i,j) (j-1)*Nx+i;

    for j=1:Ny
        for i=1:Nx
            n=index(i,j);
            %top and bottom: Dirichilet(v=0)
            if j == 1 || j == Ny
                A(n,:)=0;
                A(n,n)=1;
                b(n)=0;
                continue;
            end
            A(n,n)=coeff_center;
            if i > 1
                A(n, index(i-1,j))=coeff_x;
            else
                A(n,n)=A(n,n) - coeff_x;
            end
            if i < Nx
                A(n, index(i+1,j))=coeff_x;
            else
                A(n,n)=A(n,n) - coeff_x;
            end
            A(n, index(i,j-1))=coeff_y;
            A(n, index(i,j+1))=coeff_y;
            b(n)=v_old(j,i);
        end
    end
    v_vec=A \ b;
    v_new=reshape(v_vec, [Nx, Ny])';
end

