clc; clear; close all;

Nx=256; 
Ny=64;
Lx=2.0; 
Ly=1.0;
dx=Lx/Nx; 
dy=Ly/Ny;
x=linspace(0, Lx, Nx);
y=linspace(0, Ly, Ny);
[X, Y]=meshgrid(x, y);

rho=1.0;
nu=0.001;
U_inf=10;
dt=0.001;
steps=3000;
beta=1.5;  % SOR over-relaxation factor

cx=1.0; 
cy=0.5; 
R=0.1;
Nb=100;
theta=linspace(0, 2*pi, Nb);
Xb=cx+R*cos(theta);
Yb=cy+R*sin(theta);

% Field includes ghost cells
u=U_inf * ones(Ny+2, Nx+2);
v=zeros(Ny+2, Nx+2);
p=zeros(Ny+2, Nx+2);
ut=u; vt=v;

FxL_old=zeros(1, Nb);
FyL_old=zeros(1, Nb);
alpha=0.5;
Fx_total=zeros(1, steps);
Fy_total=zeros(1, steps);

%Time loop
for n=1:steps
    [u, v, p]=apply_bc(u, v, p, U_inf);
    ut=GS_diffuse(u, nu, dt, dx, dy);
    vt=GS_diffuse(v, nu, dt, dx, dy);

    [uL, vL]=vel_interpol(ut, vt, Xb, Yb, dx, dy);
    epsilon=1000;
    FxL=-epsilon * uL;
    FyL=-epsilon * vL;
    FxL=alpha *FxL+(1-alpha) *FxL_old;
    FyL=alpha *FyL+(1-alpha) *FyL_old;
    FxL_old=FxL; 
    FyL_old=FyL;
    Fx_total(n)=sum(FxL);  %drag_force
    Fy_total(n)=sum(FyL);  %lift force

    [fx, fy]= spreadf(FxL, FyL, Xb, Yb, Nx, Ny, dx, dy);
    ut=ut+dt * fx/rho;
    vt=vt+dt * fy/rho;
    div=((ut(2:end-1,3:end)  -ut(2:end-1,2:end-1))/dx+(vt(3:end,2:end-1) -vt(2:end-1,2:end-1))/dy);
    rhs=zeros(Ny+2, Nx+2);
    rhs(2:end-1,2:end-1)=div/dt;
    p=solve_poisson(p, rhs, dx, dy, beta);

    % Velocity correction
    u(2:end-1,2:end-1)= ut(2:end-1,2:end-1)-dt *(p(2:end-1,3:end)-p(2:end-1,2:end-1))/dx;
    v(2:end-1,2:end-1) =vt(2:end-1,2:end-1)-dt* (p(3:end,2:end-1)-p(2:end-1,2:end-1))/dy;

    % Plot every 200 steps
    if mod(n,10) == 0
        uc=0.5 *(u(2:end-1,2:end-1) +u(2:end-1 ,3:end));
        vc=  0.5 * (v(2:end-1,2:end-1)+  v(3:end,2:end-1));
        omega=(v(2:end-1,3:end) -v(2:end-1,1:end- 2 ))/ (2*dx)-  (u(3:end,2:end-1) -u(1:end-2,2:end-1))/(2*dy);

        figure(1); clf;
        subplot(2,1,1);
        contourf(X, Y, omega, 100, 'LineColor', 'none'); colorbar();
        colormap(turbo); 
        caxis([-100 100]);
        hold on; fill(cx+R*cos(theta), cy+R*sin(theta), 'c');
        title(['Vorticity at step ', num2str(n)]);
        axis equal tight;

        subplot(2,1,2);
        quiver(X(1:4:end,1:4:end), Y(1:4:end,1:4:end), uc(1:4:end,1:4:end), vc(1:4:end,1:4:end), 3);
        hold on; fill(cx+R*cos(theta), cy+R*sin(theta), 'k');
        title('Velocity Field'); axis equal tight;

        time=dt * (1:steps);
        drawnow;
    end
end

function [u,v,p]=apply_bc(u,v,p,Uinf)
    u(:,1)=Uinf; u(:,end)=u(:,end-1);
    v(:,1)=0;     v(:,end)=0;
    u(1,:)=u(2,:); u(end,:)=u(end-1,:);
    v(1,:)=0;     v(end,:)=0;
    p(:,1)=p(:,2); p(:,end)=p(:,end-1);
    p(1,:)=p(2,:); p(end,:)=p(end-1,:);
end
function u=GS_diffuse(u, nu, dt, dx, dy)
    [Ny,Nx]=size(u);
    for iter=1:50
        u_old=u;
        for j=2:Ny-1
            for i=2:Nx-1
                u(j,i)=(u_old(j,i)+dt*nu*( (u(j+1,i)+u(j-1,i)-2*u(j,i))/dy^2+(u(j,i+1)+u(j,i-1)-2*u(j,i))/dx^2));
            end
        end
    end
end
function [uL, vL]=vel_interpol(u, v, Xb, Yb, dx, dy)
    Nb=length(Xb);
    uL=zeros(1,Nb); vL=zeros(1,Nb);
    for k=1:Nb
        i=floor(Xb(k)/dx)+2; j=floor(Yb(k)/dy)+2;
        uL(k)=u(j,i);
        vL(k)=v(j,i);
    end
end


function [fx, fy]=spreadf(FxL, FyL, Xb, Yb, Nx, Ny, dx, dy)
    fx=zeros(Ny+2, Nx+2); fy=zeros(Ny+2, Nx+2);
    for k=1:length(Xb)
        i=floor(Xb(k)/dx)+2; j=floor(Yb(k)/dy)+2;
        fx(j,i)=fx(j,i)+FxL(k);
        fy(j,i)=fy(j,i)+FyL(k);
    end
end


function p=solve_poisson(p, rhs, dx, dy, beta)
    [Ny,Nx]=size(p);
    for it=1:200
        p_old=p;
        for j=2:Ny-1
            for i=2:Nx-1
                p(j,i)=(1-beta)*p(j,i)+beta*0.25 * (p(j+1,i)+p(j-1,i)+p(j,i+1)+p(j,i-1)-dx^2*rhs(j,i));
            end
        end
        if max(abs(p(:)-p_old(:))) < 1e-6
            break;
        end
    end
end