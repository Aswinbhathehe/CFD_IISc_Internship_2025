clc; clear;

Nx=256; Ny=64;
Lx=4; Ly=1;
dx=Lx/Nx; dy=Ly/Ny;
x=linspace(0, Lx-dx, Nx);
y=linspace(0, Ly-dy, Ny);
[X, Y]=meshgrid(x, y);
AR=Lx/Ly;

kappa=0.125;
a=kappa*Ly;
Reps=[1, 3];
Uw=0.1;
G=2*Uw/Ly;
rho=1;

theta=linspace(0, 2*pi, 100);
Nb=length(theta);
Xb0=a*cos(theta);  
Yb0=a*sin(theta);


ytilda_values=linspace(-0.35, 0.35, 41); % same as Fig. 3
lift_vs_y=zeros(length(Reps), length(ytilda_values));

for r=1:length(Reps)
    Rep=Reps(r);
    nu=G*a^2/Rep;
    dt=0.1;  % LBM-compatible value
    mu=rho*nu;

    for j=1:length(ytilda_values)
        y0=(ytilda_values(j)+0.5)*Ly;
        x0=2; % center in x
        ux=zeros(Ny, Nx); uy=zeros(Ny, Nx);
        fx=zeros(Ny, Nx); fy=zeros(Ny, Nx);
        

        Xb=x0+Xb0;
        Yb=y0+Yb0;

        for i=1:Ny
            uy(i,:)=0;
            ux(i,:)=Uw*(y(i)-0)/Ly;
        end


        for iter=1:1000

            Ub=interpolate_vel(ux, uy, Xb, Yb, x, y, dx, dy);
            Up=mean(Ub,1); 
            Omega=mean((Ub(:,1).*Yb0(:)-Ub(:,2).*Xb0(:))/a^2); 

            Up_local=[Up(1)+Omega*(-Yb0(:)), Up(2)+Omega*( Xb0(:))];

            alpha=1000;  % spring stiffness 
            Fb=alpha*(Up_local-Ub);

            fx(:)=0; fy(:)=0;
            [fx, fy]=spreadf(fx, fy, Fb, Xb, Yb, x, y, dx, dy);
            ux=ux+dt*fx/rho;
            uy=uy+dt*fy/rho;
        end

        lift_vs_y(r,j)=sum(Fb(:,2));
    end
end


Fl_dimless=lift_vs_y ./(rho*Uw^2*a*kappa^2); % equation taken from paper

figure;
hold on
colors=lines(length(Reps));
for i=1:length(Reps)
    plot(ytilda_values, Fl_dimless(i,:), 'LineWidth', 2, 'Color', colors(i,:));
end
xlabel('$\tilde{y}_0$', 'Interpreter', 'latex')
ylabel('$\tilde{F}_L$', 'Interpreter', 'latex')
legend(arrayfun(@(r) sprintf('Re_p=%d', r), Reps, 'UniformOutput', false))
title('Lift force vs transverse position (IBM)')
grid on
function Ub=interpolate_vel(ux, uy, Xb, Yb, x, y, dx, dy)
    [Ny, Nx]=size(ux);
    Nb=length(Xb); Ub=zeros(Nb,2);
    for k=1:Nb
        i0=floor(Xb(k)/dx); j0=floor(Yb(k)/dy);
        for ii=-1:2
            for jj=-1:2
                i=mod(i0+ii-1, Nx)+1;
                j=min(max(j0+jj,1), Ny); % non-periodic in y
                phi=delta((Xb(k)-x(i))/dx)*delta((Yb(k)-y(j))/dy);
                Ub(k,1)=Ub(k,1)+ux(j,i)*phi;
                Ub(k,2)=Ub(k,2)+uy(j,i)*phi;
            end
        end
    end
end


function [fx, fy]=spreadf(fx, fy, Fb, Xb, Yb, x, y, dx, dy)
    [Ny, Nx]=size(fx);
    Nb=length(Xb);
    for k=1:Nb
        i0=floor(Xb(k)/dx); j0=floor(Yb(k)/dy);
        for ii=-1:2
            for jj=-1:2
                i=mod(i0+ii,Nx)+1;
                j=mod(j0+jj,Ny)+1;
                phi=delta((Xb(k)-x(i))/dx)*delta((Yb(k)-y(j))/dy);
                fx(j,i)=fx(j,i)+Fb(k,1)*phi*dx*dy;
                fy(j,i)=fy(j,i)+Fb(k,2)*phi*dx*dy;
            end
        end
    end
end


function val=delta(r)
    r=abs(r);
    if r < 1
        val=0.125*(3-2*r+sqrt(1+4*r-4*r^2));
    elseif r < 2
        val=0.125*(5-2*r-sqrt(-7+12*r-4*r^2));
    else
        val=0;
    end
end
