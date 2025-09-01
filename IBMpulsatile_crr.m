clc; clear; close all;

% Parameters
Ny=48;
Ly=0.2;
dy=Ly/Ny;
visc=0.01;
rho=1;

U0=1;
f=5;
omega=2*pi*f;
T=2*pi/omega;
dt=T/2000;

tsteps=5500;
dpx_amp=dpx_calc(U0, omega, visc, rho, Ly);


y=linspace(-dy/2, Ly+dy/2, Ny+2)';  
yc=y(2:end-1);                      
u=zeros(Ny+2, 1);  % including ghost nodes
ut=zeros(Ny+2, 1);

y_ib=[0, Ly];  
x_ib=ones(size(y_ib));  % dummy x since it's 1D

% Plot setup
figure('Name','IBM Womersley Flow'); 
subplot(2,1,1);
hProf=plot(u(2:end-1), yc, 'b-', 'LineWidth', 2); hold on;
hAnal=plot(u(2:end-1), yc, 'r--', 'LineWidth', 2);
xlabel('u'); ylabel('y'); title('Velocity Profile'); grid on; legend('Numerical','Analytical','Location','south');

subplot(2,1,2);
hErr=plot(0,0); xlabel('Timestep'); ylabel('Relative L2 Error'); grid on;
err_l2_hist=zeros(tsteps,1);


u(2:end-1)=vel_womers(yc, 0, dpx_amp, omega, visc, rho, Ly);

% Time loop
for n=1:tsteps
    t=n * dt;
    dpdx=dpx_amp * sin(omega * t);
    u_ib_desired=[0; 0];
    u_ib=interpolation(u, y_ib, y, dy);
    f_ib_old=zeros(size(u_ib));
    alpha=0.3;

    f_ib_raw=(u_ib_desired-u_ib)/dt;
    f_ib=alpha * f_ib_old+(1-alpha) * f_ib_raw;
    f_ib_old=f_ib;
    F=spreadf(f_ib, y_ib, y, dy);

    for j=2:Ny+1
        diffu_n=visc * (u(j+1)-2*u(j)+u(j-1))/dy^2;
        diffu_np1=visc * (ut(j+1)-2*ut(j)+ut(j-1))/dy^2;
        f_avg=dpdx/rho+F(j);
        ut(j)=u(j)+dt * (0.5 * (diffu_n+diffu_np1)+f_avg);
    end
    u=ut;

    %analytical vel
    u_anal=vel_womers(yc, t, dpx_amp, omega, visc, rho, Ly);

    % Error
    err=u(2:end-1)-u_anal;
    err_l2_hist(n)=norm(err)/norm(u_anal);

    if mod(n,20)==0 || n==1
        set(hProf,'XData',u(2:end-1), 'YData', yc);
        set(hAnal,'XData',u_anal, 'YData', yc);
        set(hErr, 'XData', 1:n, 'YData', err_l2_hist(1:n));
        drawnow;
    end
end

fprintf('Final Relative L2 Error: %.2e\n', err_l2_hist(end));

% IBM helper funcs
function u_ib=interpolation(u, y_ib, y, dy)
    u_ib=zeros(size(y_ib));
    for k=1:length(y_ib)
        j=floor(y_ib(k)/dy)+1;
        wy=(y_ib(k)-y(j))/dy;
        u_ib(k)=(1-wy)*u(j)+wy*u(j+1);
    end
end

function F=spreadf(f_ib, y_ib, y, dy)
    F=zeros(size(y));
    for k=1:length(y_ib)
        j=floor(y_ib(k)/dy)+1;
        wy=(y_ib(k)-y(j))/dy;
        F(j)  =F(j)  +(1-wy) * f_ib(k);
        F(j+1)=F(j+1)+wy * f_ib(k);
    end
end

function dpx_amp=dpx_calc(U0, omega, visc, rho, Ly)
    i=1i;
    spatial_factor=abs(1-cosh(sqrt(i * omega/visc)*0)/cosh(sqrt(i * omega/visc)* Ly/2));
    dpx_amp=-U0 * omega * rho/spatial_factor;
end

function u=vel_womers(y, t, dpx_amp, omega, visc, rho, Ly)
    i=1i;
    y_shifted=y-Ly/2;
    denom=cosh(sqrt(i * omega/visc) * Ly/2);
    spatial_part=1-cosh(sqrt(i * omega/visc) * y_shifted)/denom;
    time_factor=exp(i * (omega * t-pi/2));
    prefactor=dpx_amp/(i * omega * rho);
    u_complex=prefactor * spatial_part * time_factor;
    u=real(u_complex);
end
