clc; clear; close all;

y=48;
Ly=0.2;
dy=Ly/y;

visc=0.01;
rho=1;

% Womersleyflow parameters
U0=1;
f=5;
omega=2*pi*f;
T=2*pi/omega;
dt=T/2000;% 2000 points per cycle

dpdx_amp=dpx_calc(U0, omega, visc, rho, Ly);  % analytical dp/dx amplitude
tsteps=5400;

yc=linspace(dy/2, Ly - dy/2, y)'; 
u=womer_vel(yc, 0, dpdx_amp, omega, visc, rho, Ly);  % t=0

% Initial acceleration consistency to reduce startup error
du_dt0=womersley_acceleration(yc, 0, dpdx_amp, omega, visc, rho, Ly);
u=u+0.5*dt*du_dt0;
ut=u;

err_l2_hist=zeros(tsteps,1);

% Plot Setup
figure('Name','Pulsating/Womersley Flow'); 
subplot(2,1,1);
hProfile=plot(u, yc, 'bo-', 'DisplayName', 'Numerical'); hold on;
hExact=plot(u, yc, 'r-', 'DisplayName', 'Analytical');
xlabel('u (m/s)'); ylabel('y'); grid on;
legend; title('Velocity Profile');

subplot(2,1,2);
hError=plot(0,0); xlabel('Time Step'); ylabel('Relative L2 Error'); grid on;
title(sprintf('L2 Error vs Time Step (dt=%.2e s)', dt));

% Time-stepping Loop
old_dpx=dpdx_amp;
for n=1:tsteps
    t=n*dt;
    dpdx=dpdx_amp*sin(omega*t);
    smooth_dpx=0.5*(dpdx+old_dpx);
    old_dpx=dpdx;

    % Crank–Nicolson time integration (semi-implicit)
    for j=2:y-1
        diffu_n=visc*(u(j+1) - 2*u(j)+u(j-1))/dy^2;
        frocin_n=old_dpx/rho;
         diffu_np1=visc*(ut(j+1) - 2*ut(j)+ut(j-1))/dy^2;
        frocin_np1=dpdx/rho;
        ut(j)=u(j)+dt/2*(diffu_n+frocin_n+diffu_np1+frocin_np1);
    end

    % No-slip BCs
    ut(1)=0;
    ut(end)=0;
    u=ut;


    u_exact=womer_vel(yc, t, dpdx_amp, omega, visc, rho, Ly);
    u_exact(1)=0; u_exact(end)=0;


    err=u - u_exact;
    err_l2=norm(err,2);
    norm_ref=norm(u_exact,2)+1e-10;
    err_l2_hist(n)=err_l2/norm_ref;

    % Plot every few steps
    if mod(n, 20) == 0 || n == 1
        set(hProfile, 'XData', u, 'YData', yc);
        set(hExact, 'XData', u_exact, 'YData', yc);
        set(hError, 'XData', 1:n, 'YData', err_l2_hist(1:n));
        drawnow;
    end
end


fprintf('Final Relative L2 Error: %.2e\n', err_l2_hist(end));

%helper functions 
function dpdx_amp=dpx_calc(U0, omega, visc, rho, Ly)
    i=1i;
    lambda=sqrt(i*omega/visc);
    spatial_factor=abs(1 - cosh(lambda*0)/cosh(lambda*Ly/2));
    dpdx_amp=-U0*omega*rho/spatial_factor;
end

function u=womer_vel(y, t, dpdx_amp, omega, visc, rho, Ly)
    i=1i;
    lambda=sqrt(i*omega/visc);
    y_shifted=y - Ly/2;
    denom=cosh(lambda*Ly/2);
    spatial_part=1 - cosh(lambda*y_shifted)/denom;
    time_factor=exp(i*(omega*t - pi/2));
    prefactor=dpdx_amp/(i*omega*rho);
    u_complex=prefactor*spatial_part*time_factor;
    u=real(u_complex);
end

function du_dt=womersley_acceleration(y, t, dpdx_amp, omega, visc, rho, Ly)
    i=1i;
    lambda=sqrt(i*omega/visc);
    y_shifted=y - Ly/2;
    denom=cosh(lambda*Ly/2);
    spatial_part=1 - cosh(lambda*y_shifted)/denom;
    time_factor=omega*exp(i*(omega*t+pi/2));
    prefactor=dpdx_amp/rho;
    du_dt_complex=prefactor*spatial_part*time_factor;
    du_dt=real(du_dt_complex);
end
