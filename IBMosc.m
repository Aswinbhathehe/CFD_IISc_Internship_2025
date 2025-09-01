clc; clear;

x=32;
y=32;
Lx=1;
Ly=.2;
dx=Lx/x;
dy=Ly/y;

visc=.05;
U0=10;                % Amplitude
f=5 ;               % Frequency (Hz)
omega=2*pi*f;     % Angular frequency

Ub=0;                 % Bottom wall velocity (stationary)

CFL=0.5;
dt1=1e6;
dt2=CFL*min(dx, dy)/abs(U0);
dt=0.9*min(dt1, dt2)


[X, Y]=meshgrid(dx/2:dx:Lx-dx/2, dy/2:dy:Ly-dy/2);

y_ib=Ly*ones(1, x);
x_ib=linspace(dx/2, Lx-dx/2, x);

u =zeros(y+2, x+2);
v =zeros(y+2, x+2);
ut=zeros(y+2, x+2);
vt=zeros(y+2, x+2);
p =zeros(y+2, x+2);

fx=zeros(y+2, x+2);
fy=zeros(y+2, x+2);

% Plot setup
figure;
subplot(1,2,1);
prof_line=plot(zeros(y,1), linspace(dy/2,Ly-dy/2,y), 'b-', 'LineWidth', 2); hold on;
anal_line=plot(zeros(y,1), linspace(dy/2,Ly-dy/2,y), 'r--', 'LineWidth', 2);
legend('Numerical','Analytical','Location','south');
xlabel('u'); ylabel('y'); title('Velocity Profile at Lx/2'); grid on;

subplot(1,2,2);
err_plot=semilogy(0,0,'k');
xlabel('Timestep'); ylabel('Relative L2 Error'); title('Convergence'); grid on;
tsteps=2000;
err_l2_hist=zeros(tsteps,1);


% 
% %Video setup
% k=VideoWriter('IBM_OscShearFlow.mp4', 'MPEG-4');
% k.FrameRate=10;
% open(k);

%Time loop
for n=1:tsteps
    time=n*dt;
    Ut=U0*sin(omega*time);
    u_desired=Ut*ones(size(x_ib));
    u_ib=interpolate(u, x_ib, y_ib, dx, dy);
    f_ib=(u_desired-u_ib)/dt;
    [fx, fy]=spreadf(f_ib, zeros(size(f_ib)), x_ib, y_ib, size(u), dx, dy);

    ut=u+dt*(visc*laplacian(u, dx, dy)+fx);
    vt=v+dt*(visc*laplacian(v, dx, dy)+fy);

    divut=(ut(2:end-1,3:end)-ut(2:end-1,2:end-1))/dx+(vt(3:end,2:end-1)-vt(2:end-1,2:end-1))/dy;
    rhs=divut/dt;
    [p,~]=sor_solver(p, rhs, Lx, Ly, x, y);

    u(2:end-1,2:end-1)=ut(2:end-1,2:end-1)-dt*(p(2:end-1,3:end)-p(2:end-1,2:end-1))/dx;
    v(2:end-1,2:end-1)=vt(2:end-1,2:end-1)-dt*(p(3:end,2:end-1)-p(2:end-1,2:end-1))/dy;

    u(end,:)=Ut;    % Top wall velocity
    u(1,:)  =Ub;    % Bottom wall stationary

    % Error analysis
    yc=linspace(dy/2, Ly-dy/2, y);
    uc=0.5*(u(2:end-1,2:end-1)+u(2:end-1,3:end));
    u_profile=mean(uc,2);
    ua=vel_anal(yc, time, U0, omega, visc);

    err=sqrt(sum((u_profile-ua').^2)/length(ua))/sqrt(sum(ua.^2)/length(ua));
    err_l2_hist(n)=err;

    % Plot
    if mod(n,10) == 0
        set(prof_line, 'XData', u_profile, 'YData', yc);
        set(anal_line, 'XData', ua, 'YData', yc);
        set(err_plot, 'XData', 1:n, 'YData', err_l2_hist(1:n));
        drawnow;


        %         % Save frame to video
        % frame=getframe(gcf);
        % writeVideo(k, frame);
    end
end


% %Close Video
% close(k);
% disp('Video saved as IBM_OscShearFlow.mp4');


function u_analytical=vel_anal(y, t, U0, omega, visc)
        alpha=sqrt(omega/(2*visc));
        Ly=y(end);
        y_from_top=Ly-y;
        u_analytical=U0*exp(-alpha*y_from_top) .* sin(omega*t-alpha*y_from_top);
end


function L=laplacian(f, dx, dy)
    L=zeros(size(f));
    L(2:end-1,2:end-1)=(f(2:end-1,3:end)-2*f(2:end-1,2:end-1)+f(2:end-1,1:end-2))/dx^2+ (f(3:end,2:end-1)-2*f(2:end-1,2:end-1)+f(1:end-2,2:end-1))/dy^2;
end

function u_ib=interpolate(u, x_ib, y_ib, dx, dy)
    u_ib=zeros(size(x_ib));
    for k=1:length(x_ib)
        i=floor(x_ib(k)/dx)+1;
        j=floor(y_ib(k)/dy)+1;
        wx=(x_ib(k)-(i-1)*dx)/dx;
        wy=(y_ib(k)-(j-1)*dy)/dy;
        u_ib(k)=(1-wx)*(1-wy)*u(j,i)+wx*(1-wy)*u(j,i+1)+(1-wx)*wy*u(j+1,i)+wx*wy*u(j+1,i+1);
    end
end

function [fx, fy]=spreadf(fx_ib, fy_ib, x_ib, y_ib, size_u, dx, dy)
    fx=zeros(size_u);
    fy=zeros(size_u);
    for k=1:length(x_ib)
        i=floor(x_ib(k)/dx)+1;
        j=floor(y_ib(k)/dy)+1;
        wx=(x_ib(k)-(i-1)*dx)/dx;
        wy=(y_ib(k)-(j-1)*dy)/dy;
        fx(j,i) =fx(j,i) +(1-wx)*(1-wy)*fx_ib(k);
        fx(j,i+1)=fx(j,i+1) +wx*(1-wy)*fx_ib(k);
        fx(j+1,i) =fx(j+1,i) +(1-wx)*wy*fx_ib(k);
        fx(j+1,i+1)=fx(j+1,i+1) +wx*wy*fx_ib(k);

        fy(j,i)=fy(j,i) +(1-wx)*(1-wy)*fy_ib(k);
        fy(j,i+1)=fy(j,i+1)+wx*(1-wy)*fy_ib(k);
        fy(j+1,i)=fy(j+1,i)+(1-wx)*wy*fy_ib(k);
        fy(j+1,i+1)  =fy(j+1,i+1)  +wx*wy*fy_ib(k);
    end
end

function [p, err]=sor_solver(p, rhs, Lx, Ly, Nx, Ny)
    dx=Lx/Nx;
    dy=Ly/Ny;
    B=1.9;
    tol=1e-10;
    maxit=8000;
    err=1e10;
    it=0;

    while err > tol && it < maxit
        p_old=p;
        for i=2:Nx+1
            for j=2:Ny+1
                if i < Nx+1 && j < Ny+1
                    p(j,i)=(1-B)*p(j,i)+B*0.5*((dy^2*(p(j,i+1)+p(j,i-1))+ dx^2*(p(j+1,i)+p(j-1,i))-dx^2*dy^2*rhs(j,i))/(2*(dx^2+dy^2)));
                end
            end
        end
        err=norm(p(:)-p_old(:), 2);
        it=it+1;
    end
end
