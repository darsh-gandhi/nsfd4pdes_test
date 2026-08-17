clear; close all;

% performing consistency check on the nsfd scheme to solve the heat
% equation by decreasing the repeatedly decreasing the mesh size and
% solving the pde

%% heat equation set up

alpha=1; L=1; T=1;
u_ex = @(x,t) sin(pi*x/L).*exp(-alpha*(pi/L)^2*t);

%% check on nsfd scheme

n_checks = 5; % # of times we refine the mesh
subints = 8; % # of initial subintervals of the spatial domain
ratio = 0.5; % ratio of delta_t/(delta_x)^2

delta_x = -1.0*ones(n_checks,1);
errors = -1.0*ones(n_checks,1);
order = -1.0*zeros(n_checks,1);

for k=1:n_checks
    N = subints*2^(k-1);
    delta_t = ratio*(L/N)^2; delta_x(k)=L/N;

    [~,~,~,err,~] = nsfd_solver(N,delta_t,T,alpha,L,u_ex);
    errors(k) = err;
    if k>1
        order(k) = log(errors(k-1)/errors(k))/log(2);
    end
end

%% compare to standard scheme

errors_ftcs = -1.0*ones(n_checks,1);
order_ftcs = -1.0*zeros(n_checks,1);

for k=1:n_checks
    N = subints*2^(k-1);
    delta_t = ratio*(L/N)^2;

    [~,~,~,err,~] = ftcs_solver(N,delta_t,T,alpha,L,u_ex);
    errors_ftcs(k) = err;
    if k>1
        order_ftcs(k) = log(errors_ftcs(k-1)/errors_ftcs(k))/log(2);
    end
end

%% error plots

figure;
subplot(1,2,1)
loglog(delta_x, errors_ftcs, 'rs-', 'LineWidth', 1.5, 'MarkerFaceColor', 'r'); hold on;
loglog(delta_x, errors_ftcs(1)*(delta_x/delta_x(1)).^2, 'k--', 'LineWidth', 1.2);
xlabel('\Delta x'); ylabel('max error');
legend('FTCS error', 'O(\Delta x^2)', 'Location', 'best');
title('Standard scheme');
grid on;

% figure;
subplot(1,2,2)
loglog(delta_x, errors, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b'); hold on;
loglog(delta_x, errors(1)*(delta_x/delta_x(1)).^2, 'k--', 'LineWidth', 1.2);
xlabel('\Delta x'); ylabel('max error');
legend('NSFD error', 'O(\Delta x^2)', 'Location', 'best');
title('NSFD scheme');
grid on;

%%

function [x,u_nsfd,uex,err,lambda] = nsfd_solver(N,dt,T,alpha,L,u_ex)
    dx=L/N;
    Nt = ceil(T/dt); dt=T/Nt;
    x = (0:N)'*dx;
    mu = 4*alpha*dt/(dx^2);
    phi = ((dx^2)/(4*alpha)) * (1-exp(-mu));
    lambda = alpha*phi/(dx^2);
    u_nsfd = u_ex(x,0);
    for n=1:Nt
        u_old = u_nsfd;
        u_nsfd(2:end-1) = u_old(2:end-1) + lambda*(u_old(3:end) - 2*u_old(2:end-1) + u_old(1:end-2));
        u_nsfd(1) = u_ex(0,n*dt);
        u_nsfd(end) = u_ex(L,n*dt);
    end
    uex = u_ex(x,Nt*dt);
    err = max(abs(u_nsfd-uex));
end

function [x,u_ftcs,uex,err,lambda] = ftcs_solver(N,dt,T,alpha,L,u_ex)
    dx=L/N;
    Nt = ceil(T/dt); dt=T/Nt;
    x = (0:N)'*dx;
    lambda = alpha*dt/(dx^2);
    u_ftcs = u_ex(x,0);
    for n=1:Nt
        u_old = u_ftcs;
        u_ftcs(2:end-1) = u_old(2:end-1) + lambda*(u_old(3:end) - 2*u_old(2:end-1) + u_old(1:end-2));
        u_ftcs(1) = u_ex(0,n*dt);
        u_ftcs(end) = u_ex(L,n*dt);
    end
    uex = u_ex(x,Nt*dt);
    err = max(abs(u_ftcs-uex));
end