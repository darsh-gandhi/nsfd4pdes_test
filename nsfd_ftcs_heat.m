clear; close all

%% problem

%compute soln to homogeneous heat eq using standard fd and nsfd scheme
%developed today; compare results when time-step and space-step do not
%follow the cfl cond (but a functional rel. of the two does)

%PDE: u_t - Ku_{xx} = 0, x\in\Omega=[0,1], u(t,0)=u(t,1)=0, u(0,x) = x(1-x);

%% initialize baseline params + functions
colors = 1/255 * [56 182 255; 255 222 89; 0 191 99]; % colors for plotting: blue, ylw, green

Omega = linspace(0,1,11); % spatial discretization
Nx = length(Omega);         % # of spatial elements
Delta_x = Omega(2)-Omega(1); % spatial step size

u0 = 4*Omega.*(1-Omega); % ic
u = u0';                % soln (row) vector (defined by ic for now)

K = 1; % thermal diffusivity
T = 1; % final time = 1 units

Delta_t = 0.9*(Delta_x)^2/(2*K); % time-step (def in accordance with CFL)
Nt = round(T/Delta_t);  % # of timesteps (large for sanity check)

function [CFL, alpha, A] = def_rhs_mat(delta_t,delta_x,K,Nx)
    CFL = delta_t/(delta_x)^2 <= 1/(2*K); % logical check for CFL cond
    
    alpha = K*delta_t/(delta_x)^2;
    A = diag(alpha*ones(1,Nx-3),-1) + diag((1-2*alpha)*ones(1,Nx-2)) + diag(alpha*ones(1,Nx-3),1); %sparse matrix rep of linear system
end

%% standard FTCS implementation (with CFL met)

[CFL, ~, A] = def_rhs_mat(Delta_t,Delta_x,K,Nx);

if(CFL)
    fprintf('CFL condition met.\n')
else
    fprintf('CFL condition not met.\n')
end

start_time = cputime;
figure(1)
for i=1:Nt
    u(2:Nx-1) = A*u(2:Nx-1);
    % fprintf('max u = %.6f\n', max(u)); failsafe

    plot(Omega,u,'LineWidth',1.5,'Color',colors(2,:))
    hold on
    plot(Omega,u0,'LineWidth',1.5,'Color',colors(1,:))
    xlabel('x')
    ylabel('u')
    ylim([0,1])
    legend('soln','ic')
    title('Standard FTCS with Satisfied CFL')
    set(gca,'FontSize',14)
    pause(0.001)
    hold off
end
comp_time_baseFTCS = cputime-start_time

%% FTCS implementation with CFL broken

u_2=u0';
delta_t = 1.5*(Delta_x)^2/(2*K);
Nt = round(T/delta_t);
[cfl, ~, A_2] = def_rhs_mat(delta_t,Delta_x,K,Nx);

if(cfl)
    fprintf('CFL condition met.\n')
else
    fprintf('CFL condition not met.\n')
end

start_time = cputime;
figure(2)
for i=1:Nt
    u_2(2:Nx-1) = A_2*u_2(2:Nx-1);
    % fprintf('max u = %.6f\n', max(u)); failsafe

    plot(Omega,u_2,'LineWidth',1.5,'Color',colors(2,:))
    hold on
    plot(Omega,u0,'LineWidth',1.5,'Color',colors(1,:))
    xlabel('x')
    ylabel('u')
    % ylim([0,1])
    legend('soln','ic')
    title('Standard FTCS with Broken CFL')
    set(gca,'FontSize',14)
    pause(0.001)
    hold off
end
comp_time_brokenFTCS = cputime-start_time

%% NSFD implementation (using params to break CFL)

%NSFD function: (1-e^(-2Kr))/2K, r := delta_t/(delta_x)^2

function phi = phi(delta_t,delta_x,K)
    phi = (1-exp(-2*K*delta_t/(delta_x)^2))/(2*K);
end

function A = def_rhs_nsfd(delta_t,delta_x,K,Nx)
    alpha = K*phi(delta_t,delta_x,K);
    A = diag(alpha*ones(1,Nx-3),-1) + diag((1-2*alpha)*ones(1,Nx-2)) + diag(alpha*ones(1,Nx-3),1); %sparse matrix rep of linear system
end

delta_t = 5*(Delta_x)^2/(2*K);
Nt = round(T/delta_t);

u_nsfd = u0';
A_nsfd = def_rhs_nsfd(delta_t,Delta_x,K,Nx);

start_time = cputime;
figure(3)
for i=1:Nt
    u_nsfd(2:Nx-1) = A_nsfd*u_nsfd(2:Nx-1);
    % fprintf('max u = %.6f\n', max(u)); failsafe

    plot(Omega,u_nsfd,'LineWidth',1.5,'Color',colors(2,:))
    hold on
    plot(Omega,u0,'LineWidth',1.5,'Color',colors(1,:))
    xlabel('x')
    ylabel('u')
    % ylim([0,1])
    legend('soln','ic')
    title('Nonstandard FTCS for Broken CFL')
    set(gca,'FontSize',14)
    pause(0.001)
    hold off
end
comp_time_nsfd = cputime-start_time