% ---------------------------------------------------%
%             Orbit-Trnasfer Problem                 %
% ---------------------------------------------------%
% Solve the following optimal control problem:       %
% Maximize t_f                                       %
% subject to the differential equation constraints   %
%   dr/dt       = v_r                                %
%   dtheta/dt   = v_theta/r                          %
%   dv_r/dt     = v_theta^2/r - mu/r^2 + T*u_1/m     %
%   dv_theta/dt = -v_r*v_theta/r + T*u_2/m           %
% the equality path constraint                       %
%   u_1^2 + u_2^2 = 1                                %
% and the boundary conditions                        %
%   r(0)         = 1                                 %
%   theta(0)     = 0                                 %
%   v_r(0)       = 0                                 %
%   v_theta(0)   = sqrt(mu/r(0))                     %
%   m(0)         = 1                                 %
%   r(t_f)       = 1.5                               %
%   v_r(t_f)     = 0                                 %
%   v_theta(t_f) = sqrt(mu/r(t_f))                   %

 close all; clear all;
% -------------------------------------------------- %
% BEGIN: DO NOT ALTER THE FOLLOWING LINES OF CODE!!! %
% -------------------------------------------------- %
global igrid CONSTANTS psStuff nstates ncontrols npaths
% -------------------------------------------------- %
% END:   DO NOT ALTER THE FOLLOWING LINES OF CODE!!! %
% -------------------------------------------------- %
path = 'C:\Users\elias\Documents\UF_Classes\EML6934\orthogonalCollocationAdigator\orthogonalCollocationAdigator';
addpath(genpath(path))

CONSTANTS.MU = 1;
CONSTANTS.m0 = 1;
% CONSTANTS.mdot = 0.0749;
% CONSTANTS.T = 0.1405;
CONSTANTS.ve = 1.8658344;
nstates      = 5;
ncontrols    = 3;
npaths       = 1;

% Bounds on State and Control
r0 = 1;   theta0 = 0; vr0 = 0; vtheta0 = 1; m0 = 1;
rf = 1.5;          vrf    = 0; vthetaf = sqrt(1/rf);  

rmin      = 0.5;  rmax      = 5;
thetamin  = 0;    thetamax  = 4*pi;
vrmin     = -10;  vrmax     = 10;
vthetamin = -10;  vthetamax = 10;
mmin      = 0.1;  mmax      = m0;
u1min     = -10;  u1max     = 10;
u2min     = -10;  u2max     = 10;
u3min     =  0;   u3max     = 0.1405;
t0min     = 0;    t0max     = 0;
tfmin     = 0;    tfmax     = 5;

% % Set polynomial degree in each mesh interval
% N = 50;
% 
% % Set locations of mesh points
% % meshPoints = linspace(-1,1,10).';
% meshPoints = [-1; 1];
% 
% % Allocate Vector of Polynomial Degrees
% polyDegrees = N*ones(length(meshPoints)-1,1);
% 
% % Compute LGR points, weights, and differentiation matrix.
% [tau,w,D] = lgrPS(meshPoints,polyDegrees);
% NLGR = length(w);
% psStuff.tau = tau;
% psStuff.w   = w;
% psStuff.D   = D;

% Set polynomial degree and number of intervals
N = 4;
numIntervals = 32;

% Create Leguandre Gauss Points
meshPoints  = linspace(-1,1,numIntervals+1).';  
polyDegrees = N*ones(numIntervals,1);
[tau,w,D]   = lgrPS(meshPoints,polyDegrees);
psStuff.tau = tau; psStuff.w = w; psStuff.D = D; NLGR = length(w);

% Set the bounds on the NLP variables.
zrmin = rmin*ones(length(tau),1);
zrmax = rmax*ones(length(tau),1);
zrmin(1)   = r0; zrmax(1)   = r0;
zrmin(end) = rf; zrmax(end) = rf;

zthetamin = thetamin*ones(length(tau),1);
zthetamax = thetamax*ones(length(tau),1);
zthetamin(1) = theta0; zthetamax(1) = theta0;

zvrmin = vrmin*ones(length(tau),1);
zvrmax = vrmax*ones(length(tau),1);
zvrmin(1) = vr0; zvrmax(1) = vr0;
zvrmin(end) = vrf; zvrmax(end) = vrf;

zvthetamin = vthetamin*ones(length(tau),1);
zvthetamax = vthetamax*ones(length(tau),1);
zvthetamin(1)   = vtheta0; zvthetamax(1)   = vtheta0;
zvthetamin(end) = vthetaf; zvthetamax(end) = vthetaf;

zmmin = mmin*ones(length(tau),1);
zmmax = mmax*ones(length(tau),1);
zmmin(1) = m0; zmmax(1) = m0;

zu1min = u1min*ones(length(tau)-1,1);
zu1max = u1max*ones(length(tau)-1,1);

zu2min = u2min*ones(length(tau)-1,1);
zu2max = u2max*ones(length(tau)-1,1);

zu3min = u3min*ones(length(tau)-1,1);
zu3max = u3max*ones(length(tau)-1,1);

zmin = [zrmin; zthetamin; zvrmin; zvthetamin; zmmin; zu1min; zu2min; zu3min; t0min; tfmin];
zmax = [zrmax; zthetamax; zvrmax; zvthetamax; zmmax; zu1max; zu2max; zu3max; t0max; tfmax];

% Set the bounds on the NLP constraints
% There are NSTATES sets of defect constraints.
defectMin = zeros(nstates*(length(tau)-1),1);
defectMax = zeros(nstates*(length(tau)-1),1);
% There is one path constraint
pathMin = ones(length(tau)-1,1); pathMax = ones(length(tau)-1,1);
% I dont believe there is nonlinear event constraint
eventMin = [];   eventMax = [];
objMin   = -inf; objMax   = inf;
Fmin = [objMin; defectMin; pathMin; eventMin];
Fmax = [objMax; defectMax; pathMax; eventMax];

% Supply an initial guess
rguess      = linspace(r0,rf,NLGR+1).';
thetaguess  = linspace(theta0,theta0,NLGR+1).';
vrguess     = linspace(vr0,vrf,NLGR+1).';
vthetaguess = linspace(vtheta0,vtheta0,NLGR+1).';
mguess      = linspace(m0,m0,NLGR+1).';
u1guess     = linspace(1,1,NLGR).';
u2guess     = linspace(0,0,NLGR).';
u3guess     = linspace(u3min,u3max,NLGR).';
t0guess     = 0;
tfguess     = 3.5;
z0 = [rguess;thetaguess;vrguess;vthetaguess;mguess;u1guess;u2guess;u3guess;t0guess;tfguess];

%-----------------------------------------------------------------%
% Generate derivatives and sparsity pattern using Adigator        %
%-----------------------------------------------------------------%
% - Constraint Function Derivatives
xsize  = size(z0);
x      = adigatorCreateDerivInput(xsize,'z0');
output = adigatorGenJacFile('orbitTransferFun',{x});
S_jac  = output.JacobianStructure;
[iGfun,jGvar] = find(S_jac);

% - Objective Function Derivatives
xsize  = size(z0);
x      = adigatorCreateDerivInput(xsize,'z0');
output = adigatorGenJacFile('orbitTransferObj',{x});
grd_structure = output.JacobianStructure;

%-----------------------------------------------------------------%
% set IPOPT callback functions
%-----------------------------------------------------------------%
funcs.objective   = @(Z)orbitTransferObj(Z);
funcs.gradient    = @(Z)orbitTransferGrd(Z);
funcs.constraints = @(Z)orbitTransferCon(Z);
funcs.jacobian    = @(Z)orbitTransferJac(Z);
funcs.jacobianstructure = @()orbitTransferJacPat(S_jac);
options.ipopt.hessian_approximation = 'limited-memory';

%-----------------------------------------------------------------%
% Set IPOPT Options %
%-----------------------------------------------------------------%
options.ipopt.tol = 1e-8;
options.ipopt.linear_solver = 'ma57'; %'mumps';
options.ipopt.max_iter = 2000;
options.ipopt.mu_strategy = 'adaptive';
options.ipopt.ma57_automatic_scaling = 'yes';
options.ipopt.print_user_options = 'yes';
options.ipopt.output_file = ['orbitTransfer','IPOPTinfo.txt']; % print output file
options.ipopt.print_level = 5; % set print level default

options.lb = zmin; % Lower bound on the variables.
options.ub = zmax; % Upper bound on the variables.
options.cl = Fmin; % Lower bounds on the constraint functions.
options.cu = Fmax; % Upper bounds on the constraint functions.

%-----------------------------------------------------------------%
% Call IPOPT
%-----------------------------------------------------------------%
[z, info] = ipopt(z0,funcs,options);

%-----------------------------------------------------------------%
% extract lagrange multipliers from ipopt output, info
%-----------------------------------------------------------------%
Fmul = info.lambda;

% Extract the state and control from the decision vector z.
% Remember that the state is approximated at the LGR points
% plus the final point, while the control is only approximated 
% at only the LGR points.
r      = z(1:NLGR+1);
theta  = z(NLGR+2:2*(NLGR+1));
vr     = z(2*(NLGR+1)+1:3*(NLGR+1));
vtheta = z(3*(NLGR+1)+1:4*(NLGR+1));
m      = z(4*(NLGR+1)+1:5*(NLGR+1));
u1     = z(5*(NLGR+1)+1:5*(NLGR+1)+NLGR);
u2     = z(5*(NLGR+1)+NLGR+1:5*(NLGR+1)+2*NLGR);
u3     = z(5*(NLGR+1)+2*NLGR+1:5*(NLGR+1)+3*NLGR);
alpha  = mod(atan2(u1,u2),2*pi)*180/pi;
t0     = z(end-1);
tf     = z(end);
t      = (tf-t0)*(tau+1)/2+t0;
tLGR   = t(1:end-1);
%-----------------------------------------------------------------%
% Extract the Lagrange multipliers corresponding                  %
% the defect constraints.                                         %
%-----------------------------------------------------------------%
multipliersDefects = Fmul(2:nstates*NLGR+1);
multipliersDefects = reshape(multipliersDefects,NLGR,nstates);
%-----------------------------------------------------------------%
% Compute the costates at the LGR points via transformation       %
%-----------------------------------------------------------------%
costateLGR = inv(diag(w))*multipliersDefects;
%-----------------------------------------------------------------%
% Compute the costate at the tau=+1 via transformation            %
%-----------------------------------------------------------------%
costateF = D(:,end).'*multipliersDefects;
%-----------------------------------------------------------------%
% Now assemble the costates into a single matrix                  %
%-----------------------------------------------------------------%
costate = [costateLGR; costateF];    
lamr = costate(:,1); lamtheta = costate(:,2);
lamvr = costate(:,3); lamvtheta = costate(:,4);

%-------------%
% Plot Results 
%-------------%
figure(1);
subplot(2,2,1);
plot(t,r,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$r(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

subplot(2,2,2);
plot(t,theta,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\theta(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

subplot(2,2,3);
plot(t,vr,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$v_r(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

subplot(2,2,4);
plot(t,vtheta,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$v_\theta(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

figure(2);
plot(tLGR,alpha,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\alpha(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

figure(3);
subplot(1,2,1);
plot(tLGR,u1,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u_1(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

subplot(1,2,2);
plot(tLGR,u2,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u_2(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

figure(4);
subplot(2,2,1);
plot(t,lamr,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\lambda_r(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

subplot(2,2,2);
plot(t,lamtheta,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\lambda_\theta(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

subplot(2,2,3);
plot(t,lamvr,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\lambda_{v_r}(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

subplot(2,2,4);
plot(t,lamvtheta,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\lambda_{v_\theta}(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

figure; hold on; grid minor
plot(t,r);
plot(t,theta);
plot(t,vr);
plot(t,vtheta);
plot(t,m);
xlabel('t','Interpreter','LaTeX')
legend('r(t)','$\theta(t)$','vr(t)','v$\theta(t)$','m(t)','Interpreter','LaTeX')
set(gcf,'color','white')
set(gca,'fontweight','bold','fontsize',10)
title('States for trajectory that minimized fuel consumption')
