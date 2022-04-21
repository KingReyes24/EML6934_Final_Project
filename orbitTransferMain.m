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
format long
%  clear all
% -------------------------------------------------- %
% BEGIN: DO NOT ALTER THE FOLLOWING LINES OF CODE!!! %
% -------------------------------------------------- %
global igrid CONSTANTS psStuff nstates ncontrols npaths path_constraint maximize_mass
% -------------------------------------------------- %
% END:   DO NOT ALTER THE FOLLOWING LINES OF CODE!!! %
% -------------------------------------------------- %
%% GENERATE PATH THAT CONTAINS IPOPT and ADIGATOR
path = 'C:\Users\elias\Documents\UF_Classes\EML6934\Final_Project\EML6934_Final_Project';
addpath(genpath(path))
addpath('C:\Users\elias\Downloads\optiMEXFiles_mexw64_2_28')
%% Save figures and latex table. 
% Do not change this if debugging
save_figs          = 0;
create_latex_table = 0;

%% Set path constrant and objective function descision
path_constraint= 0; % is there a path contraint?
maximize_mass  = 0; % maximize m(tf), else min tf

%% Set polynomial degrees and intervals to loop through
n_list = [3 ];        % polynomial degrees
k_list = [4 8 ]; % intervals

%% Set Globals
% set gloabl constants
CONSTANTS.MU = 1;
CONSTANTS.m0 = 1;
CONSTANTS.ve = 1.8658344;
% set number of states 
nstates = 5;

% set number of controls and paths
if path_constraint
    ncontrols = 3; 
    npaths    = 1;
else
    ncontrols = 2;
    npaths    = 0;
end

%% Solve Optimal Control Problem
% counter for table
count  = 1;

for nIdx = 1:numel(n_list)
    
    for kIdx = 1:numel(k_list)
        
    % Set polynomial degree and number of intervals
    N              = n_list(nIdx); % number of polynomial degree
    numIntervals   = k_list(kIdx); % number of intervals

    % Bounds on State and Control
    % thetaf and mf are free
    r0 = 1;   theta0 = 0; vr0 = 0; vtheta0 = 1; m0 = 1;
    rf = 1.5;             vrf = 0; vthetaf = sqrt(1/rf);  

    rmin      = 1;    rmax      = 1.5;
    thetamin  = 0;    thetamax  = 4*pi;
    vrmin     = -10;  vrmax     = 10;
    vthetamin = -10;  vthetamax = 10;
    mmin      = 0.1;  mmax      = m0;
    t0min     = 0;    t0max     = 0;
    tfmin     = 0;    tfmax     = 25;
    u3min     = 0;    u3max     = 0.1405; % thrust
    if path_constraint
        % set bounds for control with path contraint
        u1min = -1.5;  u1max = 1.5;     % sin(beta)
        u2min = -1.5;  u2max = 1.5;     % cos(beta)
    else
        % set bound for beta 
        u1min = -4*pi; u1max = 4*pi;  % beta
        u2min = 0;     u2max = 0;     % empty
    end

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
    zvrmin(1)   = vr0; zvrmax(1)   = vr0;
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
    
    % Contruct vector of bounds. 
    if path_constraint
        % using u1, u2, u3
        zmin = [zrmin; zthetamin; zvrmin; zvthetamin; zmmin; zu1min; zu2min; zu3min; t0min; tfmin];
        zmax = [zrmax; zthetamax; zvrmax; zvthetamax; zmmax; zu1max; zu2max; zu3max; t0max; tfmax];
    else
        % using u1, u3
        zmin = [zrmin; zthetamin; zvrmin; zvthetamin; zmmin; zu1min; zu3min; t0min; tfmin];
        zmax = [zrmax; zthetamax; zvrmax; zvthetamax; zmmax; zu1max; zu3max; t0max; tfmax];
    end
    
    % Set the bounds on the NLP constraints
    % There are NSTATES sets of defect constraints.
    defectMin = zeros(nstates*(length(tau)-1),1);
    defectMax = zeros(nstates*(length(tau)-1),1);
    if path_constraint 
        % There is a path constraint
        pathMin = ones(length(tau)-1,1); pathMax = ones(length(tau)-1,1);
    else
        % No path constraint
        pathMin = []; pathMax = [];
    end
    % I dont believe there is nonlinear event constraint
    eventMin = [];   eventMax = [];
    objMin   = -inf; objMax   = inf;
    Fmin = [objMin; defectMin; pathMin; eventMin];
    Fmax = [objMax; defectMax; pathMax; eventMax];

    % Supply an initial guess
    rguess      = linspace(r0,rf,NLGR+1).';
    thetaguess  = linspace(theta0,2.5,NLGR+1).';
    vrguess     = linspace(vr0,vrf,NLGR+1).';
    vthetaguess = linspace(vtheta0,vthetaf,NLGR+1).';
    mguess      = linspace(mmax,mmin,NLGR+1).';
    u1guess     = linspace(u1min,u1max,NLGR).';
    u3guess     = linspace(u3max,u3max,NLGR).';
    t0guess     = 0;
    tfguess     = 3;

    if path_constraint
        % contruct u1guess and u2guess for path contraint
        u1guess = linspace(1,1,NLGR).';
        u2guess = linspace(0,0,NLGR).';
        z0 = [rguess;thetaguess;vrguess;vthetaguess;mguess;u1guess;u2guess;u3guess;t0guess;tfguess];
    else  
        z0 = [rguess;thetaguess;vrguess;vthetaguess;mguess;u1guess;u3guess;t0guess;tfguess];
    end

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
    options.ipopt.tol           = 1e-8;
    options.ipopt.linear_solver = 'ma57';%'ma57'; %'mumps';
    options.ipopt.max_iter      = 8000;
    options.ipopt.mu_strategy   = 'adaptive';
    options.ipopt.ma57_automatic_scaling = 'yes';
    options.ipopt.print_user_options = 'no';
    options.ipopt.output_file        = ['orbitTransfer','IPOPTinfo.txt']; % print output file
    options.ipopt.print_level        = 5; % set print level default

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
    if path_constraint
        u2   = z(5*(NLGR+1)+NLGR+1:5*(NLGR+1)+2*NLGR);
        u3   = z(5*(NLGR+1)+2*NLGR+1:5*(NLGR+1)+3*NLGR);
        beta = atan2(u1,u2);
        beta = unwrap(beta)*180/pi;
    else
        u3   = z(5*(NLGR+1)+NLGR+1:5*(NLGR+1)+2*NLGR);
        beta = unwrap(u1)*180/pi;
%         beta = u1;
    end
    t0     = z(end-1);
    tf     = z(end);
    t      = (tf-t0)*(tau+1)/2+t0;
    tLGR   = t(1:end-1);
    %-----------------------------------------------------------------%
    % Extract the Lagrange multipliers corresponding                  %
    % the defect constraints.                                         %
    %-----------------------------------------------------------------%
%     multipliersDefects = Fmul(2:nstates*NLGR+1);
%     multipliersDefects = reshape(multipliersDefects,NLGR,nstates);
%     %-----------------------------------------------------------------%
%     % Compute the costates at the LGR points via transformation       %
%     %-----------------------------------------------------------------%
%     costateLGR = inv(diag(w))*multipliersDefects;
%     %-----------------------------------------------------------------%
%     % Compute the costate at the tau=+1 via transformation            %
%     %-----------------------------------------------------------------%
%     costateF = D(:,end).'*multipliersDefects;
%     %-----------------------------------------------------------------%
%     % Now assemble the costates into a single matrix                  %
%     %-----------------------------------------------------------------%
%     costate = [costateLGR; costateF];    
%     lamr = costate(:,1); lamtheta = costate(:,2);
%     lamvr = costate(:,3); lamvtheta = costate(:,4);

    %-----------------------------------------------------------------%
    % Get planer coordinates
    %-----------------------------------------------------------------%
    x_transfer  = r.*cos(theta);
    y_transfer  = r.*sin(theta);
    theta_orbit = linspace(0,2*pi,100);
    x_orbit_1   = r0*cos(theta_orbit);
    y_orbit_1   = r0*sin(theta_orbit);
    x_orbit_2   = rf*cos(theta_orbit);
    y_orbit_2   = rf*sin(theta_orbit);

    %-------------%
    % Plot Results 
    %-------------%
%     close all
    
    fig_cntrl = figure;
    h1 = subplot(2,1,1);
    plot(tLGR,beta,'-o');
    ylabel('$\beta(t)$','Interpreter','LaTeX');
    xlabel('$t$','Interpreter','LaTeX')
    set(gca,'FontName','Times','FontSize',14);
    set(gcf,'color','white')
    title('Thrust angle that minimized objective function')
    grid minor;
    
    h2 = subplot(2,1,2); 
    plot(tLGR,u3,'-o');
    set(gca,'FontName','Times','FontSize',14);
    set(gcf,'color','white')
    ylabel('$T(t)$','Interpreter','LaTeX');
    title('Thrust magnitude that minimized objective function')
    xlabel('$t$','Interpreter','LaTeX')
    grid minor
    linkaxes([h1 h2],'x')
    
    fig_states = figure; hold on; grid minor
    plot(t,r,'-o'); plot(t,theta,'-o'); plot(t,vr,'-o');
    plot(t,vtheta,'-o'); plot(t,m,'-o');
    xlabel('t','Interpreter','LaTeX')
    legend('r(t)','$\theta(t)$','vr(t)','v$\theta(t)$','m(t)','Interpreter','LaTeX')
    set(gcf,'color','white')
    set(gca,'FontName','Times','fontsize',14)
    title('States that minimized objective function')

    fig_orbit = figure; hold on; grid minor
    plot(x_orbit_1,y_orbit_1)
    plot(x_orbit_2,y_orbit_2)
    plot(x_transfer,y_transfer)
    set(gca,'FontName','Times','fontsize',14)
    set(gcf,'color','white')
    legend('Initial Orbit','Final Orbit','Orbit Transfer')
    axis equal
    title('Orbit Transfer Overview')
    
    if save_figs
        if maximize_mass 
            str_obj = 'mf';
        else
            str_obj = 'tf';
        end
        plot_str = sprintf('_N%d_K%d_C%d_',N,numIntervals,ncontrols);
        str_list = {'control','states','orbit'};
        for idx = 1:numel(str_list)
            name = [str_list{idx} plot_str str_obj];
            print(figure(idx),name,'-depsc')
        end
    end
    if path_constraint 
        fig_path = figure;
        subplot(1,2,1);
        plot(tLGR,u1,'-o');
        xl = xlabel('$t$','Interpreter','LaTeX');
        yl = ylabel('$u_1(t)$','Interpreter','LaTeX');
        set(xl,'FontSize',14);
        set(yl,'FontSize',14);
        set(gca,'FontName','Times','FontSize',14);
        grid on;

        subplot(1,2,2);
        plot(tLGR,u2,'-o');
        xl = xlabel('$t$','Interpreter','LaTeX');
        yl = ylabel('$u_2(t)$','Interpreter','LaTeX');
        set(xl,'FontSize',14);
        set(yl,'FontSize',14);
        set(gca,'FontName','Times','FontSize',14);
        set(gcf,'color','white')
        grid on;
        
        if save_figs
            str_path = ['path' plot_str str_obj];
            print(fig_path,str_path,'-depsc')
        end
    end
    
%     figure(1);
%     subplot(2,2,1);
%     plot(t,r,'-o');
%     xl = xlabel('$t$','Interpreter','LaTeX');
%     yl = ylabel('$r(t)$','Interpreter','LaTeX');
%     set(xl,'FontSize',14);
%     set(yl,'FontSize',14);
%     set(gca,'FontName','Times','FontSize',14);
%     set(gcf,'color','white')
%     grid on;
% 
%     subplot(2,2,2);
%     plot(t,theta,'-o');
%     xl = xlabel('$t$','Interpreter','LaTeX');
%     yl = ylabel('$\theta(t)$','Interpreter','LaTeX');
%     set(xl,'FontSize',14);
%     set(yl,'FontSize',14);
%     set(gca,'FontName','Times','FontSize',14);
%     set(gcf,'color','white')
%     grid on;
% 
%     subplot(2,2,3);
%     plot(t,vr,'-o');
%     xl = xlabel('$t$','Interpreter','LaTeX');
%     yl = ylabel('$v_r(t)$','Interpreter','LaTeX');
%     set(xl,'FontSize',14);
%     set(yl,'FontSize',14);
%     set(gca,'FontName','Times','FontSize',14);
%     set(gcf,'color','white')
%     grid on;
% 
%     subplot(2,2,4);
%     plot(t,vtheta,'-o');
%     xl = xlabel('$t$','Interpreter','LaTeX');
%     yl = ylabel('$v_\theta(t)$','Interpreter','LaTeX');
%     set(xl,'FontSize',14);
%     set(yl,'FontSize',14);
%     set(gca,'FontName','Times','FontSize',14);
%     set(gcf,'color','white')
%     grid on;

%     figure;
%     subplot(2,2,1);
%     plot(t,lamr,'-o');
%     xl = xlabel('$t$','Interpreter','LaTeX');
%     yl = ylabel('$\lambda_r(t)$','Interpreter','LaTeX');
%     set(xl,'FontSize',14);
%     set(yl,'FontSize',14);
%     set(gca,'FontName','Times','FontSize',14);
%     set(gcf,'color','white')
%     grid on;
% 
%     subplot(2,2,2);
%     plot(t,lamtheta,'-o');
%     xl = xlabel('$t$','Interpreter','LaTeX');
%     yl = ylabel('$\lambda_\theta(t)$','Interpreter','LaTeX');
%     set(xl,'FontSize',14);
%     set(yl,'FontSize',14);
%     set(gca,'FontName','Times','FontSize',14);
%     set(gcf,'color','white')
%     grid on;
% 
%     subplot(2,2,3);
%     plot(t,lamvr,'-o');
%     xl = xlabel('$t$','Interpreter','LaTeX');
%     yl = ylabel('$\lambda_{v_r}(t)$','Interpreter','LaTeX');
%     set(xl,'FontSize',14);
%     set(yl,'FontSize',14);
%     set(gca,'FontName','Times','FontSize',14);
%     set(gcf,'color','white')
%     grid on;
% 
%     subplot(2,2,4);
%     plot(t,lamvtheta,'-o');
%     xl = xlabel('$t$','Interpreter','LaTeX');
%     yl = ylabel('$\lambda_{v_\theta}(t)$','Interpreter','LaTeX');
%     set(xl,'FontSize',14);
%     set(yl,'FontSize',14);
%     set(gca,'FontName','Times','FontSize',14);
%     set(gcf,'color','white')
%     grid on;
    
    %--------------%
    % save results %
    %--------------%
    degree(count,1)      = N;
    intervals(count,1)   = numIntervals;
    iterations(count,1)  = info.iter;
    cpu_time(count,1)    = info.cpu;
    final_time(count,1)  = tf;
    final_mass(count,1)  = m(end);
    solved_info(count,1) = info.status;
    count = count + 1;
    end
end
table_mat = horzcat(degree,intervals,iterations,cpu_time,final_time,final_mass,solved_info);
T =  array2table(table_mat);
VarNames = {'Degree','Intervals','Iterations','CPU Time','tf','mf','Solved_Status'};
T.Properties.VariableNames = VarNames;
if create_latex_table
    str_table = sprintf('table_C%d_',ncontrols);
    table2latex(T,[str_table str_obj])
end

arrayfun(@(x) set( x,'windowstyle', 'docked'),findall(groot,'type','figure')) ;