% ---------------------------------------------------%
%             Orbit-Trnasfer Problem                 %
% ---------------------------------------------------%
close all; clear all;
format long
% -------------------------------------------------- %
% BEGIN: DO NOT ALTER THE FOLLOWING LINES OF CODE!!! %
% -------------------------------------------------- %
global igrid CONSTANTS psStuff nstates ncontrols npaths path_constraint maximize_mass
% -------------------------------------------------- %
% END:   DO NOT ALTER THE FOLLOWING LINES OF CODE!!! %
% -------------------------------------------------- %

%% Save figures and latex table. 
% Set to 0 when debugging
save_figs          = 0;
create_latex_table = 0;
plot_comparisons   = 1;

if save_figs
    set(0,'DefaultFigureWindowStyle','normal')
end
%% Set path constrant and objective function descision
path_constraint= 0; % is there a path contraint?
maximize_mass  = 1; % maximize m(tf), else min tf

%% Set polynomial degrees and intervals to loop through
n_list = [3 4];         % polynomial degrees
k_list = [2 4 8 16 32]; % intervals

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
    if maximize_mass
        rmax  = 1.6;
        mmin  = 0.7;
        tfmax = 100;
    end
    if path_constraint
        % set bounds for control with path contraint
%         u1min = -1.5;  u1max = 1.5;     % sin(beta)
%         u2min = -1.5;  u2max = 1.5;     % cos(beta)
        u1min = -1.5;  u1max = 1.5;     % sin(beta)
        u2min = -1.5;  u2max = 1.5;     % cos(beta)
    else
        % set bound for beta 
        u1min = -2*pi; u1max = 2*pi;  % beta
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
    zu3min = u3min*ones(length(tau)-1,1);
    zu3max = u3max*ones(length(tau)-1,1);
    
    % Contruct vector of bounds. 
    if path_constraint
        % using u1, u2, u3
        zu2min = u2min*ones(length(tau)-1,1);
        zu2max = u2max*ones(length(tau)-1,1);

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

    objMin = -inf; 
    objMax = inf;
    Fmin   = [objMin; defectMin; pathMin]; 
    Fmax   = [objMax; defectMax; pathMax];

    % Supply an initial guess
    rguess      = linspace(r0,rf,NLGR+1).';
    thetaguess  = linspace(theta0,2.5,NLGR+1).';
    vrguess     = linspace(vr0,vr0,NLGR+1).';
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
    options.ipopt.max_iter      = 4000;
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
    end
    t0     = z(end-1);
    tf     = z(end);
    t      = (tf-t0)*(tau+1)/2+t0;
    tLGR   = t(1:end-1);

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
    close all
    if maximize_mass 
        str_obj  = 'mf';
        str_obj1 = 'maximized ';
    else
        str_obj  = 'tf';
        str_obj1 = 'minimized ';
    end    

    fig_cntrl = figure;
    h1 = subplot(2,1,1);
    plot(tLGR,beta,'-o');
    ylabel('$\beta(t)$','Interpreter','LaTeX');
    xlabel('$t$','Interpreter','LaTeX')
    set(gca,'FontName','Times','FontSize',14);
    set(gcf,'color','white')
    title(['Thrust angle that ' str_obj1 str_obj])
    grid minor;
    
    h2 = subplot(2,1,2); 
    plot(tLGR,u3,'-o');
    set(gca,'FontName','Times','FontSize',14);
    set(gcf,'color','white')
    ylabel('$T(t)$','Interpreter','LaTeX');
    title(['Thrust magnitude that ' str_obj1 str_obj])
    xlabel('$t$','Interpreter','LaTeX')
    grid minor
    linkaxes([h1 h2],'x')
    
    fig_states = figure; hold on; grid minor
    plot(t,r,'-o'); plot(t,theta,'-o'); plot(t,vr,'-o');
    plot(t,vtheta,'-o'); plot(t,m,'-o');
    xlabel('t','Interpreter','LaTeX')
    legend('r(t)','$\theta(t)$','$v_r(t)$','$v_\theta(t)$','m(t)','Interpreter','LaTeX')
    set(gcf,'color','white')
    set(gca,'FontName','Times','fontsize',14)
    title(['States of trajectory that ' str_obj1 str_obj])

    fig_orbit = figure;   
    polarplot(theta_orbit,1*ones(1,numel(theta_orbit)))
    hold on;
    polarplot(theta_orbit,1.5*ones(1,numel(theta_orbit)))
    polarplot(theta,r); 
    set(gca,'FontName','Times','fontsize',14)
    set(gcf,'color','white')
    legend('Initial Orbit','Final Orbit','Orbit Transfer')
    title('Orbit transfer overview')
    
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
        suptitle(['Path constriant control that ' str_obj1 str_obj])
        
        if save_figs
            str_path = ['path' plot_str str_obj];
            print(fig_path,str_path,'-depsc')
        end
    end
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
nlistvec = repmat(n_list',1,numel(k_list));
klistvec = repmat(k_list,numel(n_list),1);
runtime  = reshape(cpu_time,[],2);
runtime  = runtime';
if maximize_mass
    mfmat  = reshape(final_mass,[],2);
    objmat = mfmat';
   
else
    tfmat  = reshape(final_time,[],2);
    objmat = tfmat';
end
if plot_comparisons
    fig_obj = figure; hold on;
    plot(k_list,objmat(1,:),'-o')
    plot(k_list,objmat(2,:),'-o')
    xlabel('Intervals')
    ylabel(['Objective function ' str_obj])
    set(gcf,'color','white')
    legend('Degree 3','Degree 4')
    set(gca,'fontweight','bold','fontsize',14, 'XMinorGrid','on','YMinorGrid','on')
    title('Objective Function Performance')

    fig_runtime = figure; hold on;
    plot(k_list,runtime(1,:),'-o')
    plot(k_list,runtime(2,:),'-o')
    legend('Degree 3','Degree 4')
    set(gcf,'color','white')
    set(gca,'fontweight','bold','fontsize',14, 'XMinorGrid','on','YMinorGrid','on')
    xlabel('Intervals')
    ylabel('Execution Time (s)')
    title('Execution Time Performance')

    if save_figs
        str_control = sprintf('_c%d_',ncontrols);
        str_path = ['runtime' str_control str_obj];
        print(fig_runtime,str_path,'-depsc')

        str_path = ['obj' str_control str_obj];
        print(fig_obj,str_path,'-depsc')

    end
end

% set(0,'DefaultFigureWindowStyle','docked')
% arrayfun(@(x) set( x,'windowstyle', 'docked'),findall(groot,'type','figure')) ;