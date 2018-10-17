classdef LC_in_square < handle
    
    % Created by YW, 2018-07-26
    % Find Linear Shape and Timing Responses of a Limit Cycle with Sliding Components 
    % under instantaneous/static perturbation
    %
    % Minimal examples to run (with default parameter and initial values):
    %   M = LC_in_square;
    %   M.solve;
    %   M.plot;
    %
    % ===========================
    % Detailed instructions:
    % ===========================
    % A class of LC_in_square enables us to create a model object and
    % change its parameter values interactively.
    %
    % 1, First step of running this solver is always to call its constructor:
    %
    %   M = LC_in_square(varargin)
    %
    % with the following optional Name-Value pairs
    %       varon       -- 'false' by default, setting to 'true' the model will solve the variational problem
    %       xinit,vinit -- the initial values of x and v (variational problem)
    %       tmax        -- the endtime of the ode integration
    %       alpha       -- model parameter alpha with default value 0.2
    %       omega       -- model parameter omega with default value 1
    %       nu          -- relative change in frequency for regions below and above the wedge defined by function inWedge
    %       eps         -- perturbation added to parameters: (alpha, omega) -> (alpha+eps, omega-eps) 
    % You can enter any number of inputs to the constructor since all of them have pre-defined default values.
    % However, inputs determine which problem will be solved
    %
    %     a) To only find solution trajectory of the model, run
    %
    %                     M = LC_in_square('xinit', xinit) 
    %
    %     b) To solve variational problem w.r.t. instantaneous perturbation, run
    %
    %                     M = LC_in_square('varOn', true, 'xinit', xinit, 'vinit', vinit)
    %
    %     c) For a uniform static perturbation: (alpha, omega) -> (alpha + eps, omega - eps),
    %        construct the model with a **SCALAR** nu as input
    %
    %          i) To find unperturbed solution and the shape response curve, run
    %
    %                     M = LC_in_square('varOn', true, 'xinit', xinit, 'vinit', vinit, 'nu', nu)
    
    %             where nu needs to be computed using iPRC
    %
    %          ii) To find perturbed solution, run 
    %
    %                     M = LC_in_square('xinit', xinit_perturb, 'alpha', alpha + eps, 'omega', omega + eps) 
    %              or
    %                     M = LC_in_square('xinit', xinit_perturb, 'eps', eps)
    %
    %              where nu takes the default value 0
    %
    %     d) For piecewise static perturbations, where the perturbation (alpha, omega) -> (alpha + eps, omega - eps)
    %        is only present in the region above the wedge, 
    %        construct the model with a **VECTOR** nu with the first element nu_below and the second element nu_above
    %
    %          i) To find unperturbed solution and the corresponding iSRC, run
    %
    %          M = LC_in_square('varOn', true, 'xinit', xinit, 'vinit', vinit, 'nu', [nu_below, nu_above])
    %
    %
    %        ii) To find perturbed solution, run
    %
    %          M = LC_in_square('xinit', xinit, 'vinit', vinit, 'nu', [0, 0], 'eps', eps)
    %
    %     e) To find iPRC,
    %                    M = LC_in_square('xinit', xinit);
    %                    M.solve;
    %                    M.find_prc(z0);   % z0 is the initial condition for iPRC found by running find_prc_monodromy
    %                    M.plot_prc
    %
    %  2, After constructing a model, typing M in command window and entering will display
    %     all current properties set. In addition to the above eight, there are also:
    %
    %           yinit   -- concatenation of xinit and vinit if varOn = true
    %          domain   -- current domain (0 interior, 1-4 walls)
    %               t   -- time array
    %            yext   -- full solution array, each new solution will be appended to new row (yext = [yext; ynew])
    %             prc   -- phase response curve
    %              t0   -- current time (scalar)
    %              y0   -- current solution (same size as yinit)
    %            Salt   -- array of Saltation matrix (v+ = Sv-)
    %         t_event   -- keeping the times of events (wall hitting)
    %    domain_event   -- domains where trajectory enters
    %          t_exit   -- keeping the times of leaving walls
    %     domain_exit   -- domains where trajectory leaves
    %       Jump_exit   -- array of Jump matrix at exit (z- = Jz+)
    %     isPiecewise   -- 0 if nu is typed in as a scalar, 1 if nu is typed in as a vector
    %
    %   It is not recommended to manually change any of these values because
    %   they are all counted as model "outputs"
    %   You can check their values any time using the dot reference, such as M.domain
    %
    % 3, Type "methods(M)" will show the available methods in the model:
    %
    %           solve   -- requires no input, solve the model with given initial values
    %            plot   -- requires no input, plot the solutions (alway solve before plot)
    %        find_prc   -- requires the initial value for iPRC/lTRC as an input,
    %                      solve the adjoint equation backward in time with given input (always solve the model before find_prc)
    %        plot_prc   -- requires no input, plot the iPRC/lTRC solutions (always find_prc before plot_prc)
    %      findPeriod   -- requires running the model with tmax big enough
    %          LC_ODE   -- ODE of the model, used internally
    %
    % ===========================
    % Some other usage examples:
    % ===========================
    % ** Estimate solution period
    %
    %  >> M=LC_in_square('xinit',[1,0.2],'vinit',[1 0],'tmax',20);
    %  >> M.solve
    %  >> M.findPeriod
    %   ans =
    %
    %        6.7662
        
    properties(Constant)
        reltol = 1e-13; % ode15s tolerance
        abstol = 1e-13; % ode15s tolerance
    end
    
    properties(SetAccess = protected)
        xinit
        vinit
        yinit
    end
    
    properties
        nu_below
        nu_above
        varOn
        tmax
        alpha
        omega
        eps
        domain
        t = [];
        reverseTspan=[]; % Full unique time after being reversed (remove duplicate)
        yext = [];
        prct=[]; % time for prc
        prc=[]; % phase response solution
        
        t_exit = []; % Record times exiting a wall or crossing the open/close bdry
        Jump_exit = {}; % Record inverse jump matrices associated with exit
    end
    
    properties(Access = protected)
        S0
        t0
        y0
    end
    
    properties(SetAccess = private)
        isPiecewise
    end
    
    methods
        
        function model = LC_in_square(varargin)    
            
            p = inputParser;
            addOptional(p, 'varOn', false, @(x)validateattributes(x,...
                {'logical'},{'nonempty'}));
            addOptional(p, 'xinit', [0.9, -0.9], @(x)validateattributes(x,...
                {'numeric'},{'nonempty'}));
            addOptional(p, 'vinit', [0.01, 0], @(x)validateattributes(x,...
                {'numeric'},{'nonempty'}));
            addOptional(p, 'tmax', 6.766182958128617, @(x)validateattributes(x,...
                {'numeric'},{'nonempty'}));
            addOptional(p, 'alpha', 0.2, @(x)validateattributes(x,...
                {'numeric'},{'nonempty'}));
            addOptional(p, 'omega', 1, @(x)validateattributes(x,...
                {'numeric'},{'nonempty'}));
            addOptional(p, 'nu', 0, @(x)validateattributes(x,...
                {'numeric'},{'nonempty'}));
            addOptional(p, 'eps', 0, @(x)validateattributes(x,...
                {'numeric'},{'nonempty'}));
            
            parse(p, varargin{:})
            
            nu = p.Results.nu;
            
            % If the input nu is a vector, then the perturbation is piecewise uniform; 
            % if nu is a scalar, then the perturbation is uniform
            if length(nu) > 1
                model.isPiecewise = true;
                model.nu_below = nu(1);
                model.nu_above = nu(2);
            else
                model.isPiecewise = false;
                model.nu_below = nu;
                model.nu_above = nu;
            end
            
            model.eps = p.Results.eps;            
            model.alpha = p.Results.alpha;
            model.omega = p.Results.omega;
            model.tmax = p.Results.tmax;
            model.vinit = p.Results.vinit;
            model.xinit = p.Results.xinit;
            model.varOn = p.Results.varOn;            
            model.xinit = reshape(model.xinit,1,length(model.xinit));
            model.vinit = reshape(model.vinit,1,length(model.vinit));
            
            model.domain = 0;
            if model.varOn
                model.yinit = [model.xinit, model.vinit];
                model.S0 = eye(2);
            else
                model.yinit = model.xinit;
            end
            model.t0 = 0;
            model.y0 = model.yinit;
        end
        
        function solve(model)
            % Initialize
            model.t0 = 0;
            model.y0 = model.yinit;
            model.domain = 0;
            model.t = []; % Full time
            model.yext = []; % Full solution
            
            model.t_exit = [];
            model.Jump_exit = {};
            
            options0=odeset('BDF','on','RelTol',model.reltol,'AbsTol',model.abstol,'Events',@model.dom0_to_wall);
            options1=odeset('BDF','on','RelTol',model.reltol,'AbsTol',model.abstol,'Events',@model.wall1_exit);
            options2=odeset('BDF','on','RelTol',model.reltol,'AbsTol',model.abstol,'Events',@model.wall2_exit);
            options3=odeset('BDF','on','RelTol',model.reltol,'AbsTol',model.abstol,'Events',@model.wall3_exit);
            options4=odeset('BDF','on','RelTol',model.reltol,'AbsTol',model.abstol,'Events',@model.wall4_exit);
            
            options1_ext=odeset('BDF','on','RelTol',model.reltol,'AbsTol',model.abstol,'Events',@model.wall1_exit_ext);
            options2_ext=odeset('BDF','on','RelTol',model.reltol,'AbsTol',model.abstol,'Events',@model.wall2_exit_ext);
            options3_ext=odeset('BDF','on','RelTol',model.reltol,'AbsTol',model.abstol,'Events',@model.wall3_exit_ext);
            options4_ext=odeset('BDF','on','RelTol',model.reltol,'AbsTol',model.abstol,'Events',@model.wall4_exit_ext);
            while model.t0 < model.tmax
                switch model.domain
                    case 0 % interior
                        if model.varOn
                            model_ode = @model.LC_ODE_ext;
                            model_opt = options0;
                        else
                            model_ode = @model.LC_ODE;
                            model_opt = options0;
                        end
                        [tnew,ynew,TE,YE,IE] = ode45(model_ode,[model.t0,model.tmax],model.y0,model_opt); % integrate forwards in time until a wall is encountered
                        model.updateSolution(tnew,ynew);
                        model.updateCurrent(tnew,ynew,TE,YE,IE);
                        
                        model.domain = IE; % new domain is the encountered wall
                    case 1 % x=1 wall
                        if model.varOn
                            model.multiplySaltation(TE,YE,'enter');
                            model_ode = @model.LC_ODE_ext;
                            model_opt = options1_ext;
                        else
                            model_ode = @model.LC_ODE;
                            model_opt = options1;
                        end
                        [tnew,ynew,TE,YE,~] = ode45(model_ode,[model.t0,model.tmax],model.y0,model_opt); % integrate forwards in time until the wall is exited
                        model.updateSolution(tnew,ynew);
                        model.updateCurrent(tnew,ynew,TE,YE,0);
                        model.storeJump(TE); % save exit time and jump matrix needed for later finding PRC
                        model.domain=0; % new domain is the interior
                    case 2 % y=1 wall
                        if model.varOn
                            model.multiplySaltation(TE,YE,'enter');
                            model_ode = @model.LC_ODE_ext;
                            model_opt = options2_ext;
                        else
                            model_ode = @model.LC_ODE;
                            model_opt = options2;
                        end
                        [tnew,ynew,TE,YE,~] = ode45(model_ode,[model.t0,model.tmax],model.y0,model_opt); % integrate forwards in time until the wall is exited
                        model.updateSolution(tnew,ynew);
                        model.updateCurrent(tnew,ynew,TE,YE,0);
                        model.storeJump(TE); % save exit time and jump matrix needed for later finding PRC
                        model.domain=0; % new domain is the interior
                    case 3 % x=-1 wall
                        if model.varOn
                            model.multiplySaltation(TE,YE,'enter');
                            model_ode = @model.LC_ODE_ext;
                            model_opt = options3_ext;
                        else
                            model_ode = @model.LC_ODE;
                            model_opt = options3;
                        end
                        [tnew,ynew,TE,YE,~] = ode45(model_ode,[model.t0,model.tmax],model.y0,model_opt); % integrate forwards in time until the wall is exited
                        model.updateSolution(tnew,ynew);
                        model.updateCurrent(tnew,ynew,TE,YE,0);
                        model.storeJump(TE); % save exit time and jump matrix needed for later finding PRC
                        model.domain=0; % new domain is the interior
                    case 4 % y=-1 wall
                        if model.varOn
                            model.multiplySaltation(TE,YE,'enter');
                            model_ode = @model.LC_ODE_ext;
                            model_opt = options4_ext;
                        else
                            model_ode = @model.LC_ODE;
                            model_opt = options4;
                        end
                        [tnew,ynew,TE,YE,~] = ode45(model_ode,[model.t0,model.tmax],model.y0,model_opt); % integrate forwards in time until the wall is exited
                        model.updateSolution(tnew,ynew);
                        model.updateCurrent(tnew,ynew,TE,YE,0);
                        model.storeJump(TE); % save exit time and jump matrix needed for later finding PRC
                        model.domain=0; % new domain is the interior
                end
            end
        end
        
        function plot_prc(model)
            figure
            subplot(2,1,1)
            plot(model.t,model.yext(:,1:2),'linewidth',2)
            legend('x','y')
            xlim([0 model.tmax])
            set(gca,'FontSize',18)
            xlabel('Time','interpreter','latex','fontsize',25)
            
            subplot(2,1,2)
            plot(model.prct,model.prc(:,1:2),'linewidth',2)
            legend('Z_x','Z_y')
            xlim([0 model.tmax])
            xlabel('Time')
            legend('x-direction','y-direction')
            title('Phase response curve')
            grid on
            
            set(gca,'FontSize',18)
            xlabel('Time','interpreter','latex','fontsize',25)
            ylabel('iPRC','interpreter','latex','fontsize',25)
        end
        
        function plot(model)
            figure
            plot(model.yext(:,1),model.yext(:,2),'k','linewidth',2)
            axis([-1.1 1.1 -1.1 1.1])
            
            if model.varOn
                figure
                plot(model.t,model.yext(:,3),'linewidth',2)
                hold on
                plot(model.t,model.yext(:,4),'linewidth',2)
                xlabel('Time')
                ylabel('v1,v2')
                legend('v1','v2')
                grid on
            end
        end
        
        function T = findPeriod(model)
            y_temp = model.yext(:,2);
            t_temp = model.t;
            k = []; % indicies of peaks
            for i = 1:length(t_temp)
                if i ~= 1 && (y_temp(i) == 1) && (y_temp(i-1) ~=1)
                    % log the index when the signal reaches 1.
                    k(end+1) = i;
                end
            end
            
            if length(k) < 2
                warning(['Unable to calculate the period. '...
                    'Only found %d peak in the solution, consider using a larger tmax!'], length(k));
            end
            
            T = mean(diff(t_temp(k)));
        end
        
        function find_prc(model, z0)
            
            if nargin < 2
                z0 = [1 0];
            end
            
            model.prct = [];
            model.prc = [];
            options_prc=odeset('BDF','on','RelTol',model.reltol,...
                'AbsTol',model.abstol,'Events',@model.exit_wall); % integration will stop when a wall is entered backwards in time
            
            [model.reverseTspan, Ind] = unique(wrev(model.t),'stable');
            if isempty(model.reverseTspan)
                error('Solve the model first before calling find_prc!');
            end
            xmat = model.yext(wrev(Ind),1:2);
            
            dom = 0;
            counter = 0;
            TE = inf;
            
            while true
                switch dom
                    case 0 % interior or on along a wall, except at exit points
                        T = model.reverseTspan(model.reverseTspan <= TE);
                        if T == 0
                            break;
                        end
                        
                        % integrate backwards in time until a wall exit point is encountered
                        [tnew,znew,TE,YE,IE]=ode15s(@model.LC_ODE_prc,T,z0,options_prc,xmat);
                        model.prct = [model.prct; tnew];
                        model.prc = [model.prc; znew];
                        dom=1;
                        
                        if ~isempty(IE)
                            IE = IE(end);
                            TE = TE(end);
                        end
                        if counter >= numel(model.t_exit)
                            break;
                        end
                        
                    case 1 % at wall exit points
                        
                        % apply the jump matrix, eq. 3.25
                        J=model.Jump_exit{IE};
                        z0 = znew(end,1:2)*J';
                        counter = counter + 1;
                        dom=0;
                end
            end
        end
    end
    
    methods(Hidden)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% ODE for Non-Variational Problem %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dydt = LC_ODE(model,~,y,domain_overload)
            if nargin > 3
                ODEdomain = domain_overload;
            else
                ODEdomain = model.domain;
            end
            
            [a, b] = model.alphaFcn(y(1),y(2));
            dydt=[a,-b;b,a]*y;
            
            switch ODEdomain
                case 1
                    dydt(1)=min(dydt(1),0); % only allow negative dx/dt or else zero
                case 2
                    dydt(2)=min(dydt(2),0); % only allow negative dy/dt or else zero
                case 3
                    dydt(1)=max(dydt(1),0); % only allow positive dx/dt or else zero
                case 4
                    dydt(2)=max(dydt(2),0); % only allow positive dy/dt or else zero
            end
        end
        
        function [value,isterminal,direction]=dom0_to_wall(~,~,y)
            value=[...
                y(1)-1;...  % when x crosses 1 from below (wall 1)
                y(2)-1;...  % when y crosses 1 from below (wall 2)
                y(1)+1;...  % when x crosses -1 from above (wall 3)
                y(2)+1];    % when y crosses -1 from below (wall 4)
            isterminal=[1;1;1;1]; % stop integration and return
            direction=[1;1;-1;-1]; % "value" should be increasing
        end
        
        function [value,isterminal,direction]=wall1_exit(model,~,y)
            % when the *unconstrained* value of dx/dt decreases through zero, return
            y(1)=1;
            [a, b] = model.alphaFcn(y(1),y(2));
            dydt=[a,-b;b,a]*y;
            value=dydt(1);
            isterminal=1;
            direction=-1;
            
        end
        
        function [value,isterminal,direction]=wall2_exit(model,~,y)
            % when the *unconstrained* value of dy/dt decreases through zero, return
            y(2)=1;
            [a, b] = model.alphaFcn(y(1),y(2));
            dydt=[a,-b;b,a]*y;
            value=dydt(2);
            isterminal=1;
            direction=-1;
            
        end
        
        function [value,isterminal,direction]=wall3_exit(model,~,y)
            % when the *unconstrained* value of dx/dt increases through zero, return
            y(1)=-1;
            [a, b] = model.alphaFcn(y(1),y(2));
            dydt=[a,-b;b,a]*y;
            value=dydt(1);
            isterminal=1;
            direction=1;
            
        end
        
        function [value,isterminal,direction]=wall4_exit(model,~,y)
            % when the *unconstrained* value of dy/dt increases through zero, return
            y(2)=-1;
            [a, b] = model.alphaFcn(y(1),y(2));
            dydt=[a,-b;b,a]*y;
            value=dydt(2);
            isterminal=1;
            direction=1;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% ODE for Variational Problem %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dydt=LC_ODE_ext(model,~,y,domain_overload)
            if nargin > 3
                ODEdomain = domain_overload;
            else
                ODEdomain = model.domain;
            end
            x=y(1:2);
            v=y(3:4);
            [a, b] = model.alphaFcn(x(1),x(2));
            dxdt=[a,-b;b,a]*x;
            switch ODEdomain
                case 0 % interior
                    dvdt=[a,-b;b,a]*v;
                case 1 % x=1 wall
                    dxdt(1)=min(dxdt(1),0); % only allow negative dx/dt or else zero
                    dvdt=[0,0;0,a]*v;
                    % dvdt(2)=dvdt(2)+y(2);
                case 2 % y=1 wall
                    dxdt(2)=min(dxdt(2),0); % only allow negative dy/dt or else zero
                    dvdt=[a,0;0,0]*v;
                    % dvdt(1)=dvdt(1)+y(1);
                case 3 % x=-1 wall
                    dxdt(1)=max(dxdt(1),0); % only allow positive dx/dt or else zero
                    dvdt=[0,0;0,a]*v;
                    % dvdt(2)=dvdt(2)+y(2);
                case 4 % y=-1 wall
                    dxdt(2)=max(dxdt(2),0); % only allow positive dy/dt or else zero
                    dvdt=[a,0;0,0]*v;
                    % dvdt(1)=dvdt(1)+y(1);
            end
            % add nonhomogeneous terms to variational problem for sustained perturbation
            % nu*F(x(t)) + DFeps/Deps; DF/Deps=[x+y; -x+y]
            dvdt = dvdt + model.addon(ODEdomain,x(1),x(2),dxdt);
            dydt=[dxdt;dvdt];
        end
        
        function dzdt=LC_ODE_prc(model,t,z,xmat)
            
            xvec = interp1(model.reverseTspan,xmat,t);
            
            [a, b] = model.alphaFcn(xvec(1),xvec(2));
            switch model.checkdomain(xvec)
                case 0 % interior
                    DF=[a,-b;b,a]; % Jacobian for interior, from eq. 5.46
                case 1 % x=1 wall
                    DF=[0,0;0,a];  % Jacobian for sliding region, from eqs. 5.46 and 4.34  [JPG: should it be [0,0;b,a]  ?]
                case 2 % y=1 wall
                    DF=[a,0;0,0];  % Jacobian for sliding region, from eqs. 5.46 and 4.34  [JPG: should it be [a,-b;0,0] ?]
                case 3 % x=-1 wall
                    DF=[0,0;0,a];  % Jacobian for sliding region, from eqs. 5.46 and 4.34  [JPG: should it be [0,0;b,a]  ?]
                case 4 % y=-1 wall
                    DF=[a,0;0,0];  % Jacobian for sliding region, from eqs. 5.46 and 4.34  [JPG: should it be [a,-b;0,0] ?]
            end
            dzdt=-DF'*[z(1); z(2)]; % adjoint equation, eq. 2.7
        end
        
        function [value,isterminal,direction]=wall1_exit_ext(model,~,y)
            % when the *unconstrained* value of dx/dt decreases through zero, return
            y(1)=1;
            [a, b] = model.alphaFcn(y(1),y(2));
            value=b*y(2)-a;
            isterminal=1;
            direction=1;
            
        end
        
        function [value,isterminal,direction]=wall2_exit_ext(model,~,y)
            % when the *unconstrained* value of dy/dt decreases through zero, return
            y(2)=1;
            [a, b] = model.alphaFcn(y(1),y(2));
            value=b*y(1)+a;
            isterminal=1;
            direction=-1;
            
        end
        
        function [value,isterminal,direction]=wall3_exit_ext(model,~,y)
            % when the *unconstrained* value of dx/dt increases through zero, return
            y(1)=-1;
            [a, b] = model.alphaFcn(y(1),y(2));
            value=b*y(2)+a;
            isterminal=1;
            direction=-1;
            
        end
        
        function [value,isterminal,direction]=wall4_exit_ext(model,~,y)
            % when the *unconstrained* value of dy/dt increases through zero, return
            y(2)=-1;
            [a, b] = model.alphaFcn(y(1),y(2));
            value=b*y(1)-a;
            isterminal=1;
            direction=1;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Exit function for PRC Problem %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [value,isterminal,direction]=exit_wall(model,t,~,~)
            value = model.t_exit - t;
            isterminal=ones(size(model.t_exit));
            direction=ones(size(model.t_exit));
        end
    end
    
    methods(Hidden)
        function updateSolution(model,tnew,ynew)
            model.t = [model.t; tnew];
            model.yext = [model.yext; ynew];
        end
        
        function updateCurrent(model,tnew,ynew,TE,YE,IE)
            if ~isempty(TE)
                model.t0 = TE(end);
                model.y0 = YE(end,:);
                model.removeError(IE(end));
            else
                model.t0 = max(tnew);
                model.y0 = ynew(end,:);
            end
        end
        
        function removeError(model, IE)
            switch IE
                case 1
                    model.y0(1)=1;
                case 2
                    model.y0(2)=1;
                case 3
                    model.y0(1)=-1;
                case 4
                    model.y0(2)=-1;
            end
        end
        
        function multiplySaltation(model,TE,YE,flag) % multiply fundMatrix by saltation matrix S [JPG: the exit flag is never used. is this because we already know the saltation matrix at liftoff is identity? should it be used anyway for validation?]
            na=[1,0]; % na is the normal vector to wall 1 and 3
            nb=[0,1]; % nb is the normal vector to wall 2 and 4
            X=YE(:,1:2)';
            if ~isempty(TE)
                dydtdom = model.LC_ODE(TE,X); % vector field of the wall
                dydt0 = model.LC_ODE(TE,X,0); % vector field of the interior
                if strcmp(flag,'enter')
                    F_plus  = dydtdom; % new vector field after entering wall
                    F_minus = dydt0;   % old vector field before entering wall
                elseif strcmp(flag,'exit')
                    F_plus  = dydt0;   % new vector field after exiting wall
                    F_minus = dydtdom; % old vector field before exiting wall
                else
                    error('flag for multiplySaltation must be "enter" or "exit"!');
                end
                switch model.domain
                    case 1
                        n = na;
                    case 2
                        n = nb;
                    case 3
                        n = na;
                    case 4
                        n = nb;
                end
                S1=eye(2)+(F_plus-F_minus)*n/(n*F_minus); % saltation matrix, eq. 3.24
                model.y0(3:4) = [YE(3),YE(4)]*S1';
            end
            
        end
        
        
        function storeJump(model, TE)
            if isempty(TE)
                return;
            end
            
            model.t_exit = [model.t_exit, TE];
            switch model.domain
                case 1
                    J0 = jumpx();
                case 2
                    J0 = jumpy();
                case 3
                    J0 = jumpx();
                case 4
                    J0 = jumpy();
            end
            model.Jump_exit{end+1}=J0;
        end
        
        function h = addon(model, ODEdomain, x, y, dxdt)
            
            % if both nu's are zero, the nonhomogeneous addon is zero
            if (model.nu_above == 0) && (model.nu_below == 0) 
                h = 0;
            end
            
            % if it is a uniformly perturbed problem and nu's are not 0
            if ~model.isPiecewise && (model.nu_below ~= 0)
                h = model.nu_below*dxdt;
                if ODEdomain == 0
                    h(1) = h(1) + x + y;
                    h(2) = h(2) + y - x;
                elseif ODEdomain == 1 || ODEdomain == 3
                    h(2) = h(2) + y - x;
                elseif ODEdomain == 2 || ODEdomain == 4
                    h(1) = h(1) + x + y;
                end
            end
            
            % if it is a piecewisely perturbed problem when the perturbation only exists in the region above the wedge,
            if model.isPiecewise && (model.nu_above ~=0) && (model.nu_below ~=0)
                if ~inWedge(x,y) % dF/deps=0 in region below wedge
                    h = model.nu_below*dxdt;
                else % in region above the wedge where the perturbation exists
                    h = model.nu_above*dxdt;                    
                    if ODEdomain == 0
                        h(1) = h(1) + x + y;
                        h(2) = h(2) + y - x;
                    elseif ODEdomain == 1 || ODEdomain == 3
                        h(2) = h(2) + y - x;
                    elseif ODEdomain == 2 || ODEdomain == 4
                        h(1) = h(1) + x + y;
                    end
                end
            end
        end
        
        
        function [a, b] = alphaFcn(model, x, y)
            % if the perturbation is piecewise
            if model.isPiecewise
                if inWedge(x, y)  % perturbation is present in region above the wedge
                    a = model.alpha + model.eps; 
                    b = model.omega - model.eps; 
                else              % perturbation is absent in region below the wedge
                    a = model.alpha;
                    b = model.omega;
                end
            % if the perturbation is uniform, perturbation is present over the full region    
            else
                a = model.alpha + model.eps;     
                b = model.omega - model.eps;     
            end
        end
        
        function domain = checkdomain(~,xinit)
            domain = 0;
            x=xinit(1);
            y=xinit(2);
            %             if ~model.varOn
            %                 x = model.xinit(1);
            %                 y = model.xinit(2);
            
            if x==1
                domain = 1;
            elseif y==1
                domain = 2;
            elseif x==-1
                domain = 3;
            elseif y==-1
                domain = 4;
            end
        end
        %         end
    end
end

%% Jump matrices
function Jy=jumpy() % Jump matrix when entering y=-1 or 1 backward in time
Jy=[1 0;0 0];
end

function Jx=jumpx() % Jump matrix when entering x=1 or -1 backward in time
Jx=[0 0;0 1];
end

function tf = inWedge(x,y)
if (x + y >=0 && y - x >=0)
    tf = true;
else
    tf = false;
end
end