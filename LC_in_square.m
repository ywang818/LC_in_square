classdef LC_in_square < handle
    
    % Modified by YW, 2018-07-26
    % Class-based model in Cube Solver [JPG: "Cube"? LC in 3D?]
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
    %   M = LC_in_square(varOn,xinit,vinit,tmax,a,nu_below,nu_above,eps)
    %
    % where 
    %       varon    -- 'false' by default, setting to 'true' the model will solve the variational problem
    %       xinit,vinit -- the initial values of x and v (variational problem) [JPG: u is used for pert throughout paper. can we use u instead of v in code?]
    %       tmax        -- the endtime of the ode integration  
    %       a        -- model parameter alpha with default value 0.2
    %       nu_below,nu_above -- relative change in frequency for regions below and above the wedge defined by function inWedge [JPG: these are the time rescaling factors?]
    %       eps      -- perturbation added to parameter alpha [JPG: persistent? in wedge only?]
    % You can enter any number of inputs to the constructor since all of them have pre-defined default values.
    % However, inputs determine which problem will be solved
    %
    %     a) To only find solution trajectory of the model, run 
    %
    %                     M = LC_in_square(false,xinit,vinit,tmax,alpha) where vinit can be any 1 by 2 vector [JPG: vinit is discarded if varOn is false?]
    %
    %     b) To solve variational problem w.r.t. instantaneous perturbation, run
    %
    %                     M = LC_in_square(true, xinit, vinit, T0, alpha) [JPG: "... where T0 is exactly one period"]
    %
    %     c) For a uniform static perturbation, 
    %          i) to find unperturbed solution and the shape response curve, construct the model with the same nonzero nu as input:
    %
    %                     M = LC_in_square(true, xinit, vinit, T0, alpha,nu1,nu1)
    %             where nu_below = nu_above = nu1 needs to be computed using iPRC [JPG: where in the paper does computing nu's from iPRC happen?]
    %
    %          ii)to find perturbed solution with alpha -> alpha+eps, run
    %
    %                     M = LC_in_square(false, xinit_perturb, vinit, tmax,alpha+eps) or 
    %                     M = LC_in_square(false, xinit_perturb,vinit,tmax, alpha,nu,nu,eps) [JPG: does "or" here mean they are equivalent? passing eps as the final arg is the same as adding it to alpha in the first case?]
    %
    %     d) For piecewise static perturbations wher the static perturbation only exists in the region above the wedge,
    %        type in two different nu and the first one is nu_below and the second one is nu_above:
    %
    %                     M = LC_in_square(varOn,xinit,vinit,tmax,alpha,nu_below,nu_above,eps),
    %        if nu_below~=nu_above, eps will only be added to alpha in the region above the wedge
    %    
    %     e) To find iPRC,
    %                    model = LC_in_square(false, xinit);
    %                    model.solve;
    %                    model.find_prc(z0);   % z0 is the initial condition for iPRC
    %                    model.plot_prc
    %
    % 2, After construct a model, typing M in command window and entering will display
    % all current properties set. In addition to the above four, there are also: 
    %
    %           yinit   -- concatenation of xinit and vinit [JPG: "... if varOn is true"]
    %          domain   -- current domain (0 interior, 1-4 walls)
    %               t   -- time array
    %            yext   -- full solution array, each new solution will be appended to new row (yext = [yext; ynew]) [JPG: "ext" means what?]
    %             prc   -- phase response curve
    %              t0   -- current time (scalar) [JPG: t0 not a public property]
    %              y0   -- current solution (same size as yinit) [JPG: y0 not a public property]
    %            Salt   -- array of Saltation matrix (v+ = Sv-) [JPG: Salt not defined. meant "S1"? also not a public property]
    %         t_event   -- keeping the times of events (wall hitting) [JPG: t_event not a public property]
    %    domain_event   -- domains where trajectory enters [JPG: domain_event not defined at all]
    %          t_exit   -- keeping the times of leaving walls
    %     domain_exit   -- domains where trajectory leaves [JPG: domain_exit not defined at all]
    %       Jump_exit   -- array of Jump matrix at exit (z- = Jz+)
    % 
    % It is not recommended to manually change any of these values because
    % they are all counted as model "outputs"
    % You can check their values any time using the dot reference, such as M.grasper [JPG: grasper not defined at all]
    %  
    % 3, Type "methods(M)" will show the available methods in the model:  
    % 
    %           solve   -- requires no input, solve the model with given initial values
    %            plot   -- requires no input, plot the solutions (alway solve before plot) [JPG: also plot_prc]
    %      findPeriod   -- requires two guesses L and R, using Bisection method to find an estimated period in [L,R] of the ODE [JPG: this method doesn't actually take any args]
    %          LC_ODE   -- ODE of the model, used internally [JPG: LC_ODE not a public M method]
    %
    % ===========================
    % Some other usage examples:
    % ===========================
    % ** Estimate solution period
    % 
    %   >> M = LC_in_square;
    %   >> M.findPeriod [JPG: this example doesn't work. must both increase tmax and then solve first.]
    %   ...................
    %   ans =
    %
    %        6.766
    
    

    
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
        S0; % [JPG: this is initialized but never used. S1 is used later. is S0 obsolete?]
        t0
        y0
    end
    
    methods
        function model = LC_in_square(varOn,xinit,vinit,tmax,a,nu_below,nu_above,eps)
            % Default values
%             model.nu_below = 0.478781045930474;
            model.nu_below = 0;
            model.nu_above = 0;
            model.eps = 0;
            model.alpha = 0.2;
            model.tmax = 6.766182958128617;
            model.vinit = [0.01, 0];
            model.xinit = [0.9, -0.9];
            model.varOn = false;
            if nargin > 0
                model.varOn = varOn;
            end
            if nargin > 1
                model.xinit = reshape(xinit,1,length(xinit));
            end
            if nargin > 2
                model.vinit = reshape(vinit,1,length(vinit));
            end
            if nargin > 3
                model.tmax = tmax;
            end
            if nargin > 4
                model.alpha = a;
            end
            if nargin > 5
                model.nu_below = nu_below;
            end
            if nargin > 6
                model.nu_above = nu_above;
            end
            if nargin > 7
                model.eps = eps;
            end
            
            model.domain = 0;
            if model.varOn
                model.yinit = [model.xinit, model.vinit];
                model.S0 = eye(2); % [JPG: never used]
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
            model.domain = 0; % [JPG: assuming we start in interior. what if start on wall?]
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
                disp(['entering domain ', num2str(model.domain), ' at time t = ', num2str(model.t0)])
                switch model.domain
                    case 0 % interior
                        if model.varOn
                            [tnew,ynew,TE,YE,IE] = ode45(@model.LC_ODE_ext,[model.t0,model.tmax],model.y0,options0);
                        else
                            [tnew,ynew,TE,YE,IE] = ode45(@model.LC_ODE,[model.t0,model.tmax],model.y0,options0);
                        end
                        model.updateSolution(tnew,ynew);
                        model.updateCurrent(tnew,ynew,TE,YE,IE);
                        
                        model.domain = IE;
                    case 1 % x=1 wall
                        if model.varOn
                            model.multiplySaltation(TE,YE,'enter');
                            [tnew,ynew,TE,YE,~] = ode45(@model.LC_ODE_ext,[model.t0,model.tmax],model.y0,options1_ext);
                        else
                            [tnew,ynew,TE,YE,~] = ode45(@model.LC_ODE,[model.t0,model.tmax],model.y0,options1);
                        end
                        
                        model.updateSolution(tnew,ynew);
                        model.updateCurrent(tnew,ynew,TE,YE,0);
                        model.storeJump(TE);
                        model.domain=0;
                    case 2 % y=1 wall
                        if model.varOn
                            model.multiplySaltation(TE,YE,'enter');
                            [tnew,ynew,TE,YE,~] = ode45(@model.LC_ODE_ext,[model.t0,model.tmax],model.y0,options2_ext);
                        else
                            [tnew,ynew,TE,YE,~] = ode45(@model.LC_ODE,[model.t0,model.tmax],model.y0,options2);
                        end
                        model.updateSolution(tnew,ynew);
                        model.updateCurrent(tnew,ynew,TE,YE,0);
                        model.storeJump(TE);
                        model.domain=0;
                    case 3 % x=-1 wall
                        if model.varOn
                            model.multiplySaltation(TE,YE,'enter');
                            [tnew,ynew,TE,YE,~] = ode45(@model.LC_ODE_ext,[model.t0,model.tmax],model.y0,options3_ext);
                        else
                            [tnew,ynew,TE,YE,~] = ode45(@model.LC_ODE,[model.t0,model.tmax],model.y0,options3);
                        end
                        model.updateSolution(tnew,ynew);
                        model.updateCurrent(tnew,ynew,TE,YE,0);
                        model.storeJump(TE);
                        model.domain=0;
                    case 4 % y=-1 wall
                        if model.varOn
                            model.multiplySaltation(TE,YE,'enter');
                            [tnew,ynew,TE,YE,~] = ode45(@model.LC_ODE_ext,[model.t0,model.tmax],model.y0,options4_ext);
                        else
                            [tnew,ynew,TE,YE,~] = ode45(@model.LC_ODE,[model.t0,model.tmax],model.y0,options4);
                        end
                        model.updateSolution(tnew,ynew);
                        model.updateCurrent(tnew,ynew,TE,YE,0);
                        model.storeJump(TE);
                        model.domain=0;
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
            [~,ind] = findpeaks(model.yext(:,2),'MinPeakProminence',1);
            if length(ind) < 10
                warning(['!!! Caution, the period found might be inaccurate, '...
                'only found %d peaks (target: 10), consider using longer tmax'], length(ind));
            end
            T = mean(diff(model.t(ind(1:end))));
        end
        
        function find_prc(model, z0)
            
            if nargin < 2
                z0 = [1 0];
            end
            
            model.prct = [];
            model.prc = [];
            options_prc=odeset('BDF','on','RelTol',model.reltol,...
                'AbsTol',model.abstol,'Events',@model.exit_wall);
           
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
                    case 0
                        T = model.reverseTspan(model.reverseTspan <= TE);
                        if T == 0
                            break;
                        end
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
                    case 1
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
            
            a = model.alphaFcn(y(1),y(2));
            dydt=[a,-1;1,a]*y;
            
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
                y(2)+1];    % when y crosses -1 from below (wall 4) [JPG: should say "from above"]
            isterminal=[1;1;1;1]; % stop integration and return
            direction=[1;1;-1;-1]; % "value" should be increasing
        end
        
        function [value,isterminal,direction]=wall1_exit(model,~,y)
            % when the *unconstrained* value of dx/dt decreases through zero, return
            y(1)=1;
            a = model.alphaFcn(y(1),y(2));
            dydt=[a,-1;1,a]*y;
            value=dydt(1);
            isterminal=1;
            direction=-1;
            
        end
        
        function [value,isterminal,direction]=wall2_exit(model,~,y)
            % when the *unconstrained* value of dy/dt decreases through zero, return
            y(2)=1;
            a = model.alphaFcn(y(1),y(2));
            dydt=[a,-1;1,a]*y;
            value=dydt(2);
            isterminal=1;
            direction=-1;
            
        end
        
        function [value,isterminal,direction]=wall3_exit(model,~,y)
            % when the *unconstrained* value of dx/dt increases through zero, return
            y(1)=-1;
            a = model.alphaFcn(y(1),y(2));
            dydt=[a,-1;1,a]*y;
            value=dydt(1);
            isterminal=1;
            direction=1;
            
        end
        
        function [value,isterminal,direction]=wall4_exit(model,~,y)
            % when the *unconstrained* value of dy/dt increases through zero, return
            y(2)=-1;
            a = model.alphaFcn(y(1),y(2));
            dydt=[a,-1;1,a]*y;
            value=dydt(2);
            isterminal=1;
            direction=1;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% ODE for Variational Problem %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dydt=LC_ODE_ext(model,~,y,domain_overload) % flow on the interior
            if nargin > 3
                ODEdomain = domain_overload;
            else
                ODEdomain = model.domain;
            end
            x=y(1:2);
            v=y(3:4);
            a = model.alphaFcn(x(1),x(2));
            dxdt=[a,-1;1,a]*x;
            switch ODEdomain
                case 0
                    dvdt=[a,-1;1,a]*v;
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
            % nu * F(x(t)) + DF/Deps; DF/Deps=[x; y] 
            dvdt = dvdt + model.addon(ODEdomain,x(1),x(2),dxdt);
            dydt=[dxdt;dvdt];
        end
        
        function dzdt=LC_ODE_prc(model,t,z,xmat) % flow on the interior
            
            xvec = interp1(model.reverseTspan,xmat,t);
            
            a = model.alphaFcn(xvec(1),xvec(2));
            DF=[a,-1;1,a];
            switch model.checkdomain(xvec)
                case 1 % x=1 wall
                    %dxdt(1)=min(dxdt(1),0); % only allow negative dx/dt or else zero
                    DF=[0,0;0,a];
                case 2 % y=1 wall
                    %dxdt(2)=min(dxdt(2),0); % only allow negative dy/dt or else zero
                    DF=[a,0;0,0];
                case 3 % x=-1 wall
                    %dxdt(1)=max(dxdt(1),0); % only allow positive dx/dt or else zero
                    DF=[0,0;0,a];
                case 4 % y=-1 wall
                    %dxdt(2)=max(dxdt(2),0); % only allow positive dy/dt or else zero
                    DF=[a,0;0,0];
            end
            dzdt=-DF'*[z(1); z(2)];
        end
        
        function [value,isterminal,direction]=wall1_exit_ext(model,~,y)
            % when the *unconstrained* value of dx/dt decreases through zero, return
            y(1)=1;
            a = model.alphaFcn(y(1),y(2));
            value=y(2)-a;
            isterminal=1;
            direction=1;
            
        end
        
        function [value,isterminal,direction]=wall2_exit_ext(model,~,y)
            % when the *unconstrained* value of dy/dt decreases through zero, return
            y(2)=1;
            a = model.alphaFcn(y(1),y(2));
            value=y(1)+a;
            isterminal=1;
            direction=-1;
            
        end
        
        function [value,isterminal,direction]=wall3_exit_ext(model,~,y)
            % when the *unconstrained* value of dx/dt increases through zero, return
            y(1)=-1;
            a = model.alphaFcn(y(1),y(2));
            value=y(2)+a;
            isterminal=1;
            direction=-1;
            
        end
        
        function [value,isterminal,direction]=wall4_exit_ext(model,~,y)
            % when the *unconstrained* value of dy/dt increases through zero, return
            y(2)=-1;
            a = model.alphaFcn(y(1),y(2));
            value=y(1)-a;
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
        
        function removeError(model, IE) % [JPG: does this function correct the trajectory when a wall is encountered? is this necessary because the solver would otherwise just barely overstep and put us outside the square?]
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
            S1=[];
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
                    error("flag for multiplySaltation must be 'enter' or 'exit'!");
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
                S1=eye(2)+(F_plus-F_minus)*n/(n*F_minus); % saltation matrix, eq. 3.22
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
            
            if (model.nu_above == 0) && (model.nu_below == 0) % if both nu's are zero, the perturbation is instantaneous, so the nonhomogeneous addon is zero
                h = 0;
            end
            
            if (model.nu_above == model.nu_below) && (model.nu_below ~= 0) % if the input nu's are the same and nonzero
                h = model.nu_below*dxdt;
                if ODEdomain == 0 
                    h(1) = h(1) + x;
                    h(2) = h(2) + y;
                elseif ODEdomain == 1 || ODEdomain == 3
                    h(2) = h(2) + y;
                elseif ODEdomain == 2 || ODEdomain == 4
                    h(1) = h(1) + x;
                end
            end
            
            if (model.nu_above ~= model.nu_below) && (model.nu_above ~=0)&& (model.nu_below ~=0) % if the nu's are different and nonzero in different subregions
                if ~inWedge(x,y) % in region below wedge
                    h = model.nu_below*dxdt;
                else % in region above the wedge
                    h = model.nu_above*dxdt;
                    if ODEdomain == 0
                        h(1) = h(1) + x;
                        h(2) = h(2) + y;
                    elseif ODEdomain == 1 || ODEdomain == 3
                        h(2) = h(2) + y;
                    elseif ODEdomain == 2 || ODEdomain == 4
                        h(1) = h(1) + x;
                    end
                end
            end
        end
        
        
        function a = alphaFcn(model, x, y)
            if model.nu_above ~= model.nu_below % if there are regions with different nu's
                if inWedge(x, y)
                    a = model.alpha + model.eps; % perturbation is only present in the given region
                else
                    a = model.alpha;             
                end
            elseif model.nu_above == model.nu_below % if nu's are the same 
                a = model.alpha + model.eps;     % perturbation is present all the time
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