% Solve the nonhomogeneous variational equation for the iSRC using piecewise uniform rescaling 
% under a nonuniform static perturbation that is only present above the wedge
%        (alpha, omega) -> (alpha + eps, omega - eps)
% !! Need to run local_TRC_plot to find T1_above and T1_below, linear shifts in time in regions above and below the wedge

%  Region I: above wedge, Region II: below wedge

%  For eps over (0.00006:0.01:0.1), compute the norms of the actual (LC_eps - LC) and
%  approximated (eps*iSRC) displacements between perturbed and unperturbed LC solutions

T0=6.766182958186305; % intrinsic oscillator period
xinit=[1,0];  % initial condition for LC

% see 'SRC_nonuniform_perturbation_piecewise_nu_plot' for how to find the following values
x_in=[0.811100985386121   0.811100985460158]; % coordinate of the unperturbed entry point into the region 1
T0_above=1.691545739499871; % unperturbed total time spent in region I
T0_below=5.074637218686434; % unperturbed total time spent in region II

% Compute rescaling factors in different regions
T1=2.694391001334606; % obtained from prc_plot
T1_above=2.167992616327751; % obtained from local_TRC_plot
nu_above=T1_above/T0_above; % obtained from local_TRC_plot
nu_below=(T1-T1_above)/T0_below; % compute the nu in the region below the wedge

epsvec = (0.00006:0.01:0.1);     % perturbation size
normu       = zeros(size(epsvec));     % norm of approximated displacements
normNumDisp = zeros(size(epsvec));     % norm of actual displacements
normdiff    = zeros(size(epsvec));     % relative difference between the actual norm and approximated norm

% Initial condition for iSRC at the entry point into region I,
% see 'SRC_nonuniform_perturbation_piecewise_nu_plot' for how to find it
vinit = [3.074e-10 -4.7892e-10];

for i=1:length(epsvec)
    
    eps=epsvec(i);   % perturbation size
    x_in_pert = vinit*eps + x_in; % approximated initial condition for the perturbed LC
    
    % Find the period Teps of the LC under perturbation
    model_pert = LC_in_square('xinit', xinit, 'vinit', [0 0], ...
        'tmax', 20*T0, 'nu', [0, 0], 'eps', eps);
    model_pert.solve;
    Teps=model_pert.findPeriod;   % find perturbed period
    
    % Take x_in_pert as initial cond, compute the perturbed solution on [0 Teps]
    model_pert = LC_in_square('xinit', x_in_pert, 'vinit', [0 0], 'tmax', Teps, ...
        'nu', [0, 0], 'eps', eps);
    model_pert.solve;
    
    % Separate the solution into two segments, one in region I and the other in region II
    ind_above_wedge_pert=(model_pert.yext(:,1) + model_pert.yext(:,2) >=0) & (model_pert.yext(:,2) - model_pert.yext(:,1) >=0) & (model_pert.t<6); % index for above wedge; % index for above wedge
    time_above_wedge_pert=model_pert.t(ind_above_wedge_pert);  % time above wedge
    x_above_wedge_pert=model_pert.yext(ind_above_wedge_pert,1:2);
    T0_above_pert=time_above_wedge_pert(end)-time_above_wedge_pert(1); % time elapsed in the wedge
    T0_below_pert=Teps-T0_above_pert;    
    ind_below_wedge_pert=(model_pert.t >= T0_above_pert); % index for above wedge
    time_below_wedge_pert=[time_above_wedge_pert(end); model_pert.t(ind_below_wedge_pert)];  % time above wedge
    x_below_wedge_pert=[x_above_wedge_pert(end,:); model_pert.yext(ind_below_wedge_pert,1:2)];
    
    % Compute the unperturbed solution and the iSRC with piecewise nu's
    model = LC_in_square('varOn', true, 'xinit', x_in, 'vinit', vinit, ...
        'tmax', T0, 'nu', [nu_below,nu_above]);
    model.solve
    
    % Separate the solution into two segments, one in region I and the other in region II
    ind_above_wedge=(model.yext(:,1) + model.yext(:,2) >=0) & (model.yext(:,2) - model.yext(:,1) >=0) & (model.t<6); % index for above wedge; % index for above wedge
    time_above_wedge=model.t(ind_above_wedge);  % time above wedge
    x_above_wedge=model.yext(ind_above_wedge,1:2);
    u_above_wedge=model.yext(ind_above_wedge,3:4);
    
    ind_below_wedge=(model.t >= T0_above); % index for above wedge
    time_below_wedge=[time_above_wedge(end); model.t(ind_below_wedge)];  % time above wedge
    x_below_wedge=[x_above_wedge(end,:); model.yext(ind_below_wedge,1:2)];
    u_below_wedge=[u_above_wedge(end,:); model.yext(ind_below_wedge,3:4)];
    
    
    % rescale the perturbed time to be the same as the unperturbed time in region I and II, respectively
    [tspan_above_wedge, Ind_above_wedge] = unique((time_above_wedge_pert./T0_above_pert).*T0_above,'stable'); % rescale time such that it has time [0 T0_above]
    x_unique_above_wedge = x_above_wedge_pert(Ind_above_wedge,:); % obtain xr value after time is rescaled
    x_pert_interp_above_wedge = interp1(tspan_above_wedge,x_unique_above_wedge,time_above_wedge); % interpolate using unperturbed solution
    
    [tspan_below_wedge, Ind_below_wedge] = unique(time_above_wedge(end)+ ((time_below_wedge_pert-time_above_wedge_pert(end))./(time_below_wedge_pert(end)-time_above_wedge_pert(end))).*T0_below,'stable'); % rescale time such that it has time [0 T0_above]
    x_unique_below_wedge = x_below_wedge_pert(Ind_below_wedge,:); % obtain xr value after time is rescaled
    x_pert_interp_below_wedge = interp1(tspan_below_wedge,x_unique_below_wedge,time_below_wedge); % interpolate using unperturbed solution
    
    % actual displacement obtained from subtracting unperturbed LC from perturbed LC
    displacement=[x_pert_interp_above_wedge; x_pert_interp_below_wedge ]- [x_above_wedge; x_below_wedge];
    
    % evaluate the norm of the actual displacement vector
    displacement(isnan(displacement))=0;
    normNumDisp(i)=norm(displacement);
    
    % approximated displacement using iSRC
    u = [u_above_wedge; u_below_wedge];

    % evaluate the norm of the approximated displacement
    normu(i)=norm(eps*u);
    
    % compute the relative difference between the actual norm and approximated norm
    normdiff(i)=(normNumDisp(i)-normu(i))/normNumDisp(i); 
        
end

%% Plot the norm of the actual displacements and the approximated displacements (eps*iSRC) vs eps
figure
plot(epsvec,normNumDisp,'b','linewidth',2)
hold on
plot(epsvec,normu,'b:','linewidth',2)
xlabel('$\varepsilon$','interpreter','latex','fontsize',30)
ylabel('$\rm norm(displacement)$','interpreter','latex','fontsize',30)
set(gca,'FontSize',18)

% Plot the relative difference between the norms vs eps
figure
plot(epsvec,normdiff,'b','linewidth',2)
xlabel('$\varepsilon$','interpreter','latex','fontsize',30)
ylabel('relative difference in norm','interpreter','latex','fontsize',30)
set(gca,'FontSize',18)

