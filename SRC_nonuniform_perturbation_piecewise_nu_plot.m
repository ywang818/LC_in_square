% Solve the nonhomogeneous variational equation for the iSRC using piecewise uniform rescaling 
% under a nonuniform static perturbation that is only present above the wedge
%        (alpha, omega) -> (alpha + eps, omega - eps) 
% !! Need to run local_TRC_plot to find T1_above and T1_below, linear shifts in time in regions above and below the wedge

%  Region I: above wedge, Region II: below wedge
%  Consider eps = 0.1


T0=6.766182958186305; % intrinsic oscillator period

eps = 0.1;   % perturbation size
varOn = true; % run variational problem
xinit=[1,0];  % initial condition for LC

% Compute x_in (coordinate of the entry point into region above the wedge), 
% T0_above and T0_below (time that LC spends in regions above and below the wedge)

% model = LC_in_square(false,xinit,[0 0],T0,alpha,0,0); 
% model.solve
% ind_above_wedge=(model.yext(:,1) + model.yext(:,2) >=0) & (model.yext(:,2) - model.yext(:,1) >=0); % index for above wedge
% time_above_wedge=model.t(ind_above_wedge);  % time above wedge
% x_above_wedge=model.yext(ind_above_wedge,:);
% x_in=x_above_wedge(1,:);
% T0_above=time_above_wedge(end)-time_above_wedge(1); % time elapsed in the wedge
% T0_below=T0-T0_above;

x_in=[0.811100985386121   0.811100985460158]; % coordinate of the unperturbed entry point into the region 1
T0_above=1.691545739499871; % time that LC spends in region I
T0_below=5.074637218686434; % time that LC spends in region II

% Find the period Teps of the LC under perturbation 
model_pert = LC_in_square('xinit', xinit, 'vinit', [0 0], ...
    'tmax', 20*T0, 'nu', [0, 0], 'eps', eps); 
model_pert.solve;
Teps=model_pert.findPeriod;   % find perturbed period

% Find x_in_pert, the coordinate of the perturbed entry point into region I 
model_pert = LC_in_square('xinit', xinit, 'vinit', [0 0], 'tmax', Teps, ...
    'nu', [0, 0], 'eps', eps);
model_pert.solve;
ind_above_wedge_pert=(model_pert.yext(:,1) + model_pert.yext(:,2) >=0) & (model_pert.yext(:,2) - model_pert.yext(:,1) >=0); % index for above wedge
x_above_wedge_pert=model_pert.yext(ind_above_wedge_pert,:);
x_in_pert=x_above_wedge_pert(1,:); 

% Take x_in_pert as the initial condition, compute the perturbed LC over [0 Teps]
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

% find the initial value for the iSRC
vinit = (x_in_pert-x_in)/eps;  

% Compute rescaling factors in different regions
T1=2.694391001334606;            % obtained from prc_plot
T1_above=2.167992616327751;      % obtained from local_TRC_plot
nu_above=T1_above/T0_above;      % obtained from local_TRC_plot
nu_below=(T1-T1_above)/T0_below; % compute the nu in the region below the wedge

% Compute the unperturbed solution and the iSRC with piecewise nu's found from local timing response curve
model = LC_in_square('varOn', true, 'xinit', x_in, 'vinit', vinit, ...
    'tmax', T0, 'nu', [nu_below,nu_above]);  
model.solve

% Separate the solution into two segments, one in region I and the other in region II
ind_above_wedge=(model.yext(:,1) + model.yext(:,2) >=0) & (model.yext(:,2) - model.yext(:,1) >=0) & (model.t<6); % index for above wedge; % index for above wedge
time_above_wedge=model.t(ind_above_wedge);  % time above wedge
x_above_wedge=model.yext(ind_above_wedge,1:2);
ind_below_wedge=(model.t >= T0_above); % index for above wedge
time_below_wedge=[time_above_wedge(end); model.t(ind_below_wedge)];  % time above wedge
x_below_wedge=[x_above_wedge(end,:); model.yext(ind_below_wedge,1:2)];

% Rescale the perturbed time in region I to be the same as the unperturbed time for LC in region I
[tspan_above_wedge, Ind_above_wedge] = unique((time_above_wedge_pert./T0_above_pert).*T0_above,'stable'); % rescale time such that it has time [0 T0_above]
x_unique_above_wedge = x_above_wedge_pert(Ind_above_wedge,:); % obtain xr value after time is rescaled
x_pert_interp_above_wedge = interp1(tspan_above_wedge,x_unique_above_wedge,time_above_wedge); % interpolate using unperturbed solution

% Rescale the perturbed time in region II to be the same as the unperturbed time for LC in region II
[tspan_below_wedge, Ind_below_wedge] = unique(time_above_wedge(end)+ ((time_below_wedge_pert-time_above_wedge_pert(end))./(time_below_wedge_pert(end)-time_above_wedge_pert(end))).*T0_below,'stable'); % rescale time such that it has time [0 T0_above]
x_unique_below_wedge = x_below_wedge_pert(Ind_below_wedge,:); % obtain xr value after time is rescaled
x_pert_interp_below_wedge = interp1(tspan_below_wedge,x_unique_below_wedge,time_below_wedge); % interpolate using unperturbed solution

% After rescaling, the time LC under perturbation spends in region I is the  same as T0_above
T0_above_pert_rescale=tspan_above_wedge(end)-tspan_above_wedge(1); % T0_above_pert after rescaling

% actual displacement obtained from subtracting unperturbed LC from perturbed LC
displacement=[x_pert_interp_above_wedge; x_pert_interp_below_wedge ]- [x_above_wedge; x_below_wedge];

%% plot actual displacement obtained from numerical computation vs approximated disp using iSRC
figure
set(gcf,'Position',[0 0 720 520])
subplot(2,1,1)
plot([time_above_wedge; time_below_wedge],displacement(:,1),'k','linewidth',3)
hold on
plot(model.t,eps*model.yext(:,3),'r:','linewidth',2)
plot([T0_above_pert_rescale, T0_above_pert_rescale], [-0.2 0.2],'m','linewidth',2)
plot([T0_above, T0_above], [-0.2 0.2],'b:','linewidth',2)
xlim([0 T0])
ylim([-0.2 0.2])
xlabel('\rm time (ms)','interpreter','latex','fontsize',25)
ylabel('$x_{\varepsilon}(t)-x(t)$','interpreter','latex','fontsize',25)
legend({'actual','approximation'},'Interpreter','latex')
set(gca,'FontSize',18)
model.draw_wall_contact_rectangles
title('Piecewise uniform rescaling','Interpreter','latex','FontWeight','normal','Fontsize',20)

subplot(2,1,2)
plot([time_above_wedge; time_below_wedge],displacement(:,2),'k','linewidth',3)
hold on
plot(model.t,eps*model.yext(:,4),'r:','linewidth',2)
plot([T0_above_pert_rescale, T0_above_pert_rescale], [-0.2 0.2],'m','linewidth',2)
plot([T0_above, T0_above], [-0.2 0.2],'b:','linewidth',2)
xlim([0 T0])
ylim([-0.2 0.2])
xlabel('\rm time (ms)','interpreter','latex','fontsize',25)
ylabel('$y_{\varepsilon}(t)-y(t)$','interpreter','latex','fontsize',25)
legend({'actual','approximation'},'Interpreter','latex')
set(gca,'FontSize',18)
model.draw_wall_contact_rectangles
