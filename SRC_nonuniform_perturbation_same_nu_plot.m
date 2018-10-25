% Solve the nonhomogeneous variational equation for the iSRC using uniform rescaling 
% under a nonuniform static perturbation that is only present above the wedge
%        (alpha, omega) -> (alpha + eps, omega - eps) 
% !! Need to run prc_plot to find T1, linear shift in period

%  Region I: above wedge, Region II: below wedge
%  Consider eps = 0.1

T0=6.766182958186305; % intrinsic oscillator period

eps = 0.1;   % perturbation on the parameter
varOn = true; % run variational problem
xinit=[1,0];  % initial condition for LC

% see 'SRC_nonuniform_perturbation_piecewise_nu_plot' for how to find the following values
x_in=[0.811100985386121   0.811100985460158]; % coordinate of the unperturbed entry point into the region 1
T0_above=1.691545739499871; % unperturbed total time spent in region I

% Find the period Teps of the LC under perturbation 
model_pert = LC_in_square('xinit', xinit, 'vinit',[0 0],...
    'tmax', 20*T0, 'nu', [0,0], 'eps', eps); 
model_pert.solve;
Teps=model_pert.findPeriod;   

% Find x_in_pert, the coordinate of the perturbed entry point into region I 
model_pert = LC_in_square('xinit', xinit, 'vinit',[0 0],...
    'tmax', Teps, 'nu', [0,0], 'eps', eps);
model_pert.solve;
ind_above_wedge_pert=(model_pert.yext(:,1) + model_pert.yext(:,2) >=0) & (model_pert.yext(:,2) - model_pert.yext(:,1) >=0); % index for above wedge
x_above_wedge_pert=model_pert.yext(ind_above_wedge_pert,:);
x_in_pert=x_above_wedge_pert(1,:); 

% Take x_in_pert as the new initial point, compute the perturbed solution over [0 Teps]
model_pert = LC_in_square('xinit', x_in_pert, 'vinit',[0 0],...
    'tmax', Teps, 'nu', [0,0], 'eps', eps);
model_pert.solve;

% Compute the initial value for the iSRC
vinit = (x_in_pert-x_in)/eps;  

% Compute the global uniform rescaling
T1=2.694391001334606; % obtained from prc_plot % the relative linear shift in period, computed from prc_plot
nu=T1/T0; 

% Compute the unperturbed solution and the iSRC with uniform rescaing
model = LC_in_square('varOn', true, 'xinit', x_in, 'vinit', vinit, ...
    'tmax', T0, 'nu', [nu,nu]);  
model.solve

% Rescale the time of the perturbed LC to be the same as model.t, the time for the unperturbed LC
[tspan1, Ind1] = unique((model_pert.t./Teps).*T0,'stable'); % rescale time such that it has time [0 T0]
x_unique = model_pert.yext(Ind1,1:2); % obtain x,y values after time is rescaled
x_pert_interp = interp1(tspan1,x_unique,model.t); % interpolate using unperturbed solution

% find the time that LC under perturbation spends in the region above the wedge
ind_above_wedge_pert=(x_pert_interp(:,1) + x_pert_interp(:,2) >=0) & (x_pert_interp(:,2) - x_pert_interp(:,1) >=0) & (model.t<6); % index for above wedge; % index for above wedge
time_above_wedge_pert=model.t(ind_above_wedge_pert);  % time above wedge after rescaling
%x_above_wedge_pert=x_pert_interp(ind_above_wedge_pert,1:2);
T0_above_pert=time_above_wedge_pert(end)-time_above_wedge_pert(1); % time elapsed in the wedge

% actual displacement obtained from subtracting unperturbed LC from perturbed LC
displacement=x_pert_interp - model.yext(:,1:2);
 
%% plot actual displacement obtained from numerical computation vs approximated disp using iSRC
figure
set(gcf,'Position',[0 0 720 520])
subplot(2,1,1)
plot(model.t,displacement(:,1),'k','linewidth',3)
hold on
plot(model.t,eps*model.yext(:,3),'r:','linewidth',2)
plot([T0_above_pert, T0_above_pert], [-0.2 0.2],'m','linewidth',2)
plot([T0_above, T0_above], [-0.2 0.2],'b:','linewidth',2)
xlim([0 T0])
ylim([-0.2 0.2])
xlabel('\rm time (ms)','interpreter','latex','fontsize',25)
ylabel('$x_{\varepsilon}(t)-x(t)$','interpreter','latex','fontsize',25)
legend({'actual','approximation'},'Interpreter','latex')
set(gca,'FontSize',18)
model.draw_wall_contact_rectangles

title('Uniform rescaling','Interpreter','latex','FontWeight','normal','Fontsize',20)

subplot(2,1,2)
plot(model.t,displacement(:,2),'k','linewidth',3)
hold on
plot(model.t,eps*model.yext(:,4),'r:','linewidth',2)
plot([T0_above_pert, T0_above_pert], [-0.2 0.2],'m','linewidth',2)
plot([T0_above, T0_above], [-0.2 0.2],'b:','linewidth',2)
xlim([0 T0])
ylim([-0.2 0.2])
xlabel('\rm time (ms)','interpreter','latex','fontsize',25)
ylabel('$y_{\varepsilon}(t)-y(t)$','interpreter','latex','fontsize',25)
legend({'actual','approximation'},'Interpreter','latex')
set(gca,'FontSize',18)
model.draw_wall_contact_rectangles

