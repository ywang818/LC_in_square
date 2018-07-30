% Solve the nonhomogeneous variational equation for the shape response curve 
% under static perturbation: alpha -> alpha+eps over region I

% Need to run local timing response to find nu_above and nu_below, associated with region I and II;

T0=6.766182958186305; % intrinsic oscillator period

alpha = 0.2;  % model parameter
eps = 0.01;   % perturbation on the parameter
varOn = true; % run variational problem
xinit=[1,0];  % initial condition for both perturbed and unperturbed limit cycle solution

% model = LC_in_square(false,xinit,[0 0],T0,alpha,0,0); 
% model.solve
% 
% ind_above_wedge=(model.yext(:,1) + model.yext(:,2) >=0) & (model.yext(:,2) - model.yext(:,1) >=0); % index for above wedge
% time_above_wedge=model.t(ind_above_wedge);  % time above wedge
% x_above_wedge=model.yext(ind_above_wedge,:);
% yinit=x_above_wedge(1,:);
% 
% T0_above=time_above_wedge(end)-time_above_wedge(1); % time elapsed in the wedge
% T0_below=T0-T0_above;

yinit=[0.811100985386121   0.811100985460158]; % coordinate of the unperturbed entry point into the region 1
T0_above=1.691545739499871; % unperturbed total time spent in region I
T0_below=5.074637218686434; % unperturbed total time spent in region II

model_pert = LC_in_square(false, xinit, [0 0], 20*T0, alpha,0.1,0,eps); 
model_pert.solve;
Teps=model_pert.findPeriod;   % find perturbed period

model_pert = LC_in_square(false, xinit, [0 0], Teps, alpha,0.1,0,eps);
model_pert.solve;
ind_above_wedge_pert=(model_pert.yext(:,1) + model_pert.yext(:,2) >=0) & (model_pert.yext(:,2) - model_pert.yext(:,1) >=0); % index for above wedge
x_above_wedge_pert=model_pert.yext(ind_above_wedge_pert,:);
yinit_pert=x_above_wedge_pert(1,:); % find the coordinate of the perturbed entry point into the region 1

% Take yinit_pert as initial cond, compute the perturbed solution on [0 Teps]
model_pert = LC_in_square(false, yinit_pert, [0 0], Teps, alpha,0.1,0,eps);
model_pert.solve;

ind_above_wedge_pert=(model_pert.yext(:,1) + model_pert.yext(:,2) >=0) & (model_pert.yext(:,2) - model_pert.yext(:,1) >=0) & (model_pert.t<6); % index for above wedge; % index for above wedge
time_above_wedge_pert=model_pert.t(ind_above_wedge_pert);  % time above wedge
x_above_wedge_pert=model_pert.yext(ind_above_wedge_pert,1:2);

T0_above_pert=time_above_wedge_pert(end)-time_above_wedge_pert(1); % time elapsed in the wedge
T0_below_pert=Teps-T0_above_pert;

ind_below_wedge_pert=(model_pert.t >= T0_above_pert); % index for above wedge
time_below_wedge_pert=[time_above_wedge_pert(end); model_pert.t(ind_below_wedge_pert)];  % time above wedge
x_below_wedge_pert=[x_above_wedge_pert(end,:); model_pert.yext(ind_below_wedge_pert,1:2)];

vinit = (yinit_pert-yinit)/eps;  % find the initial value for the SRC
disp(['vinit is ' num2str(eps*vinit)])

T1=0.840776293874193; % obtained from prc_plot
T1_above=0.400602124648554; % obtained from local_TRC_plot
nu_above=0.236935852187003; % obtained from local_TRC_plot
nu_below=(T1-T1_above)/T0_below; % compute the nu in the region below the wedge
% nu=T1/T0; nu_above=nu; nu_below=nu; % if rescale using the same relative
% change of the full period, then we wouldn't get a nice result

% Compute the unperturbed solution and the SRC with piecewise nu's
model = LC_in_square(true,yinit,vinit,T0,alpha,nu_below,nu_above);  
model.solve

ind_above_wedge=(model.yext(:,1) + model.yext(:,2) >=0) & (model.yext(:,2) - model.yext(:,1) >=0) & (model.t<6); % index for above wedge; % index for above wedge
time_above_wedge=model.t(ind_above_wedge);  % time above wedge
x_above_wedge=model.yext(ind_above_wedge,1:2);

ind_below_wedge=(model.t >= T0_above); % index for above wedge
time_below_wedge=[time_above_wedge(end); model.t(ind_below_wedge)];  % time above wedge
x_below_wedge=[x_above_wedge(end,:); model.yext(ind_below_wedge,1:2)];


% rescale the new time to be the same as unperturbed time
[tspan_above_wedge, Ind_above_wedge] = unique((time_above_wedge_pert./T0_above_pert).*T0_above,'stable'); % rescale time such that it has time [0 T0_above]
x_unique_above_wedge = x_above_wedge_pert(Ind_above_wedge,:); % obtain xr value after time is rescaled
x_pert_interp_above_wedge = interp1(tspan_above_wedge,x_unique_above_wedge,time_above_wedge); % interpolate using unperturbed solution

[tspan_below_wedge, Ind_below_wedge] = unique(time_above_wedge(end)+ ((time_below_wedge_pert-time_above_wedge_pert(end))./(time_below_wedge_pert(end)-time_above_wedge_pert(end))).*T0_below,'stable'); % rescale time such that it has time [0 T0_above]
x_unique_below_wedge = x_below_wedge_pert(Ind_below_wedge,:); % obtain xr value after time is rescaled
x_pert_interp_below_wedge = interp1(tspan_below_wedge,x_unique_below_wedge,time_below_wedge); % interpolate using unperturbed solution

displacement=[x_pert_interp_above_wedge; x_pert_interp_below_wedge ]- [x_above_wedge; x_below_wedge];
 
%%
figure
subplot(2,1,1)
plot(time_above_wedge,x_above_wedge(:,1),'b','linewidth',3)
hold on
plot(time_above_wedge,x_pert_interp_above_wedge(:,1),'r:','linewidth',2)
plot(time_below_wedge,x_below_wedge(:,1),'b','linewidth',3)
plot([T0_above, T0_above], [-1 1],'g:','linewidth',2)

plot(time_below_wedge,x_pert_interp_below_wedge(:,1),'r:','linewidth',2)
plot([T0_above_pert, T0_above_pert], [-1 1],'g:','linewidth',2)
xlim([0 T0])
ylim([-1.1 1.1])
xlabel('time','interpreter','latex','fontsize',25)
ylabel('$x$','interpreter','latex','fontsize',25,'rot',0)
legend('unperturbed','perturbed')
set(gca,'FontSize',18)


subplot(2,1,2)
plot([time_above_wedge; time_below_wedge],displacement(:,1),'b','linewidth',3)
hold on
plot(model.t,eps*model.yext(:,3),'r:','linewidth',2)
plot([T0_above, T0_above], [-eps eps],'g:','linewidth',2)

xlim([0 T0])
xlabel('time','interpreter','latex','fontsize',25)
ylabel('displacement','interpreter','latex','fontsize',25)
legend('numerics','linear response approximation')
set(gca,'FontSize',18)