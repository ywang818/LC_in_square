% Solve the nonhomogeneous variational equation for the shape response curve under static perturbation: alpha -> alpha+eps

% Need to run prc_plot.m first to find nu1


T0=6.766182958128617; % intrinsic period of the oscillator
nu1=0.498678192193239; % relative change in frequency computed from prc_plot

alpha = 0.2;
eps = 0.001; % static perturbation on model parameter alpha
alpha_pert = alpha + eps;
varOn = true;
yinit=[1, alpha];   % liftoff point for unperturbed trajectory
yinit_pert=[1, alpha_pert]; % liftoff point for perturbed trajectory
vinit = yinit_pert-yinit;
tmax = 40;

% find the unperturbed solutions and the shape response curves
model = LC_in_square(varOn, yinit, vinit/eps, T0, alpha,nu1,nu1);
model.solve;

% find the new period after perturbation using function findPeriod 
model_pert = LC_in_square(false, yinit_pert, [0 0], tmax, alpha_pert);
model_pert.solve;
Teps=model_pert.findPeriod;

% find perturbed solutions over [0 Teps]
model_pert = LC_in_square(false, yinit_pert, [0 0], Teps, alpha_pert,0);
model_pert.solve;

% rescale time to be the same as unperturbed time
[tspan1, Ind1] = unique((model_pert.t./Teps).*T0,'stable'); % rescale time such that it has time [0 T0]
x_unique = model_pert.yext(Ind1,1:2); % obtain x,y values after time is rescaled
x_pert_interp = interp1(tspan1,x_unique,model.t); % interpolate using unperturbed solution

%%
figure
subplot(2,1,1)
plot(model.t,model.yext(:,1),'linewidth',2)
hold on
plot(model.t,x_pert_interp(:,1),'linewidth',2)

xlim([0 T0])
ylim([-1.1 1.1])
ylabel('$x$','interpreter','latex','fontsize',30,'rot',0)
legend('UnPerturbed','Perturbed')
set(gca,'FontSize',18,'linewidth',0.5)

subplot(2,1,2)
plot(model.t,x_pert_interp(:,1)-model.yext(:,1),'linewidth',2)
hold on
plot(model.t,eps*model.yext(:,3),'linewidth',2)
xlim([0 T0])
xlabel('$\rm time (ms)$','interpreter','latex','fontsize',30)
% ylim([-1.1 1.1])
ylabel('$x$','interpreter','latex','fontsize',30,'rot',0)
legend('numerical displacement','linear response approximation')
set(gca,'FontSize',18,'linewidth',0.5)

title('x-component of SRC','FontWeight','normal','fontsize',26)


figure
subplot(2,1,1)
plot(model.t,model.yext(:,2),'linewidth',2)
hold on
plot(model.t,x_pert_interp(:,2),'linewidth',2)

xlim([0 T0])
ylim([-1.1 1.1])
ylabel('$x$','interpreter','latex','fontsize',30,'rot',0)
legend('UnPerturbed','Perturbed')
set(gca,'FontSize',18,'linewidth',0.5)


subplot(2,1,2)
plot(model.t,x_pert_interp(:,2)-model.yext(:,2),'linewidth',2)
hold on
plot(model.t,eps*model.yext(:,4),'linewidth',2)
xlim([0 T0])
% ylim([-1.1 1.1])
legend('numerical displacement','linear response approximation')
set(gca,'FontSize',18,'linewidth',0.5)
ylabel('$\mathbf{x}_{1,2}$','interpreter','latex','fontsize',30,'rot',0)

title('y-component of SRC','FontWeight','normal','fontsize',26)

