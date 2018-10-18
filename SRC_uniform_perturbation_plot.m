% Solve the nonhomogeneous variational equation for the iSRC 
% under a uniform static perturbation: (alpha, omega) -> (alpha + eps, omega - eps)
% !! Need to run prc_plot.m first to find T1, the linear shift in period

% Consider eps = 0.01


T0=6.766182958128617;  % intrinsic period of the oscillator
T1=10.802054306767772; % relative linear shift in period, computed from prc_plot
nu1=T1/T0;             % rescaling factor used in the equation for iSRC

alpha = 0.2;           % default value for alpha
omega = 1;             % default value for omega
eps = 0.01;            % size of the static perturbation 
alpha_pert = alpha + eps;              % perturbed value for alpha
omega_pert = omega - eps;              % perturbed value for omega
yinit=[1, alpha/omega];                % liftoff point at x=1 for unperturbed trajectory
yinit_pert=[1, alpha_pert/omega_pert]; % liftoff point at x=1 for perturbed trajectory
vinit = (yinit_pert-yinit)/eps;        % initial value for the iSRC
tmax = 80; 

% Find the unperturbed solution and the shape response curve
model = LC_in_square('varOn', true, 'xinit', yinit, 'vinit', vinit, 'tmax', T0, 'nu', nu1);
model.solve;

% Find the period Teps of the LC under perturbation 
model_pert = LC_in_square('xinit', yinit_pert, 'vinit', [0 0], 'tmax', tmax, 'alpha', alpha_pert, 'omega', omega_pert);
model_pert.solve;
Teps=model_pert.findPeriod;

% Find the perturbed solution over [0 Teps]
model_pert = LC_in_square('xinit', yinit_pert, 'vinit', [0 0], 'tmax', Teps, 'alpha', alpha_pert, 'omega', omega_pert);
model_pert.solve;

% Rescale the time of the perturbed LC to be the same as model.t, the time for the unperturbed LC
[tspan1, Ind1] = unique((model_pert.t./Teps).*T0,'stable'); % rescale time such that it has time [0 T0]
x_unique = model_pert.yext(Ind1,1:2); % obtain x,y values after time is rescaled
x_pert_interp = interp1(tspan1,x_unique,model.t); % interpolate using unperturbed solution

% actual displacement obtained from subtracting unperturbed LC from perturbed LC 
displacement=x_pert_interp - model.yext(:,1:2);

%% Plot the actual displacement and the approximated displacement (eps*iSRC) together

figure
subplot(2,1,1)
plot(model.t,displacement(:,1),'k','linewidth',3)
hold on
plot(model.t,eps*model.yext(:,3),'r:','linewidth',2)
xlim([0 T0])
xlabel('\rm time (ms)','interpreter','latex','fontsize',25)
ylabel('$x_{\varepsilon}(t)-x(t)$ ','interpreter','latex','fontsize',25)
legend({'actual','approximation'},'Interpreter','latex')
set(gca,'FontSize',18)
model.draw_wall_contact_rectangles

subplot(2,1,2)
plot(model.t,displacement(:,2),'k','linewidth',3)
hold on
plot(model.t,eps*model.yext(:,4),'r:','linewidth',2)
xlim([0 T0])
xlabel('\rm time (ms)','interpreter','latex','fontsize',25)
ylabel('$y_{\varepsilon}(t)-y(t)$','interpreter','latex','fontsize',25)
legend({'actual','approximation'},'Interpreter','latex')
set(gca,'FontSize',18)
model.draw_wall_contact_rectangles


% figure
% subplot(2,1,1)
% plot(model.t,model.yext(:,1),'linewidth',2)
% hold on
% plot(model.t,x_pert_interp(:,1),'linewidth',2)
% 
% xlim([0 T0])
% ylim([-1.1 1.1])
% ylabel('$x$','interpreter','latex','fontsize',30,'rot',0)
% legend('UnPerturbed','Perturbed')
% set(gca,'FontSize',18,'linewidth',0.5)
% 
% subplot(2,1,2)
% plot(model.t,x_pert_interp(:,1)-model.yext(:,1),'linewidth',2)
% hold on
% plot(model.t,eps*model.yext(:,3),'linewidth',2)
% xlim([0 T0])
% xlabel('$\rm time (ms)$','interpreter','latex','fontsize',30)
% % ylim([-1.1 1.1])
% ylabel('$x$','interpreter','latex','fontsize',30,'rot',0)
% legend('numerical displacement','linear response approximation')
% set(gca,'FontSize',18,'linewidth',0.5)
% 
% title('x-component of SRC','FontWeight','normal','fontsize',26)
% 
% 
% figure
% subplot(2,1,1)
% plot(model.t,model.yext(:,2),'linewidth',2)
% hold on
% plot(model.t,x_pert_interp(:,2),'linewidth',2)
% 
% xlim([0 T0])
% ylim([-1.1 1.1])
% ylabel('$x$','interpreter','latex','fontsize',30,'rot',0)
% legend('UnPerturbed','Perturbed')
% set(gca,'FontSize',18,'linewidth',0.5)
% 
% 
% subplot(2,1,2)
% plot(model.t,x_pert_interp(:,2)-model.yext(:,2),'linewidth',2)
% hold on
% plot(model.t,eps*model.yext(:,4),'linewidth',2)
% xlim([0 T0])
% % ylim([-1.1 1.1])
% legend('numerical displacement','linear response approximation')
% set(gca,'FontSize',18,'linewidth',0.5)
% ylabel('$\mathbf{x}_{1,2}$','interpreter','latex','fontsize',30,'rot',0)
% 
% title('y-component of SRC','FontWeight','normal','fontsize',26)
% 
