% Generate solutions to the homogeneous variatinal equation 
% with ICs [-1, 0] and [0, 1] at the liftoff point on the wall x=1
alpha=0.2;
S=[];
T0=6.766182958128606; % intrinsic period of the oscillator
yinit=[1, alpha]; % liftoff point
vinit=[0, 0.1];

modele2 = LC_in_square('varOn', true, 'xinit', yinit,'vinit', vinit,'tmax', T0);
modele2.solve

% find the trajectory when the initial condition at the liftoff point is perturbed by [0 0.1] 
model_pert = LC_in_square('varOn', true, 'xinit', yinit + vinit,'tmax', T0);
model_pert.solve;

% Rescale the time of the perturbed LC to be the same as model.t, the time for the unperturbed LC
[tspan1, Ind1] = unique(model_pert.t,'stable'); % get rid of repetitions
x_unique = model_pert.yext(Ind1,1:2); % obtain x,y values after repetitions are removed
x_pert_interp = interp1(tspan1,x_unique,modele2.t); % interpolate the perturbed soltuion using unperturbed time vector

% actual displacement obtained from subtracting unperturbed LC from perturbed LC 
displacement=x_pert_interp - modele2.yext(:,1:2);

%% 
figure
set(gcf,'Position',[0 0 1200 830])
subplot(2,2,1)
plot(modele2.t, modele2.yext(:,1),'k', 'linewidth',2)
hold on
plot(model_pert.t, model_pert.yext(:,1),'r:','linewidth',3)
xlim([0 T0])
ylim([-1.1 1.1])
xlabel('$\rm time$','interpreter','latex','fontsize',25)
legend('original','perturbed','Location','north')
set(gca,'FontSize',18)
text(-2,0.9, '$\textbf{(A)}$','Interpreter','latex','FontSize',28,'Color','k')
ylabel('$x$','interpreter','latex','fontsize',30,'rot',0)
modele2.draw_wall_contact_rectangles

subplot(2,2,2)
plot(modele2.t,displacement(:,1),'k', 'linewidth',2)
hold on
plot(modele2.t, modele2.yext(:,3),'r:', 'linewidth',3)
xlim([0 T0])
ylim([-0.15 0.15])
xlabel('$\rm time$','interpreter','latex','fontsize',25)
legend('actual','approximation','Location','northwest')
set(gca,'FontSize',18)
ylabel('$\textbf{u}_x$','interpreter','latex','fontsize',25,'rot',0)
text(-2,0.13, '$\textbf{(B)}$','Interpreter','latex','FontSize',28,'Color','k')
modele2.draw_wall_contact_rectangles

subplot(2,2,3)
plot(modele2.t, modele2.yext(:,2),'k', 'linewidth',2)
hold on
plot(model_pert.t, model_pert.yext(:,2),'r:','linewidth',3)
xlim([0 T0])
ylim([-1.1 1.1])
xlabel('$\rm time$','interpreter','latex','fontsize',25)
legend('original','perturbed','Location','northeast')
set(gca,'FontSize',18)
text(-2,0.9, '$\textbf{(C)}$','Interpreter','latex','FontSize',28,'Color','k')
ylabel('$y$','interpreter','latex','fontsize',30,'rot',0)
modele2.draw_wall_contact_rectangles

subplot(2,2,4)
plot(modele2.t,displacement(:,2), 'k', 'linewidth',2)
hold on
plot(modele2.t, modele2.yext(:,4),'r:', 'linewidth',3)
xlim([0 T0])
ylim([-0.15 0.15])
xlabel('$\rm time$','interpreter','latex','fontsize',25)
legend('actual','approximation','Location','north')
set(gca,'FontSize',18)
ylabel('$\textbf{u}_y$','interpreter','latex','fontsize',25,'rot',0)
modele2.draw_wall_contact_rectangles

text(-2,0.13, '$\textbf{(D)}$','Interpreter','latex','FontSize',28,'Color','k')