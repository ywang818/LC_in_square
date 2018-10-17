% Generate solutions to the homogeneous variatinal equation 
% with ICs [-1, 0] and [0, 1] at the liftoff point on the wall x=1

S=[];
alpha=0.2;
T0=6.766182958128606; % intrinsic period of the oscillator for alpha=0.2
yinit=[1, alpha]; % liftoff point

modele1 = LC_in_square('varOn', true, 'xinit', yinit, 'vinit', [-1 0], 'tmax', T0, 'alpha', alpha);
modele1.solve

modele2 = LC_in_square('varOn', true, 'xinit', yinit, 'vinit', [ 0 1], 'tmax', T0, 'alpha', alpha);
modele2.solve

%% 
figure
set(gcf,'Position',[50 800 800 420])
subplot(1,2,1)
plot(modele1.t, modele1.yext(:,3:4),'linewidth',2)
xlim([0 T0])
ylim([-1.1 1.1])
xlabel('$\rm time (ms)$','interpreter','latex','fontsize',25)
ylabel('$\textbf{u}$','interpreter','latex','fontsize',25,'rot',0)
legend('x-direction','y-direction')
set(gca,'FontSize',18)

subplot(1,2,2)
plot(modele2.t, modele2.yext(:,3:4),'linewidth',2)
xlim([0 T0])
ylim([-1.1 1.1])
xlabel('$\rm time (ms)$','interpreter','latex','fontsize',25)
ylabel('$\textbf{u}$','interpreter','latex','fontsize',25,'rot',0)
legend('x-direction','y-direction')
set(gca,'FontSize',18)