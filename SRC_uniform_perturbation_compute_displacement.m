% Solve the nonhomogeneous variational equation for the iSRC 
% under a uniform static perturbation: (alpha, omega) -> (alpha + eps, omega - eps)
% !! Need to run prc_plot.m first to find nu1

% For eps over (0:0.001:0.01), compute the norms of the actual (LC_eps - LC) and
% approximated (eps*iSRC) displacements between perturbed and unperturbed LC solutions

T0=6.766182958128617;  % intrinsic period of the oscillator
T1=10.802054306767772; % obtained from prc_plot % the relative linear shift in period, computed from prc_plot
nu1=T1/T0;

alpha = 0.2; omega = 1;   % default value for alpha, omega
yinit=[1, alpha/omega];   % liftoff point at x=1 for unperturbed trajectory
tmax = 80;

% Find an initial condition for the iSRC when eps = 0.0001
eps=0.0001;
alpha_pert = alpha + eps;
omega_pert = omega - eps;
yinit_pert=[1, alpha_pert/omega_pert]; % liftoff point for perturbed trajectory
vinit = (yinit_pert-yinit)/eps;        % initial value for SRC

epsvec = (0:0.001:0.01);
normdiff = zeros(size(epsvec));  % norm of actual displacements 
normu = zeros(size(epsvec));     % norm of approximated displacements 

for i=1:length(epsvec)
eps = epsvec(i); % static perturbation size
yinit_pert = yinit + eps*vinit; % approximated initial condition for the perturbed LC

% find the unperturbed solution and the shape response curve
model = LC_in_square('varOn', true, 'xinit', yinit, 'vinit', vinit, 'tmax', T0, 'nu', nu1);
model.solve;

% Find the period Teps of the LC under perturbation 
model_pert = LC_in_square('xinit', yinit_pert, 'vinit', [0 0], 'tmax', tmax, 'alpha', alpha + eps, 'omega', omega - eps);
model_pert.solve;
Teps=model_pert.findPeriod;

% find perturbed solutions over [0 Teps]
model_pert = LC_in_square('xinit', yinit_pert, 'vinit', [0 0], 'tmax', Teps, 'alpha', alpha + eps, 'omega', omega - eps);
model_pert.solve;

% Rescale the time of the perturbed LC to be the same as model.t, the time for the unperturbed LC
[tspan1, Ind1] = unique((model_pert.t./Teps).*T0,'stable'); % rescale time such that it has time [0 T0]
x_unique = model_pert.yext(Ind1,1:2);                       % obtain x,y values after time is rescaled
x_pert_interp = interp1(tspan1,x_unique,model.t);           % interpolate using unperturbed solution

% actual displacement obtained from subtracting unperturbed LC from perturbed LC 
displacement=x_pert_interp - model.yext(:,1:2);
displacement(isnan(displacement))=0;

% evaluate the norm of the actual displacement vector
normdiff(i)=norm(displacement);

% evaluate the norm of the displacement vector estimated from eps*iSRC
normu(i)=norm(eps*model.yext(:,3:4));
end

%% Plot the norm of the actual displacements and the approximated displacements (eps*iSRC) vs eps
figure
set(gcf,'Position',[0 0 720 520])
plot(epsvec,normdiff,'k','linewidth',3)
hold on
plot(epsvec,normu,'r:','linewidth',3)
%plot(log(epsvec),log(normdiff),'linewidth',2)
xlim([0 0.01])
ylim([0 0.515])
xlabel('$\varepsilon$','interpreter','latex','fontsize',30)
ylabel('norm $\rm (\gamma_{\varepsilon}(t)-\gamma(t))$','interpreter','latex','fontsize',30)
legend({' actual',' approximation'},'Interpreter','latex','Location','northwest')
set(gca,'FontSize',18)


