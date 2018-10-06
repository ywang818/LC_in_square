% Solve the nonhomogeneous variational equation for the iSRC using uniform rescaling
% under a nonuniform static perturbation that is only present above the wedge
%        (alpha, omega) -> (alpha + eps, omega - eps)
% !! Need to run prc_plot to find T1, linear shift in period

%  Region I: above wedge, Region II: below wedge

%  For eps over (0.00006:0.01:0.1), compute the norms of the actual (LC_eps - LC) and
%  approximated (eps*iSRC) displacements between perturbed and unperturbed LC solutions

T0=6.766182958186305; % intrinsic oscillator period
T1=2.694391001334606; % obtained from prc_plot % the relative linear shift in period, computed from prc_plot
nu=T1/T0; %  the uniform scaling.

xinit=[1,0];  % initial condition for both perturbed and unperturbed limit cycle solution

% see 'SRC_nonuniform_perturbation_piecewise_nu_plot' for how to find the following values
x_in=[0.811100985386121   0.811100985460158]; % coordinate of the unperturbed entry point into the region 1
T0_above=1.691545739499871; % unperturbed total time spent in region I

epsvec = (0.00006:0.01:0.1);
normu = zeros(size(epsvec));       % norm of approximated displacements
normNumDisp = zeros(size(epsvec)); % norm of actual displacements
normdiff = zeros(size(epsvec));    % relative difference between the actual norm and approximated norm

% Initial condition for iSRC at the entry point into region I,
% see 'SRC_nonuniform_perturbation_piecewise_nu_plot' for how to find it
vinit = [3.074e-10 -4.7892e-10];

for i=1:length(epsvec)
    
    eps=epsvec(i);   % perturbation size
    x_in_pert = vinit*eps + x_in; % approximated initial condition for the perturbed LC
    
    % Find the period Teps of the LC under perturbation
    model_pert = LC_in_square('xinit', xinit, 'vinit',[0 0],...
        'tmax', 20*T0, 'nu', [0,0], 'eps', eps);
    model_pert.solve;
    Teps=model_pert.findPeriod;   % find perturbed period
            
    % Take x_in_pert as initial cond, compute the perturbed solution on [0 Teps]
    model_pert = LC_in_square('xinit', x_in_pert, 'vinit',[0 0],...
        'tmax', Teps, 'nu', [0,0], 'eps', eps);
    model_pert.solve;
    
    % Compute the unperturbed solution and the iSRC with uniform nu
    model = LC_in_square('varOn', true, 'xinit', x_in, 'vinit', vinit, ...
        'tmax', T0, 'nu', [nu,nu]);
    model.solve
    
    % Rescale the time of the perturbed LC to be the same as model.t, the time for the unperturbed LC
    [tspan1, Ind1] = unique((model_pert.t./Teps).*T0,'stable'); % rescale time such that it has time [0 T0]
    x_unique = model_pert.yext(Ind1,1:2); % obtain x,y values after time is rescaled
    x_pert_interp = interp1(tspan1,x_unique,model.t); % interpolate using unperturbed solution
    
    % actual displacement obtained from subtracting unperturbed LC from perturbed LC
    displacement=x_pert_interp - model.yext(:,1:2);
    
    % evaluate the norm of the actual displacement vector
    displacement(isnan(displacement))=0;    
    normNumDisp(i)=norm(displacement);
    
    % evaluate the norm of the approximated displacement
    normu(i)=norm(eps*model.yext(:,3:4));
    
    % compute the relative difference between the actual norm and approximated norm
    normdiff(i)=(normNumDisp(i)-normu(i))/normNumDisp(i);     
end

%% Plot the norm of the actual displacements and the approximated displacements (eps*iSRC) vs eps
figure
plot(epsvec,normNumDisp,'r','linewidth',2)
hold on
plot(epsvec,normu,'r:','linewidth',2)
%plot(log(epsvec),log(normdiff),'linewidth',2)
xlabel('$\varepsilon$','interpreter','latex','fontsize',30)
ylabel('$\rm norm(displacement)$','interpreter','latex','fontsize',30)
set(gca,'FontSize',18)

% Plot the relative difference between the norms vs eps
figure
plot(epsvec,normdiff,'r','linewidth',2)
% hold on
% plot(epsvec,normdiff1,'r','linewidth',2)
xlabel('$\varepsilon$','interpreter','latex','fontsize',30)
ylabel('$\rm relative difference$','interpreter','latex','fontsize',30)
set(gca,'FontSize',18)


