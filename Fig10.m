%% Under nonuniform perturbation.
% Plot the norm of the actual displacements and the approximated displacements (eps*iSRC) vs eps

epsvec = (0.00006:0.002:0.1);     % perturbation size

% piecewise uniform rescaling 
[normNumDisp_pw,normu_pw,normdiff_pw] = SRC_nonuniform_perturbation_piecewise_nu_compute_displacement(epsvec);
[normNumDisp_same,normu_same,normdiff_same] = SRC_nonuniform_perturbation_same_nu_compute_displacement(epsvec);

%%
figure
% one nu
plot(epsvec,normNumDisp_same,'r','linewidth',2)
hold on
plot(epsvec,normu_same,'r:','linewidth',2)

% piecewise nu
plot(epsvec,normNumDisp_pw,'b','linewidth',2)
plot(epsvec,normu_pw,'b:','linewidth',2)


xlabel('$\varepsilon$','interpreter','latex','fontsize',30)
ylabel('$\rm norm(displacement)$','interpreter','latex','fontsize',30)
set(gca,'FontSize',18)
axis([0 0.1 0 6])
legend('uniform rescaling - actual','uniform rescaling - approximation','piecewise uniform rescaling - actual','piecewise uniform rescaling - approximation')

% Plot the relative difference between the norms vs eps
figure
% one nu
plot(epsvec,normdiff_same,'r','linewidth',2)
hold on
% piecewise nu
plot(epsvec,normdiff_pw,'b','linewidth',2)

xlabel('$\varepsilon$','interpreter','latex','fontsize',30)
ylabel('relative difference in norm','interpreter','latex','fontsize',30)
set(gca,'FontSize',18)
axis([0 0.1 -0.01 0.07])
legend('uniform rescaling','piecewise uniform rescaling');

