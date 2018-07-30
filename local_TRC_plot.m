% Consider two regions, above and below the wedge
% Generate local timing response curve for the region I that is above the wedge, usng the function find_prc in LC_in_square
% and compute the piecewise relative change in frequency nu_above_wedge and nu_below_wedge


alpha=0.2;
T0=6.766182958128617;

% xinit=[0.6547,1];
% model = LC_in_square(false, xinit);
% model.solve;
% ind_close_end=(abs(model.yext(:,1)+model.yext(:,2))<5e-3) & (model.yext(:,2)>=0);
% ind_close_ini=(abs(model.yext(:,1)-model.yext(:,2))<5e-3) & (model.yext(:,2)>=0);
% x_exit=model.yext(ind_close_end,:); % the exit point from region 1
% x_entry=model.yext(ind_close_ini,:); % the entry point into region 1

x_entry=[0.811047481339979   0.811181234361245];  % the entry point into region 1
x_exit=[-0.813550109151894   0.809463778128291]; % the exit point from region 1
model = LC_in_square(false, x_exit);  % Compute the solution trajectory  
model.solve;

% Using the solution computed from above to compute the total time spent in the region above the wedge
ind_above_wedge=(model.yext(:,1) + model.yext(:,2) >=0) & (model.yext(:,2) - model.yext(:,1) >=0);
time_above_wedge=model.t(ind_above_wedge);
T0_above_wedge=time_above_wedge(end)-time_above_wedge(1); % total time spent in region I

% Compute the solution trajectory that begins from the entry point and ends at the exit point
model = LC_in_square(false, x_entry,[0,0],T0_above_wedge);
model.solve;

% Compute the value of lTRC at the exit point, trcinit,
% This will be the initial condition for the lTRC since we will integrate the adjoint equation backward in time
zinit=[1,1]'; % the normal vector to the exit boundary of the wedge (y=-x)
dummy=0;
f10=model.LC_ODE(dummy,x_exit',model.checkdomain(x_exit)); % vector field at the exit point
rescale=f10'*zinit;
trcinit=-zinit'/(rescale); % the value for the lTRC at the exit point 

% compute the lTRC using the initial cond trcinit
model.find_prc(trcinit); 
% model.plot_prc;          % plot the result

% integral of the lTRC over the region above the wedge
int=-trapz(model.prct,wrev(model.yext(:,1)).*model.prc(:,1)+wrev(model.yext(:,2)).*model.prc(:,2));% time is backward,so the integral has the opposite sign 

T1_above_wedge=int; 
nu_above_wedge=T1_above_wedge/T0_above_wedge; % relative change in frequency in the region above the wedge

disp('nu_above_wedge is')
disp(nu_above_wedge)


figure
plot(model.prct,model.prc(:,1:2),'linewidth',2)
legend('Z_x','Z_y')
xlim([0 model.tmax])
xlabel('$\rm time (ms)$','interpreter','latex','fontsize',30)
legend('x-direction','y-direction')
title('Local timing response curve')
% grid on
set(gca,'FontSize',18)
hold on
plot([T0_above_wedge T0_above_wedge], [-2 2],'g-.','linewidth',2)
plot([0 0], [-2 2],'b-.','linewidth',2)
