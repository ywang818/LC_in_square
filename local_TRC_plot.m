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
% x_out=model.yext(ind_close_end,:); % the exit point from region 1
% x_in=model.yext(ind_close_ini,:); % the entry point into region 1

x_in=[0.811047481339979   0.811181234361245];  % the entry point into region 1
x_out=[-0.813550109151894   0.809463778128291]; % the exit point from region 1
model = LC_in_square('xinit', x_out);  % Compute the solution trajectory from x_out
model.solve;

% Compute the total time spent in the region I, that is above the wedge
ind_above_wedge=(model.yext(:,1) + model.yext(:,2) >=0) & (model.yext(:,2) - model.yext(:,1) >=0);
time_above_wedge=model.t(ind_above_wedge);
T0_above_wedge=time_above_wedge(end)-time_above_wedge(1); % total time spent in region I

% To only compute the lTRC in region I, compute the solution from x_in to x_out
model = LC_in_square('xinit', x_in, 'vinit', [0,0], 'tmax', T0_above_wedge);
model.solve;

% To compute the lTRC over the full cycle, compute the solution with IC at x_out on [0,T0]
% model = LC_in_square('xinit', x_out, [0,0],T0);
% model.solve;


% Compute the boundary value of lTRC at the exit point, lTRCinit,
% This will be the initial condition for the lTRC since we will integrate the adjoint equation backward in time
lTRCinit0=[1,1]'; % the normal vector to the exit boundary of the wedge (y=-x)
dummy=0;
f10=model.LC_ODE(dummy,x_out',model.checkdomain(x_out)); % vector field at the exit point
rescale=f10'*lTRCinit0;
lTRCinit=-lTRCinit0'/(rescale); % the value for the lTRC at the exit point 

model.find_prc(lTRCinit); 
% model.plot_prc;          % visualize the lTRC result

% Compute lTRC(xin)* (xin_pert-xin)/eps, the first term in T1_above_wedge
vinit = [3.074e-10 -4.7892e-10]; % (xin_pert-xin)/eps, the IC for iSRC at the entry point into region I, obtained from 'shape_response_curve_piecewise_nu_plot'
T1_above_wedge_1=model.prc(end,:)*vinit'; 

% integral of the lTRC over the region above the wedge,
% there is a minus sign in front of the integral since time is backward
int=-trapz(model.prct,(wrev(model.yext(:,1))+wrev(model.yext(:,2))).*model.prc(:,1)+(-wrev(model.yext(:,1))+wrev(model.yext(:,2))).*model.prc(:,2));% time is backward,so the integral has the opposite sign 

T1_above_wedge_int=int; 
T1_above_wedge = T1_above_wedge_1 + T1_above_wedge_int; % relative shift in time in region I
nu_above_wedge=T1_above_wedge/T0_above_wedge; % relative change in frequency in the region above the wedge

disp('nu_above_wedge is')
disp(nu_above_wedge)

%%
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
plot([T0_above_wedge T0_above_wedge], [-2 2],'b-.','linewidth',2)
plot([0 0], [-2 2],'g-.','linewidth',2)
text(0.03,-1.5,'$t_A$','Interpreter','latex','FontSize',30,'Color','g')
text(1.55,-1.5,'$t_B$','Interpreter','latex','FontSize',30,'Color','b')
