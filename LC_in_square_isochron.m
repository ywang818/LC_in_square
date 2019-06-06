% Compute isochrons (level curves of phase) for LC_in_square

% For a grid of initial locations in [-1 1 -1 1], find the time (modulo the period T0)
% to reach the liftoff point on the east wall. This is equivalent to the
% asymptotic phase function. Isochron (the level curve of phase function)
% is computed using pcolor and contour. 

T0=6.766182958128606; % intrinsic oscillator period

% Take 600*600 many initial locations
k=600;  
xgrid=linspace(-0.99999999,0.99999999,k);
ygrid=linspace(-0.99999999,0.99999999,k);
[xmesh,ymesh]=meshgrid(xgrid,ygrid);
isochron=nan(size(xmesh));

% For each initial location, find the time (mod T0) to reach the liftoff
% point (1, alpha) (when 'isochronOn' is true, integration in LC_in_squre will stop when
% the liftoff point is reached.)
for i=1:k
    for j=1:k
        x0=xmesh(i,j); y0=ymesh(i,j);
        model=LC_in_square('isochronOn', true, 'xinit',[x0,y0], 'tmax',50);
        model.solve;
        isochron(i,j)=mod(model.t(end),T0);
    end
end
% Save data. 
save isochron_full600_20190505.mat 

%% find trajectories
% Find the LC trajectory
alpha=0.2;
model_lc = LC_in_square('xinit',[1,alpha],'tmax',T0);
model_lc.solve;

% Find the trajectory for osculating trajectory (below is its analytic expression)
t=(0:-0.01:-20);
xgamma_1=sqrt(alpha^2+1).*exp(alpha*t).*cos(atan(alpha)+t);
ygamma_1=sqrt(alpha^2+1).*exp(alpha*t).*sin(atan(alpha)+t);

%% plot isochron curves (gray) together with limit cycle (black solid)
% and the oscullating curve originated from the liftoff point (black dashed)
 
load isochron_full600_20190505.mat  %load data generated from above

figure
contour(xmesh,ymesh,isochron,50,'linewidth',2)
colormap(hsv)
hold on
plot(model_lc.yext(:,1), model_lc.yext(:,2),'k','linewidth',2)
plot(xgamma_1,ygamma_1,'k-.','linewidth',2)
axis([-1.1 1.1 -1.1 1.1])
axis square
% set(gca,'Position',[.06 .14 .85 .80]) % x,y, width, length
set(gca,'FontSize',18)
xlabel('$x$','interpreter','latex','fontsize',25)
ylabel('$y$','interpreter','latex','fontsize',25,'rot',0)
plot(1,0.2,'kp','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',25)
text(1.1,0.2,'$\theta=0$','Interpreter','latex','FontSize',18,'Color','k')
box off

% another way to plot:
% figure
% pcolor(xmesh,ymesh,isochron)
% shading flat
% colormap(hsv)
% hold on
% contour(xmesh,ymesh,isochron,50,'linewidth',2,'Color',[0.45 0.45 0.45]);
% plot(model_lc.yext(:,1), model_lc.yext(:,2),'k','linewidth',2)
% plot(xgamma_1,ygamma_1,'k-.','linewidth',2)
% axis([-1.1 1.1 -1.1 1.1])
% axis square
% % set(gca,'Position',[.06 .14 .85 .80]) % x,y, width, length
% set(gca,'FontSize',18)
% xlabel('$x$','interpreter','latex','fontsize',25)
% ylabel('$y$','interpreter','latex','fontsize',25,'rot',0)
% plot(1,0.2,'kp','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',25)
% text(1.1,0.2,'$\theta=0$','Interpreter','latex','FontSize',18,'Color','k')
% box off
% % grid on
% 
% figure
% pcolor(xmesh,ymesh,isochron)
% shading flat
% colormap(hsv)
% hold on
% contour(xmesh,ymesh,isochron,50,'linewidth',2,'Color',[0.45 0.45 0.45]);
% plot(model_lc.yext(:,1), model_lc.yext(:,2),'k','linewidth',2)
% plot(xgamma_1,ygamma_1,'k-.','linewidth',2)
% axis([0 1.1 -1.1 0.3])
% % set(gca,'Position',[.06 .14 .85 .80]) % x,y, width, length
% set(gca,'FontSize',18)
% xlabel('$x$','interpreter','latex','fontsize',25)
% ylabel('$y$','interpreter','latex','fontsize',25,'rot',0)
% plot(1,0.2,'kp','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',25)
% text(1.05,0.2,'$\theta=0$','Interpreter','latex','FontSize',18,'Color','k')
% box off