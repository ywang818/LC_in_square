% Fig4.m: plot fig4

% Find the LC trajectory
alpha=0.2;
model_lc = LC_in_square('xinit',[1,alpha],'tmax',T0);
model_lc.solve;

% Find the osculating trajectory (below is its analytic expression)
t=(0:-0.01:-20);
xgamma_1=sqrt(alpha^2+1).*exp(alpha*t).*cos(atan(alpha)+t);
ygamma_1=sqrt(alpha^2+1).*exp(alpha*t).*sin(atan(alpha)+t);

% Compute vector field
[x,y]=meshgrid(-0.9:0.2:0.9,-0.9:0.2:0.9);
dxdt=alpha.*x-y; dydt=x+alpha.*y;

% plot trajectory
model_lc.plot
subplot(1,2,2)
hold on
quiver(x,y,dxdt,dydt,'r')
plot(xgamma_1,ygamma_1,'k-.','linewidth',2)
plot(1,0.2,'kp','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',25)
axis square
% set(gca,'Position',[.06 .14 .85 .80]) % x,y, width, length
set(gca,'FontSize',18)
xlabel('$x$','interpreter','latex','fontsize',25)
ylabel('$y$','interpreter','latex','fontsize',25,'rot',0)
plot(1,0.2,'kp','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',25)



