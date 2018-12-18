% Compute isochrons (level curves of phase) for LC_in_square

T0=6.766182958186305; % intrinsic oscillator period

k=300;
xgrid=linspace(-0.99999999,0.99999999,k);
ygrid=linspace(-0.99999999,0.99999999,k);
[xmesh,ymesh]=meshgrid(xgrid,ygrid);
isochron=zeros(k.^2,3);


r=1;
for i=1:k
    for j=1:k
        x0=xmesh(i,j); y0=ymesh(i,j);
        [t,y]=LC_in_square('xinit', [x0,y0],'tmax',50);
        isochron(r,1:2)=[x0,y0];
        isochron(r,3)=mod(t(end),T);
        r=r+1;
    end
end
save ini_20170408.mat 

% to use,
% load ini_20170408.mat

% xlabel('x coord')
% ylabel('y coord')
% grid on
% print -dpdf LC_in_square_plot1.pdf

