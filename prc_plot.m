% generate iPRC (Figure 7) using the initial cond calculated from find_INI_prc.m 

alpha=0.2;
T0=6.766182958128617;  % intrinsic period of the oscillator
xinit=[1,0.200000001]; % initial condition for the trajectory
model = LC_in_square(false, xinit);
model.solve;
M=find_prc_INI(model);
disp('The eigenvalues of the Monodromy matrix are')
disp(eig(M));

[V, D]=eig(M');
ind_e1=(abs(diag(D))==max(abs(diag(D))));%find index where the eigenvalue is 1
% if ind_e1(1)==1, zinit=V(:,1);else zinit=V(:,2);end
zinit=V(:,ind_e1);
dummy=0;
f10=model.LC_ODE(dummy,model.yext(end,1:2)',model.checkdomain(model.yext(end,1:2)));

rescale=f10'*zinit;
z0=zinit/(rescale);
disp('The initial condtion for iPRC is')
disp(z0);

model.find_prc(z0)
model.plot_prc

%% find the relative change in frequency (nu) for shape response curve

% if the perturbation is present throughout the whole period

int=-trapz(model.prct,wrev(model.yext(:,1)).*model.prc(:,1)+wrev(model.yext(:,2)).*model.prc(:,2));% time is backward,so the integral has the opposite sign 
T1=-int; 
nu1=T1/T0;

% if the perturbation is only present over the triangle, 
% T1 is therefore only the integral over the triangle
ind_above_triangle=(model.yext(:,1) + model.yext(:,2) >=0) & (model.yext(:,2) - model.yext(:,1) >=0);
time_above_triangle=model.t(ind_above_triangle);
x_above_triangle=model.yext(ind_above_triangle,:);

iprc_time_above_triangle=model.prct(flipud(ind_above_triangle));
iprc_above_triangle=model.prc(flipud(ind_above_triangle),:);

int_pert=-trapz(iprc_time_above_triangle,wrev(x_above_triangle(:,1)).*iprc_above_triangle(:,1)+wrev(x_above_triangle(:,2)).*iprc_above_triangle(:,2));% time is backward,so the integral has the opposite sign 

T1_pert=-int_pert; 
nu_pert=T1_pert/T0;


