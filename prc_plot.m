% generate iPRC (Figure 7) and compute the global relative change in frequency nu1

alpha=0.2;
T0=6.766182958128617;  % intrinsic period of the oscillator
xinit=[1,0.200000001]; % initial condition for the trajectory
model = LC_in_square(false, xinit);
model.solve;
M=find_prc_monodromy(model);  % find the Monodromy matrix for the iPRC
disp('The eigenvalues of the Monodromy matrix are')
disp(eig(M));                 

[V, D]=eig(M');              
ind_e1=(abs(diag(D))==max(abs(diag(D)))); %find index where the eigenvalue is 1
% if ind_e1(1)==1, zinit=V(:,1);else zinit=V(:,2);end
zinit=V(:,ind_e1);                        %find the eigenvector associated with eigenvalue 1
dummy=0;
f10=model.LC_ODE(dummy,model.yext(end,1:2)',model.checkdomain(model.yext(end,1:2))); % compute the initial vector field

rescale=f10'*zinit;   % compute the scaling factor
z0=zinit/(rescale);   % compute the initial value for iPRC
disp('The initial condtion for iPRC is')
disp(z0);

model.find_prc(z0)    % Solve the adjoint equation with initial condition for the iPRC
model.plot_prc        % Plot iPRC

%% find the relative change in frequency (nu) for computing shape response curve

% if the perturbation is always present 
int=-trapz(model.prct,wrev(model.yext(:,1)).*model.prc(:,1)+wrev(model.yext(:,2)).*model.prc(:,2));% time is backward,so the integral has the opposite sign 
T1=-int; 
nu1=T1/T0;

% if the perturbation is only present over part of the region, say, over the wedge 
ind_above_wedge=(model.yext(:,1) + model.yext(:,2) >=0) & (model.yext(:,2) - model.yext(:,1) >=0);
time_above_wedge=model.t(ind_above_wedge);
x_above_wedge=model.yext(ind_above_wedge,:);

iprc_time_above_wedge=model.prct(flipud(ind_above_wedge));
iprc_above_wedge=model.prc(flipud(ind_above_wedge),:);

int_above_wedge=-trapz(iprc_time_above_wedge,wrev(x_above_wedge(:,1)).*iprc_above_wedge(:,1)+wrev(x_above_wedge(:,2)).*iprc_above_wedge(:,2));% time is backward,so the integral has the opposite sign 

T1_above_wedge=-int_above_wedge; 
nu1_above_wedge=T1_above_wedge/T0;


