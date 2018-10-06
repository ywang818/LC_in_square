% Generate iPRC by integrating the adjoint equation backward in time
% and compute the global linear change in period T1

T0=6.766182958128617;  % intrinsic period of the oscillator

% Fix an initial condition x0 on the LC, compute the LC over [0, T0]
% Notice that the initial condition is made slightly off the liftoff point 
% to trigger the exit-wall event in backwards time, at which iPRC will be updated
xinit=[1,0.2000000001];                        
model = LC_in_square('xinit', xinit);
model.solve;

% Find the Monodromy (fundamental) matrix at one period for the iPRC,
% which should have a single eigenvector with unit eigenvalue 
M=find_prc_monodromy(model);  
disp('The eigenvalues of the Monodromy matrix are')
disp(eig(M));                 

[V, D]=eig(M');              
ind_e1=(abs(diag(D))==max(abs(diag(D)))); %find index where the eigenvalue is 1
zinit=V(:,ind_e1);                        %find the eigenvector associated with eigenvalue 1
dummy=0;
f10=model.LC_ODE(dummy,model.yext(end,1:2)',model.checkdomain(model.yext(end,1:2))); % compute the initial vector field
rescale=f10'*zinit;   % compute the scaling factor

% compute the initial value z0 for iPRC
z0=zinit/(rescale);   
disp('The initial condtion for iPRC is')
disp(z0);

% Solve the adjoint equation with the initial condition z0 for the complete iPRC over [0 T0]
model.find_prc(z0)    
model.plot_prc        % Plot iPRC

%% find the relative change in time T1 in response to the static perturbation (alpha, omega) -> (alpha + eps, omega - eps)
% This is needed for computing the iSRC 

% if the perturbation is present over the full cycle, T1 = int_0^T0 (z(s) * dF/deps) ds
% Time for iPRC is reversed compared with the time for LC,so the integral has a negative sign
int=-trapz(model.prct,(wrev(model.yext(:,1))+wrev(model.yext(:,2))).*model.prc(:,1)+(-wrev(model.yext(:,1))+wrev(model.yext(:,2))).*model.prc(:,2)); 
T1=-int; 
nu1=T1/T0;

% if the perturbation is only present in the region above the wedge, T1 = int_{time above the wedge} (z(s) * dF/deps) ds
ind_above_wedge=(model.yext(:,1) + model.yext(:,2) >=0) & (model.yext(:,2) - model.yext(:,1) >=0);
time_above_wedge=model.t(ind_above_wedge);
x_above_wedge=model.yext(ind_above_wedge,:);

iprc_time_above_wedge=model.prct(flipud(ind_above_wedge));
iprc_above_wedge=model.prc(flipud(ind_above_wedge),:);

int_above_wedge=-trapz(iprc_time_above_wedge,(wrev(x_above_wedge(:,1))+wrev(x_above_wedge(:,2))).*iprc_above_wedge(:,1)+(-wrev(x_above_wedge(:,1))+wrev(x_above_wedge(:,2))).*iprc_above_wedge(:,2));% time is backward,so the integral has the opposite sign 

T1_above_wedge=-int_above_wedge; 
nu1_above_wedge=T1_above_wedge/T0;


