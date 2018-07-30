% Generate Figure 5: solution to the homogeneous variatinal equation

S=[];alpha=0.2;
T=6.766182958128606; % intrinsic period of the oscillator
varOn = true;
yinit=[1, alpha];

modele1 = LC_in_square(varOn,yinit,[-1 0],T,alpha);
modele1.solve
modele1.plot

modele2 = LC_in_square(varOn,yinit,[0 1],T,alpha);
modele2.solve
modele2.plot