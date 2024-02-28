%Code for the forward LTK model using the results of the inverse problem run 
%in a single point

clear all
close all

%Load original data and VIF for comparison
load('RR_GBM_Persistent_Smoothdata_15.mat')
load('Data_PersistentCOH.mat','VIF_mod')
VIF=VIF_mod;
ydata=CA_XY_time_Smo;
time=xdata;

Tmax=time(end);  
t_VIF=time;

%Load the results of the inverse problem
load('RR_GBM_Inverse_Persistent_Smooth15.mat')
ktr=x(1);
ve=x(2);
vp=x(3);
leakage=x(4);

%Parameter vector from the optimization
par=[ktr,ve,vp,leakage];

%Initial conditions;
Ce_0=0; 
Cl_0=0;
Init=[Ce_0;Cl_0];

%Forward LTK model
[t,y] = ode45(@(t,y) LTK_ODE(t,y,par,VIF,t_VIF),t_VIF,Init,[]);
C_for=y(:,1)+y(:,2)+vp*VIF';
t_for=t;

%Plot of the forward model with the optimized parameter values and the
%original data
plot(t_VIF,CA_XY_time_Smo,'r',t_for,C_for,'k')

%MLE evaluation
Sbar=1/(size(CA_XY_time_Smo,1))*sum(CA_XY_time_Smo);
MLE=sum ((CA_XY_time_Smo'-C_for).^2);
Stot=sum((CA_XY_time_Smo-Sbar).^2); 
R2=1-MLE/Stot;

%Saving option
%save('RR_GBM_Forward_with_inverse_Persistent_Smooth15.mat','CA_XY_time_Smo','t_VIF','t_for','C_for')
