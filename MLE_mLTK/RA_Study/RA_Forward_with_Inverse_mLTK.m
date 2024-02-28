%Code for the forward mLTK model using the results of the inverse problem run 
%in a single point

clear all
close all

%Load original insilico data obtained with real VIF for comparison
load('RA_mLTK_Persistent.mat')

%Initial conditions
C_0=0;  
Init=C_0;
time=t_for;

%Load the results of the inverse problem
load('RA_Inverse_Persistent.mat')
ktr=x(1);
ve=x(2);
vp=x(3);
leakage=x(4);

%Parameter vector from the optimization
par=[ktr,ve,vp,leakage];

%Forward mLTK model
[t,y] = ode45(@(t,y) mLTK_ODE(t,y,par,VIF,time),time,Init,[]);
C=y+vp*VIF;

%Plot of the forward model with the optimized parameter values and the
%original insilico data
plot(t,C,'k','LineWidth',0.5)
hold on
plot(t_for,C_for,'r')

%MLE evaluation
Sbar=1/(size(C_for,1))*sum(C_for);
MLE=sum ((C_for-C).^2);
Stot=sum((C_for-Sbar).^2);
R2=1-MLE/Stot;

%Saving option
%save('RA_Inverse_Persistent_Sol.mat','C','t','x')


