%Code for the forward LTK model using the results of the inverse problem run 
%in a single point

clear all
close all

%Load original data and VIF for comparison
load('QIN-BREAST02-01.mat')
CA_XY_time=squeeze(CONCS(X,Y,Z,:));
for k=1:size(CA_XY_time,1)
   if isnan(CA_XY_time(k))==1 
      CA_XY_time(k)=0;
   end
end
time=[1:6.69:6.69*55]; 
t_VIF=time;

%Load the results of the inverse problem
load('RR_Breast_Inverse_Persistent.mat')
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
C_for=y(:,1)+y(:,2)+vp*VIF;
t_for=t;

%Plot of the forward model with the optimized parameter values and the
%original data
plot(t_VIF,CA_XY_time,'r',t_for,C_for,'k')

%MLE evaluation 
Sbar=1/(size(CA_XY_time,1))*sum(CA_XY_time);
MLE=sum ((CA_XY_time-C_for).^2);
Stot=sum((CA_XY_time-Sbar).^2);  
R2=1-MLE/Stot;

%Saving option
%save('RR_Breast_Forward_with_inverse_Persistent.mat','CA_XY_time','t_VIF','t_for','C_for')
