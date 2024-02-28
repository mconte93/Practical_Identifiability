%Code for the forward LTK model to generate insilico data with insilico VIF

clear all
close all

ktr=0.0025;
ve=0.0;
vp=0.05;
lambda=0.001;

%Parameter vector
par=[ktr,ve,vp,lambda];

%Initial conditions
Ce_0=0;
Cl_0=0;
Init=[Ce_0;Cl_0];

%Load VIF produced with the code 'ArtificialVIF.m'
load('ArtificialVIF.mat')
time=time_VIF;

%Forward LTK model
[t,y] = ode45(@(t,y) LTK_ODE(t,y,par,VIF,time),time,Init,[]);
C_for=y(:,1)+y(:,2)+vp*VIF';
t_for=t;

%Saving data option 
%save ('AA_LTK_Persistent.mat','t_for','C_for','par','VIF')

