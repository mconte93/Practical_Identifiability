%Code for the forward mLTK model to generate insilico data with real VIF
clear all
close all

ktr=0.008;
ve=0.001;
vp=0.9;
lambda=0;

%Parameter vector
par=[ktr,ve,vp,lambda];

%Initial conditions
C_0=0;  
Init=C_0;

%Load real VIF estimated from the data
load('VIF.mat')
time=time_VIF;

%Forward mLTK model
[t,y] = ode45(@(t,y) mLTK_ODE(t,y,par,VIF,time),time,Init,[]);
C_for=y+vp*VIF;
t_for=t;

%Saving data option 
%save ('RA_mLTK_Persistent.mat','t_for','C_for','par','VIF')
