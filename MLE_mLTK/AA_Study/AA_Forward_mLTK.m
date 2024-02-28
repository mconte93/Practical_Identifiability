%Code for the forward mLTK model to generate insilico data with insilico VIF

clear all
close all

ktr=0.0015;
ve=0.1;
vp=0.2;
lambda=0.01;

%Parameter vector
par=[ktr,ve,vp,lambda];

%Initial conditions
C_0=0;  
Init=C_0;

%Load VIF produced with the code 'ArtificialVIF.m'
load('ArtificialVIF.mat')
time=t;

%Forward mLTK model
[t,y] = ode45(@(t,y) mLTKM_ODE(t,y,par,VIF,time),time,Init,[]);
C_for=y+vp*VIF;
t_for=t;

%Saving data option 
%save ('AA_mLTK_Persistent.mat','t_for','C_for','par','VIF')

