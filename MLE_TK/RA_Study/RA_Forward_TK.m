%Code for the forward TK model to generate insilico data with real VIF
clear all
close all

ktr=0.0015;
ve=0.3;

%Parameter vector
par=[ktr,ve];

%Initial conditions
C_0=0;  
Init=C_0;

%Load real VIF estimated from the data
load('VIF.mat')
time=time_VIF;
VIF=VIF';

%Forward TK model
[t,y] = ode45(@(t,y) TK_ODE(t,y,par,VIF,time),time,Init,[]);
C_for=y;
t_for=t;

%Saving data option 
%save ('RA_TK_Persistent.mat','t_for','C_for','par','VIF')