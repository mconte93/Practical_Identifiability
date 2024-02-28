%Code for the forward TK model to generate insilico data with insilico VIF

clear all
close all

ktr=0.007;
ve=0.2;

%Parameter vector
par=[ktr,ve];

%Initial condition
C_0=0;  
Init=C_0;

%Load VIF produced with the code 'ArtificialVIF.m'
load('ArtificialVIF.mat')
time=t;

%Forward TK model
[t,y] = ode45(@(t,y) TK_ODE(t,y,par,VIF,time),time,Init,[]);
C_for=y;
t_for=t;

%Saving data option 
%save ('AA_TK_Persistent.mat','t_for','C_for','par','VIF')