%Code to generate data for the RR-GBM case from the originial data with a 
% smoothing process at different levels. We consider the persistent type
% curve

clear all
close all

%Load original data with no smoth
load('Data_PersistentCOH.mat','CA_mod','time','VIF_mod')
xdata=time;
ydata=CA_XY_time;

%Smoothing 
B = smoothdata(ydata,"SamplePoints",xdata,'SmoothingFactor',.05);
CA_XY_time_Smo=B;

%Saving option
%save('RR_GBM_Persistent_Smoothdata_05.mat','CA_XY_time_Smo','xdata')


