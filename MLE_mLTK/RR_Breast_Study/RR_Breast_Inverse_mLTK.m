%Code to run the inverse problem on the mLTK model with real GBM data and
%real VIF to obtain the optimal estimation of the model parameters

clear all
close all

%% Visualizations for the three types of curves related to the voxels of 
% coordinates (X1,Y1,Z1), (X2,Y2,Z2), and (X3,Y3,Z3).

global X1 Y1 Z1 X2 Y2 Z2 X3 Y3 Z3 

load('QIN-BREAST02-01.mat')

X1=80;
Y1=112;    
Z1=5;    
time=[1:6.69:6.69*55];

CA_XY_time=squeeze(CONCS(:,:,Z1,:));

for i=1:size(CA_XY_time,1)
    for j=1:size(CA_XY_time,2)
        for k=1:size(CA_XY_time,3)
            if isnan(CA_XY_time(i,j,k))==1 
            CA_XY_time(i,j,k)=0;
            end
        end
    end
end

subplot(3,2,1)
imagesc(squeeze(DCE(:,:,Z1,end)))
colormap gray
hold on
plot3(X1,Y1,10,'r*')
hold off
title('Slice ',num2str(Z1),'interpreter','latex')
axis square

subplot(3,2,2)
plot(time,squeeze(CA_XY_time(X1,Y1,:)),'k')
title('Type I','interpreter','latex')
xlabel('Time [s]','interpreter','latex')
ylabel('CA concentration','interpreter','latex')
axis square


X2=78;
Y2=118;   
Z2=5;    

CA_XY_time=squeeze(CONCS(:,:,Z2,:));

for i=1:size(CA_XY_time,1)
    for j=1:size(CA_XY_time,2)
        for k=1:size(CA_XY_time,3)
            if isnan(CA_XY_time(i,j,k))==1 
            CA_XY_time(i,j,k)=0;
            end
        end
    end
end

subplot(3,2,3)
imagesc(squeeze(DCE(:,:,Z2,end)))
hold on
plot3(X2,Y2,10,'r*')
hold off
title('Slice ',num2str(Z2),'interpreter','latex')
axis square

subplot(3,2,4)
plot(time,squeeze(CA_XY_time(X2,Y2,:)),'k')
title('Type II','interpreter','latex')
xlabel('Time [s]','interpreter','latex')
ylabel('CA concentration','interpreter','latex')
axis square


X3=81; 
Y3=117; 
Z3=5;   

CA_XY_time=squeeze(CONCS(:,:,Z3,:));

for i=1:size(CA_XY_time,1)
    for j=1:size(CA_XY_time,2)
        for k=1:size(CA_XY_time,3)
            if isnan(CA_XY_time(i,j,k))==1 
            CA_XY_time(i,j,k)=0;
            end
        end
    end
end

subplot(3,2,5)
imagesc(squeeze(DCE(:,:,Z3,end)))
hold on
plot3(X3,Y3,10,'r*')
title('Slice ',num2str(Z3),'interpreter','latex')
hold off
axis square

subplot(3,2,6)
plot(time,squeeze(CA_XY_time(X3,Y3,:)),'k')
title('Type III','interpreter','latex')
xlabel('Time [s]','interpreter','latex')
ylabel('CA concentration','interpreter','latex')
axis square  

%% Inverse problem in each of the three indicated voxel

global X1 Y1 Z1 X2 Y2 Z2 X3 Y3 Z3 

ktr_0=.005;
ve_0=.1;
vp_0=.1;
leakage_0=.005;

%Optimization of the parameter values with PS algorithm
x0=[ktr_0,ve_0,vp_0,leakage_0];
fun=@fun_RR_mLTK_PS;
lb = [0 1.e-3 1.e-3 0];
ub = [1 1 1 1];

options=optimoptions('particleswarm','SwarmSize',200,'HybridFcn',@fmincon,'MaxIterations',500,'ObjectiveLimit',1e-16,'InitialSwarmMatrix',x0,'PlotFcn',@pswplotbestf);
nvar=4;
[x,fval,exitflag,output] = particleswarm(fun,nvar,lb,ub,options);
formatstring = 'particleswarm reached the value %f using %d function evaluations.\n';
fprintf(formatstring,fval,output.funccount)

%Saving option
%save('RR_Breast_Inverse_Persistent.mat','x','X','Y','Z')

function [MLE] = fun_RR_mLTK_PS(x)

global X1 Y1 Z1 X2 Y2 Z2 X3 Y3 Z3 

    X=X1;
    Y=Y1;
    Z=Z1;
    
    %Load the data 
    load('QIN-BREAST02-01.mat')
    time=[1:6.69:6.69*55];
    CA_XY_time=squeeze(CONCS(X,Y,Z,:));
    for k=1:size(CA_XY_time,1)
       if isnan(CA_XY_time(k))==1 
          CA_XY_time(k)=0;
       end
    end
    
    xdata=time;
    ydata=CA_XY_time';
    
    %Model parameter
    ktr=x(1);
    ve=x(2);
    vp=x(3);
    leakage=x(4);
       
    par=[ktr,ve,vp,leakage];
    
    %Initial conditions
    C_0=0;   
    time=xdata;
    Init=C_0;
    
    %Forward mLTK model
    [t,y] = ode45(@(t,y) mLTK_ODE(t,y,par,VIF,time),time,Init,[]);
    C_approx=y+vp*VIF;
    
    %MLE evaluation
    MLE=0;
    for i=1:size(C_approx,1)
        res=ydata(i)-C_approx(i);
        if ~isnan(res)
           MLE=MLE+res.^2;
        end
    end
    
        
end