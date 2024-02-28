%Code for obtaining the results used in the profile likelihood plot for
%persistent curve in the RR-Breast case

clear all
close all

%Load results of the inverse problem
load('RR_Breast_Inverse_Persistent.mat')
ktr_opt=x(1);
ve_opt=x(2);
vp_opt=x(3);
leakage_opt=x(4); 

%Load forward model data obtained from the results of the inverse problem in a single point 
load('RR_Breast_Forward_with_inverse_Persistent.mat','CA_XY_time','t_VIF')
ydata=CA_XY_time;
xdata=t_VIF;
load('QIN-BREAST02-01.mat')

%Initial conditions
C_0=ydata(1);   
Init=C_0;

%Forward mLTK model
par_opt=[ktr_opt,ve_opt,vp_opt,leakage_opt];
[t,y] = ode45(@(t,y) mLTK_ODE(t,y,par_opt,VIF,xdata),xdata,Init,[]);
C_opt=y+vp_opt*VIF;

%Evaluate MLE function in the optimum parameter estimation
Sbar=(1/size(ydata,1))*sum(ydata);
Stot=sum((ydata-Sbar).^2); 
MLE_opt=sum((ydata-C_opt).^2);
R2_opt=1-MLE_opt/Stot;

%Paramenters vectors for the optimization
ktr_vec=linspace(ktr_opt/2,2*ktr_opt,10);
ve_vec=linspace(ve_opt/2,2*ve_opt,10);
vp_vec=linspace(vp_opt/2,2*vp_opt,10);
leakage_vec=linspace(leakage_opt/2,2*leakage_opt,10);

%%  ktr fix, optimization of the other parameter
global ktr

for i=1:size(ktr_vec,2)

ve_0=.1;
vp_0=.1;
leakage_0=.005;
ktr=ktr_vec(i);

%Optimization of the parameter values with PS algorithm
x0=[ve_0,vp_0,leakage_0];
lb = [1.e-3 1.e-3 0];
ub = [1 1 1];

options=optimoptions('particleswarm','SwarmSize',200,'HybridFcn',@fmincon,'MaxIterations',500,'ObjectiveLimit',1e-16,'InitialSwarmMatrix',x0,'PlotFcn',@pswplotbestf);
fun=@fun_mLTK_PS_RR_Breast_ktrPe;
nvar=3;
[x,fval,exitflag,output] = particleswarm(fun,nvar,lb,ub,options);

%Forward mLTK model
par_vec=[ktr_vec(i),x(1),x(2),x(3)];
[t,y] = ode45(@(t,y) mLTK_ODE(t,y,par_vec,VIF,xdata),xdata,Init,[]);
C=y+x(2)*VIF;

%MLE evaluation
res=ydata-C;
MLE(i)=sum(res.^2);
R2(i)=1-MLE(i)/Stot; 
xsave(:,i)=x;

end 

%Plot of MLE
figure(1)     
plot(ktr_vec,MLE,ktr_opt,MLE_opt,'r*')
xlabel('k_{tr}')
ylabel('MLE')
title('Type III')
axis square

%Saving options
%save('MLE_ktr_Persistent','MLE','R2','ktr_vec','MLE_opt','ktr_opt','R2_opt','xsave')

 %%  ve fix, optimization of the other parameter
global ve

for i=1:size(ve_vec,2)
  
ktr_0=.005;
vp_0=.1;
leakage_0=.005;
ve=ve_vec(i);

%Optimization of the parameter values with PS algorithm
x0=[ktr_0,vp_0,leakage_0];
lb = [0 1.e-3 0];
ub = [1 1 1];

options=optimoptions('particleswarm','SwarmSize',200,'HybridFcn',@fmincon,'MaxIterations',500,'ObjectiveLimit',1e-16,'InitialSwarmMatrix',x0,'PlotFcn',@pswplotbestf);
fun=@fun_mLTK_PS_RR_Breast_vePe;
nvar=3;
[x,fval,exitflag,output] = particleswarm(fun,nvar,lb,ub,options);

%Forward mLTK model
par_vec=[x(1),ve_vec(i),x(2),x(3)];
[t,y] = ode45(@(t,y) mLTK_ODE(t,y,par_vec,VIF,xdata),xdata,Init,[]);
C=y+x(2)*VIF;

%MLE evaluation
res=ydata-C;
MLE(i)=sum(res.^2);
R2(i)=1-MLE(i)/Stot;  
 xsave(:,i)=x;

end 

%Plot of MLE
figure(2)    
plot(ve_vec,MLE,ve_opt,MLE_opt,'r*')
xlabel('v_e')
ylabel('MLE')
title('Type II')
axis square

%Saving options
%save('MLE_ve_Persistent','MLE','R2','ve_vec','MLE_opt','ve_opt','R2_opt','xsave')

%%  vp fix, optimization of the other parameter
global vp

for i=1:size(vp_vec,2)
  
ktr_0=.005;
ve_0=.1;
leakage_0=.005;
vp=vp_vec(i);

%Optimization of the parameter values with PS algorithm
x0=[ktr_0,ve_0,leakage_0];
lb = [0 1.e-3 0];
ub = [1 1 1];

options=optimoptions('particleswarm','SwarmSize',200,'HybridFcn',@fmincon,'MaxIterations',500,'ObjectiveLimit',1e-16,'InitialSwarmMatrix',x0,'PlotFcn',@pswplotbestf);
fun=@fun_mLTK_PS_RR_Breast_vpPe;
nvar=3;
[x,fval,exitflag,output] = particleswarm(fun,nvar,lb,ub,options);

%Forward mLTK model
par_vec=[x(1),x(2),vp_vec(i),x(3)];
[t,y] = ode45(@(t,y) mLTK_ODE(t,y,par_vec,VIF,xdata),xdata,Init,[]);
C=y+x(2)*VIF;

%MLE evaluation
res=ydata-C;
MLE(i)=sum(res.^2);
R2(i)=1-MLE(i)/Stot; 
xsave(:,i)=x; 

end 

%Plot of MLE
figure(3)    
plot(vp_vec,MLE,vp_opt,MLE_opt,'r*')
xlabel('v_p')
ylabel('MLE')
title('Type II')
axis square 

%Saving options
%save('MLE_vp_Persistent.mat','MLE','R2','vp_vec','MLE_opt','vp_opt','R2_opt','xsave')

%%  leakage fix, optimization of the other parameter
global leakage

for i=1:size(leakage_vec,2)
  
ktr_0=.005;
ve_0=.1;
vp_0=.1;
leakage=leakage_vec(i);

%Optimization of the parameter values with PS algorithm
x0=[ktr_0,ve_0,vp_0];
lb = [0 1.e-3 1.e-3];
ub = [1 1 1];

options=optimoptions('particleswarm','SwarmSize',200,'HybridFcn',@fmincon,'MaxIterations',500,'ObjectiveLimit',1e-16,'InitialSwarmMatrix',x0,'PlotFcn',@pswplotbestf);
fun=@fun_mLTK_PS_RR_Breast_leakagePe;
nvar=3;
[x,fval,exitflag,output]  = particleswarm(fun,nvar,lb,ub,options);

%Forward mLTK model
par_vec=[x(1),x(2),x(3),leakage_vec(i)];
[t,y] = ode45(@(t,y) mLTK_ODE(t,y,par_vec,VIF,xdata),xdata,Init,[]);
C=y+x(2)*VIF;

%MLE evaluation
res=ydata-C;
MLE(i)=sum(res.^2);
R2(i)=1-MLE(i)/Stot; 
xsave(:,i)=x;

end

%Plot of MLE
figure(4)
plot(leakage_vec,MLE,leakage_opt,MLE_opt,'r*')
xlabel('\lambda')
ylabel('MLE')
title('Type I')
axis square

%Saving options
%save('MLE_leakage_Persistent.mat','MLE','R2','leakage_vec','MLE_opt','leakage_opt','R2_opt','xsave')



