%Code for obtaining the results used in the profile likelihood plot for
%plateau curve in the RR-GBM case

clear all
close all

%Load results of the inverse problem
load('RR_GBM_Inverse_Plateau.mat')
ktr_opt=x(1);
ve_opt=x(2);

%Load forward model data obtained from the results of the inverse problem in a single point 
load('RR_GBM_Forward_with_inverse_Plateau.mat','CA_XY_time','t_VIF')
ydata=CA_XY_time;
xdata=t_VIF;
 
%Load real VIF data
load('QIN-GBM-TR-10.mat','VIF')

%Initial conditions
C_0=ydata(1);   
Init=C_0;

%Forward TK model
par_opt=[ktr_opt,ve_opt];
[t,y] = ode45(@(t,y) TK_ODE(t,y,par_opt,VIF,xdata),xdata,Init,[]);
C_opt=y;

%Evaluate MLE function in the optimum parameter estimation
Sbar=(1/size(ydata,1))*sum(ydata);
Stot=sum((ydata-Sbar).^2); 
MLE_opt=sum((ydata-C_opt).^2);
R2_opt=1-MLE_opt/Stot;

%Paramenters vectors for the optimization
ktr_vec=linspace(ktr_opt/2,2*ktr_opt,10);
ve_vec=linspace(ve_opt/2,2*ve_opt,10);

%%  ktr fix, optimization of the other parameter
global ktr

for i=1:size(ktr_vec,2)

ve_0=.1;
ktr=ktr_vec(i);

%Optimization of the parameter values with PS algorithm
x0=[ve_0];
lb = [1.e-3];
ub = [1];

options=optimoptions('particleswarm','SwarmSize',200,'HybridFcn',@fmincon,'MaxIterations',500,'ObjectiveLimit',1e-16,'InitialSwarmMatrix',x0,'PlotFcn',@pswplotbestf);
fun=@fun_TK_PS_RR_GBM_ktrPl;
nvar=1;
[x,fval,exitflag,output] = particleswarm(fun,nvar,lb,ub,options);

%Forward TK model
par_vec=[ktr_vec(i),x(1)];
[t,y] = ode45(@(t,y) TK_ODE(t,y,par_vec,VIF,xdata),xdata,Init,[]);
C=y;

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
%save('SEE_ktr_Plateau','MLE','R2','ktr_vec','MLE_opt','ktr_opt','R2_opt','xsave')

%%  ve fix, optimization of the other parameter
global ve

for i=1:size(ve_vec,2)
  
ktr_0=.005;
ve=ve_vec(i);

%Optimization of the parameter values with PS algorithm
x0=[ktr_0];
lb = [0];
ub = [1];

options=optimoptions('particleswarm','SwarmSize',200,'HybridFcn',@fmincon,'MaxIterations',500,'ObjectiveLimit',1e-16,'InitialSwarmMatrix',x0,'PlotFcn',@pswplotbestf);
fun=@fun_TK_PS_RR_GBM_vePl;
nvar=1;
[x,fval,exitflag,output] = particleswarm(fun,nvar,lb,ub,options);
 
%Forward TK model
par_vec=[x(1),ve_vec(i)];
[t,y] = ode45(@(t,y) TK_ODE(t,y,par_vec,VIF,xdata),xdata,Init,[]);
C=y;

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
%save('SEE_ve_Plateau','MLE','R2','ve_vec','MLE_opt','ve_opt','R2_opt','xsave')
