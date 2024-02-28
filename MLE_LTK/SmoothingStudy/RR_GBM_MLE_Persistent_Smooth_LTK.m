%Code for obtaining the results used in the profile likelihood plot for
%persistent curve in the RR-GBM case with different levels of smoothing in 
%the data

%%  ktr fix, optimization of the other parameter, .05% of smoothing
clear all
close all

%Load results of the inverse problem
load('RR_GBM_Inverse_Persistent_Smooth05.mat')
ktr_opt=x(1);
ve_opt=x(2);
vp_opt=x(3);
leakage_opt=x(4); 

%Load forward model data obtained from the results of the inverse problem in a single point 
load('RR_GBM_Forward_with_inverse_Persistent_Smooth05.mat')
ydata=CA_XY_time_Smo';
xdata=t_VIF;

%Load real VIF data
load('Data_PersistentCOH.mat','VIF_mod')
VIF=VIF_mod;

%Initial conditions
Ce_0=0;
Cl_0=0;
Init=[Ce_0;Cl_0];

%Forward LTK model
par_opt=[ktr_opt,ve_opt,vp_opt,leakage_opt];
[t,y] = ode45(@(t,y) LTK_ODE(t,y,par_opt,VIF,xdata),xdata,Init,[]);
C_opt=y(:,1)+y(:,2)+vp_opt*VIF';

%Evaluate MLE function in the optimum parameter estimation
Sbar=(1/size(ydata,1))*sum(ydata);
Stot=sum((ydata-Sbar).^2); 
MLE_opt=sum((ydata-C_opt').^2);
R2_opt=1-MLE_opt/Stot;

%Paramenters vectors for the optimization
ktr_vec=linspace(ktr_opt/2,2*ktr_opt,10);

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
fun=@fun_LTK_PS_RR_Smooth05_ktr;
nvar=3;
[x,fval,exitflag,output] = particleswarm(fun,nvar,lb,ub,options);

%Forward LTK model
par_vec=[ktr_vec(i),x(1),x(2),x(3)];
[t,y] = ode45(@(t,y) LTK_ODE(t,y,par_vec,VIF,xdata),xdata,Init,[]);
C=y(:,1)+y(:,2)+x(2)*VIF';

%MLE evaluation
res=ydata-C';
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

%Sabing option
save('MLE_ktr_Persistent_Smo05','MLE','R2','ktr_vec','MLE_opt','ktr_opt','R2_opt','xsave')


%% ktr fix, optimization of the other parameter, .10% of smoothing
clear all
close all

%Load results of the inverse problem
load('RR_GBM_Inverse_Persistent_Smooth10.mat')
ktr_opt=x(1);
ve_opt=x(2);
vp_opt=x(3);
leakage_opt=x(4); 

%Load forward model data obtained from the results of the inverse problem in a single point 
load('RR_GBM_Forward_with_inverse_Persistent_Smooth10.mat')
ydata=CA_XY_time_Smo';
xdata=t_VIF;

%Load real VIF data
load('Data_PersistentCOH.mat','VIF_mod')
VIF=VIF_mod;

%Initial conditions
Ce_0=0;
Cl_0=0;
Init=[Ce_0;Cl_0];

%Forward LTK model
par_opt=[ktr_opt,ve_opt,vp_opt,leakage_opt];
[t,y] = ode45(@(t,y) LTK_ODE(t,y,par_opt,VIF,xdata),xdata,Init,[]);
C_opt=y(:,1)+y(:,2)+vp_opt*VIF';

%Evaluate MLE function in the optimum parameter estimation
Sbar=(1/size(ydata,1))*sum(ydata);
Stot=sum((ydata-Sbar).^2); 
MLE_opt=sum((ydata-C_opt').^2);
R2_opt=1-MLE_opt/Stot;

%Paramenters vectors for the optimization
ktr_vec=linspace(ktr_opt/2,2*ktr_opt,10);

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
fun=@fun_LTK_PS_RR_Smooth10_ktr;
nvar=3;
[x,fval,exitflag,output] = particleswarm(fun,nvar,lb,ub,options);

%Forward LTK model
par_vec=[ktr_vec(i),x(1),x(2),x(3)];
[t,y] = ode45(@(t,y) LTK_ODE(t,y,par_vec,VIF,xdata),xdata,Init,[]);
C=y(:,1)+y(:,2)+x(2)*VIF';

%MLE evaluation
res=ydata-C';
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

%Sabing option
%save('MLE_ktr_Persistent_Smooth10','MLE','R2','ktr_vec','MLE_opt','ktr_opt','R2_opt','xsave')

 %% ktr fix, optimization of the other parameter, .15% of smoothing
clear all
close all

%Load results of the inverse problem
load('RR_GBM_Inverse_Persistent_Smooth15.mat')
ktr_opt=x(1);
ve_opt=x(2);
vp_opt=x(3);
leakage_opt=x(4); 

%Load forward model data obtained from the results of the inverse problem in a single point 
load('RR_GBM_Forward_with_inverse_Persistent_Smooth15.mat')
ydata=CA_XY_time_Smo';
xdata=t_VIF;

%Load real VIF data
load('Data_PersistentCOH.mat','VIF_mod')
VIF=VIF_mod;

%Initial conditions
Ce_0=0;
Cl_0=0;
Init=[Ce_0;Cl_0];

%Forward LTK model
par_opt=[ktr_opt,ve_opt,vp_opt,leakage_opt];
[t,y] = ode45(@(t,y) LTK_ODE(t,y,par_opt,VIF,xdata),xdata,Init,[]);
C_opt=y(:,1)+y(:,2)+vp_opt*VIF';

%Evaluate MLE function in the optimum parameter estimation
Sbar=(1/size(ydata,1))*sum(ydata);
Stot=sum((ydata-Sbar).^2); 
MLE_opt=sum((ydata-C_opt').^2);
R2_opt=1-MLE_opt/Stot;

%Paramenters vectors for the optimization
ktr_vec=linspace(ktr_opt/2,2*ktr_opt,10);

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
fun=@fun_LTK_PS_RR_Smooth15_ktr;
nvar=3;
[x,fval,exitflag,output] = particleswarm(fun,nvar,lb,ub,options);

%Forward LTK model
par_vec=[ktr_vec(i),x(1),x(2),x(3)];
[t,y] = ode45(@(t,y) LTK_ODE(t,y,par_vec,VIF,xdata),xdata,Init,[]);
C=y(:,1)+y(:,2)+x(2)*VIF';

%MLE evaluation
res=ydata-C';
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

%Sabing option
%save('MLE_ktr_Persistent_Smo15','MLE','R2','ktr_vec','MLE_opt','ktr_opt','R2_opt','xsave')