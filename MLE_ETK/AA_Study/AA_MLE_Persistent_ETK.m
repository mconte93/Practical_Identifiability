%Code for obtaining the results used in the profile likelihood plot for
%persistent curve in the AA case

clear all
close all

%Load results of the inverse problem
load('AA_Inverse_PS_Persistent.mat')
ktr_opt=x(1);
ve_opt=x(2);
vp_opt=x(3); 

%Load the insilico data
load('AA_ETK_Persistent.mat')
ydata=C_for;
xdata=t_for;

%Initial conditions
C_0=ydata(1);   
Init=C_0;

%Forward ETK model
par_opt=[ktr_opt,ve_opt,vp_opt];
[t,y] = ode45(@(t,y) ETK_ODE(t,y,par_opt,VIF,xdata),xdata,Init,[]);
C_opt=y+vp_opt*VIF';

%Evaluate MLE function in the optimum parameter estimation
Sbar=(1/size(ydata,1))*sum(ydata);
Stot=sum((ydata-Sbar).^2); 
MLE_opt=sum((ydata-C_opt).^2);
R2_opt=1-MLE_opt/Stot;

%Paramenters vectors for the optimization
ktr_vec=linspace(ktr_opt/2,2*ktr_opt,10);
ve_vec=linspace(ve_opt/2,2*ve_opt,10);
vp_vec=linspace(vp_opt/2,2*vp_opt,10);

%%  ktr fix, optimization of the other parameter
global ktr     

for i=1:size(ktr_vec,2)

ve_0=.1;
vp_0=.1;
ktr=ktr_vec(i);

x0=[ve_0,vp_0];
lb = [1.e-3 1.e-3];
ub = [1 1];

%Optimization of the parameter values with PS algorithm
options=optimoptions('particleswarm','SwarmSize',200,'HybridFcn',@fmincon,'MaxIterations',500,'ObjectiveLimit',1e-16,'InitialSwarmMatrix',x0,'PlotFcn',@pswplotbestf);
fun=@fun_ETK_PS_AA_ktrPe;
nvar=2;
[x,fval,exitflag,output] = particleswarm(fun,nvar,lb,ub,options);

%Forward ETK model
par_vec=[ktr_vec(i),x(1),x(2)];
[t,y] = ode45(@(t,y) ETK_ODEArtificial(t,y,par_vec,VIF,xdata),xdata,Init,[]);
C=y+x(2)*VIF';

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
ylabel('SSE')
title('Type I')
axis square

%Saving option
%save('MLE_ktr_Persistent.mat','MLE','R2','ktr_vec','MLE_opt','ktr_opt','R2_opt','xsave')

 %%  ve fix, optimization of the other parameter
global ve

for i=1:size(ve_vec,2)
  
ktr_0=.005;
vp_0=.1;
ve=ve_vec(i);

%Optimization of the parameter values with PS algorithm
x0=[ktr_0,vp_0];
lb = [0 1.e-3];
ub = [1 1];

options=optimoptions('particleswarm','SwarmSize',200,'HybridFcn',@fmincon,'MaxIterations',500,'ObjectiveLimit',1e-16,'InitialSwarmMatrix',x0,'PlotFcn',@pswplotbestf);
fun=@fun_ETK_PS_AA_vePe;
nvar=2;
[x,fval,exitflag,output] = particleswarm(fun,nvar,lb,ub,options);
 
%Forward ETK model
par_vec=[x(1),ve_vec(i),x(2)];
[t,y] = ode45(@(t,y) ETK_ODE(t,y,par_vec,VIF,xdata),xdata,Init,[]);
C=y+x(2)*VIF';

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
title('Type I')
axis square

%Saving option
%save('MLE_ve_Persistent.mat','MLE','R2','ve_vec','MLE_opt','ve_opt','R2_opt','xsave')

%%  vp fix, optimization of the other parameter
global vp

for i=1:size(vp_vec,2)
  
ktr_0=.005;
ve_0=.1;
vp=vp_vec(i);

%Optimization of the parameter values with PS algorithm
x0=[ktr_0,ve_0];
lb = [0 1.e-3];
ub = [1 1];

options=optimoptions('particleswarm','SwarmSize',200,'HybridFcn',@fmincon,'MaxIterations',500,'ObjectiveLimit',1e-16,'InitialSwarmMatrix',x0,'PlotFcn',@pswplotbestf);
fun=@fun_ETK_PS_AA_vpPe;
nvar=2;
[x,fval,exitflag,output] = particleswarm(fun,nvar,lb,ub,options);

%Forward ETK model
par_vec=[x(1),x(2),vp_vec(i)];
[t,y] = ode45(@(t,y) ETK_ODE(t,y,par_vec,VIF,xdata),xdata,Init,[]);
C=y+vp_vec(i)*VIF';

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
ylabel('SSE')
title('Type I')
axis square

%Saving option
%save('MLE_vp_Persistent.mat','MLE','R2','vp_vec','MLE_opt','vp_opt','R2_opt','xsave')
