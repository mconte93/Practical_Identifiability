%Code for obtaining the results used in the profile likelihood plot for
%plateau curve in the AA case

clear all
close all

%Load results of the inverse problem
load('AA_Inverse_Plateau.mat')
ktr_opt=x(1);
vp_opt=x(2); 

%Load the insilico data
load('AA_PM_Plateau.mat')
ydata=C_for;
xdata=t_for;

%Initial conditions
C_0=ydata(1);   
Init=C_0;

%Forward PM model
par_opt=[ktr_opt,vp_opt];
[t,y] = ode45(@(t,y) PM_ODE(t,y,par_opt,VIF,xdata),xdata,Init,[]);
C_opt=y+vp_opt*VIF';

%Evaluate MLE function in the optimum parameter estimation
Sbar=(1/size(ydata,1))*sum(ydata);
Stot=sum((ydata-Sbar).^2); 
MLE_opt=sum((ydata-C_opt).^2);
R2_opt=1-MLE_opt/Stot;

%Paramenters vectors for the optimization
ktr_vec=linspace(ktr_opt/2,2*ktr_opt,10);
vp_vec=linspace(vp_opt/2,2*vp_opt,10);

%%  ktr fix, optimization of the other parameter
global ktr   

for i=1:size(ktr_vec,2)

vp_0=.1;
ktr=ktr_vec(i);

%Optimization of the parameter values with PS algorithm
x0=[vp_0];
lb = [1.e-3];
ub = [1];

options=optimoptions('particleswarm','SwarmSize',200,'HybridFcn',@fmincon,'MaxIterations',500,'ObjectiveLimit',1e-16,'InitialSwarmMatrix',x0,'PlotFcn',@pswplotbestf);
fun=@fun_PM_PS_AA_ktrPl;
nvar=1;
[x,fval,exitflag,output] = particleswarm(fun,nvar,lb,ub,options);

%Forward PM model
par_vec=[ktr_vec(i),x(1)];
[t,y] = ode45(@(t,y) PM_ODE(t,y,par_vec,VIF,xdata),xdata,Init,[]);
C=y+x(1)*VIF';

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
title('Type I')
axis square

%Saving option
%save('MLE_ktr_Plateau.mat','MLE','R2','ktr_vec','MLE_opt','ktr_opt','R2_opt','xsave')

%%  vp fix, optimization of the other parameter
global vp

for i=1:size(vp_vec,2)
  
ktr_0=.005;
vp=vp_vec(i);

%Optimization of the parameter values with PS algorithm
x0=[ktr_0];
lb = [0];
ub = [1];

options=optimoptions('particleswarm','SwarmSize',200,'HybridFcn',@fmincon,'MaxIterations',500,'ObjectiveLimit',1e-16,'InitialSwarmMatrix',x0,'PlotFcn',@pswplotbestf);
fun=@fun_PM_PS_AA_vpPl;
nvar=1;
[x,fval,exitflag,output] = particleswarm(fun,nvar,lb,ub,options);

%Forward PM model
par_vec=[x(1),vp_vec(i)];
[t,y] = ode45(@(t,y) PM_ODE(t,y,par_vec,VIF,xdata),xdata,Init,[]);
C=y+vp_vec(i)*VIF';

%MLE evaluation
res=ydata-C;
MLE(i)=sum(res.^2);
R2(i)=1-MLE(i)/Stot;
xsave(:,i)=x;

end 

%Plot of MLE
figure(2)    
plot(vp_vec,MLE,vp_opt,MLE_opt,'r*')
xlabel('v_p')
ylabel('MLE')
title('Type I')
axis square

%Saving option
%save('MLE_vp_Plateau.mat','MLE','R2','vp_vec','MLE_opt','vp_opt','R2_opt','xsave')