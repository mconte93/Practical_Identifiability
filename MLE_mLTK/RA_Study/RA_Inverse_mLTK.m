%Code to run the inverse problem on the mLTK model with insilico data and
%VIF estimated from the data to obtain the optimal estimation of the model 
%parameters

clear all
close all

ktr_0=.005;
ve_0=.1;
lambda_0=.005;
vp_0=.1;

%Optimization of the parameter values with PS algorithm
x0=[ktr_0 ve_0 vp_0 lambda_0];
fun=@fun_mLTK_PS_RA;
lb = [0 1.e-3 1.e-3 0];
ub = [1 1 1 1];

options=optimoptions('particleswarm','SwarmSize',200,'HybridFcn',@fmincon,'MaxIterations',500,'ObjectiveLimit',1e-16,'InitialSwarmMatrix',x0,'PlotFcn',@pswplotbestf);
nvar=4;
[x,fval,exitflag,output] = particleswarm(fun,nvar,lb,ub,options);
formatstring = 'particleswarm reached the value %f using %d function evaluations.\n';
fprintf(formatstring,fval,output.funccount)

%Saving option
%save('RA_Inverse_Persistent.mat','x')


function [MLE] = fun_mLTK_PS_RA(x)

%Function to evaluate the MLE function 

%Load the data   
load('RA_mLTK_Persistent.mat');
ydata=C_for; 
xdata=t_for;

%Model parameter
   ktr=x(1);
   ve=x(2);
   vp=x(3);
   lambda=x(4);
   
   par_temp=[ktr,ve,vp,lambda];

   %Initial conditions
   C_0=0;  
   Init=C_0;

   %Forward mLTK model
   [t,y] = ode45(@(t,y) mLTK_ODE(t,y,par_temp,VIF,xdata),xdata,Init,[]);
   C_approx=y+vp*VIF';
   
   %MLE evaluation
   MLE=0;
   for i=1:size(C_approx,1)
       res=ydata(i)-C_approx(i);
       if ~isnan(res)
          MLE=MLE+res.^2;
       end
   end

    
end