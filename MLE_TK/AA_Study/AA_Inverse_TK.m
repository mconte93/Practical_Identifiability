%Code to run the inverse problem on the TK model with insilico data and
%insilico VIF to obtain the optimal estimation of the model parameters

clear all
close all

ktr_0=.005;
ve_0=.1;

%Optimization of the parameter values with PS algorithm
x0=[ktr_0,ve_0];
fun=@fun_TK_PS_AA;
lb = [0 1.e-3];
ub = [1 1];

options=optimoptions('particleswarm','SwarmSize',200,'HybridFcn',@fmincon,'MaxIterations',500,'ObjectiveLimit',1e-16,'InitialSwarmMatrix',x0,'PlotFcn',@pswplotbestf);
nvar=2;
[x,fval,exitflag,output] = particleswarm(fun,nvar,lb,ub,options);
formatstring = 'particleswarm reached the value %f using %d function evaluations.\n';
fprintf(formatstring,fval,output.funccount)

%Saving option
%save('AA_Inverse_Persistent.mat','x')


function [MLE] = fun_TK_PS_AA(x)

%Function to evaluate the MLE function 
    
    %Load the data 
    load('ArtificialVIF_TK_WashOut.mat');
    ydata=C_for; 
    xdata=t_for;

    %Model parameter
    ktr=x(1);
    ve=x(2);
    par_temp=[ktr,ve];

    %Initial conditions
    C_0=0;   
    Init=C_0;

    %Forward TK model
   [t,y] = ode45(@(t,y) TK_ODE(t,y,par_temp,VIF,t_for),t_for,Init,[]);
    C_approx=y;
   
   %MLE evaluation
   MLE=0;
   for i=1:size(C_approx,1)
       res=ydata(i)-C_approx(i);
       if ~isnan(res)
          MLE=MLE+res.^2;
       end
   end

    
end