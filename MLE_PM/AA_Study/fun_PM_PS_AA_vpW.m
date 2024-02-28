function [MLE] = fun_PM_PS_AA_vpW(x)

%Function to evaluate the MLE function in the optimization process with
%vp parameter fix and wash-out curve profile in the AA case

load('AA_PM_WashOut.mat');
ydata=C_for; 
xdata=t_for;

global vp
   
   ktr=x(1);

   %Vector parameter and initial conditions
   par_temp=[ktr,vp];
   C_0=0;   
   Init=C_0;

   %Forward PM model
   [t,y] = ode45(@(t,y) PM_ODE(t,y,par_temp,VIF,xdata),xdata,Init,[]);
   C_approx=y+vp*VIF';
   
   %MLE evalutation
   MLE=0;
   for i=1:size(C_approx,1)
       res=ydata(i)-C_approx(i);
       if ~isnan(res)
          MLE=MLE+res.^2;
       end
   end

    
end