function [MLE] = fun_PM_PS_RR_Breast_vpPl(x)

%Function to evaluate the MLE function in the optimization process with
%ve parameter fix and plateau curve profile in the RR case

load('RR_Breast_Forward_with_inverse_Plateau.mat','CA_XY_time','t_VIF')
ydata=CA_XY_time;
xdata=t_VIF;
load('QIN-BREAST02-01.mat','VIF')

global vp
   
   ktr=x(1);

    %Vector parameter and initial conditions
   par_temp=[ktr,vp];
   C_0=0;   
   Init=C_0;

   %Forward PM model
   [t,y] = ode45(@(t,y) PM_ODE(t,y,par_temp,VIF,xdata),xdata,Init,[]);
   C_approx=y+vp*VIF;

   %MLE evalutation
   MLE=0;
   for i=1:size(C_approx,1)
       res=ydata(i)-C_approx(i);
       if ~isnan(res)
          MLE=MLE+res.^2;
       end
   end

    
end