function [MLE] = fun_PM_PS_AA_ktrPl(x)

%Function to evaluate the MLE function in the optimization process with
%ktr parameter fix and plateau curve profile in the AA case

load('AA_PM_Plateau.mat');
ydata=C_for; 
xdata=t_for;

global ktr
   
   vp=x(1);

   %Vector parameter and initial conditions
   par_temp=[ktr,vp];
   C_0=0;   
   Init=C_0;

   %Forward PM model
   [t,y] = ode45(@(t,y) PM_ODE(t,y,par_temp,VIF,xdata),xdata,Init,[]);
   C_approx=y+vp*VIF';
   
   %MLE evalutaion
   MLE=0;
   for i=1:size(C_approx,1)
       res=ydata(i)-C_approx(i);
       if ~isnan(res)
          MLE=MLE+res.^2;
       end
   end

    
end