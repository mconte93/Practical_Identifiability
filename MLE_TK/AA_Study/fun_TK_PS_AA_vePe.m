function [MLE] = fun_TK_PS_AA_vePe(x)

%Function to evaluate the MLE function in the optimization process with
%ve parameter fix and persistent curve profile in the AA case

load('AA_TK_Persistent.mat');
ydata=C_for; 
xdata=t_for;

global ve
   
   ktr=x(1);

   %Vector parameter and initial conditions
   par_temp=[ktr,ve];
   C_0=0;   
   Init=C_0;

   %Forward TK model
   [t,y] = ode45(@(t,y) TK_ODE(t,y,par_temp,VIF,xdata),xdata,Init,[]);
   C_approx=y;
   
   %MLE evalutation
   MLE=0;
   for i=1:size(C_approx,1)
       res=ydata(i)-C_approx(i);
       if ~isnan(res)
          MLE=MLE+res.^2;
       end
   end

    
end