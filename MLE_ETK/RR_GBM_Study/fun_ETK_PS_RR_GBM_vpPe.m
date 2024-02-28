function [MLE] = fun_ETK_PS_RR_GBM_vpPe(x)

%Function to evaluate the MLE function in the optimization process with
%vp parameter fix and persistent curve profile in the RR case

load('Forward_with_inverse_GBM_Persistent.mat','CA_XY_time','t_VIF')
ydata=CA_XY_time;
xdata=t_VIF;

load('Data_PersistentCOH.mat','VIF_mod')
VIF=VIF_mod;

global vp
   
   ktr=x(1);
   ve=x(2);

   %Vector parameter and initial conditions
   par_temp=[ktr,ve,vp];
   C_0=0;   
   Init=C_0;

   %Forward ETK model
   [t,y] = ode45(@(t,y) ETK_ODE(t,y,par_temp,VIF,xdata),xdata,Init,[]);
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