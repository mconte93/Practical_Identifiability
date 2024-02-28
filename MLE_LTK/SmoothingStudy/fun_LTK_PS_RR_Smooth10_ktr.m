function [MLE] = fun_LTK_PS_RR_Smooth10_ktr(x)

%Function to evaluate the MLE function in the optimization process with
%ktr parameter fix and persistent curve profile in the RR-GBM case with .10% of
%smooth in the original data

load('RR_GBM_Forward_with_inverse_Persistent_Smooth10.mat','CA_XY_time_Smo','t_VIF')
ydata=CA_XY_time_Smo;
xdata=t_VIF;

global ktr
   
   ve=x(1);
   vp=x(2);
   leakage=x(3);

   %Vector parameter and initial conditions
   par_temp=[ktr,ve,vp,leakage];
   Ce_0=0;   
   Cl_0=0;
   Init=[Ce_0;Cl_0];

   %Forward LTK model
   [t,y] = ode45(@(t,y) LTK_ODE(t,y,par_temp,VIF,xdata),xdata,Init,[]);
   C_approx=y(:,1)+y(:,2)+vp*VIF';

   %MLE evalutaion
   MLE=0;
   for i=1:size(C_approx,1)
       res=ydata(i)-C_approx(i);
       if ~isnan(res)
          MLE=MLE+res.^2;
       end
   end

    
end