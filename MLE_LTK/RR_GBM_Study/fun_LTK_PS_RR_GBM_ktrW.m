function [MLE] = fun_LTK_PS_RR_GBM_ktrW(x)

%Function to evaluate the MLE function in the optimization process with
%ktr parameter fix and wash-out curve profile in the RR case

load('RR_GBM_Forward_with_inverse_WashOut.mat','CA_XY_time','t_VIF')
ydata=CA_XY_time;
xdata=t_VIF;
load('QIN-GBM-TR-10.mat','VIF')

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
   
   %MLE evalutation
   MLE=0;
   for i=1:size(C_approx,1)
       res=ydata(i)-C_approx(i);
       if ~isnan(res)
          MLE=MLE+res.^2;
       end
   end

    
end