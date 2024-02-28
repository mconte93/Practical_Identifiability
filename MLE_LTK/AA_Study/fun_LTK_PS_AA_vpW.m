function [MLE] = fun_LTK_PS_AA_vpW(x)

%Function to evaluate the MLE function in the optimization process with
%vp parameter fix and wash-out curve profile in the AA case

load('AA_LTK_WashOut.mat');
ydata=C_for; 
xdata=t_for;

global vp
   
   ktr=x(1);
   ve=x(2);
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