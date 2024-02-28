function [MLE] = fun_mLTK_PS_RA_ktrPe(x)

%Function to evaluate the MLE function in the optimization process with
%ktr parameter fix and persistent curve profile in the RA case

load('RA_mLTK_Persistent.mat');
ydata=C_for; 
xdata=t_for;


global ktr
   
   ve=x(1);
   vp=x(2);
   leakage=x(3);

   %Vector parameter and initial conditions
   par_temp=[ktr,ve,vp,leakage];
   C_0=0;   
   Init=C_0;

    %Forward mLTK model
   [t,y] = ode45(@(t,y) mLTK_ODE(t,y,par_temp,VIF,xdata),xdata,Init,[]);
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