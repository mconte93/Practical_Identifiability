function [MLE] = fun_LTK_PS_AA_10Noise_Ktr(x)

%Function to evaluate the MLE function in the optimization process with
%ktr parameter fix and persistent curve profile in the AA case with 5% of
%noise in the VIF

load('AA_Inverse_Persistent_10Noise.mat');
ydata=C_for; 
xdata=t_for;

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