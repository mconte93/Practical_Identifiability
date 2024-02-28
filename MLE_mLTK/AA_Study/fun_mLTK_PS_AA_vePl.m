function [MLE] = fun_mLTK_PS_AA_vePl(x)

%Function to evaluate the MLE function in the optimization process with
%ve parameter fix and plateau curve profile in the AA case

load('AA_mLTK_Plateau.mat');
ydata=C_for; 
xdata=t_for;


global ve
   
   ktr=x(1);
   vp=x(2);
   leakage=x(3);

   %Vector parameter and initial conditions
   par_temp=[ktr,ve,vp,leakage];
   C_0=0;   
   Init=C_0;

   %Forward mLTK model
   [t,y] = ode45(@(t,y) mLTK_ODE(t,y,par_temp,VIF,xdata),xdata,Init,[]);
   C_approx=y+vp*VIF;
   
  %MLE evalutaion
   MLE=0;
   for i=1:size(C_approx,1)
       res=ydata(i)-C_approx(i);
       if ~isnan(res)
          MLE=MLE+res.^2;
       end
   end

    
end