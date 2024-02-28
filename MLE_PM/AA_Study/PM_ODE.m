function dydt = PM_ODE(t,y,par,VIF,t_VIF)

%States
C=y; 

%Parameters
   ktr=par(1);

%Equations
VIF=interp1(t_VIF,VIF,t,'linear')'; 

dydt =ktr.*VIF; 

end