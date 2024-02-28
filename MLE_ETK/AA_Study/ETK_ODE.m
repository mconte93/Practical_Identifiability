function dydt = ETK_ODE(t,y,par,VIF,t_VIF)

%States
C=y; 

%Parameters
   ktr=par(1);
   ve=par(2);

%Equations
VIF=interp1(t_VIF,VIF,t,'linear')'; 

dydt =ktr.*(VIF-(C./ve)); 

end