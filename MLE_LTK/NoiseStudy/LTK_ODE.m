function dydt = LTK_ODE(t,y,par,VIF,t_VIF)

%Function for the forward LTK model

%Parameters
   ktr=par(1);
   ve=par(2);
   vp=par(3);
   lambda=par(4);

%Equations

VIF=interp1(t_VIF,VIF,t,'linear')'; 
Ce=y(1);
Cl=y(2);

dCe=ktr.*(VIF-(Ce./ve));
dCl=lambda*VIF;

dydt=[dCe;dCl];

end