% This MATLAB program solve the equation of motion and create plots for various values of excitation frequency
t = 0:0.01:10;  	 % Simulation time

% Bridge parameters
m=100000;               % Bridge mass [kg]
zeta=0.01;     	% Damping factor [-]
k=4000000;              % Bridge stiffness [N/m]
Wn=sqrt(k/m);         % Angular natural frequency [rad/s]
fn=Wn/(2*pi)           % Natural frequency [Hz]. Note that f[Hz]=W/2*pi

% Excitation parameters
Fo=20000;                % Excitation force acting on the bridge [N]
W=[0.5*Wn 0.945*Wn 1.19*Wn]; 	% Excitation frequency of the force acting on the bridge [rad/s]
 
  for jzeta=1:length(zeta)
    for j=1:length(W)
      % Steady state solution - Displacement
      x{jzeta}(j,:)=(((Wn^2-W(j)^2)*Fo/m)/((Wn^2-W(j)^2)^2+(2*zeta(jzeta)*W(j)*Wn).^2))*cos(W(j)*t)+...
        (((2*zeta(jzeta)*W(j)*Wn)*Fo/m)/((Wn^2-W(j)^2)^2+(2*zeta(jzeta)*W(j)*Wn)^2))*sin(W(j)*t);
        
      dx{jzeta}(j,:)=(((Wn^2-W(j)^2)*Fo/m)/((Wn^2-W(j)^2)^2+(2*zeta(jzeta)*W(j)*Wn).^2))*(-W(j))*sin(W(j)*t)+...
        (((2*zeta(jzeta)*W(j)*Wn)*Fo/m)/((Wn^2-W(j)^2)^2+(2*zeta(jzeta)*W(j)*Wn)^2))*(W(j))*cos(W(j)*t);

      ddx{jzeta}(j,:)=(((Wn^2-W(j)^2)*Fo/m)/((Wn^2-W(j)^2)^2+(2*zeta(jzeta)*W(j)*Wn).^2))*(-W(j)^2)*cos(W(j)*t)+...
        (((2*zeta(jzeta)*W(j)*Wn)*Fo/m)/((Wn^2-W(j)^2)^2+(2*zeta(jzeta)*W(j)*Wn)^2))*(-W(j)^2)*sin(W(j)*t);
        
      legend_W{j}    = sprintf('f=%2.2f Hz; Fo=20 kN\n', W(j)/(2*pi));  
      xm{jzeta}(j,:) = max(x{jzeta}(j,:));
    end
      legend_zeta{jzeta}=sprintf('zeta=%2.2f\n', zeta(jzeta));
  end

 % Plot displacement,velocity, acceleration vs time for one value of omega, zeta...
  figure
  subplot(3,1,1)
    plot(t,1000*x{1}(1,:))
    grid on
    xlabel('Time, [sec]',"FontSize",12)
    ylabel('Displacement, [mm]',"FontSize",12)
    %h=legend(legend_W)
    %set (h, "fontsize", 12);
  subplot(3,1,2)
    plot(t,1000*dx{1}(1,:))    
    title('Bridge Vertical Displacement')

    grid on
    xlabel('Time, [sec]',"FontSize",12)
    ylabel('Velocity, [m/s]',"FontSize",12)
    title('Bridge Vertical Velocity',"FontSize",14)
    %h=legend(legend_W{})
    %set (h, "fontsize", 12); 
  subplot(3,1,3)
    plot(t,1000*ddx{1}(1,:))
    grid on
    xlabel('Time, [sec]',"FontSize",12)
    ylabel('Acceleration, [m/s^2]',"FontSize",12)
    title('Bridge Vertical Acceleration',"FontSize",14)
    %h=legend(legend_W{})
    %set (h, "fontsize",8); 
    filenam = sprintf('Disp_Velo_Acc_Vs_Time zeta=%1.2f Fo=%6.0f k=%9.0f.png',zeta(jzeta), Fo, k)
    saveas(gcf,filenam,'png')  
  
 % Plot displacement vs time for 3 values of omega
  figure
  plot(t,1000*x{1}(),'-x')
  grid on
  xlabel('Time, [sec]',"FontSize",16)
  ylabel('Displacement, [mm]',"FontSize",16)
  ylim([-50 50])
  title('Bridge Vertical Displacement',"FontSize",16)
  %h=legend(legend_W{})
  %set (h, "fontsize", 12);
  filenam = sprintf('Disp_Vs_Time_various_omega zeta=%1.2f Fo=%6.0f k=%9.0f.png',zeta(jzeta), Fo, k)
  saveas(gcf,filenam,'png')  

   % Plot acceleration vs time for 3 values of omega
  figure
  plot(t,1000*ddx{1}(),'-x')
  grid on
  xlabel('Time, [sec]',"FontSize",16)
  ylabel('Acceleration, [m/s^2]',"FontSize",16)
  title('Bridge Vertical Acceleration',"FontSize",16)
  %h=legend(legend_W{})
  %set (h, "fontsize", 12);
  filenam = sprintf('Acc_Vs_Time_various_omega zeta=%1.2f Fo=%6.0f k=%9.0f.png',zeta(jzeta), Fo, k)
  saveas(gcf,filenam,'png') 

 
