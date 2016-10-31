% Project 3
% Clair Cunningham

clear all; close all; clc
fid = fopen('Project_3.txt','w+');

% %Sets the units of your root object (screen) to pixels
 set(0,'units','pixels');
% 
% %Obtains this pixel information
 Pix_SS = get(0,'screensize');
% 
% %Sets the units of your root object (screen) to inches
% set(0,'units','inches');
% 
% %Obtains this inch information
% Inch_SS = get(0,'screensize');
% 
% %Calculates the resolution (pixels per inch)
% Res = Pix_SS./Inch_SS;


% initial trim conditions
u_0 = 670.360471; % ft/sec
v_0 = 0; % ft/sec
w_0 = 40.362171; % ft/sec
p_0 = 0; % ft/sec
q_0 = 0; % ft/sec
r_0 = 0; % radians
phi_0 = 0; % radians
theta_0 = atan(w_0/u_0); % radians
psi_0 = 0; % radians

% Initial Conditions and Time Step
IC = [u_0, v_0, w_0, p_0, q_0, r_0, phi_0, theta_0, psi_0];
time = [0:0.01:10]';

% get solution for trim
[time_out,x_out] = ode45(@Project_3_f,time,IC); 
u = x_out(:,1);
v = x_out(:,2);
w = x_out(:,3);
p = x_out(:,4)*180/pi;
q = x_out(:,5)*180/pi;
r = x_out(:,6)*180/pi;
phi = x_out(:,7)*180/pi;
theta = x_out(:,8)*180/pi;
psi = x_out(:,9)*180/pi;
fig = figure('OuterPosition',[0 0 Pix_SS(3)*.90 Pix_SS(4)*.90],'PaperPositionMode','auto');
subplot(3,3,1)
plot(time_out,u,'LineWidth',3.0); title({'Axial Velocity','at Trim'});xlabel('Time (sec)');ylabel('Axial Velocity (ft/sec)'); grid on
%set(gca, 'YTickLabel', num2str(get(gca,'YTick')','%d'))
subplot(3,3,2)
plot(time_out,v,'LineWidth',3.0); title({'Side Velocity','at Trim'});xlabel('Time (sec)');ylabel('Side Velocity (ft/sec)'); grid on
subplot(3,3,3)
plot(time_out,w,'LineWidth',3.0); title({'Normal Velocity','at Trim'});xlabel('Time (sec)');ylabel('Normal Velocity (ft/sec)'); grid on
subplot(3,3,4)
plot(time_out,p,'LineWidth',3.0); title({'Roll Rate','at Trim'});xlabel('Time (sec)');ylabel('Roll Rate (deg/sec)'); grid on
subplot(3,3,5)
plot(time_out,q,'LineWidth',3.0); title({'Pitch Rate','at Trim'});xlabel('Time (sec)');ylabel('Pitch Rate (deg/sec)'); grid on
subplot(3,3,6)
plot(time_out,r,'LineWidth',3.0); title({'Yaw Rate','at Trim'});xlabel('Time (sec)');ylabel('Yaw Rate (deg/sec)'); grid on
subplot(3,3,7)
plot(time_out,phi,'LineWidth',3.0); title({'Bank Angle','at Trim'});xlabel('Time (sec)');ylabel('Axial Velocity (deg)'); grid on
subplot(3,3,8)
plot(time_out,theta,'LineWidth',3.0); title({'Pitch Angle','at Trim'});xlabel('Time (sec)');ylabel('Pitch Angel (deg/sec)'); grid on
subplot(3,3,9)
plot(time_out,psi,'LineWidth',3.0); title({'Heading Angle','at Trim'});xlabel('Time (sec)');ylabel('Heading Angle (deg)'); grid on

name = ['Figure_' num2str(fig.Number)];
print(fig,'-depsc','-noui','-painters',name);
fprintf(fid,'\\sectionmark{Project \\# 3\\hspace*{\\fill} Clair Cunningham \\hspace*{\\fill} Problem %d}\n',fig.Number);
fprintf(fid,'The figure below while not the exact same as the solution given has variance in the 10$^{-6}$ place only, a small number. Even if the model output is not exactly the same it still remains within a very close trim condition given no input. The rest of the solutions for the different problems follow with no variance from the solution given.');
fprintf(fid,'\\vspace*{\\fill}\\begin{figure}[H]\\centering\\includegraphics[keepaspectratio=true,height=1\\textheight,width=1\\textwidth,angle=90]{%s.eps}\n \\caption{Problem %d Graphical Solution}\\end{figure}\\vspace*{\\fill}\n\\newpage\n',name,fig.Number);

% get solution for step input of -0.5 deg to the elevator
[time_out,x_out] = ode45(@Project_3_de_f,time,IC); 
u = x_out(:,1);
v = x_out(:,2);
w = x_out(:,3);
p = x_out(:,4)*180/pi;
q = x_out(:,5)*180/pi;
r = x_out(:,6)*180/pi;
phi = x_out(:,7)*180/pi;
theta = x_out(:,8)*180/pi;
psi = x_out(:,9)*180/pi;
fig = figure('OuterPosition',[0 0 Pix_SS(3)*.90 Pix_SS(4)*.90],'PaperPositionMode','auto');
subplot(3,3,1)
plot(time_out,u,'LineWidth',3.0); title({'Axial Velocity','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Axial Velocity (ft/sec)'); grid on
subplot(3,3,2)
plot(time_out,v,'LineWidth',3.0); title({'Side Velocity','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Side Velocity (ft/sec)'); grid on
subplot(3,3,3)
plot(time_out,w,'LineWidth',3.0); title({'Normal Velocity','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Normal Velocity (ft/sec)'); grid on
subplot(3,3,4)
plot(time_out,p,'LineWidth',3.0); title({'Roll Rate','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Roll Rate (deg/sec)'); grid on
subplot(3,3,5)
plot(time_out,q,'LineWidth',3.0); title({'Pitch Rate','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Pitch Rate (deg/sec)'); grid on
subplot(3,3,6)
plot(time_out,r,'LineWidth',3.0); title({'Yaw Rate','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Yaw Rate (deg/sec)'); grid on
subplot(3,3,7)
plot(time_out,phi,'LineWidth',3.0); title({'Bank Angle','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Axial Velocity (deg)'); grid on
subplot(3,3,8)
plot(time_out,theta,'LineWidth',3.0); title({'Pitch Angle','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Pitch Angel (deg/sec)'); grid on
subplot(3,3,9)
plot(time_out,psi,'LineWidth',3.0); title({'Heading Angle','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Heading Angle (deg)'); grid on

name = ['Figure_' num2str(fig.Number)];
print(fig,'-depsc','-noui','-painters',name);
fprintf(fid,'\\sectionmark{Project \\# 3\\hspace*{\\fill} Clair Cunningham \\hspace*{\\fill} Problem %d}\n',fig.Number);
fprintf(fid,'\\vspace*{\\fill}\\begin{figure}[H]\\centering\\includegraphics[keepaspectratio=true,height=1\\textheight,width=1\\textwidth,angle=90]{%s.eps}\n \\caption{Problem %d Graphical Solution}\\end{figure}\\vspace*{\\fill}\n\\newpage\n',name,fig.Number);

% get solution for step input of -0.5 deg to the aileron
[time_out,x_out] = ode45(@Project_3_da_f,time,IC); 
u = x_out(:,1);
v = x_out(:,2);
w = x_out(:,3);
p = x_out(:,4)*180/pi;
q = x_out(:,5)*180/pi;
r = x_out(:,6)*180/pi;
phi = x_out(:,7)*180/pi;
theta = x_out(:,8)*180/pi;
psi = x_out(:,9)*180/pi;
fig = figure('OuterPosition',[0 0 Pix_SS(3)*.90 Pix_SS(4)*.90],'PaperPositionMode','auto');
subplot(3,3,1)
plot(time_out,u,'LineWidth',3.0); title({'Axial Velocity','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Axial Velocity (ft/sec)'); grid on
subplot(3,3,2)
plot(time_out,v,'LineWidth',3.0); title({'Side Velocity','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Side Velocity (ft/sec)'); grid on
subplot(3,3,3)
plot(time_out,w,'LineWidth',3.0); title({'Normal Velocity','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Normal Velocity (ft/sec)'); grid on
subplot(3,3,4)
plot(time_out,p,'LineWidth',3.0); title({'Roll Rate','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Roll Rate (deg/sec)'); grid on
subplot(3,3,5)
plot(time_out,q,'LineWidth',3.0); title({'Pitch Rate','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Pitch Rate (deg/sec)'); grid on
subplot(3,3,6)
plot(time_out,r,'LineWidth',3.0); title({'Yaw Rate','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Yaw Rate (deg/sec)'); grid on
subplot(3,3,7)
plot(time_out,phi,'LineWidth',3.0); title({'Bank Angle','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Axial Velocity (deg)'); grid on
subplot(3,3,8)
plot(time_out,theta,'LineWidth',3.0); title({'Pitch Angle','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Pitch Angel (deg/sec)'); grid on
subplot(3,3,9)
plot(time_out,psi,'LineWidth',3.0); title({'Heading Angle','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Heading Angle (deg)'); grid on

name = ['Figure_' num2str(fig.Number)];
print(fig,'-depsc','-noui','-painters',name);
fprintf(fid,'\\sectionmark{Project \\# 3\\hspace*{\\fill} Clair Cunningham \\hspace*{\\fill} Problem %d}\n',fig.Number);
fprintf(fid,'\\vspace*{\\fill}\\begin{figure}[H]\\centering\\includegraphics[keepaspectratio=true,height=1\\textheight,width=1\\textwidth,angle=90]{%s.eps}\n \\caption{Problem %d Graphical Solution}\\end{figure}\\vspace*{\\fill}\n\\newpage\n',name,fig.Number);

% get solution for step input of -2.0 deg to the rudder
[time_out,x_out] = ode45(@Project_3_dr_f,time,IC); 
u = x_out(:,1);
v = x_out(:,2);
w = x_out(:,3);
p = x_out(:,4)*180/pi;
q = x_out(:,5)*180/pi;
r = x_out(:,6)*180/pi;
phi = x_out(:,7)*180/pi;
theta = x_out(:,8)*180/pi;
psi = x_out(:,9)*180/pi;
fig = figure('OuterPosition',[0 0 Pix_SS(3)*.90 Pix_SS(4)*.90],'PaperPositionMode','auto');
subplot(3,3,1)
plot(time_out,u,'LineWidth',3.0); title({'Axial Velocity','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Axial Velocity (ft/sec)'); grid on
subplot(3,3,2)
plot(time_out,v,'LineWidth',3.0); title({'Side Velocity','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Side Velocity (ft/sec)'); grid on
subplot(3,3,3)
plot(time_out,w,'LineWidth',3.0); title({'Normal Velocity','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Normal Velocity (ft/sec)'); grid on
subplot(3,3,4)
plot(time_out,p,'LineWidth',3.0); title({'Roll Rate','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Roll Rate (deg/sec)'); grid on
subplot(3,3,5)
plot(time_out,q,'LineWidth',3.0); title({'Pitch Rate','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Pitch Rate (deg/sec)'); grid on
subplot(3,3,6)
plot(time_out,r,'LineWidth',3.0); title({'Yaw Rate','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Yaw Rate (deg/sec)'); grid on
subplot(3,3,7)
plot(time_out,phi,'LineWidth',3.0); title({'Bank Angle','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Axial Velocity (deg)'); grid on
subplot(3,3,8)
plot(time_out,theta,'LineWidth',3.0); title({'Pitch Angle','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Pitch Angel (deg/sec)'); grid on
subplot(3,3,9)
plot(time_out,psi,'LineWidth',3.0); title({'Heading Angle','At Elevator Step -0.5 Deg'});xlabel('Time (sec)');ylabel('Heading Angle (deg)'); grid on

name = ['Figure_' num2str(fig.Number)];
print(fig,'-depsc','-noui','-painters',name);
fprintf(fid,'\\sectionmark{Project \\# 3\\hspace*{\\fill} Clair Cunningham \\hspace*{\\fill} Problem %d}\n',fig.Number);
fprintf(fid,'\\vspace*{\\fill}\\begin{figure}[H]\\centering\\includegraphics[keepaspectratio=true,height=1\\textheight,width=1\\textwidth,angle=90]{%s.eps}\n \\caption{Problem %d Graphical Solution}\\end{figure}\\vspace*{\\fill}\n\\newpage\n',name,fig.Number);
%print('-P\\meprint2\gle-2120-pr02c',figure(1))
fclose(fid);