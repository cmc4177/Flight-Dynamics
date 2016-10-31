% Project_4_m.m
% Clair Cunningham

clear all; close all; clc
fid = fopen('Project_4.txt','w+');

set(0,'Units','pixels')
Pix_SS = get(0,'screensize');
%Position Center At 90 percent screen in vector form :[left bottom width height]
width = Pix_SS(3)*.90;
height = Pix_SS(4)*.90;
taskbar = 40;
left = (Pix_SS(3)-width-Pix_SS(1))/2;
bottom = (Pix_SS(4)-height-Pix_SS(2))/2+taskbar;
pos_size = [left bottom width height-taskbar];

%% Numbers
% Conversion factor from radians to degrees
r2d = 180/pi;

%Conversion factor from degrees to radians
d2r = pi/180;
pd2pr = 180/pi;

% Aircraft Properties
i_xx = 8890.63; %slug-ft^2
i_yy = 71973.5; %slug-ft^2
i_zz = 77141.1; %slug-ft^2
i_xz = 181.119; %slug-ft^2
i_xy = 0.0; %slug-ft^2
i_yz = 0.0; %slug-ft^2
I = [i_xx -i_xy -i_xz; -i_xy i_yy -i_yz; -i_xz -i_yz i_zz];
Iinv = inv(I);
S = 300; % ft^2
cbar = 11.31; % ft
b = 30; % ft
mass = 762.8447; % slugs


% Properties of Air
rho = 0.0014962376; % slugs/ft^2
g = 32.17561865; % ft/sec^2
% qbar = 0.5*rho*V_t^2; % lbf/ft^2

%+++++++++++++++++++++++++ Degree Format +++++++++++++++++++++++++
% Initial "Trim" Conditions
alpha_0_d = 3.6102915; % deg
beta_0_d = 0; % deg
p_0_d = 0; % deg
q_0_d = 0; % deg
r_0_d = 0; % deg
phi_0_d = 0; % deg
theta_0_d = 3.6102915; % deg
psi_0_d = 0; % deg
de_0_d = -3.03804303; % deg
da_0_d = 0; % deg
dr_0_d = 0; % deg
df_0_d = 1.5; % deg

% Coefficient of Lift
c_L_0 = 0.004608463; %      confirmed
c_L_alpha_d = 0.0794655; %1/deg
c_L_q_d = 0.0508476; %1/deg
c_L_alphadot_d = 0.0; %1/deg
c_L_de_d = 0.0121988; %1/deg
c_L_df_d = 0.0144389; %1/deg

% Coefficient of Drag
c_D_0 = 0.01192128;
c_D_alpha_d = 0.00550063; %1/deg
c_D_q_d = 0.00315057; %1/deg
c_D_alphadot_d = 0.0; %1/deg
c_D_de_d = -0.000587647; %1/deg
c_d_df_d = 0.00136385; %1/deg

% Coefficient of Side Force
c_y_0 = 0.0;
c_y_beta_d = -0.0219309; %1/deg
c_y_p_d = 0.00133787; %1/deg
c_y_r_d = 0.0094053; %1/deg
c_y_da_d = 0.00049355; %1/deg
c_y_dr_d = 0.00293048; %1/deg

% Coefficient of Rolling Moment
c_l_0 = 0.0;
c_l_beta_d = -0.00173748; %1/deg
c_l_p_d = -0.00739342; %1/deg
c_l_r_d = 0.0000699792; %1/deg
c_l_da_d = -0.00213984; %1/deg
c_l_dr_d = 0.000479021; %1/deg

% Coefficient of Pitching Moment
c_m_0 = -0.02092347;
c_m_alpha_d = -0.0041873; %1/deg
c_m_q_d = -0.110661; %1/deg
c_m_alphadot_d = 0.0; %1/deg
c_m_de_d = -0.0115767; %1/deg
c_m_df_d = 0.000580220; %1/deg

% Coefficient of Yawing Moment
c_n_0 = 0.0;
c_n_beta_d = 0.00320831; %1/deg
c_n_p_d = -0.000432575; %1/deg
c_n_r_d = -0.00886783; %1/deg
c_n_da_d = -0.000206591; %1/deg
c_n_dr_d = -0.00144865; %1/deg

%+++++++++++++++++++++++++ Radian Format +++++++++++++++++++++++++
% Initial "Trim" Conditions
V_t_0 = 626.81863; % ft/sec
T_0 = 3146.482666; % lb

% Conversion factor from radians to degrees
% r2d = 180/pi;

%Conversion factor from degrees to radians
% d2r = pi/180;
% pd2pr = 180/pi;

alpha_0 = alpha_0_d*d2r;% 
beta_0 = beta_0_d*d2r; % 
p_0 = p_0_d*d2r; % 
q_0 = q_0_d*d2r; % 
r_0 = r_0_d*d2r; % 
phi_0 = phi_0_d*d2r; % 
theta_0 = theta_0_d*d2r; % 
psi_0 = psi_0_d*d2r; % 
de_0 = de_0_d*d2r; % 
da_0 = da_0_d*d2r; %
dr_0 = dr_0_d*d2r; %
df_0 = df_0_d*d2r; % 

% Coefficient of Life
c_L_0 = 0.004608463;
c_L_alpha = c_L_alpha_d*pd2pr; %1/rad
c_L_q = c_L_q_d*pd2pr; %1/rad
c_L_alphadot = c_L_alphadot_d*pd2pr; %1/rad
c_L_de = c_L_de_d*pd2pr; %1/rad
c_L_df = c_L_df_d*pd2pr; %1/rad

% Coefficient of Drag
c_D_0 = 0.01192128;
c_D_alpha = c_D_alpha_d*pd2pr; %1/rad
c_D_q = c_D_q_d*pd2pr; %1/rad
c_D_alphadot = c_D_alphadot_d*pd2pr; %1/rad
c_D_de = c_D_de_d*pd2pr; %1/rad
c_D_df = c_d_df_d*pd2pr; %1/rad

% Coefficient of Side Force
c_y_0 = 0.0;
c_y_beta = c_y_beta_d*pd2pr; %1/rad
c_y_p = c_y_p_d*pd2pr; %1/rad
c_y_r = c_y_r_d*pd2pr; %1/rad
c_y_da = c_y_da_d*pd2pr; %1/rad
c_y_dr = c_y_dr_d*pd2pr; %1/rad

% Coefficient of Rolling Moment
c_l_0 = 0.0;
c_l_beta = c_l_beta_d*pd2pr; %1/rad
c_l_p = c_l_p_d*pd2pr; %1/rad
c_l_r = c_l_r_d*pd2pr; %1/rad
c_l_da = c_l_da_d*pd2pr; %1/rad
c_l_dr = c_l_dr_d*pd2pr; %1/rad

% Coefficient of Pitching Moment
c_m_0 = -0.02092347;
c_m_alpha = c_m_alpha_d*pd2pr; %1/rad
c_m_q = c_m_q_d*pd2pr; %1/rad
c_m_alphadot = c_m_alphadot_d*pd2pr; %1/rad
c_m_de = c_m_de_d*pd2pr; %1/rad
c_m_df = c_m_df_d*pd2pr; %1/rad

% Coefficient of Yawing Moment
c_n_0 = 0.0;
c_n_beta = c_n_beta_d*pd2pr; %1/rad
c_n_p = c_n_p_d*pd2pr; %1/rad
c_n_r = c_n_r_d*pd2pr; %1/rad
c_n_da = c_n_da_d*pd2pr; %1/rad
c_n_dr = c_n_dr_d*pd2pr; %1/rad

% Quaternion Initial Conditions
e1_0 = cos(psi_0/2)*cos(theta_0/2)*cos(phi_0/2)+sin(psi_0/2)*sin(theta_0/2)*sin(phi_0/2);
e2_0 = cos(psi_0/2)*cos(theta_0/2)*sin(phi_0/2)-sin(psi_0/2)*sin(theta_0/2)*cos(phi_0/2);
e3_0 = cos(psi_0/2)*sin(theta_0/2)*cos(phi_0/2)+sin(psi_0/2)*cos(theta_0/2)*sin(phi_0/2);
e4_0 = sin(psi_0/2)*cos(theta_0/2)*cos(phi_0/2)-cos(psi_0/2)*sin(theta_0/2)*sin(phi_0/2);

%% Trim Simulation
num = 1;
%Inputs to Elevator, Aileron, or Rudder
    % Degs
    de_d = 0;
    da_d = 0;
    df_d = 0;
    dr_d = 0;
    % Radians
    de = de_d*d2r;
    da = da_d*d2r;
    df = df_d*d2r;
    dr = dr_d*d2r;
%Run Simulation    
sim('Project_4_s')
%Call up next figure
fig = figure('OuterPosition',pos_size,'PaperPositionMode','auto');
fig.Name = 'Trim Simulation';
%Changes Paper Print Orientation to landscape on a per figure basic
orient landscape
%Output plots to called figure
subplot(3,3,1)
plot(time,V_t,'LineWidth',3.0); title({'Total Velocity','At Trim'});xlabel('Time (sec)');ylabel('Axial Velocity (ft/sec)'); grid on
subplot(3,3,2)
plot(time,alpha,'LineWidth',3.0); title({'Angle-of-Attack','At Trim'});xlabel('Time (sec)');ylabel('Angle-of-Attack (deg)'); grid on
subplot(3,3,3)
plot(time,beta,'LineWidth',3.0); title({'Sideslip Angle','At Trim'});xlabel('Time (sec)');ylabel('Sideslip Angle (deg)'); grid on
subplot(3,3,4)
plot(time,p,'LineWidth',3.0); title({'Roll Rate','At Trim'});xlabel('Time (sec)');ylabel('Roll Rate (deg/sec)'); grid on
subplot(3,3,5)
plot(time,q,'LineWidth',3.0); title({'Pitch Rate','At Trim'});xlabel('Time (sec)');ylabel('Pitch Rate (deg/sec)'); grid on
subplot(3,3,6)
plot(time,r,'LineWidth',3.0); title({'Yaw Rate','At Trim'});xlabel('Time (sec)');ylabel('Yaw Rate (deg/sec)'); grid on
subplot(3,3,7)
plot(time,phi,'LineWidth',3.0); title({'Bank Angle','At Trim'});xlabel('Time (sec)');ylabel('Axial Velocity (deg)'); grid on
subplot(3,3,8)
plot(time,theta,'LineWidth',3.0); title({'Pitch Angle','At Trim'});xlabel('Time (sec)');ylabel('Pitch Angle (deg)'); grid on
subplot(3,3,9)
plot(time,psi,'LineWidth',3.0); title({'Heading Angle','At Trim'});xlabel('Time (sec)');ylabel('Heading Angle (deg)'); grid on
% %'Paperposition', [left bottom width height]
  set(fig,'PaperPositionMode', 'manual', 'PaperUnits','Inches', 'Paperposition',[0.0 0.0 11 8.5]) 
%  print('-P\\meprint2\gle-2120-pr02c',fig)
name = [strrep(fig.Name,' ','_')];
print(fig,'-depsc','-noui','-painters',name);
fprintf(fid,'\\sectionmark{Project \\# 4\\hspace*{\\fill} Clair Cunningham \\hspace*{\\fill} Problem %d}\n',num);
fprintf(fid,'\\vspace*{\\fill}\\begin{figure}[H]\\centering\\includegraphics[keepaspectratio=true,height=1\\textheight,width=1\\textwidth,angle=90]{%s.eps}\n \\caption{%s Graphical Solution}\\end{figure}\\vspace*{\\fill}\n\\newpage\n',name,fig.Name);

%Derivative Figure and Plot
fig = figure('OuterPosition',pos_size,'PaperPositionMode','auto');
fig.Name = 'Trim Simulation Derivative';
%Changes Paper Print Orientation to landscape on a per figure basic
orient landscape
%Output derivative plots to called figure
subplot(3,3,1)
plot(time,V_tdot,'LineWidth',3.0); title({'Total Acceleration','At Trim'});xlabel('Time (sec)');ylabel('Acceleration (ft^2/sec)'); grid on
subplot(3,3,2)
plot(time,alphadot,'LineWidth',3.0); title({'Angle-of-Attack Rate','At Trim'});xlabel('Time (sec)');ylabel('Angle-of-Attack Rate (rad/sec)'); grid on
subplot(3,3,3)
plot(time,betadot,'LineWidth',3.0); title({'Sideslip Angle Rate','At Trim'});xlabel('Time (sec)');ylabel('Sideslip Angle Rate (rad/sec)'); grid on
subplot(3,3,4)
plot(time,pdot,'LineWidth',3.0); title({'Roll Acceleration','At Trim'});xlabel('Time (sec)');ylabel('Roll Acceleration (rad^2/sec)'); grid on
subplot(3,3,5)
plot(time,qdot,'LineWidth',3.0); title({'Pitch Acceleration','At Trim'});xlabel('Time (sec)');ylabel('Pitch Acceleration (rad^2/sec)'); grid on
subplot(3,3,6)
plot(time,rdot,'LineWidth',3.0); title({'Yaw Acceleration','At Trim'});xlabel('Time (sec)');ylabel('Yaw Acceleration (rad^2/sec)'); grid on
subplot(3,3,7)
plot(time,phidot,'LineWidth',3.0); title({'Bank Angle Rate','At Trim'});xlabel('Time (sec)');ylabel('Bank Angle Rate (rad/sec)'); grid on
subplot(3,3,8)
plot(time,thetadot,'LineWidth',3.0); title({'Pitch Angle Rate','At Trim'});xlabel('Time (sec)');ylabel('Pitch Angle Rate (rad/sec)'); grid on
subplot(3,3,9)
plot(time,psidot,'LineWidth',3.0); title({'Heading Angle Rate','At Trim'});xlabel('Time (sec)');ylabel('Heading Angle Rate (rad/sec)'); grid on
% %'Paperposition', [left bottom width height]
  set(fig,'PaperPositionMode', 'manual', 'PaperUnits','Inches', 'Paperposition',[0.0 0.0 11 8.5]) 
%  print('-P\\meprint2\gle-2120-pr02c',fig)
% name = [strrep(fig.Name,' ','_')];
% print(fig,'-depsc','-noui','-painters',name);
% fprintf(fid,'\\sectionmark{Project \\# 4\\hspace*{\\fill} Clair Cunningham \\hspace*{\\fill} Problem %d}\n',num);
% fprintf(fid,'\\vspace*{\\fill}\\begin{figure}[H]\\centering\\includegraphics[keepaspectratio=true,height=1\\textheight,width=1\\textwidth,angle=90]{%s.eps}\n \\caption{%s Graphical Solution}\\end{figure}\\vspace*{\\fill}\n\\newpage\n',name,fig.Name);

%Quaternion Figure and Plot
fig = figure('OuterPosition',pos_size,'PaperPositionMode','auto');
fig.Name = 'Trim Simulation Quaternion Model';
%Changes Paper Print Orientation to landscape on a per figure basic
orient landscape
%Output Euler Angles Derived from Quaternion Model
subplot(2,2,1)
plot(time,phi_quat,'LineWidth',3.0); title({'Bank Angle Quaternion','At Trim'});xlabel('Time (sec)');ylabel('Angle (deg)'); grid on
subplot(2,2,2)
plot(time,theta_quat,'LineWidth',3.0); title({'Pitch Angle Quaternion','At Trim'});xlabel('Time (sec)');ylabel('Angle (deg'); grid on
subplot(2,2,3)
plot(time,psi_quat,'LineWidth',3.0); title({'Heading Angle Quaternion','At Trim'});xlabel('Time (sec)');ylabel('Angle (deg)'); grid on
% %'Paperposition', [left bottom width height]
  set(fig,'PaperPositionMode', 'manual', 'PaperUnits','Inches', 'Paperposition',[0.0 0.0 11 8.5]) 
%  print('-P\\meprint2\gle-2120-pr02c',fig)
name = [strrep(fig.Name,' ','_')];
print(fig,'-depsc','-noui','-painters',name);
fprintf(fid,'\\sectionmark{Project \\# 4\\hspace*{\\fill} Clair Cunningham \\hspace*{\\fill} Problem %d}\n',num);
fprintf(fid,'\\vspace*{\\fill}\\begin{figure}[H]\\centering\\includegraphics[keepaspectratio=true,height=1\\textheight,width=1\\textwidth,angle=90]{%s.eps}\n \\caption{%s Graphical Solution}\\end{figure}\\vspace*{\\fill}\n\\newpage\n',name,fig.Name);
%% Elevator Input Simulation
num = num + 1;
%Inputs to Elevator, Aileron, or Rudder
    % Degs
    de_d = -0.5;
    da_d = 0;
    df_d = 0;
    dr_d = 0;
    % Radians
    de = de_d*d2r;
    da = da_d*d2r;
    df = df_d*d2r;
    dr = dr_d*d2r;

%Run Simulation    
sim('Project_4_s')
%Call up next figure
fig = figure('OuterPosition',pos_size,'PaperPositionMode','auto');
fig.Name = 'Elevator Input Simulation';
%Changes Paper Print Orientation to landscape on a per figure basic
orient landscape
%Output plots to called figure
subplot(3,3,1)
plot(time,V_t,'LineWidth',3.0); title({'Total Velocity','At Elevator Input'});xlabel('Time (sec)');ylabel('Axial Velocity (ft/sec)'); grid on
subplot(3,3,2)
plot(time,alpha,'LineWidth',3.0); title({'Angle-of-Attack','At Elevator Input'});xlabel('Time (sec)');ylabel('Angle-of-Attack (deg)'); grid on
subplot(3,3,3)
plot(time,beta,'LineWidth',3.0); title({'Sideslip Angle','At Elevator Input'});xlabel('Time (sec)');ylabel('Sideslip Angle (deg)'); grid on
subplot(3,3,4)
plot(time,p,'LineWidth',3.0); title({'Roll Rate','At Elevator Input'});xlabel('Time (sec)');ylabel('Roll Rate (deg/sec)'); grid on
subplot(3,3,5)
plot(time,q,'LineWidth',3.0); title({'Pitch Rate','At Elevator Input'});xlabel('Time (sec)');ylabel('Pitch Rate (deg/sec)'); grid on
subplot(3,3,6)
plot(time,r,'LineWidth',3.0); title({'Yaw Rate','At Elevator Input'});xlabel('Time (sec)');ylabel('Yaw Rate (deg/sec)'); grid on
subplot(3,3,7)
plot(time,phi,'LineWidth',3.0); title({'Bank Angle','At Elevator Input'});xlabel('Time (sec)');ylabel('Axial Velocity (deg)'); grid on
subplot(3,3,8)
plot(time,theta,'LineWidth',3.0); title({'Pitch Angle','At Elevator Input'});xlabel('Time (sec)');ylabel('Pitch Angle (deg)'); grid on
subplot(3,3,9)
plot(time,psi,'LineWidth',3.0); title({'Heading Angle','At Elevator Input'});xlabel('Time (sec)');ylabel('Heading Angle (deg)'); grid on
% %'Paperposition', [left bottom width height]
  set(fig,'PaperPositionMode', 'manual', 'PaperUnits','Inches', 'Paperposition',[0.0 0.0 11 8.5]) 
%  print('-P\\meprint2\gle-2120-pr02c',fig)
name = [strrep(fig.Name,' ','_')];
print(fig,'-depsc','-noui','-painters',name);
fprintf(fid,'\\sectionmark{Project \\# 4\\hspace*{\\fill} Clair Cunningham \\hspace*{\\fill} Problem %d}\n',num);
fprintf(fid,'\\vspace*{\\fill}\\begin{figure}[H]\\centering\\includegraphics[keepaspectratio=true,height=1\\textheight,width=1\\textwidth,angle=90]{%s.eps}\n \\caption{%s Graphical Solution}\\end{figure}\\vspace*{\\fill}\n\\newpage\n',name,fig.Name);

%Derivative Figure and Plot
fig = figure('OuterPosition',pos_size,'PaperPositionMode','auto');
fig.Name = 'Elevator Input Simulation Derivative';
%Changes Paper Print Orientation to landscape on a per figure basic
orient landscape
%Output derivative plots to called figure
subplot(3,3,1)
plot(time,V_tdot,'LineWidth',3.0); title({'Total Acceleration','At Elevator Input'});xlabel('Time (sec)');ylabel('Acceleration (ft^2/sec)'); grid on
subplot(3,3,2)
plot(time,alphadot,'LineWidth',3.0); title({'Angle-of-Attack Rate','At Elevator Input'});xlabel('Time (sec)');ylabel('Angle-of-Attack Rate (rad^2/sec)'); grid on
subplot(3,3,3)
plot(time,betadot,'LineWidth',3.0); title({'Sideslip Angle Rate','At Elevator Input'});xlabel('Time (sec)');ylabel('Sideslip Angle Rate (rad/sec)'); grid on
subplot(3,3,4)
plot(time,pdot,'LineWidth',3.0); title({'Roll Acceleration','At Elevator Input'});xlabel('Time (sec)');ylabel('Roll Acceleration (rad^2/sec)'); grid on
subplot(3,3,5)
plot(time,qdot,'LineWidth',3.0); title({'Pitch Acceleration','At Elevator Input'});xlabel('Time (sec)');ylabel('Pitch Acceleration (rad^2/sec)'); grid on
subplot(3,3,6)
plot(time,rdot,'LineWidth',3.0); title({'Yaw Acceleration','At Elevator Input'});xlabel('Time (sec)');ylabel('Yaw Acceleration (rad^2/sec)'); grid on
subplot(3,3,7)
plot(time,phidot,'LineWidth',3.0); title({'Bank Angle Rate','At Elevator Input'});xlabel('Time (sec)');ylabel('Bank Angle Acceleration (rad/sec)'); grid on
subplot(3,3,8)
plot(time,thetadot,'LineWidth',3.0); title({'Pitch Angle Rate','At Elevator Input'});xlabel('Time (sec)');ylabel('Pitch Angle Rate (rad/sec)'); grid on
subplot(3,3,9)
plot(time,psidot,'LineWidth',3.0); title({'Heading Angle Rate','At Elevator Input'});xlabel('Time (sec)');ylabel('Heading Angle Rate (rad/sec)'); grid on
% %'Paperposition', [left bottom width height]
  set(fig,'PaperPositionMode', 'manual', 'PaperUnits','Inches', 'Paperposition',[0.0 0.0 11 8.5]) 
%  print('-P\\meprint2\gle-2120-pr02c',fig)
% name = [strrep(fig.Name,' ','_')];
% print(fig,'-depsc','-noui','-painters',name);
% fprintf(fid,'\\sectionmark{Project \\# 4\\hspace*{\\fill} Clair Cunningham \\hspace*{\\fill} Problem %d}\n',num);
% fprintf(fid,'\\vspace*{\\fill}\\begin{figure}[H]\\centering\\includegraphics[keepaspectratio=true,height=1\\textheight,width=1\\textwidth,angle=90]{%s.eps}\n \\caption{%s Graphical Solution}\\end{figure}\\vspace*{\\fill}\n\\newpage\n',name,fig.Name);

%Quaternion Figure and Plot
fig = figure('OuterPosition',pos_size,'PaperPositionMode','auto');
fig.Name = 'Elevator Input Simulation Quaternion Model';
%Changes Paper Print Orientation to landscape on a per figure basic
orient landscape
%Output Euler Angles Derived from Quaternion Model
subplot(2,2,1)
plot(time,phi_quat,'LineWidth',3.0); title({'Bank Angle Quaternion','At Elevator Input'});xlabel('Time (sec)');ylabel('Angle (deg)'); grid on
subplot(2,2,2)
plot(time,theta_quat,'LineWidth',3.0); title({'Pitch Angle Quaternion','At Elevator Input'});xlabel('Time (sec)');ylabel('Angle (deg)'); grid on
subplot(2,2,3)
plot(time,psi_quat,'LineWidth',3.0); title({'Heading Angle Quaternion','At Elevator Input'});xlabel('Time (sec)');ylabel('Angle (deg)'); grid on
% %'Paperposition', [left bottom width height]
  set(fig,'PaperPositionMode', 'manual', 'PaperUnits','Inches', 'Paperposition',[0.0 0.0 11 8.5]) 
%  print('-P\\meprint2\gle-2120-pr02c',fig)
name = [strrep(fig.Name,' ','_')];
print(fig,'-depsc','-noui','-painters',name);
fprintf(fid,'\\sectionmark{Project \\# 4\\hspace*{\\fill} Clair Cunningham \\hspace*{\\fill} Problem %d}\n',num);
fprintf(fid,'\\vspace*{\\fill}\\begin{figure}[H]\\centering\\includegraphics[keepaspectratio=true,height=1\\textheight,width=1\\textwidth,angle=90]{%s.eps}\n \\caption{%s Graphical Solution}\\end{figure}\\vspace*{\\fill}\n\\newpage\n',name,fig.Name);
%% Aileron Input Simulation
num = num + 1;
%Inputs to Elevator, Aileron, or Rudder
    % Degs
    de_d = 0;
    da_d = -0.5;
    df_d = 0;
    dr_d = 0; % deg
    % Radians
    de = de_d*d2r;
    da = da_d*d2r;
    df = df_d*d2r;
    dr = dr_d*d2r;

%Run Simulation    
sim('Project_4_s')
%Call up next figure
fig = figure('OuterPosition',pos_size,'PaperPositionMode','auto');
fig.Name = 'Aileron Input Simulation';
%Changes Paper Print Orientation to landscape on a per figure basic
orient landscape
%Output plots to called figure
subplot(3,3,1)
plot(time,V_t,'LineWidth',3.0); title({'Total Velocity','At Aileron Input'});xlabel('Time (sec)');ylabel('Axial Velocity (ft/sec)'); grid on
subplot(3,3,2)
plot(time,alpha,'LineWidth',3.0); title({'Angle-of-Attack','At Aileron Input'});xlabel('Time (sec)');ylabel('Angle-of-Attack (deg)'); grid on
subplot(3,3,3)
plot(time,beta,'LineWidth',3.0); title({'Sideslip Angle','At Aileron Input'});xlabel('Time (sec)');ylabel('Sideslip Angle (deg)'); grid on
subplot(3,3,4)
plot(time,p,'LineWidth',3.0); title({'Roll Rate','At Aileron Input'});xlabel('Time (sec)');ylabel('Roll Rate (deg/sec)'); grid on
subplot(3,3,5)
plot(time,q,'LineWidth',3.0); title({'Pitch Rate','At Aileron Input'});xlabel('Time (sec)');ylabel('Pitch Rate (deg/sec)'); grid on
subplot(3,3,6)
plot(time,r,'LineWidth',3.0); title({'Yaw Rate','At Aileron Input'});xlabel('Time (sec)');ylabel('Yaw Rate (deg/sec)'); grid on
subplot(3,3,7)
plot(time,phi,'LineWidth',3.0); title({'Bank Angle','At Aileron Input'});xlabel('Time (sec)');ylabel('Axial Velocity (deg)'); grid on
subplot(3,3,8)
plot(time,theta,'LineWidth',3.0); title({'Pitch Angle','At Aileron Input'});xlabel('Time (sec)');ylabel('Pitch Angle (deg)'); grid on
subplot(3,3,9)
plot(time,psi,'LineWidth',3.0); title({'Heading Angle','At Aileron Input'});xlabel('Time (sec)');ylabel('Heading Angle (deg)'); grid on
% %'Paperposition', [left bottom width height]
  set(fig,'PaperPositionMode', 'manual', 'PaperUnits','Inches', 'Paperposition',[0.0 0.0 11 8.5]) 
%  print('-P\\meprint2\gle-2120-pr02c',fig)
name = [strrep(fig.Name,' ','_')];
print(fig,'-depsc','-noui','-painters',name);
fprintf(fid,'\\sectionmark{Project \\# 4\\hspace*{\\fill} Clair Cunningham \\hspace*{\\fill} Problem %d}\n',num);
fprintf(fid,'\\vspace*{\\fill}\\begin{figure}[H]\\centering\\includegraphics[keepaspectratio=true,height=1\\textheight,width=1\\textwidth,angle=90]{%s.eps}\n \\caption{%s Graphical Solution}\\end{figure}\\vspace*{\\fill}\n\\newpage\n',name,fig.Name);

%Derivative Figure and Plot
fig = figure('OuterPosition',pos_size,'PaperPositionMode','auto');
fig.Name = 'Aileron Input Simulation Derivative';
%Changes Paper Print Orientation to landscape on a per figure basic
orient landscape
%Output derivative plots to called figure
subplot(3,3,1)
plot(time,V_tdot,'LineWidth',3.0); title({'Total Acceleration','At Aileron Input'});xlabel('Time (sec)');ylabel('Acceleration (ft^2/sec)'); grid on
subplot(3,3,2)
plot(time,alphadot,'LineWidth',3.0); title({'Angle-of-Attack Rate','At Aileron Input'});xlabel('Time (sec)');ylabel('Angle-of-Attack Rate (rad^2/sec)'); grid on
subplot(3,3,3)
plot(time,betadot,'LineWidth',3.0); title({'Sideslip Angle Rate','At Aileron Input'});xlabel('Time (sec)');ylabel('Sideslip Angle Rate (rad/sec)'); grid on
subplot(3,3,4)
plot(time,pdot,'LineWidth',3.0); title({'Roll Acceleration','At Aileron Input'});xlabel('Time (sec)');ylabel('Roll Acceleration (rad^2/sec)'); grid on
subplot(3,3,5)
plot(time,qdot,'LineWidth',3.0); title({'Pitch Acceleration','At Aileron Input'});xlabel('Time (sec)');ylabel('Pitch Acceleration (rad^2/sec)'); grid on
subplot(3,3,6)
plot(time,rdot,'LineWidth',3.0); title({'Yaw Acceleration','At Aileron Input'});xlabel('Time (sec)');ylabel('Yaw Acceleration (rad^2/sec)'); grid on
subplot(3,3,7)
plot(time,phidot,'LineWidth',3.0); title({'Bank Angle Rate','At Aileron Input'});xlabel('Time (sec)');ylabel('Bank Angle Acceleration (rad/sec)'); grid on
subplot(3,3,8)
plot(time,thetadot,'LineWidth',3.0); title({'Pitch Angle Rate','At Aileron Input'});xlabel('Time (sec)');ylabel('Pitch Angle Rate (rad/sec)'); grid on
subplot(3,3,9)
plot(time,psidot,'LineWidth',3.0); title({'Heading Angle Rate','At Aileron Input'});xlabel('Time (sec)');ylabel('Heading Angle Rate (rad/sec)'); grid on
% %'Paperposition', [left bottom width height]
  set(fig,'PaperPositionMode', 'manual', 'PaperUnits','Inches', 'Paperposition',[0.0 0.0 11 8.5]) 
%  print('-P\\meprint2\gle-2120-pr02c',fig)
name = [strrep(fig.Name,' ','_')];
print(fig,'-depsc','-noui','-painters',name);
fprintf(fid,'\\sectionmark{Project \\# 4\\hspace*{\\fill} Clair Cunningham \\hspace*{\\fill} Problem %d}\n',num);
fprintf(fid,'\\vspace*{\\fill}\\begin{figure}[H]\\centering\\includegraphics[keepaspectratio=true,height=1\\textheight,width=1\\textwidth,angle=90]{%s.eps}\n \\caption{%s Graphical Solution}\\end{figure}\\vspace*{\\fill}\n\\newpage\n',name,fig.Name);

%Quaternion Figure and Plot
fig = figure('OuterPosition',pos_size,'PaperPositionMode','auto');
fig.Name = 'Aileron Input Simulation Quaternion Model';
%Changes Paper Print Orientation to landscape on a per figure basic
orient landscape
%Output Euler Angles Derived from Quaternion Model
subplot(2,2,1)
plot(time,phi_quat,'LineWidth',3.0); title({'Bank Angle Quaternion','At Aileron Input'});xlabel('Time (sec)');ylabel('Angle (deg)'); grid on
subplot(2,2,2)
plot(time,theta_quat,'LineWidth',3.0); title({'Pitch Angle Quaternion','At Aileron Input'});xlabel('Time (sec)');ylabel('Angle (deg)'); grid on
subplot(2,2,3)
plot(time,psi_quat,'LineWidth',3.0); title({'Heading Angle Quaternion','At Aileron Input'});xlabel('Time (sec)');ylabel('Angle (deg)'); grid on
% %'Paperposition', [left bottom width height]
  set(fig,'PaperPositionMode', 'manual', 'PaperUnits','Inches', 'Paperposition',[0.0 0.0 11 8.5]) 
%  print('-P\\meprint2\gle-2120-pr02c',fig)
name = [strrep(fig.Name,' ','_')];
print(fig,'-depsc','-noui','-painters',name);
fprintf(fid,'\\sectionmark{Project \\# 4\\hspace*{\\fill} Clair Cunningham \\hspace*{\\fill} Problem %d}\n',num);
fprintf(fid,'\\vspace*{\\fill}\\begin{figure}[H]\\centering\\includegraphics[keepaspectratio=true,height=1\\textheight,width=1\\textwidth,angle=90]{%s.eps}\n \\caption{%s Graphical Solution}\\end{figure}\\vspace*{\\fill}\n\\newpage\n',name,fig.Name);

% %% Flap Input Simulation
% num = num + 1;
% %Inputs to Elevator, Aileron, or Rudder
%     % Degs
%     de_d = 0;
%     da_d = 0;
%     df_d = 0;
%     dr_d = 0;
%     % Radians
%     de = de_d*d2r;
%     da = da_d*d2r;
%     df = df_d*d2r;
%     dr = dr_d*d2r;
% 
% %Run Simulation    
% sim('Project_4_s')
% %Call up next figure
% fig = figure('OuterPosition',pos_size,'PaperPositionMode','auto');
% fig.Name = 'Flap Input Simulation';
% %Changes Paper Print Orientation to landscape on a per figure basic
% orient landscape
% %Output plots to called figure
% subplot(3,3,1)
% plot(time,V_t,'LineWidth',3.0); title({'Total Velocity','At Flap Input'});xlabel('Time (sec)');ylabel('Axial Velocity (ft/sec)'); grid on
% subplot(3,3,2)
% plot(time,alpha,'LineWidth',3.0); title({'Angle-of-Attack','At Flap Input'});xlabel('Time (sec)');ylabel('Angle-of-Attack (deg)'); grid on
% subplot(3,3,3)
% plot(time,beta,'LineWidth',3.0); title({'Sideslip Angle','At Flap Input'});xlabel('Time (sec)');ylabel('Sideslip Angle (deg)'); grid on
% subplot(3,3,4)
% plot(time,p,'LineWidth',3.0); title({'Roll Rate','At Flap Input'});xlabel('Time (sec)');ylabel('Roll Rate (deg/sec)'); grid on
% subplot(3,3,5)
% plot(time,q,'LineWidth',3.0); title({'Pitch Rate','At Flap Input'});xlabel('Time (sec)');ylabel('Pitch Rate (deg/sec)'); grid on
% subplot(3,3,6)
% plot(time,r,'LineWidth',3.0); title({'Yaw Rate','At Flap Input'});xlabel('Time (sec)');ylabel('Yaw Rate (deg/sec)'); grid on
% subplot(3,3,7)
% plot(time,phi,'LineWidth',3.0); title({'Bank Angle','At Flap Input'});xlabel('Time (sec)');ylabel('Axial Velocity (deg)'); grid on
% subplot(3,3,8)
% plot(time,theta,'LineWidth',3.0); title({'Pitch Angle','At Flap Input'});xlabel('Time (sec)');ylabel('Pitch Angle (deg)'); grid on
% subplot(3,3,9)
% plot(time,psi,'LineWidth',3.0); title({'Heading Angle','At Flap Input'});xlabel('Time (sec)');ylabel('Heading Angle (deg)'); grid on
% % %'Paperposition', [left bottom width height]
% set(fig,'PaperPositionMode', 'manual', 'PaperUnits','Inches', 'Paperposition',[0.0 0.0 11 8.5]) 
% %  print('-P\\meprint2\gle-2120-pr02c',fig)
% name = [strrep(fig.Name,' ','_')];
% print(fig,'-depsc','-noui','-painters',name);
% fprintf(fid,'\\sectionmark{Project \\# 4\\hspace*{\\fill} Clair Cunningham \\hspace*{\\fill} Problem %d}\n',num);
% fprintf(fid,'\\vspace*{\\fill}\\begin{figure}[H]\\centering\\includegraphics[keepaspectratio=true,height=1\\textheight,width=1\\textwidth,angle=90]{%s.eps}\n \\caption{%s Graphical Solution}\\end{figure}\\vspace*{\\fill}\n\\newpage\n',name,fig.Name);
% 
% %Derivative Figure and Plot
% fig = figure('OuterPosition',pos_size,'PaperPositionMode','auto');
% fig.Name = 'Flap Input Simulation Derivative';
% %Changes Paper Print Orientation to landscape on a per figure basic
% orient landscape
% %Output derivative plots to called figure
% subplot(3,3,1)
% plot(time,V_tdot,'LineWidth',3.0); title({'Total Acceleration','At Flap Input'});xlabel('Time (sec)');ylabel('Acceleration (ft^2/sec)'); grid on
% subplot(3,3,2)
% plot(time,alphadot,'LineWidth',3.0); title({'Angle-of-Attack Rate','At Flap Input'});xlabel('Time (sec)');ylabel('Angle-of-Attack Rate (rad^2/sec)'); grid on
% subplot(3,3,3)
% plot(time,betadot,'LineWidth',3.0); title({'Sideslip Angle Rate','At Flap Input'});xlabel('Time (sec)');ylabel('Sideslip Angle Rate (rad/sec)'); grid on
% subplot(3,3,4)
% plot(time,pdot,'LineWidth',3.0); title({'Roll Acceleration','At Flap Input'});xlabel('Time (sec)');ylabel('Roll Acceleration (rad^2/sec)'); grid on
% subplot(3,3,5)
% plot(time,qdot,'LineWidth',3.0); title({'Pitch Acceleration','At Flap Input'});xlabel('Time (sec)');ylabel('Pitch Acceleration (rad^2/sec)'); grid on
% subplot(3,3,6)
% plot(time,rdot,'LineWidth',3.0); title({'Yaw Acceleration','At Flap Input'});xlabel('Time (sec)');ylabel('Yaw Acceleration (rad^2/sec)'); grid on
% subplot(3,3,7)
% plot(time,phidot,'LineWidth',3.0); title({'Bank Angle Rate','At Flap Input'});xlabel('Time (sec)');ylabel('Bank Angle Acceleration (rad/sec)'); grid on
% subplot(3,3,8)
% plot(time,thetadot,'LineWidth',3.0); title({'Pitch Angle Rate','At Flap Input'});xlabel('Time (sec)');ylabel('Pitch Angle Rate (rad/sec)'); grid on
% subplot(3,3,9)
% plot(time,psidot,'LineWidth',3.0); title({'Heading Angle Rate','At Flap Input'});xlabel('Time (sec)');ylabel('Heading Angle Rate (rad/sec)'); grid on
% %'Paperposition', [left bottom width height]
%  set(fig,'PaperPositionMode', 'manual', 'PaperUnits','Inches', 'Paperposition',[0.0 0.0 11 8.5]) 
% %  print('-P\\meprint2\gle-2120-pr02c',fig)
% name = [strrep(fig.Name,' ','_')];
% print(fig,'-depsc','-noui','-painters',name);
% fprintf(fid,'\\sectionmark{Project \\# 4\\hspace*{\\fill} Clair Cunningham \\hspace*{\\fill} Problem %d}\n',num);
% fprintf(fid,'\\vspace*{\\fill}\\begin{figure}[H]\\centering\\includegraphics[keepaspectratio=true,height=1\\textheight,width=1\\textwidth,angle=90]{%s.eps}\n \\caption{%s Graphical Solution}\\end{figure}\\vspace*{\\fill}\n\\newpage\n',name,fig.Name);
%  
% %Quaternion Figure and Plot
% fig = figure('OuterPosition',pos_size,'PaperPositionMode','auto');
% fig.Name = 'Flap Input Simulation Quaternion Model';
% %Changes Paper Print Orientation to landscape on a per figure basic
% orient landscape
% %Output Euler Angles Derived from Quaternion Model
% subplot(2,2,1)
% plot(time,phi_quat,'LineWidth',3.0); title({'Bank Angle Quaternion','At Flap Input'});xlabel('Time (sec)');ylabel('Angle (deg)'); grid on
% subplot(2,2,2)
% plot(time,theta_quat,'LineWidth',3.0); title({'Pitch Angle Quaternion','At Flap Input'});xlabel('Time (sec)');ylabel('Angle (deg)'); grid on
% subplot(2,2,3)
% plot(time,psi_quat,'LineWidth',3.0); title({'Heading Angle Quaternion','At Flap Input'});xlabel('Time (sec)');ylabel('Angle (deg)'); grid on
% % %'Paperposition', [left bottom width height]
% %  set(fig,'PaperPositionMode', 'manual', 'PaperUnits','Inches', 'Paperposition',[0.0 0.0 11 8.5]) 
% %  print('-P\\meprint2\gle-2120-pr02c',fig)
% name = [strrep(fig.Name,' ','_')];
% print(fig,'-depsc','-noui','-painters',name);
% fprintf(fid,'\\sectionmark{Project \\# 4\\hspace*{\\fill} Clair Cunningham \\hspace*{\\fill} Problem %d}\n',num);
% fprintf(fid,'\\vspace*{\\fill}\\begin{figure}[H]\\centering\\includegraphics[keepaspectratio=true,height=1\\textheight,width=1\\textwidth,angle=90]{%s.eps}\n \\caption{%s Graphical Solution}\\end{figure}\\vspace*{\\fill}\n\\newpage\n',name,fig.Name);
%% Rudder Input Simulation
num = num +1;
%Inputs to Elevator, Aileron, or Rudder
    % Degs
    de_d = 0;
    da_d = 0;
    df_d = 0;
    dr_d = -2.0; % deg
    % Radians
    de = de_d*d2r;
    da = da_d*d2r;
    df = df_d*d2r;
    dr = dr_d*d2r;

%Run Simulation    
sim('Project_4_s')
%Call up next figure
fig = figure('OuterPosition',pos_size,'PaperPositionMode','auto');
fig.Name = 'Rudder Input Simulation';
%Changes Paper Print Orientation to landscape on a per figure basic
orient landscape
%Output plots to called figure
subplot(3,3,1)
plot(time,V_t,'LineWidth',3.0); title({'Total Velocity','At Rudder Input'});xlabel('Time (sec)');ylabel('Axial Velocity (ft/sec)'); grid on
subplot(3,3,2)
plot(time,alpha,'LineWidth',3.0); title({'Angle-of-Attack','At Rudder Input'});xlabel('Time (sec)');ylabel('Angle-of-Attack (deg)'); grid on
subplot(3,3,3)
plot(time,beta,'LineWidth',3.0); title({'Sideslip Angle','At Rudder Input'});xlabel('Time (sec)');ylabel('Sideslip Angle (deg)'); grid on
subplot(3,3,4)
plot(time,p,'LineWidth',3.0); title({'Roll Rate','At Rudder Input'});xlabel('Time (sec)');ylabel('Roll Rate (deg/sec)'); grid on
subplot(3,3,5)
plot(time,q,'LineWidth',3.0); title({'Pitch Rate','At Rudder Input'});xlabel('Time (sec)');ylabel('Pitch Rate (deg/sec)'); grid on
subplot(3,3,6)
plot(time,r,'LineWidth',3.0); title({'Yaw Rate','At Rudder Input'});xlabel('Time (sec)');ylabel('Yaw Rate (deg/sec)'); grid on
subplot(3,3,7)
plot(time,phi,'LineWidth',3.0); title({'Bank Angle','At Rudder Input'});xlabel('Time (sec)');ylabel('Axial Velocity (deg)'); grid on
subplot(3,3,8)
plot(time,theta,'LineWidth',3.0); title({'Pitch Angle','At Rudder Input'});xlabel('Time (sec)');ylabel('Pitch Angle (deg)'); grid on
subplot(3,3,9)
plot(time,psi,'LineWidth',3.0); title({'Heading Angle','At Rudder Input'});xlabel('Time (sec)');ylabel('Heading Angle (deg)'); grid on
% %'Paperposition', [left bottom width height]
  set(fig,'PaperPositionMode', 'manual', 'PaperUnits','Inches', 'Paperposition',[0.0 0.0 11 8.5]) 
%  print('-P\\meprint2\gle-2120-pr02c',fig)
name = [strrep(fig.Name,' ','_')];
print(fig,'-depsc','-noui','-painters',name);
fprintf(fid,'\\sectionmark{Project \\# 4\\hspace*{\\fill} Clair Cunningham \\hspace*{\\fill} Problem %d}\n',num);
fprintf(fid,'\\vspace*{\\fill}\\begin{figure}[H]\\centering\\includegraphics[keepaspectratio=true,height=1\\textheight,width=1\\textwidth,angle=90]{%s.eps}\n \\caption{%s Graphical Solution}\\end{figure}\\vspace*{\\fill}\n\\newpage\n',name,fig.Name);

%Derivative Figure and Plot
fig = figure('OuterPosition',pos_size,'PaperPositionMode','auto');
fig.Name = 'Rudder Input Simulation Derivative';
%Changes Paper Print Orientation to landscape on a per figure basic
orient landscape
%Output derivative plots to called figure
subplot(3,3,1)
plot(time,V_tdot,'LineWidth',3.0); title({'Total Acceleration','At Rudder Input'});xlabel('Time (sec)');ylabel('Acceleration (ft^2/sec)'); grid on
subplot(3,3,2)
plot(time,alphadot,'LineWidth',3.0); title({'Angle-of-Attack Rate','At Rudder Input'});xlabel('Time (sec)');ylabel('Angle-of-Attack Rate (rad^2/sec)'); grid on
subplot(3,3,3)
plot(time,betadot,'LineWidth',3.0); title({'Sideslip Angle Rate','At Rudder Input'});xlabel('Time (sec)');ylabel('Sideslip Angle Rate (rad/sec)'); grid on
subplot(3,3,4)
plot(time,pdot,'LineWidth',3.0); title({'Roll Acceleration','At Rudder Input'});xlabel('Time (sec)');ylabel('Roll Acceleration (rad^2/sec)'); grid on
subplot(3,3,5)
plot(time,qdot,'LineWidth',3.0); title({'Pitch Acceleration','At Rudder Input'});xlabel('Time (sec)');ylabel('Pitch Acceleration (rad^2/sec)'); grid on
subplot(3,3,6)
plot(time,rdot,'LineWidth',3.0); title({'Yaw Acceleration','At Rudder Input'});xlabel('Time (sec)');ylabel('Yaw Acceleration (rad^2/sec)'); grid on
subplot(3,3,7)
plot(time,phidot,'LineWidth',3.0); title({'Bank Angle Rate','At Rudder Input'});xlabel('Time (sec)');ylabel('Bank Angle Acceleration (rad/sec)'); grid on
subplot(3,3,8)
plot(time,thetadot,'LineWidth',3.0); title({'Pitch Angle Rate','At Rudder Input'});xlabel('Time (sec)');ylabel('Pitch Angle Rate (rad/sec)'); grid on
subplot(3,3,9)
plot(time,psidot,'LineWidth',3.0); title({'Heading Angle Rate','At Rudder Input'});xlabel('Time (sec)');ylabel('Heading Angle Rate (rad/sec)'); grid on
% %'Paperposition', [left bottom width height]
  set(fig,'PaperPositionMode', 'manual', 'PaperUnits','Inches', 'Paperposition',[0.0 0.0 11 8.5]) 
%  print('-P\\meprint2\gle-2120-pr02c',fig)
% name = [strrep(fig.Name,' ','_')];
% print(fig,'-depsc','-noui','-painters',name);
% fprintf(fid,'\\sectionmark{Project \\# 4\\hspace*{\\fill} Clair Cunningham \\hspace*{\\fill} Problem %d}\n',num);
% fprintf(fid,'\\vspace*{\\fill}\\begin{figure}[H]\\centering\\includegraphics[keepaspectratio=true,height=1\\textheight,width=1\\textwidth,angle=90]{%s.eps}\n \\caption{%s Graphical Solution}\\end{figure}\\vspace*{\\fill}\n\\newpage\n',name,fig.Name);

%Quaternion Figure and Plot
fig = figure('OuterPosition',pos_size,'PaperPositionMode','auto');
fig.Name = 'Rudder Input Simulation Quaternion Model';
%Changes Paper Print Orientation to landscape on a per figure basic
orient landscape
%Output Euler Angles Derived from Quaternion Model
subplot(2,2,1)
plot(time,phi_quat,'LineWidth',3.0); title({'Bank Angle Quaternion','At Rudder Input'});xlabel('Time (sec)');ylabel('Angle (deg)'); grid on
subplot(2,2,2)
plot(time,theta_quat,'LineWidth',3.0); title({'Pitch Angle Quaternion','At Rudder Input'});xlabel('Time (sec)');ylabel('Angle (deg)'); grid on
subplot(2,2,3)
plot(time,psi_quat,'LineWidth',3.0); title({'Heading Angle Quaternion','At Rudder Input'});xlabel('Time (sec)');ylabel('Angle (deg)'); grid on

% %'Paperposition', [left bottom width height]
  set(fig,'PaperPositionMode', 'manual', 'PaperUnits','Inches', 'Paperposition',[0.0 0.0 11 8.5]) 
%  print('-P\\meprint2\gle-2120-pr02c',fig)
 name = [strrep(fig.Name,' ','_')];
print(fig,'-depsc','-noui','-painters',name);
fprintf(fid,'\\sectionmark{Project \\# 4\\hspace*{\\fill} Clair Cunningham \\hspace*{\\fill} Problem %d}\n',num);
fprintf(fid,'\\vspace*{\\fill}\\begin{figure}[H]\\centering\\includegraphics[keepaspectratio=true,height=1\\textheight,width=1\\textwidth,angle=90]{%s.eps}\n \\caption{%s Graphical Solution}\\end{figure}\\vspace*{\\fill}\n\\newpage\n',name,fig.Name);

%% End of File Commands
name = 'Project_4_s';
load_system(name);
%modelhandle = get_param('name', 'Handle')
handles = find_system(name, 'FindAll', 'On', 'SearchDepth', 10, ...
    'regexp', 'on', 'blocktype', 'port');
list = get(handles,'Path');
if ~iscell(list)
    list = {list};
end
list = unique(list);
if ~iscell(handles)
    handle = {handles};
end
% add main model
list{end+1} = name;
% GEt only last part of path, that is, after the last /
[r1, r2] = regexpi(list, '[^/]+$', 'tokens', 'match');

% Convert to usable format
names = [r2{:}]';

% Cells of printNames.
% Just rename every non-alphanumeric char to _, all space to ''
printNames = regexprep(names', {'\s', '\W'}, {'', '_'});

for i = 1 : length(list)
    item = char(list(i));
modelhandle = get_param(item, 'Handle');
set(modelhandle,'PaperPositionMode', 'manual', 'PaperUnits','Inches', 'Paperposition',[0.0 0.0 11 8.5])
print(['-s' item],'-depsc','-noui','-painters',printNames{i});
fprintf(fid,'\\sectionmark{Project \\# 4\\hspace*{\\fill} Clair Cunningham \\hspace*{\\fill} %s}\n',strrep(printNames{i},'_','\_'));
fprintf(fid,'\\vspace*{\\fill}\\begin{figure}[H]\\centering\\includegraphics[keepaspectratio=true,height=1\\textheight,width=1\\textwidth,angle=90]{%s.eps}\n \\caption{Simulink Diagram}\\end{figure}\\vspace*{\\fill}\n\\newpage\n',printNames{i});
end
fclose(fid);