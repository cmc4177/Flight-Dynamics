%Project_3_dr_f.m
function dx = Project_3_dr_f(t,x)
% set up variables to incoming conditions
u = x(1); % ft/sec
v = x(2); % ft/sec
w = x(3); % ft/sec
p = x(4); % ft/sec
q = x(5); % ft/sec
r = x(6); % radians
phi = x(7); % radians
theta = x(8); % radians
psi = x(9); % radians

% % command troubleshooting
% u = u_0;
% v = v_0;
% w = w_0;
% p = p_0;
% q = q_0;
% r = r_0;
% phi = phi_0;
% theta = theta_0;
% psi = psi_0;

% Initial Trim Conditions
T = 3767.207337; % thrust in lbs
de = -2.9846046; % deg
da = 0; % deg
dr = 0; % deg

if t >= 1
    dr = dr + -2.0; % deg
end

% aircraft properties
i_xx = 8691.46164;
i_yy = 70668.585;
i_zz = 70418.67355;
i_xz = 151.43836;
i_xy = 0;
i_yz = 0;

s = 300; % ft^2
cbar = 11.32; % ft
b = 30; % ft
m = 756.5262463; % mass in slugs
rho = 0.0012669984; % slugs/ft^3
g = 32.17561865; % ft/sec^2

% other variables
qbar = 0.5*rho*u^2;

% intitial aerodynamic conditions; where except for "not" units (/deg);
c_x_0 = -2.13536e-02;
c_x_u = 1.289018e-04;
c_x_w = -2.17775e-03;
c_x_q = 2.1928052e-04;
c_x_de = 1.386632e-03;
c_y_0 = 0;
c_y_v = -6.4490425e-02;
c_y_p = 1.33481e-03;
c_y_r = 9.401418e-03;
c_y_da = 4.618436e-04;
c_y_dr = 2.991717e-03;
c_z_0 = 5.092263e-02;
c_z_u = -4.3444023e-04;
c_z_w = -1.9946051e-03;
c_z_q = -5.3473522e-02;
c_z_de = -1.2167892e-02;
c_l_0 = 0;
c_l_v = -1.75539e-03;
c_l_p = -7.392626e-03;
c_l_r = 5.910111e-05;
c_l_da = -2.089358e-03;
c_l_dr = 4.7651867e-04;
c_m_0 = -1.39985e-02;
c_m_u = -1.15335756e-04;
c_m_w = -1.16313463e-03;
c_m_q = -6.08086182e-01;
c_m_de = -4.632495451e-02;
c_n_0 = 0;
c_n_v = 5.1988574e-03;
c_n_p = -4.294548e-04;
c_n_r = -8.6047784e-03;
c_n_da = -1.95539e-04;
c_n_dr = -5.50282873e-03;

% aerodynamic model;  where u,v,w have units ft/sec;
c_x = c_x_0 + c_x_u*u + c_x_w*w + (cbar/(2*u))*c_x_q*q*(180/pi) + c_x_de*de;
c_y = c_y_0 + c_y_v*v + (b/(2*u))*(c_y_p*p + c_y_r*r)*(180/pi) + c_y_da*da + c_y_dr*dr;
c_z = c_z_0 + c_z_u*u + c_z_w*w + (cbar/(2*u))*c_z_q*q*(180/pi) + c_z_de*de;
c_l = c_l_0 + c_l_v*v + (b/(2*u))*(c_l_p*p + c_l_r*r)*(180/pi) + c_l_da*da + c_l_dr*dr;
c_m = c_m_0 + c_m_u*u + c_m_w*w + (cbar/(2*u))*c_m_q*q*(180/pi) + c_m_de*de;
c_n = c_n_0 + c_n_v*v + (b/(2*u))*(c_n_p*p + c_n_r*r)*(180/pi) + c_n_da*da + c_n_dr*dr;

% external forces
f_a_x = qbar*s*c_x;
f_a_y = qbar*s*c_y;
f_a_z = qbar*s*c_z;
f_T_x = T;
f_T_y = 0;
f_T_z = 0;

% external moments
m_e_x = qbar*s*b*c_l;
m_e_y = qbar*s*cbar*c_m;
m_e_z = qbar*s*b*c_n;

% Solves for the rate change of velocity from force equations
udot = -g*sin(theta)+(f_a_x+f_T_x)/m+v*r-w*q;
vdot = g*sin(phi)*cos(theta)+(f_a_y+f_T_y)/m+w*p-u*r;
wdot = g*cos(phi)*cos(theta)+(f_a_z+f_T_z)/m-v*p+u*q;
UVW = [udot;vdot;wdot];
% Solves for rate change of moments
I = [i_xx -i_xy -i_xz; -i_xy i_yy -i_yz; -i_xz -i_yz i_zz]; % The "inertia" matrix
B = [q*r*(i_yy-i_zz)+(q^2-r^2)*i_yz-p*r*i_xy+p*q*i_xz+m_e_x; p*r*(i_zz-i_xx)+(r^2-p^2)*i_xz-p*q*i_yz+q*r*i_xy+m_e_y; p*q*(i_xx-i_yy)+(p^2-q^2)*i_xy-q*r*i_xz+p*r*i_yz+m_e_z];
A =  inv(I)*B; % rate of change matrix
% Solves for rate change of angles
phidot = p+q*sin(phi)*tan(theta)+r*cos(phi)*tan(theta);
thetadot = q*cos(phi)-r*sin(phi);
psidot = (q*sin(phi)+r*cos(phi))*sec(theta);
POW = [phidot;thetadot;psidot];

% Spits out derivatives
dx = [UVW;A;POW]; % don't forget to add the semi-colon back after checking to make sure all rate changes are zero or close enough to the millionth
% pause
end