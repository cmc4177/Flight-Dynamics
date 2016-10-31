%odex_ex1_m.m

clear all, close all, clc

% get solution
time = [0:0.01:10]'; % ' means transpose
IC = [0,0];
[time_out,x_out] = ode45(@ode45_ex1_f,time,IC);

% plot the results
figure(1), plot(time_out,x_out(:,1)),grid,xlabel('time (sec)'),ylabel('x1 (---)')
figure(2), plot(time_out,x_out(:,2)),grid,xlabel('time (sec)'),ylabel('x2 (---)')