##Nina Haas
##solves system of ordinary differential equations defined in F.m
##plots resulting approximations

function [y,t] = solve_ODE_system()
  NSTEP = 1000000;
  IOSTEP = 50;
  DT = 1e-3;
  y0 = [1 2 3]';
  T = DT*NSTEP;
  [y,t] = AB3(@F,y0,T,DT,IOSTEP);
  
  figure;
  plot(t,y(1,:));
  hold on;
  plot(t,y(2,:));
  plot(t,y(3,:));
  
  figure;
  plot3(y(1,:),y(2,:),y(3,:));