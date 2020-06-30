% ParticleSimulation.m by Yuanzheng Wen,wenyuanzheng@stu.cdut.edu.cn
% Simulation of charged particles motion in Earth magnetic field
% This is a matlab file
% but you can also change it into a function file

% 2020-06-01

close all; clear; clc; 
global c m rel Bfield
global B0 Re;
rel=0;
% Parameters for dipole (Earth) and double dipole (Earth) B field
e = 1.602176565e-19; % Elementary charge (Coulomb)
m_pr = 1.672621777e-27; % Proton mass (kg)
m_el = 9.10938291e-31; % Electron mass (kg)
c = 299792458; % speed of light (m/s)
Re = 6378137; % meter (Earth radius)

% Defintions for diople B field
Bdipole=@(x,y,z,t)[3*x*z*(-B0*Re^3/power(x*x+y*y+z*z,5.0/2.0)),...
3*y*z*(-B0*Re^3/power(x*x+y*y+z*z,5.0/2.0)),...
(2*z*z-x*x-y*y)*(-B0*Re^3/power(x*x+y*y+z*z,5.0/2.0))];
q=e;m=m_pr;
B0=3.07e-5; % Tesla
Bfield=@(x,y,z,t)Bdipole(x,y,z,t);

% Initial and Calculating
K=1e7; % kinetic energy in eV
x0=4*Re; y0=0*Re; z0=0*Re;
pitch_angle=30.0; % initial angle between velocity and mag.field (degrees)
K=K*e  ; % convert to Joule
v=c/sqrt(1+(m_pr*c^2)/K); % replace m_pr with m_el for electron  
vz0=v*cos(pitch_angle*pi/180);
vy0=v*sin(pitch_angle*pi/180); vx0=0;
gamma=1.0/sqrt(1-(vx0^2+vy0^2+vz0^2)/c^2);
T=2*pi*gamma*m/(abs(q)*sqrt(sum(Bfield(x0,y0,z0,0).^2)));
dt=T/10;
tend=1000*T; % end time

% Solving equation
yy0=[x0,y0,z0,vx0,vy0,vz0];
options = odeset('RelTol',1e-4); % 17-04-07 13:38
[t,y]=ode45('SolveNewtonLorenz',[0:dt:tend],yy0,options);

% Plotting
plot3(x0,y0,z0,'o'); hold on; 
plot3(y(:,1),y(:,2),y(:,3),'r');
xlabel('x');ylabel('y');zlabel('z');
colormap(jet)
cla;
[X,Y,Z]=sphere(20);surf(X,Y,Z);hold on; grid on;
plot3(y(:,1)./Re,y(:,2)./Re,y(:,3)./Re,'color','blue','linewidth',2.5); title('Dipole (Earth)');
xlabel('x[R_e]');ylabel('y[R_e]');zlabel('z[R_e]');
amax=2*sqrt(x0^2+y0^2+z0^2)/Re;amin=-amax;
axis equal;
