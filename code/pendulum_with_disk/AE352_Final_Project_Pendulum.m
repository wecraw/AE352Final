clear;clc;close all;
l=20;%pendulum length in meters
r=12;%disk radius in meters
g=9.8;%m/s^2
m=2000;%mass in kilograms of disk
J=m*l^2;%moment of inertia
F_tire=10000;%Newtons
L_tire=.5;%contact area of tire in meters
T_disk=1000;%torque on disk
J_disk=.5*m*r^2;%moment of inertia of disk in Newton-meters
oscillations=15;

t0=0;
theta0=0;%radians
thetadot0=.05;%radians/second (give a kick to start the simulation
phi0=0;
phidot0=1;%radians/second
tstep=0.05;
%output vectors (appended on each run to show total run of pendulum)
t_out=[];
y_out=[];
for i=1:oscillations
    if i==1
    thetamax=acos(1-thetadot0^2*l/(2*g)); %find max theta value (see notes for derivation)
    T=2*pi*sqrt(l/g)*(1+1/16*thetamax^2+11/3072*thetamax^4);%period

    infovec=[theta0 thetadot0 phi0 phidot0];
    [t, y]= ode45(@thetafunc, t0:tstep:(T/2),infovec);%go from t=0 to when pendulum completes...
    %half an oscillation (back at the bottom)
    t_out=t;
    y_out=y;
    tf=T/2;
    else
    %load in new information
    thetadot0=y(size(y,1),2);
    phi0=y(size(y,1),3);
    phidot0=y(size(y,1),4);
    
     v=thetadot0*l;%velocity at bottom
    t_contact=L_tire/v; %approximate time of contact between ride and tire
    torque_tire=F_tire*l;%torque exerted by the tire
    deltaH=torque_tire*t_contact;%change in angular momentum
    thetadotf=thetadot0+deltaH/J;%angular velocity after tire
    
    thetamax=acos(1-thetadotf^2*l/(2*g));
    T=2*pi*sqrt(l/g)*(1+1/16*thetamax^2+11/3072*thetamax^4);

    infovec=[theta0 thetadotf phi0 phidot0];
    [t, y]= ode45(@thetafunc, tf:tstep:(tf+T/2),infovec);%go from end of last run...
    %to when pendulum completes half an oscillation (back at the bottom)
    tf=tf+T/2;
    t_out=[t_out;t]
    y_out=[y_out;y]
    end
end
    

%outputs
thetaout=y_out(:,1).';
thetadotout=y_out(:,2).';
phiout=y_out(:,3).';
phidotout=y_out(:,4).';
R0to1=[ones(1,length(thetaout)) zeros(1,length(thetaout)) zeros(1,length(thetaout));...
    zeros(1,length(thetaout)) cos(thetaout) -sin(thetaout);...
    zeros(1,length(thetaout)) sin(thetaout) cos(thetaout)];
P1=[r*cos(phiout); r*sin(phiout); zeros(1,length(thetaout))];
P0=[];
P1=[];
P2=[];
Origin=[];
for i=1:length(thetaout)
    O_0to1=[0;l*sin(thetaout(i));-l*cos(thetaout(i))];
    R0to1=[1 0 0;...
    0 cos(thetaout(i)) -sin(thetaout(i));...
    0 sin(thetaout(i)) cos(thetaout(i))];
    P1=[r*cos(phiout(i)); r*sin(phiout(i)); 0];
    Origin=[Origin,O_0to1];
    P0=[P0,O_0to1+R0to1*P1];
    %^above outputs x,y,z in rows 1, 2, 3
end

%plot outputs
figure;
plot(t_out,y_out(:,1));
xlabel('Time, t, seconds')
ylabel('Angle, \theta, radians')
title('\Theta vs. t, Paul DeTrempe, AE 352 Pirate Ship Model')

%plot 3-Space
figure;
xlabel('x');
ylabel('y');
zlabel('z');
for j=1:size(P0,2)
person=plot3(P0(1,j),P0(2,j),P0(3,j),'ro');%point on disk
hold on;
O=[0 0 0];
pendulum=plot3([0, Origin(1,j)],[0,Origin(2,j)],[0,Origin(3,j)],'g');%line for pendulum
pause(.1);
hold on;
delete(pendulum);
end
