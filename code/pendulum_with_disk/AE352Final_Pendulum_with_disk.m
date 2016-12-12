clear;clc;close all;

Numpeople=40;


for k=1:Numpeople
    
phi0=(k-1)*2*pi/Numpeople;%disperse people evenly along disk
    
l=20.11;%pendulum length in meters
r=4.11;%disk radius in meters
g=9.8;%m/s^2
m=2000;%mass in kilograms of disk
J=m*l^2;%moment of inertia
F=5000;%Newtons
L_tire=.5;%contact area of tire in meters
T_disk=1000;%torque on disk
J_disk=.5*m*r^2;%moment of inertia of disk in Newton-meters
oscillations=20;


t0=0;
theta0=0;%radians
thetadot0=.05;%radians/second (give a kick to start the simulation)
phidot0=1;%radians/second (spin rate of disk)
tstep=0.1;

%output vectors (appended on each run to show total run of pendulum)
t_out=[];
y_out=[];

plot_outputs=zeros(3,706,Numpeople);%initialize plot outputs

for i=1:oscillations
    if i<oscillations/2
        F_tire=i*F;
    else
        F_tire=-F*(i-oscillations/2);
    end
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
    t_out=[t_out;t];
    y_out=[y_out;y];
    end
end
    

%outputs
thetaout=y_out(:,1).';
thetadotout=y_out(:,2).';
phiout=y_out(:,3).';
phidotout=y_out(:,4).';
%initialize rotation matrix
R0to1=[ones(1,length(thetaout)) zeros(1,length(thetaout)) zeros(1,length(thetaout));...
    zeros(1,length(thetaout)) cos(thetaout) -sin(thetaout);...
    zeros(1,length(thetaout)) sin(thetaout) cos(thetaout)];
%initialize storage vectors for 2 frames of reference
P1=[r*cos(phiout); r*sin(phiout); zeros(1,length(thetaout))];
P0=[];
P1=[];
P1dot=[];
omega=[];
P2=[];
Origin=[];
P0dot=[];
P0doubledot=[];
%for each output (timestep), use rotation matrices and formulas from class
%to calculate position, velocity, and acceleration for each point in time
for i=1:length(thetaout)
    O_0to1=[0;l*sin(thetaout(i));-l*cos(thetaout(i))];%end of pendulum position in base frame
    O_0to1dot=[0;...
        thetadotout(i)*l*cos(thetaout(i));...
        thetadotout(i)*l*sin(thetaout(i))];%movement of pendulum in base frame
    %acceleration of pendulum in base frame
    O_0to1doubledot=[0;...
        -g/l*sin(thetaout(i))*l*cos(thetaout(i))+thetadotout(i)*l*(-sin(thetaout(i)));...
        -g/l*sin(thetaout(i))*l*sin(thetaout(i))+thetadotout(i)*l*cos(thetaout(i))];%acceleration of pendulum in base frame
    %rotation matrix to base frame
    R0to1=[1 0 0;...
        0 cos(thetaout(i)) -sin(thetaout(i));...
        0 sin(thetaout(i)) cos(thetaout(i))];
    %movement in disk frame
    P1=[r*cos(phiout(i)); r*sin(phiout(i)); 0];
    omega=[0;0;phidotout(i)];%disk spin rate
    omegahat=[0 -omega(3) omega(2); omega(3) 0 -omega(1);-omega(2) omega(1) 0];%wedge matrix for disk spinning
    P1dot=[-r*omega(3)*sin(phiout(i));r*omega(3)*cos(phiout(i));0];%velocity in disk frame
    P1doubledot=[-r*omega(3)^2*cos(phiout(i));-r*omega(3)^2*sin(phiout(i));0];%acceleration in disk frame
    
    
    Origin=[Origin,O_0to1];%store different origin values (used from plotting pendulum motion)
    
    P0=[P0,O_0to1+R0to1*P1];%store position in base frame
    P0dot=[P0dot,O_0to1dot+R0to1*omegahat*P1+R0to1*P1dot];%store velocity in base frame
    P0doubledot=[P0doubledot,O_0to1doubledot+R0to1*omegahat^2*P1+2*R0to1*omegahat*P1dot+R0to1*P1doubledot];%store acceleration in base frame
    %^above outputs x,y,z in rows 1, 2, 3
end
if k==1
    %initialize position, velocity and acceleration arrays
   plouts=zeros(3,size(P0,2),Numpeople);
   vouts=zeros(3,size(P0dot,2),Numpeople);
   aouts=zeros(3,size(P0doubledot,2),Numpeople);
end
%store position, velocity, acceleration, and pendulum position
plouts(:,:,k)=P0;
vouts(:,:,k)=P0dot;
aouts(:,:,k)=P0doubledot;
OriginArray(:,:,k)=Origin;

if k<Numpeople
    clearvars -except k plouts vouts aouts Numpeople
end

end

%plot outputs
figure;
plot(t_out,y_out(:,1));
xlabel('Time, t, seconds')
ylabel('Angle, \theta, radians')
title('\Theta vs. t, Paul DeTrempe, AE 352 Pirate Ship Model')

%calculate velocity and acceleration for each rider
figure;
for ridernumber=1:Numpeople
    vplot=[];
    aplot=[];
    for ridetime=1:size(vouts,2)
        vplot=[vplot,norm(vouts(:,ridetime,ridernumber))];
        aplot=[aplot,norm(aouts(:,ridetime,ridernumber))];
    end
    plot(t_out,aplot);
    hold on;
end
xlabel('time, t, seconds');
ylabel('acceleration magnitude, a, m/s^2');

% plot 3-Space


figure;
xlabel('x');
ylabel('y');
zlabel('z');

for j=1:size(plouts,2)
axis([-25 25 -25 25 -25 25])    
for k=1:Numpeople
    people(k)=plot3(plouts(1,j,k),plouts(2,j,k),plouts(3,j,k),'ro');
end
pendulum=plot3([0, Origin(1,j)],[0,Origin(2,j)],[0,Origin(3,j)],'g');
pause(.01);
hold on;
delete(people)
delete(pendulum)
end
