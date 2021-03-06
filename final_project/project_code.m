%% Introduction
% Spacecraft Mission Design 
% Final Project 
% Donald Barnickel, Steven Calalpa, Samuel Chernov, Daniella Chung 

%% Part 0: Constant Initialization
mu=3.986e14; % Earth's gravitational parameter 
rE=6.378e6; % m, Earth's Radius
Wlb=1000; % lb, mass of spacecraft
W=Wlb/2.205; %kg, mass of spacecraft
hp=600*1.609*1000; %m, height at perigee
rp=(600*1.609)*1000+rE; %m, radius at perigee
eps=0.8; % eccentricity
hm=(200*1.609)*1000; %m, Hohmann Height to moon 
rmoon=1738e3; % m, radius of moon 
rEM=356794e3; %m, distance from earth to moon 
vEsc=sqrt((2*mu)/rE); % m/s 

%% Part 1: Thrust and Flight Path Angle 
% Perigee Parameters
a=rp/(1-eps); %m, radius of apogee
E=-mu/(2*a); %m^2/s^2, energy of the orbit 
H=sqrt(mu*a*(1-eps^2)); % kg*m^2/sec, specific angular momentum 
V=sqrt(2*(E+mu/rp)); %m/s, velocity at perigee

% Flight path angle 
Rreq=60000*1000; % m, radius @ enterance
Vreq=10000; % m/s, velocity @ enterance
phi=acos(H/(Rreq*Vreq));

% Thrust & Propellant Mass
mPay=W+796.7; % kg, mass of payload- payload & booster fuel 
Tmp=thrust(mPay); % separate function 
T=Tmp(1); % N, Falcon 9 Thrust
mp1=Tmp(2); % kg, Falcon 9 Stage 1 Prop. Mass
mp2=Tmp(3); % kg, Falcon 9 Stage 2 Prop. Mass

%% Part 2: Hohmann Transfer 
% Transfer from 600 mile perigee to 
% 200 miles away from the moon 

% Properties of leaving orbit 
ra=a*(1+eps); %m, transfer pt. from current orbit
va=H/ra; % m/s, circularized orbit velocity 

% Radii from Earth 
% ra- already achieved prior -> dist. of first orbit to earth 
rmp=rEM-rmoon; %m, Earth's center to moon surface
rmo=rmp+hm+2*rmoon;%(ADD for 12 o'clock) %m, 200 mi. away from the moons orbit, closer to earth 

% Transfer Orbit Properties 
at=(ra+rmo)/2; %m, semi-major axis of the Hohmann transfer orbit 
Et=-mu/(2*at); % energy of the transfer orbit 
vDep=sqrt(2*(Et+(mu/ra))); % m/s, velocity of leaving orbit 
Ht=ra*vDep; % specific angular momentum of the transfer orbit 
vArr=Ht/rmo; % m/s, arrival velocity 
epsT=(rmo-ra)/(rmo+ra); % epsilon of the transfer orbit

% Final Orbit properties 
muM=4.90339e12; % gravitational constant of the moon 
vf=sqrt(muM/(rmoon+hm)); % velocity of the circular orbit 

% Velocity Change 
dv1=vDep-va; % Velocity change between circular and elliptical transfer
dv2=vf-vArr; % Vel. chng btwn elliptical and final circular
dvt=abs(dv1)+abs(dv2); % total velocity change 

%% Orbit Plot 
% Sketch this orbit 
figure 
hold on

% First, plot the Earth
center=[0 0]; % center point of the graph 
colorRGB=[0 0 1]; % color to fill in the Earth
p1e=viscircles(center,rE,'Color',colorRGB,'LineWidth',3); % Create circles with Earth's radius

% Plot moon
center=[rEM 0];
colorRGB=[0 1 0];
p1m=viscircles(center,rmoon,'Color',colorRGB,'LineWidth',3); % create circles with moon's radius 

% Global theta
theta=0:0.05:2*pi; % define theta range in small increments

% Plot initial orbit 
orbXi=a*(cos(theta)-eps); % take the x component of orbit, Kepler's Laws
orbYi=a*sqrt(1-eps^2)*sin(theta); % y component of orbit, Kepler's Laws
p1i=plot(orbXi,orbYi,'c--','LineWidth',2); % plot the initial orbit

% Plot transfer orbit
orbXt=-at*(cos(theta)-epsT); % take the x component of orbit, Kepler's Laws
orbYt=-at*sqrt(1-epsT^2)*sin(theta); % y component of orbit, Kepler's Laws
p1t=plot(orbXt,orbYt,'m--','LineWidth',2); % plot the initial orbit

% Plot final circular orbit 
colorRGB=[1 0 0];
cirR=rmoon+hm; 
p1f=viscircles(center,cirR,'Color',colorRGB,'LineWidth',2,'LineStyle','--');

% Set the Earth to be the origin 
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

% Plot labels 
title('Orbit Visualization');
ylabel('Orbit Distance (m)');
xlabel('Orbit Distance (m)');
legend([p1e p1m p1i p1t p1f],'Earth','Moon','Initial Orbit','Transfer Orbit','Final Orbit');
hold off