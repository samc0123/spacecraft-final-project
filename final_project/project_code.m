%% Introduction
% Spacecraft Mission Design 
% Final Project 
% Donald Barnickel, Steven Calalpa, Samuel Chernov, Daniella Chung 

%% Part 0: Constant Initialization
mu=3.986e14;
rE=6.378e6; % m, Earth's Radius
W=1000/2.205; %kg, mass of spacecraft
hp=600*1.609*1000; %m, height at perigee
rp=(600*1.609)*1000+rE; %m, radius at perigee
eps=0.8; % eccentricity
hH=(200*1.609)*1000; %m, Hohmann Height to moon 
vEsc=sqrt((2*mu)/rE); % m/s 

%% Part 1: Thrust and Flight Path Angle 
% Perigee Parameters
a=rp/(1-eps); %m, radius of apogee
E=-mu/(2*a); %m^2/s^2, energy of the orbit 
H=sqrt(mu*a*(1-eps^2)); % kg*m^2/sec, specific angular momentum 
V=sqrt(2*(E+mu/rp)); %m/s, velocity at perigee

% Flight path angle 
phi=acos(H/(rp*V));

% Energy on Earth
aE=rE; % radius of Earth is A
Ee=-mu/(2*aE); 

% Energy Difference 
ET=abs(E-Ee);
T=((ET)/(hp))*W;



