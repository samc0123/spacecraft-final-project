function [T] = thrust(mPay)
%% Introduction Information
%{
* Calculate the thrust produced by the rocket 
* Taking parameters of the SpaceX Falcon 9 Launcher
     - Comparable payload weights 
     - Comparable orbit lift off 
* Thurst Eqn: F = q × Ve + (Pe - Pa) × Ae
     - F= Thrust
     - q= Propellant Mass Flow Rate
     - Ve= Velocity of exhaust gas
     - Pe= Pressure @ exit Nozzle 
     - Pa= Ambient Pressure 
     - Ae= Area of nozzle exit 
* Alternatively 
     - Isp=F/mdot*g0
     - F=Isp*mdot*g0
     - Only need the Isp's of the rocket 

%}

% Alternate Way
%% Constants 
g0=9.81; %m/sec^2, gravitational constant

tb1=180; %s, first burn time
tb2=346; %s, second burn time

%% Rocket Properties 
tkg=907.185; % ton to kg 
ms1=28*tkg; % stage 1 F9 mass
ms2=4.7*tkg; % stage 2 F9 mass
% Payload mass inputted- mp
IspF91=282; % sec, sea level
IspF92=340; % sec, vaccuum 

%% Propellant Properties
% Stage 1
mp1max=411*tkg; %kg, max mass of Merlin 1D
mp1min=0.4*mp1max; % 40% of maxmum mass
mp1=linspace(mp1min,mp1max,40);

% Stage 2
mp2max=73.4*tkg; %kg, max mass of Merlin 1D Vac.
mp2min=0.39*mp2max; % 39% of the max as specified by props.
mp2=linspace(mp2min,mp2max,40); 

% Totals
Vt=zeros(1,40);
Vtcrit=10000; % Velocity req'd to obtain
VtcritE=Vtcrit*0.005; % allowable error
mp1f=0;
mp2f=0;

%% Velocity
for i=1:length(mp1)
    tempMP1=mp1(i);
    tempMP2=mp2(i);


% Stage 1: Merlin 1D
    lnFrac1=(tempMP1+ms1+tempMP2+ms2+mPay)/(ms1+tempMP2+ms2+mPay);
    Vb1=g0*IspF91*log(lnFrac1); % Obtain the velocity of first stg.

% Stage 2: Merlin 1D- Vaccuum
    lnFrac2=(tempMP2+ms2+mPay)/(ms2+mPay);
    Vb2=g0*IspF92*log(lnFrac2); % velocity of scnd stg

    Vt(i)=Vb1+Vb2; % total rocket velocity
    
% Check Critical cond, select those values for mass
    
    if (Vt(i)-Vtcrit)<VtcritE
        mp1f=tempMP1;
        mp2f=tempMP2;
        ic=i; 
    end
end
% Plot of possible Propellant mass configs
figure 
yyaxis right 
plot(Vt,mp2,'LineWidth',2);
hold on
plot(Vt(ic),mp2f, 'c*', 'LineWidth',3);
hold on
ylabel('Second Stage Mass (kg)');
yyaxis left
plot(Vt,mp1,'LineWidth',2);
hold on
plot(Vt(ic),mp1f,'m*','LineWidth',3);
hold on
ylabel('First Stage Mass (kg)');
xlabel('Velocity (m/s)');
title('Propellant Mass vs. Velocity');
legend('Prop. M1','Cri. Pt. M1','Prop. M2','Cri. Pt. M2',...
    'Location','northwest');
hold off

%% Total Thurst 

% Mass Flow Rate 
mdot1=mp1f/tb1;
mdot2=mp2f/tb2;

% Thrust 
T1=mdot1*IspF91*g0; % First Thrust 
T2=mdot2*IspF92*g0; % Second Thrust 

T=[T1+T2 mp1f mp2f]; % total req'd thrust

return 



end

