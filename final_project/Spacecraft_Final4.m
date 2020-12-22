mu1 = 3.986*10^14;
rE = 6.378*10^6;
Vesc = sqrt(2*mu1/rE);
EEarth = -mu1/(2*rE);

mass = 1000/2.205; %Mass in pounds to kg
altitude = 600*1.609*1000;
rP1 = rE+altitude;

epsilon1 = .8;
rA1 = (rP1*(-1-epsilon1))/(epsilon1-1);
check = rA1/9;
a1 = (rP1+rA1)/2;
H1 = sqrt(a1*mu1*(1-(epsilon1^2)));
E1 = -mu1/(2*a1);
VP1 = H1/rP1;
VA1 = H1/rA1;
Isp = 318; %Found from Beresheet rocket parameters

MoonToEarth = 356794*1000; %Meters from Earth
MoonAltitude = 200*1.609*1000; %Meters above Moon
r4 = (MoonAltitude+(1.738*(10^6)));

% deltaE = abs(EEarth) - abs(E1);
% LaunchThrustSpecific = deltaE/altitude;
% LaunchThrust = LaunchThrustSpecific*mass;
rA3 = MoonToEarth-r4;

%Create an array to vary Apogee positions
ArbitraryApogee = (rA1):(100*1000):(rA3);

BurnTime1 = 1:1:3600; %Seperated into Burn 1 and 2, in case you want to 
BurnTime2 = 1:1:3600; %burn the 3rd impulse for longer, seeing as it's a much larger deltaV

VDiff12 = [];
VDiff23 = [];
VDiff34 = [];
VTotal = [];
Mass0 = [];
MassP1 = [];
MassP2 = [];
MassP3 = [];
Compare = [];
T1 = [];
T2 = [];
T3 = [];



for i=1:length(ArbitraryApogee)
    
    %Depending on where you are leaving from, Perigee or Apogee, change rP2 =
    %                                            rP1 or rA1, respectively.
    rP2 = rP1;
    rA2 = ArbitraryApogee(i);
    a2 = (rP2+rA2)/2;
    epsilon2 = (rA2-rP2)/(rA2+rP2);
    E2 = -mu1/(2*a2);
    H2 = sqrt(a2*mu1*(1-(epsilon2^2)));
    VP2 = H2/rP2;
    VA2 = H2/rA2;

    %Depending on where you are leaving from, Perigee or Apogee, change rP2 =
    %                                            VP1 or VA1, respectively.
    VDiff12(i) = VP2 - VP1;

    rP3 = rP2;
    
    %If we assume we're entering at 6 'oclock, the do rA2 = MoonToEarth-r3
    %If we assume we're entering at 12 'oclock, the do rA2 = MoonToEarth+r3
    rA3 = MoonToEarth-r4;
    a3 = (rP3+rA3)/2;
    epsilon3 = (rA3-rP3)/(rA3+rP3);
    E3 = -mu1/(2*a3);
    H3 = sqrt(a3*mu1*(1-(epsilon3^2)));
    VP3 = H3/rP3;
    VA3 = H3/rA3;


    %Depending on where you are leaving from, Perigee or Apogee, change rP2 =
    %                                            VP2 or VA2, respectively.
    VDiff23(i) = VP3 - VP2;

    a4 = (r4+r4)/2;
    epsilon4 = (r4-r4)/(r4+r4);
    mu2 = 4.903*10^12;
    E4 = -mu2/(2*a4);
    H4 = sqrt(a4*mu2*(1-(epsilon4^2)));
    V4 = H4/r4;

    VDiff34(i) = V4 - VA3;
    VTotal(i) = VDiff12(i) + VDiff23(i) + VDiff34(i);

    Mass3 = mass;
    Mass0(i) = Mass3/((1-(1 - (exp((-VDiff12(i))/(9.81*Isp)))))*...
        (1-(1 - (exp((-VDiff23(i))/(9.81*Isp)))))*...
        (1-(1 - (exp((-VDiff34(i))/(9.81*Isp))))));


    % Mass0 = 2500/2.205;
    MassP1(i) = Mass0(i)*(1 - (exp((-VDiff12(i))/(9.81*Isp))));
    Mass1 = Mass0(i)-MassP1(i);
    MassP2(i) = Mass1*(1 - (exp((-VDiff23(i))/(9.81*Isp))));
    Mass2 = Mass1-MassP2(i);
    MassP3(i) = Mass2*(1 - (exp((-VDiff34(i))/(9.81*Isp))));
    Mass3 = Mass2-MassP3(i);
    
    Compare(i) = VDiff12(i)-VDiff23(i);
    
    if Compare(i)<=1 && Compare(i)>=-1
        MiddlePoint = ArbitraryApogee(i)/1000;
        MiddlePointIndex = i;
    end
    
    for j=1:length(BurnTime1)
        Mdot1 = MassP1(i)/BurnTime1(j);
        Mdot2 = MassP2(i)/BurnTime1(j);
        T1(j,i) = Mdot1*9.81*Isp;
        T2(j,i) = Mdot2*9.81*Isp;
    end
    for k=1:length(BurnTime2)
        Mdot3 = MassP3(i)/BurnTime2(k);
        T3(k,i) = Mdot3*9.81*Isp;
    end
%     if i==MiddlePointIndex
%         break
%     end
    
end
txt=['Epsilon is: ;',num2str(epsilon2)];
fprintf(txt)
% Plot 
    % Row- constant apogee, different time
    % Cols.- constant time, different apogee
figure 
hold on 
% Constant Burn Time, Different Apogee L
for i=1:5
    loc=i*500;
    y=T1(loc,:)'; 
    txt=['Burn Time: ',num2str(BurnTime1(loc)),' s'];
    plot(ArbitraryApogee./1000,y,'DisplayName',txt);
end
hold off
legend show
legend('Location','northwest');
ylabel('Thrust Force (N)','FontWeight','bold');
xlabel('Apogee Distance (km)','FontWeight','bold');
title('First Impulse Thrust Forces',...
    'FontWeight','bold');





