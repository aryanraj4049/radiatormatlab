data1 = xlsread('motor.xlsx','D:D');%power
data2 = xlsread('motor.xlsx','G:G');%elapsed time
data3 = xlsread('motor.xlsx','E:E');%efficiency
data4 = xlsread('motor.xlsx','A:A');%speed
%heatloss=xlsread('motor.xlsx','F:F');%heat loss
torque=xlsread('motor.xlsx','C:C');
speed = (data4*5)/18; % in m/s

%%Ploss= ((data1*746)/2).*(1-data3)+(0.06*746*data1)/2+3000;   % in W

%%Ploss = zeros(1, 5108);
%Ploss(i)=7810;
for i=1:5106
    Ploss(i)=(1-(data3(i)/100))*(data1(i)*746)+306+254;
end
Vwater= (10/60000) ;  % m^3/sec
Rlen = 0.3000;     % Radiator length
Rwid = 0.3000;     % Radiator width
InAirT = 30 ;  % in deg celsius
Dwater = 997 ; 
Twi = 60; %degreeC
Two = Twi + Ploss/(Vwater*Dwater*4.2*1000) ;
Co = 262.074;
C1 = 877.79;
Kr = 10;
Da = 1.225;   % Density of air
Vair = (data4*5)/18;   % velocity of air w.r.t. car (m/s)
A1 = Rlen*Rwid;
A4 = 0.0397;
N = (C1)^2 + (4*((Da*0.5)*((Kr/(A1^2))+(1/(A4^2))*(Co+(0.5*Da.*Vair.*Vair)))));
D = 0.5*Da*((Kr/(A1^2))+(1/(A4^2))) ;
Q = (-C1 + sqrt(N))./D ;
V1 = Q/A1 ;
V4 = Q/A4 ;
deltap = abs(Co-(C1*Q)) ;  %static pressure in fan
flowair = (262.075-deltap)/877.797; % from linear regression model of fan
Power = deltap.*flowair;    % power used by fan
%% Calculations
cpw = 4200 ; 
cpair = 1007 ;           % J/(kg-K)
Dwater = 997 ; 
Dair= 1.225;
uair= 18.6*10^(-6);      %% viscosity air
uwater = 8.9*10^(-4);    %% viscosity of water
kwater = 0.64;     %% at 60C
kair = 26.24*10^(-3) ;          %% conductivity of air (W/m K)
kalu = 205 ;
%%hi= 2000;    %% convective HTC of water (W/m2-k)
Rdep = 0.016000;     %% Radiator depth
FinH = 8.0*10^(-3) ;     %% Fin Height
FinSpc = 1.5*10^(-3);    %% Fin Spacing
FinTh = 0.4*10^(-3);    %% Fin thickness
FinL = sqrt((FinH)^2+(FinSpc/2)^2);
twall = 3*10^(-4);      %% thickness of tube wall
TubeTh = 2*10^(-3);     %% thickness of tube
TubeWid = 16*10^(-3) ;  %% Tube Width
TubeL = Rlen ;          %% Tube Length
NTubes = Rwid/(FinH+TubeTh) ;
NRows = NTubes + 1 ;
NFSR = (2*Rlen)/(FinSpc); %% No. of fins in a single row
TArea = Rlen*Rwid ;   %%Total area 
FinArea = FinL*FinTh ;
TFinArea = FinArea*NFSR*NRows ;      %% Total fin area;
TPipeArea = TubeTh*NTubes*Rlen ;         %% Total pipe area %%ee
RemArea = TArea - (TFinArea + TPipeArea) ;    %% Remaining Area
V2 = (V1 * TArea)/RemArea ;          %% Countinuity equation
Mwater = Dwater*Vwater ;
Mair = Dair*V1*TArea;       %% Mass flow rate of air
%% Internal Flow of water
NuWater = 8.24;                                   %% depends on a/b ratio & shape consult table
DhTube = (4*TubeWid*TubeTh)/(2*(TubeWid+TubeTh));
VelTube = Vwater/(NTubes*TubeWid*TubeTh) ;            %% velocity of water in the tube
ReTube =  (Dwater*VelTube*DhTube)/(uwater);           %% Re of water in the tube
hi = (NuWater*kwater)/DhTube ;                        %% inside heat transfer coeff.
%% External Flow of air
PrAir = (uair*cpair)/kair ;        %% Prandlt No. air
Dh = (4*0.5*FinH*FinSpc)/(FinSpc +(2*FinL)) ;    %% Hydraulic diameter
ReAir = (V2*Dh*Dair)/uair ;
NuAir = 0.664* (power(ReAir,0.5))*(power(PrAir,1/3)) ; 
ho = (NuAir*kair)/Dh ;                           %% air convective heat transfer coeff.
%% Fin efficiency
Af = 2*(FinH*(TubeWid+FinTh)) ;
Abase = Rlen*TubeWid ;
Lc = FinH ;
m = sqrt((2*ho)/(kalu*FinTh)) ;
Atot = (NFSR*Af)+ Abase ;
nfin = (tanh(m*Lc))./(m*Lc);
no = 1 - ((NFSR*Af/(Atot))*(1-nfin)) ;
%% Resistances
Rwall = twall/(TubeWid*TubeTh*kalu) ; 
Ri = 1/hi ;
Res = no.*ho;
Ro = 1./Res ;
%%battery cooling 
% Constants 
MaxMotorVoltage = 430 ;
maxRPM = 6500 ;
maxMotorCurrent = 200 ;
contMotorCurrent = 100 ;

torquePercurrent = 0.6;
controllerEfficiency = 0.97;

totalMotorCurrent = torque / (torquePercurrent * controllerEfficiency);
%plot(time,totalMotorCurrent)
%ylim([295,315])
Z =0.0012 * 79 / 3; % total impedance
heatGenerated = 0;  % during one lap
for i = 1:5105
    intervals(i) = data2(i+1) - data2(i) ;
    heat(i) = totalMotorCurrent(i) * totalMotorCurrent(i) * Z * intervals(i);
    heatGenerated = heatGenerated + heat(i);
end
heat(i+1) = heat(i);

%plot(time,heat)

%heatGenerated   during one lap         0.19 MJ 

totalHeatGenerated=heatGenerated * 22  ;  % 22 laps 

% total heat generated 4.27 MJ 
% which is just tbd % 
%time_onelap = time(end) ;   % 67.28 s
%radPower=heatGenerated / time_onelap;  % 2.88 kW 
%% NTU Method    %% Done for single layer( MIN no of layer = 2)
total = Ri+Ro;
UA = 1./total ;      %% Rwall missing
Cair = Mair*cpair ;
Cwater = Mwater*cpw ;
Cmin = min(Cair,Cwater) ;
Cmax = max(Cair,Cwater) ;
C = Cmin./Cmax ;
NTU = UA./Cmin ;
E = 1- exp((power(NTU,0.22)).*(1./C).*(exp(-C.*(power(NTU,0.78)))-1)) ;  %%effectiveness
Qmax = Cmin.*(Two-InAirT) ;
Qpredic = E.*Qmax ;
Tout = Two - Qpredic/(Vwater*Dwater*4.2*1000) ;
%plot(data2,Power)
%title('Power consumption of fan in a lap')
%xlabel('Time(sec)')
%ylabel('Power(w)') 
subplot(1,2,1)
plot(data2,Ploss)
title('Heat generated in a lap')
xlabel('Time(sec)')
ylabel('Heat generated(J)')
subplot(1, 2, 2)
plot(data2,Qpredic)
title('heat dissipated in a lap')
xlabel('Time(sec)')
ylabel('Heat dissipated(J)')
%subplot(1,2,1)
%plot(data2,InWaterT)
%title('Temperature of coolant going into radiator vs time')
%xlabel('Time(sec)')
%ylabel('Temperature(c)')
%subplot(1, 2, 2)
%plot(data2,Tout)
%title('Temperature of coolant coming out of radiator vs time')
%xlabel('Time(sec)')
%ylabel('Temperature (C)')
%int1 = trapz(data2,Ploss)/1000;
%int2 = trapz(data2,Qpredic)/1000;