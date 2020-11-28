load('data.mat')
format long

% 1 Heure 
% 2 Minute 
% 3 Seconde
% 4 Pression atmos mbar
% 5 Température  ambiante °C
% 6 Humidité ambiante %
% 7 T108Ap : Température de l’air à l’entrée de la chaudière °C
% 8 T200Gb :   Température   dans   la   chambre   de   combustion   au   niveau   du   1er
% thermocouple °C
% 9 T201Gb :   Température   dans   la   chambre   de   combustion   au   niveau   du   2ème
% thermocouple °C
% 10 T202Gb :   Température   dans   la   chambre   de   combustion   au   niveau   du   3ème
% thermocouple °C
% 11 T203GB : Température des gaz à la sorite de l’échangeur de chaleur eau­fumée
% °C
% 12 T300Wa : Température de l’eau  à la sortie de l’échangeur de chaleur eaufumée °C
% 13 T304Wa : Température de l’eau à l’entrée de l’échangeur de chaleur eau­fumée
% °C
% 14 T400Wa : Température de l’eau à la sortie de l’échangeur de chaleur eau­eau
% °C
% 15 T405Wa : Température de l’eau à l’entrée de l’échangeur de chaleur eau­eau
% °C
% 16      T406Wa :   Température de l’eau à l’entrée de l’échangeur autour de la chambre
% de combustion  °C
% 17 Débit d'air m³/h
% 18 Débit gaz m³/h
% 19 F303Wa : Débit dans le circuit de l’échangeur eau­fumée l/min
% 20 F401Wa : Débit d’alimentation d’eau extérieur l/min
% 21 P107Ap : Chute de pression au travers du diaphragme Pa
% 22 CO ppm 
% 23 CO2 % 
% 24 NO ppm 
% 25 O2 % 

d20 = mean(d20_raw); d24 = mean(d24_raw); d28 = mean(d28_raw);

figure(1)
plot(1:length(d28_raw(:,8)),d28_raw(:,8:10))
title('Temperature in combustion chamber at flowrate = 28m^3/s')
xlabel('Time [s]')
ylabel('Temperature [°C]')
grid on

figure(2)
plot(1:length(d24_raw(:,8)),d24_raw(:,8:10))
title('Temperature in combustion chamber at flowrate = 24m^3/s')
xlabel('Time [s]')
ylabel('Temperature [°C]')
grid on

figure(3)
plot(1:length(d20_raw(:,8)),d20_raw(:,8:10))
title('Temperature in combustion chamber at flowrate = 20m^3/s')
xlabel('Time [s]')
ylabel('Temperature [°C]')
grid on
