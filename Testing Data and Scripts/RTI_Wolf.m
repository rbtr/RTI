% Evan Baker
% WolfRTI
% 18 April 2014
function RTI_Wolf
clc
close all
% This function will calculate RTI inhibiting vibration parameters for two 
% different density fluids layered in a cylindrical vessel

% The stability condition is given by:
% bmax ~> (0.54*g*D)/(a*At)
% where
% bmax = a*w^2 [m/s^2]
% a = oscillation amplitude [cm]
% w = oscillation frequency [Hz]
% g = gravity [m/s^2]
% D = vessel diameter [cm]
% At = the Atwood number: []
%      (ph - pl)/(ph + pl)
%      where ph = heavier fluid density [g/cm^3]
%            pl = lighter fluid density [g/cm^3]

% Since bmax = a*w^2
% substitution gives:
% bmax^2 ~> (0.54*g*D*w^2)/(At)

f = 1:300; % generate frequency list
w = 2*pi.*f;
g = 980; % define gravity [m/s^2]
D = 1.27; % define diameter (.5 inches) [cm]
Atwood = @ (ph,pl) (ph - pl)/(ph + pl); % make the Atwood number

options = 0; % for list options
% build a decision list
fprintf('Use built in values or define new ones?\n1) Define new values\n2) SAE140 and air\n3) Water and air\n4) Potassium Jodide and SAE140\n') 
options = input('(1234): '); % get user decision
if options == 1;
    pl = input('pl = '); % fluid 1 density
    ph = input('ph = '); % fluid 2 density
elseif options == 2;
    pl = 0; % g/cm^3
    ph = 0.9; % air is negligible comparitively
elseif options == 3;
    pl = 0;
    ph = 1;
elseif options == 4;
    pl = 0.9; % SAE140 g/cm^3
    ph = 1.6; % potassium jodide g/cm^3
elseif options == 5;
    pl = 0; % [g/cm^3] Density - air
    ph = 1.54; % [g/cm^3] Density - honey
end

At = Atwood(ph,pl); % calculate the Atwood number

b = sqrt((0.54*g*D*w.^2)/(At))./g; % get relative b
a = sqrt((0.54*D*g)./(At.*w.^2));
figure(1)
plot(f,b);xlabel('Frequency [Hz]');ylabel('Acceleration [Gs]');%legend('Absolute acceleration','Relative acceleration (Gs)');
%figure(2)
%plot(w./(100*pi),a);xlabel('Frequency 100\pi [Hz]*pi');ylabel('Amplitude [cm]');
%figure(3)
%b2 = a.*w.^2./g;
%plot(w./(100*pi),b2);xlabel('Frequency 100\pi [Hz]');ylabel('Acceleration [Gs]');

fbtable = [f',b'];
fprintf('Frequency [Hz] |  Acceleration [Gs]')
fbtable

end