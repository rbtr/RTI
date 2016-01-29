function RTI_UCLM

% Define constants and anonymous functions
% p is used as rho
% u is used as mu
% o is used as sigma

g = 981; % [cm/s^2] gravity
diameter = 1.27; %1.6; % [cm] Using UCLM diameter (ours is 1.27)
km = 3.68/diameter; % [dimensionless] minimum wave number for a finite container of diameter

atwood = @ (pl,ph) (ph-pl)/(ph+pl); % [dimensionless] calculate the Atwood number
kZero = @ (at,ph,ul,uh) (((2*at^(1/2))/(1+at))*((ph*g^(1/2))/(ul+uh)))^(2/3); % calculate the kzero
omegaZero = @ (at,ph,ul,uh) (((2*at^2)/(1+at))*((ph*g^2)/(ul+uh)))^(1/3); % calculate the omegazero
ampAccelMod = @ (w,a) (w^2)*(a); % calculate b, the amplitude of acceleration modulation
defD = @ (at,ul,uh,ph,w,k) ((1+at)*(ul+uh)*(k^2))/(2*ph*w); % calculate D 
defK = @ (at,o,k,ph,w) (((1+at)*o*(k^3))/(2*ph*(w^2)))-((at*k*g)/(w^2)); % calculate K 
defkc = @ (at,ph,o) (((2*at)/(1+at))*((ph*g/o)))^(2/3); % calculate kc

% Begin script proper
fprintf('Use built in values or define new ones?\n1) Define new values\n2) SAE140 and air\n3) WH01 and air\n4) Water and air\n') % build a decision list
options = input('(1234): '); % get user decision
if options == 1;
    fprintf('The denser fluid is designated h\n');
    pl = input('pl = '); % lighter fluid density
    ph = input('ph = '); % heavier fluid density
    ul = input('ul = '); % lighter fluid viscosity
    uh = input('uh = '); % heavier fluid viscosity
    o = input('sigma = '); % inteface surface tension - don't know how to calculate this
elseif options == 2;
    pl = 0; % [g/cm^3] Density - air
    ph = 0.9; % [g/cm^3] Density - SAE140
    ul = 0; % [g/cm s] Viscosity - air
    uh = 9; % [g/cm s] Viscosity - SAE140
    o = 32; % [dyn/cm] interface surface tension of SAE140
elseif options == 3;
    pl = 0; % [g/cm^3] Density - air
    ph = 1.2; % [g/cm^3] Density - WH01
    ul = 0; % [g/cm s] Viscosity - air
    uh = 12; % [g/cm s] Viscosity - WH01
    o = 62; % [dyn/cm] interface surface tension of WH01
elseif options == 4;
    pl = 0; % [g/cm^3] Density - air
    ph = 1; % [g/cm^3] Density - water
    ul = 0; % [g/cm s] Viscosity - air
    uh = .0089; % [g/cm s] Viscosity - water
    o = 72.8; % [dyn/cm] interface surface tension of water
elseif options == 5;
    pl = 0; % [g/cm^3] Density - air
    ph = 1.54; % [g/cm^3] Density - honey
    ul = 0; % [g/cm s] Viscosity - air
    uh = 20; % [g/cm s] Viscosity - honey
    o = 62; % [dyn/cm] interface surface tension of honey?
end

At = atwood(pl,ph); % Get the Atwood number for chosen fluids
k0 = kZero(At,ph,ul,uh); % call the kzero 
omega0 = omegaZero(At,ph,ul,uh); % call the omegazero
kc = defkc(At,ph,o);

Kappac = kc/k0; % Get KappaC
Kappam = km/k0; % Get KappaM

omegaBarT = ((2*pi)/(3*sqrt(3)))*Kappac*Kappam; % Get critical omegabar


%%%%%%%%
f = 1:300; % Frequency range
omega = 2*pi.*f; % Calculate omegas for the frequency range
omegaBar = omega./omega0; % Calculate omegaBar from omega0 and the omega list
b_over_g = zeros(1,length(omegaBar)); % Preallocate the vector for b/g

greaterThanCritCount = 0;%These Var's keep track of the number of points in
lessThanCritCount = 0;   %the linear/non-linear regions. They are are output.


for i = 1:length(omegaBar)      % SEE EQUATION 19 IN UCLM PRIMARY PAPER. THIS                          
    if omegaBar(i) > omegaBarT  % LOOP DETERMINES b/g VALUES BASED ON OMEGABAR
        b_over_g(i) = (sqrt(2)/sqrt(Kappam))*omegaBar(i);
        greaterThanCritCount = greaterThanCritCount + 1;
    elseif omegaBar(i) < omegaBarT
        b_over_g(i) = sqrt((8*pi/sqrt(27)))*sqrt(Kappac)*sqrt(omegaBar(i));
        lessThanCritCount = lessThanCritCount + 1;
    end
end
greaterThanCritCount %Output of counts for above/below critical omegaBar
lessThanCritCount

% THESE GRAPHS ARE AN ATTEMPT TO DUPLICATE FIG 7 IN UCLM
x = omegaBar/sqrt(Kappam);
y = b_over_g;
figure(3)
plot(x,y);xlabel('omegaBar/sqrt(Kappam)');ylabel('b/g');title('UCLM primary (valid tesing points only)')

b_over_g2 = sqrt((8*pi/sqrt(27)))*sqrt(Kappac).*sqrt(omegaBar);%b/g determined entirely using omegaBar < omegaBarT
b_over_g3 = (sqrt(2)/sqrt(Kappam)).*omegaBar;%b/g determined entirely using omegaBar > omegaBarT

%figure(4)
%plot(x,b_over_g2,'*',x,b_over_g3,'r');xlabel('omegaBar/sqrt(Kappam)');ylabel('b/g');title('UCLM primary (comaprative)')
%legend('function of sqrt(omega)','Linear Region')
%%%%%%%%%

fbtable = [f',b_over_g'];
fprintf('Frequency [Hz] |  Acceleration [Gs]')
fbtable
figure(5)
plot(f,b_over_g);xlabel('Frequency [Hz]');ylabel('Acceleration [Gs]');
end