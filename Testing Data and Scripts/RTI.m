function [omega, A, b] = RTI

% Define constants and anonymous functions
% omega = input('omega = ');


g = 9.81; % [m/s^2] | gravity
km = 3.68/1.27; % minimum wave number for a finite container of diameter D = 1.27 cm (.5 in)

atwood = @ (p1,p2) (p2-p1)/(p2+p1); % calculate the atwood number
kZero = @ (at,p2,u1,u2) (((2*at^(1/2))/(1+at))*((p2*g^(1/2))/(u1+u2)))^(2/3); % calculate the k zero
omegaZero = @ (at,p2,u1,u2) (((2*at^2)/(1+at))*((p2*g^2)/(u1+u2)))^(1/3); % calcate the omega zero
ampAccelMod = @ (w,a) (w^2)*(a); % calculate b, the amplitude of acceleration modulation
defD = @ (at,u1,u2,p2,w,k) ((1+at)*(u1+u2)*(k^2))/(2*p2*w); % define a D for later use
defK = @ (at,o,k,p2,w) (((1+at)*o*(k^3))/(2*p2*(w^2)))-((at*k*g)/(w^2)); % define K for later use
defkc = @ (at,p2,o) (((2*at)/(1+at))*((p2*g/o)))^(2/3);

% Begin script proper
% Make an option list
fprintf('Use built in values or define new ones?\n1) Define new values\n2) SAE140 and air\n3) Water and air\n')%4) SAE140 and KI (potassium iodide)\n')
options = input('(123): '); % get user decision
if options == 1;
    fprintf('The denser fluid is #2\n');
    rho1 = input('rho1 = '); % fluid 1 density
    rho2 = input('rho2 = '); % fluid 2 density
    mu1 = input('mu1 = '); % fluid 1 viscosity
    mu2 = input('mu2 = '); % fluid 2 viscosity
    sig = input('sigma = '); % inteface surface tension
elseif options == 2;
    rho1 = 0; % g/cm^3
    rho2 = 0.9; % air is negligible comparitively
    mu1 = 0; % g/cm s
    mu2 = 9; % air is negligible comparitively
    sig = 32; % dyn/cm
elseif options == 3;
    rho1 = 0;
    rho2 = 1;
    mu1 = 0; % g/cm s
    mu2 = 0.0089;
    sig = 72.8; % dyn/cm
% elseif options == 4;
    % rho1 = 0.9; % SAE140 g/cm^3
    % rho2 = 1.6; % potassium jodide g/cm^3
	% mu1 = 
	% mu2 = 
	% sig = 
end

At = atwood(rho1,rho2); % call and assign Atwood number
ko = kZero(At,rho2,mu1,mu2); % call and assign the kzero 
omegao = omegaZero(At,rho2,mu1,mu2); % call and assign the omegazero
kc = defkc(At,rho2,sig); % call and assign the kc

lambda = 100; % arbitrary lambda?
k = 2*pi/lambda;
Kappa = k/ko;
Kappao = kc/ko;
Kappam = km/ko;

omega = 1:300;
omegaBar = omega/omegao;
b = 2*g.*omega./(sqrt(At*km*g));
A = b./(omega.^2);
%fprintf('omega = %i\n A = %i\n',omega,A)
plot(omega,A);xlabel('\omega');ylabel('Amplitude');
figure(2)
plot(omegaBar/sqrt(Kappam),b/g,'r');xlabel('\omega Bar/ sqrt(\kappa M)');ylabel('b/g');

end