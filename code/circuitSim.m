clearvars
clearvars -GLOBAL
close all
set(0,'DefaultFigureWindowStyle', 'docked')

sizex = 6;
sizey = 6;

% Voltage Range

Vmin = 0.1;
Vmax = 10;

% Components

Cap = 0.25;
R1 = 1;
R2 = 2;
L = 0.2;
% R3 = R3finder(Vmin,Vmax,20);
R3 = 10;
alpha = 100;
R4 = 0.1;
Ro = 1000;
omega = 10;

% C Matrix
C = zeros(sizex,sizey);
C(2,1) = -Cap;
C(2,2) = Cap;
C(6,6) = L;

% G Matrix
G = zeros (sizex, sizey);
G(1,1) = 1;
G(2,1) = -1/R1;
G(2,2) = (1/R1) + (1/R2);
G(2,6) = -1;
G(3,3) = 1/R3;
G(3,6) = 1;
G(4,3) = -alpha/R3;
G(4,4) = 1;
G(5,4) = -R4;
G(5,5) = R4 - (1/Ro);
G(6,2) = 1;
G(6,3) = -1;

% F Vector
F = zeros(1,sizey);
stepsize = 21;

VoutVect = zeros(1,stepsize);
V3Vect = zeros(1,stepsize);

% DC Sweep

for i = -10:10
    F(1) = i;
    
    % V vector
    V = (G + omega.*C)\F';
    VoutVect(i+11) = V(5);
    V3Vect(i+11) = V(3);
end
subplot(4,3,1);

plot (linspace(-10,10,stepsize),VoutVect);
title('-10V to 10V');
hold on
plot (linspace(-10,10,stepsize), V3Vect);
legend('Vo', 'V3');
xlabel('Vin');
ylabel('V');



% AC Sweep 

F = zeros(1,sizey);
F(1) = 1;
stepsize = 100;

VoutVect = zeros(1,stepsize);
V3Vect = zeros(1,stepsize);

omega = linspace(1,100,stepsize);

for i = 1:stepsize
    V = (G + 1j*omega(i).*C)\F';
    VoutVect(i) = V(5);
    V3Vect(i) = V(3);
end    

subplot(4,3,2);

plot (omega,abs(VoutVect));
title(' 0 to 100Hz');
hold on
plot (omega, abs(V3Vect));
legend('Vo', 'V3');
xlabel('w');
ylabel('V');

gain = 20 * log(abs(VoutVect./F(1)));

subplot(4,3,3);
plot(omega, gain);
title('Gain Vo/V1 in dB');
xlabel('w');
ylabel('Vo/V1 in dB');


% Histogram Cap and gain
omega = pi;
Crand = 0.05*randn(1,stepsize) + Cap;
for i = 1:stepsize
    V = (G + 1j*omega.*Crand(i))\F';
    VoutVect(i) = V(5);
end    

gain = 20 * log(abs(VoutVect./F(1)));

subplot(4,3,4);
histogram(Crand);
xlabel('C');
ylabel('Number');

subplot(4,3,5);
histogram(gain);
xlabel('Gain (dB)');
ylabel('Number');

% TIME SIMULATION

deltaT = 1e-6;

% A vector
A = (C./deltaT) + G;

timesteps = 1000;
Vp = zeros(sizey,1);


% F vector
F = zeros(1,sizey);

% Vin and Vout Vectors
VinVect = zeros(timesteps,1);
VoutVect = zeros(timesteps,1);

% Time simulation - step function

timeVector = linspace(1,timesteps,timesteps);

for i = 2:timesteps
    
    % F vector
    
    if (i == 30)
        F(1) = 1;
    end
    
    V = A\(((C * Vp)./deltaT) + F');
  
    subplot(4,3,6)
    plot([timeVector(i-1) timeVector(i)],[Vp(1) V(1)],'-r');
    
    hold on
    
    plot([timeVector(i-1) timeVector(i)],[Vp(5) V(5)],'-b');

    pause(0.01);
    
    VinVect(i) = V(1);
    VoutVect(i) = V(5);
    
    Vp = V;
end

legend('Vin', 'Vout');
title('Time Simulation - Step function');

xlim([0 1000]);
ylim([0 12]);
xlabel('Time (ms)');
ylabel('Voltage');
subplot(4,3,7)
plot(linspace(1,1000,1000),fftshift(20*log(fft(VoutVect))));

% Time simulation - sinusoidal function

% Vin and Vout Vectors
VinVect = zeros(timesteps,1);
VoutVect = zeros(timesteps,1);
Vp = zeros(sizey,1);

for i = 2:timesteps
    
    % F vector
    
    F(1) = sin(2 * pi * (1/0.03) * timeVector(i) * deltaT);
    
    V = A\(((C * Vp)./deltaT) + F');
  
    subplot(4,3,8)
    plot([timeVector(i-1) timeVector(i)],[Vp(1) V(1)],'-r');
    
    hold on
    
    plot([timeVector(i-1) timeVector(i)],[Vp(5) V(5)],'-b');
    pause(0.01);
    
    VinVect(i) = V(1);
    VoutVect(i) = V(5);
    
    Vp = V;
end

legend('Vin', 'Vout');
title('Time Simulation - Step function');

% xlim([0 1000]);
% ylim([0 12]);
xlabel('Time (ms)');
ylabel('Voltage');
subplot(4,3,9)
plot(linspace(1,1000,1000),fftshift(20*log(fft(VoutVect))));
% Time simulation - gaussian pulse

% Vin and Vout Vectors
VinVect = zeros(timesteps,1);
VoutVect = zeros(timesteps,1);
Vp = zeros(sizey,1);

pulsepos = 30 * randi(10);
% pulsepos = 5;
delayCnt = 0;
deltaT = 0.06;

for i = 2:timesteps
    
    % F vector
    
     if (i >= pulsepos)
        delayCnt = delayCnt + 1;
        if(delayCnt == 60)
          F(1) = 1;
          delayCnt = 0;
        end
         
     end
    
    V = A\(((C * Vp)./deltaT) + F');
  
    subplot(4,3,10)
    plot([timeVector(i-1) timeVector(i)],[Vp(1) V(1)],'-r');
    
    hold on
    
    plot([timeVector(i-1) timeVector(i)],[Vp(5) V(5)],'-b');

    pause(0.01);
    
    VinVect(i) = V(1);
    VoutVect(i) = V(5);
    
    Vp = V;
end
legend('Vin', 'Vout');
title('Time Simulation - Step function');

 xlim([0 1000]);
% ylim([0 12]);
xlabel('Time (ms)');
ylabel('Voltage');
subplot(4,3,11)
plot(linspace(1,1000,1000),fftshift(20*log(fft(VoutVect))));


