%{
Author : Maharshi Gurjar
Elec 4700 Assignment 4 -  Circuit Modeling
Part 1 - No noise simulation
%}
clc; close all; clear;
set(0, 'DefaultFigureWindowStyle', 'docked')
%% Find value for R3 using FDM from assignment 2
%{
Author : Maharshi Gurjar
Elec 4700 Assignment 2 -  Finite Difference Method
Part 2
%}
%% Define G,C,b matrix and fill in
global G C b
NrNodes = 5;  % The total number of nodes in the circuit
% Define G, C, b, for a circuit (do not include additional variables).
G = zeros(NrNodes,NrNodes);
C = zeros(NrNodes,NrNodes);
b = zeros(NrNodes,1);
%set up component values
vol(1,0,1);
res(1,2,1)
cap(1,2,0.25)
res(2,0,2)
ind(2,3,0.2)
res(3,0,186)
vcvs(4,0,3,0,(100/186))
res(4,5,0.1)
res(5,0,1000)


%% DC Sweep of Circuit
VNode3 =[];
Vout =[];
for n= -10:1:10
    b(6) = n; %Mapping sweep for the input
    Voltage = G\b; % Solve for Voltage
    VNode3 = [VNode3 Voltage(3)]; %Voltage at node 3
    Vout = [Vout Voltage(5)]; %Voltage at output node
    %gain = Vout(n)/Vin(n); % Vout/Vin gain
end
figure('name','DC Sweep');
plot(-10:1:10,VNode3, -10:1:10, Vout);
title('DC Sweep','interpreter','latex')
xlabel('Input Voltage (V)')
ylabel('Ouput Voltage (V)')
legend(['V_{node3}';'V_{Ouput}'])
axis tight;
grid on;
%saveas(gcf,fullfile('D:\School Work\ELEC 4700\My 4700 Code\Assignment 4\Simulation Results','[Part1]DC Sweep.png'),'png')
%% Frequency Response of circuit
Freq = logspace(0, 5, 5000); %Define the frequency of the simulation
for n=1:length(Freq)
    w = 2*pi*Freq(n); %define omega at frequency point
    s = 1i*Freq(n); %define s at frequency point
    A = G + (s.*C); %define full G C matrix
    Voltage = A\b; %solve for voltage at frequency
    Vout(n) = abs(Voltage(5));%output at frequency
    gain(n) = 20*log(abs(Voltage(5)));%gain at frequency
end
figure('name','reseponse');
subplot(2,1,1)
semilogx(Freq, Vout);
xlabel('Frequency (Hz)');
ylabel('V_O (V)');
title('Frequncy Response','interpreter','latex');
axis auto;
grid on;
subplot(2,1,2)
semilogx(Freq, gain);
xlabel('Frequency (Hz)');
ylabel('Gain (dB)');
title('Frequency Response (dB)','interpreter','latex');
axis auto;
grid on;
%saveas(gcf,fullfile('D:\School Work\ELEC 4700\My 4700 Code\Assignment 4\Simulation Results','[Part1]FrequencySweep.png'),'png')
%% Create histogram with capacitance pertubations
std_dev = 0.05;
range = 10e3;
Rand_Cap = 0.25*std_dev.*randn(range,1);
for n=1:range
    C(1,1) = Rand_Cap(n);
    C(2,2) = Rand_Cap(n);
    C(1,2) = -Rand_Cap(n);
    C(2,1) = -Rand_Cap(n);
    s = 2*pi;
    A = G + (s.*C);
    V = A\b;
    gain(n) = 20*log10(abs(abs(V(5)))/abs(V(1)));
end
figure('name','Histogram of gain across a range of capacitors')
subplot(2,1,1)
histogram(Rand_Cap+1)
xlabel('Gain (dB)')
ylabel('Bin count')
title('Gain across a range of capacitors','interpreter','latex')
grid on;
axis auto;
subplot(2,1,2)
histogram(gain)
xlabel('Gain (dB)')
ylabel('Bin count')
title('Gain across a range of capacitors','interpreter','latex')
grid on;
axis auto;
%saveas(gcf,fullfile('D:\School Work\ELEC 4700\My 4700 Code\Assignment 4\Simulation Results','[Part1]Histogram.png'),'png')
%% Transient Analysis
clc; clear;
% Define G,C,b matrix and fill in
global G C b
NrNodes = 5;  % The total number of nodes in the circuit
% Define and create simulation parameters
Sim_Time = 1;
Steps = 1000;
delta = Sim_Time/Steps;
Time = linspace(0,1,Steps);
% Step the frequency of the input
Frequency = 1/0.03;
Frequency2 = 1/0.3;
% Define the three input signals
Input_Step = zeros(1,Steps);
Input_Sin = zeros(1,Steps);
Input_Sin2 = zeros(1,Steps);
Input_Guassian = zeros(1,Steps);

% To make the step input an actual step input
Transition_Time = 0.03;
Input_Step = double(Time>=Transition_Time);

% Fill in sin input with sin values
Input_Sin = sin(2*pi*Frequency*Time);
Input_Sin2 = sin(2*pi*Frequency2*Time);
% Create Guassian signal values
GaussianDist = makedist('Normal', 'mu', 0.06, 'sigma', 0.03);
GuassianPulse = pdf(GaussianDist, Time);
Input_Pulse = (GuassianPulse.*sin(pi*Time))/max(GuassianPulse.*sin(pi*Time));

% Do trans. sim for step input
Vin = Input_Step;
for i = 1:Steps
    % Define G, C, b, for a circuit (do not include additional variables).
    G = zeros(NrNodes,NrNodes);
    C = zeros(NrNodes,NrNodes);
    b = zeros(NrNodes,1);
    %set up component values
    
    res(1,2,1)
    cap(1,2,0.25)
    res(2,0,2)
    ind(2,3,0.2)
    res(3,0,185)
    vcvs(4,0,3,0,(100/185))
    res(4,5,0.1)
    res(5,0,1000)
    vol(1,0,Vin(i));
    if i == 1
        NewVoltage = (G + (C./delta)) \ b;
    else
        NewVoltage = (G + (C./delta)) \ (b  + (C./delta)*OldVoltage);
    end
    Voutput(i) = NewVoltage(5);
    OldVoltage = NewVoltage;
end
Xaxis = (-length(Vin)/2:length(Vin)/2-1);
figure('Name','Trans Sim Step')
subplot(2,1,1)
plot(Time,Vin,Time,Voutput)
grid on;
axis auto;
legend('Vin','Vout')
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Step response of circuit','interpreter','latex')
subplot(2,1,2)
plot(Xaxis, mag2db(abs(fftshift(fft(Vin)))),Xaxis, mag2db(abs(fftshift(fft(Voutput)))));
title('Frequency Spectrum response','interpreter','latex')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
grid on;
axis auto;
legend('Vin','Vout')
%saveas(gcf,fullfile('D:\School Work\ELEC 4700\My 4700 Code\Assignment 4\Simulation Results','[Part1]Step input response.png'),'png')

% Do trans. sim for Sin input
Vin = Input_Sin;
for i = 1:Steps
    % Define G, C, b, for a circuit (do not include additional variables).
    G = zeros(NrNodes,NrNodes);
    C = zeros(NrNodes,NrNodes);
    b = zeros(NrNodes,1);
    %set up component values
    vol(1,0,Vin(i));
    res(1,2,1)
    cap(1,2,0.25)
    res(2,0,2)
    ind(2,3,0.2)
    res(3,0,185)
    vcvs(4,0,3,0,(100/185))
    res(4,5,0.1)
    res(5,0,1000)
    
    if i == 1
        NewVoltage = (G + (C./delta)) \ b;
    else
        NewVoltage = (G + (C./delta)) \ (b  + (C./delta)*OldVoltage);
    end
    Voutput(i) = NewVoltage(5);
    OldVoltage = NewVoltage;
end
Xaxis = (-length(Vin)/2:length(Vin)/2-1);
figure('Name','Trans Sim Sin')
subplot(2,1,1)
plot(Time,Vin,Time,Voutput)
grid on;
axis auto;
legend('Vin','Vout')
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Sin response of circuit','interpreter','latex')
subplot(2,1,2)
plot(Xaxis, mag2db(abs(fftshift(fft(Vin)))),Xaxis, mag2db(abs(fftshift(fft(Voutput)))));
title('Frequency Spectrum response','interpreter','latex')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
grid on;
axis auto;
legend('Vin','Vout')
%saveas(gcf,fullfile('D:\School Work\ELEC 4700\My 4700 Code\Assignment 4\Simulation Results','[Part1]Sin1.png'),'png')

%Impulse response
Vin = Input_Pulse;
for i = 1:Steps
    % Define G, C, b, for a circuit (do not include additional variables).
    G = zeros(NrNodes,NrNodes);
    C = zeros(NrNodes,NrNodes);
    b = zeros(NrNodes,1);
    %set up component values
    
    res(1,2,1)
    cap(1,2,0.25)
    res(2,0,2)
    ind(2,3,0.2)
    res(3,0,185)
    vcvs(4,0,3,0,(100/185))
    res(4,5,0.1)
    res(5,0,1000)
    vol(1,0,Vin(i));
    if i == 1
        NewVoltage = (G + (C./delta)) \ b;
    else
        NewVoltage = (G + (C./delta)) \ (b  + (C./delta)*OldVoltage);
    end
    Voutput(i) = NewVoltage(5);
    OldVoltage = NewVoltage;
end
Xaxis = (-length(Vin)/2:length(Vin)/2-1);
figure('Name','Trans Sim Impulse')
subplot(2,1,1)
plot(Time,Vin,Time,Voutput)
grid on;
axis auto;
legend('Vin','Vout')
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Pulse response of circuit','interpreter','latex')
subplot(2,1,2)
plot(Xaxis, mag2db(abs(fftshift(fft(Vin)))),Xaxis, mag2db(abs(fftshift(fft(Voutput)))));
title('Frequency Spectrum response','interpreter','latex')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
grid on;
axis auto;
legend('Vin','Vout')
%saveas(gcf,fullfile('D:\School Work\ELEC 4700\My 4700 Code\Assignment 4\Simulation Results','[Part1]Impluse.png'),'png')

% Do trans. sim for Sin input`
Vin = Input_Sin2;
for i = 1:Steps
    % Define G, C, b, for a circuit (do not include additional variables).
    G = zeros(NrNodes,NrNodes);
    C = zeros(NrNodes,NrNodes);
    b = zeros(NrNodes,1);
    %set up component values
    res(1,2,1)
    cap(1,2,0.25)
    res(2,0,2)
    ind(2,3,0.2)
    res(3,0,185)
    vcvs(4,0,3,0,(100/185))
    res(4,5,0.1)
    res(5,0,1000)
    vol(1,0,Vin(i));
    if i == 1
        NewVoltage = (G + (C./delta)) \ b;
    else
        NewVoltage = (G + (C./delta)) \ (b  + (C./delta)*OldVoltage);
    end
    Voutput(i) = NewVoltage(5);
    OldVoltage = NewVoltage;
end
Xaxis = (-length(Vin)/2:length(Vin)/2-1);
figure('Name','Trans Sim Sin diff Freq')
subplot(2,1,1)
plot(Time,Vin,Time,Voutput)
grid on;
axis auto;
legend('Vin','Vout')
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Sin response of circuit','interpreter','latex')
subplot(2,1,2)
plot(Xaxis, mag2db(abs(fftshift(fft(Vin)))),Xaxis, mag2db(abs(fftshift(fft(Voutput)))));
title('Frequency Spectrum response','interpreter','latex')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
grid on;
axis auto;
legend('Vin','Vout')
%saveas(gcf,fullfile('D:\School Work\ELEC 4700\My 4700 Code\Assignment 4\Simulation Results','[Part1]Sin2.png'),'png')
%% Circuit with noise
clc; clear;
% Define G,C,b matrix and fill in
global G C b
NrNodes = 5;  % The total number of nodes in the circuit
% Define and create simulation parameters
Sim_Time = 1;
Steps = 1000;
delta = Sim_Time/Steps;
Time = linspace(0,1,Steps);
% Step the frequency of the input
Frequency = 1/0.03;

Input_Guassian = zeros(1,Steps);

GaussianDist = makedist('Normal', 'mu', 0.06, 'sigma', 0.03);
GuassianPulse = pdf(GaussianDist, Time);
Input_Pulse = (GuassianPulse.*sin(pi*Time))/max(GuassianPulse.*sin(pi*Time));

Vin = Input_Pulse;
Vout = zeros(1,Steps);
Cn = 0.00001;
In = 0.001;
h = 1/10e3;
W = 0:h:1;
CurrNoise = In*randn(1,numel(W));
for i = 1:Steps
    % Define G, C, b, for a circuit (do not include additional variables).
    G = zeros(NrNodes,NrNodes);
    C = zeros(NrNodes,NrNodes);
    b = zeros(NrNodes,1);
    %set up component values
    vol(1,0,Vin(i));
    cur(3,0,CurrNoise(i));
    res(1,2,1)
    cap(1,2,0.25)
    cap(3,0,Cn)
    res(2,0,2)
    ind(2,3,0.2)
    res(3,0,185)
    res(4,5,0.1)
    res(5,0,1000)
    vcvs(4,0,3,0,100/185);
    if i == 1
        NewVoltage = (G + (C./delta)) \ b;
    else
        NewVoltage = (G + (C./delta)) \ (b  + (C./delta)*OldVoltage);
    end
    Voutput(i) = NewVoltage(5);
    OldVoltage = NewVoltage;
end

C

Xaxis = (-length(Vin)/2:length(Vin)/2-1);
figure('Name','Trans Sim Noise')
subplot(2,1,1)
plot(Time,Vin,Time,Voutput)
grid on;
axis auto;
legend('Vin','Vout')
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Sin response of circuit','interpreter','latex')
subplot(2,1,2)
plot(Xaxis, mag2db(abs(fftshift(fft(Vin)))),Xaxis, mag2db(abs(fftshift(fft(Voutput)))));
title('Frequency Spectrum response','interpreter','latex')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
grid on;
axis auto;
legend('Vin','Vout')
%saveas(gcf,fullfile('D:\School Work\ELEC 4700\My 4700 Code\Assignment 4\Simulation Results','[Part2]NoisySignal.png'),'png')
% Varrying Cn

Cn = [0.00001 0.0001 0.001 0.01 0.1];
for i = 1:numel(Cn)
    for j = 1:Steps
        % Define G, C, b, for a circuit (do not include additional variables).
        G = zeros(NrNodes,NrNodes);
        C = zeros(NrNodes,NrNodes);
        b = zeros(NrNodes,1);
        %set up component values
        vol(1,0,Vin(j));
        cur(3,0,CurrNoise(j));
        res(1,2,1)
        cap(1,2,0.25)
        cap(3,0,Cn(i))
        res(2,0,2)
        ind(2,3,0.2)
        res(3,0,185)
        res(4,5,0.1)
        res(5,0,1000)
        vcvs(4,0,3,0,100/185);
        if j == 1
            NewVoltage = (G + (C./delta)) \ b;
        else
            NewVoltage = (G + (C./delta)) \ (b  + (C./delta)*OldVoltage);
        end
        Voutput(j) = NewVoltage(5);
        OldVoltage = NewVoltage;
    end
    Xaxis = (-length(Vin)/2:length(Vin)/2-1);
    sgtitle('Delta Cn')
    if i == 1
        figure(21)
        subplot(5,2,i)
        plot(Time,Vin,Time,Voutput)
        grid on;
        axis auto;
        legend('Vin','Vout')
        xlabel('Time (s)')
        ylabel('Voltage (V)')
        title(sprintf('Frequency Specetrum Response, C_n = %f',Cn(i)))
        hold on;
        subplot(5,2,i+1)
        plot(Xaxis, mag2db(abs(fftshift(fft(Vin)))),Xaxis, mag2db(abs(fftshift(fft(Voutput)))));
        title('Frequency Spectrum response','interpreter','latex')
        xlabel('Frequency (Hz)')
        ylabel('Magnitude (dB)')
        grid on;
        axis auto;
        legend('Vin','Vout')
        hold on;
    else
        figure(21)
        subplot(5,2,(i*2-1))
        plot(Time,Vin,Time,Voutput)
        grid on;
        axis auto;
        legend('Vin','Vout')
        xlabel('Time (s)')
        ylabel('Voltage (V)')
        title(sprintf('Frequency Specetrum Response, C_n = %.4f',Cn(i)))
        hold on;
        subplot(5,2,i*2)
        plot(Xaxis, mag2db(abs(fftshift(fft(Vin)))),Xaxis, mag2db(abs(fftshift(fft(Voutput)))));
        title('Frequency Spectrum response','interpreter','latex')
        xlabel('Frequency (Hz)')
        ylabel('Magnitude (dB)')
        grid on;
        axis auto;
        legend('Vin','Vout')
        hold on;
    end
end
%saveas(gcf,fullfile('D:\School Work\ELEC 4700\My 4700 Code\Assignment 4\Simulation Results','[Part2]VarryingCn.png'),'png')

%Varrying Time

Steps = [10 100 1e3 1e4];

Cn = 0.00001;
for i = 1: numel(Steps)
    delta = 1/Steps(i);
    Time = linspace(0,1,Steps(i));
    GaussianDist = makedist('Normal', 'mu', 0.06, 'sigma', 0.03);
    GuassianPulse = pdf(GaussianDist, Time);
    Input_Pulse = (GuassianPulse.*sin(pi*Time))/max(GuassianPulse.*sin(pi*Time));
    Vin = Input_Pulse;
    Voutput=zeros(1,Steps(i));
    for j = 1:Steps(i)
        % Define G, C, b, for a circuit (do not include additional variables).
        G = zeros(NrNodes,NrNodes);
        C = zeros(NrNodes,NrNodes);
        b = zeros(NrNodes,1);
        %set up component values
        vol(1,0,Vin(j));
        cur(3,0,CurrNoise(j));
        res(1,2,1)
        cap(1,2,0.25)
        cap(3,0,Cn)
        res(2,0,2)
        ind(2,3,0.2)
        res(3,0,185)
        res(4,5,0.1)
        res(5,0,1000)
        vcvs(4,0,3,0,100/185);
        if j == 1
            NewVoltage = (G + (C./delta)) \ b;
        else
            NewVoltage = (G + (C./delta)) \ (b  + (C./delta)*OldVoltage);
        end
        Voutput(j) = NewVoltage(5);
        OldVoltage = NewVoltage;
    end
    Xaxis = (-length(Vin)/2:length(Vin)/2-1);
    sgtitle('Delta t')
    if i == 1
        figure(22)
        subplot(4,2,i)
        plot(Time,Vin,Time,Voutput)
        grid on;
        axis auto;
        legend('Vin','Vout')
        xlabel('Time (s)')
        ylabel('Voltage (V)')
        title(sprintf('Frequency Specetrum Response, t = %f',Steps(i)))
        hold on;
        subplot(4,2,i+1)
        plot(Xaxis, mag2db(abs(fftshift(fft(Vin)))),Xaxis, mag2db(abs(fftshift(fft(Voutput)))));
        title('Frequency Spectrum response','interpreter','latex')
        xlabel('Frequency (Hz)')
        ylabel('Magnitude (dB)')
        grid on;
        axis auto;
        legend('Vin','Vout')
        hold on;
    else
        figure(22)
        subplot(4,2,(i*2-1))
        plot(Time,Vin,Time,Voutput)
        grid on;
        axis auto;
        legend('Vin','Vout')
        xlabel('Time (s)')
        ylabel('Voltage (V)')
        title(sprintf('Frequency Specetrum Response, t = %.4f',Steps(i)))
        hold on;
        subplot(4,2,i*2)
        plot(Xaxis, mag2db(abs(fftshift(fft(Vin)))),Xaxis, mag2db(abs(fftshift(fft(Voutput)))));
        title('Frequency Spectrum response','interpreter','latex')
        xlabel('Frequency (Hz)')
        ylabel('Magnitude (dB)')
        grid on;
        axis auto;
        legend('Vin','Vout')
        hold on;
    end
end
%saveas(gcf,fullfile('D:\School Work\ELEC 4700\My 4700 Code\Assignment 4\Simulation Results','[Part2]VarryingTime.png'),'png')
    