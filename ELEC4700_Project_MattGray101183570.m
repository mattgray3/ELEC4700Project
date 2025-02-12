% Milestone 4 Matt Gray 101183570 ELEC 4700 Project 

%Default setup for plotting
set(0,'defaultaxesfontsize', 20)
set(0, 'DefaultFigureWindowStyle', 'docked')
set(0, 'DefaultLineLineWidth', 2)
set(0, 'Defaultaxeslinewidth', 2)

set(0, 'DefaultFigureWindowStyle', 'docked')

% Constants used in Simulation
c_c = 299792458; % Speed of light
c_eps_0 = 8.8542149e-12; % eo
c_eps_0_cm = c_eps_0/100; % eo in cm
c_mu_0 = 1/c_eps_0/c_c^2; % uo 
c_q = 1.60217653e-19; % electron charge
c_hb = 1.05457266913e-34; %plancks constant
c_h = c_hb*2*pi; %reduced plancks constant
kappa0 = 100; % Milestone3 --> Kappa constant imaginary part of effective refractive index
kappaStart = 1/3;
kappaStop = 2/3;


%Input parameters for L and R constants
InputParasL.E0 = 1e5; %Adjusts the amplitude of the wave
InputParasL.we = 0; %Milestone 2: Modulates the waveform **Remember in PicoSecond Scale Needs to be like e13 to modulate
InputParasL.t0 = 2e-12; % inital time of the wave
InputParasL.wg = 5e-13; %pusle width
InputParasL.phi = 0;
InputParasR = 0; % right source equal to zero since only reflections come from right side

%Material and Simulation Parameters
n_g = 3.5; % medium index
vg = c_c/n_g * 1e2; % medium velocity 
Lambda = 1550e-9; %wavelength

RL = 0.9i; % Reflectoin coefficients for L and R sides
RR = 0.9i;

plotN = 10; % frequency of plotting

L =  1000e-6 * 1e2; %Len of sim
XL = [0, L]; % X and Y Axis limits
YL = [0, InputParasL.E0];

Nz = 500; % # spatical grid pts
dz = L / (Nz-1); % spacital step size
dt =  dz/vg; % Time step size
fsync = dt * vg /dz; % sync factor

Nt =  floor(2*Nz); % total # of time steps
tmax = Nt * dt; % Max sim time
t_L = dt * Nz;

% Field Setup
z = linspace(0, L, Nz); % Linear Spaced Vector
time = nan(1,Nt);
InputL = nan(1, Nt);
InputR = nan(1, Nt);
OutputL = nan(1, Nt);
OutputR = nan(1, Nt);

g_fwhm = 3.53e+012/10; % Milestone 4
LGamma = g_fwhm * 2 * pi; % Milestone 4
Lw0 = 0.0; %Milestone 4
LGain = 0.001; %Milestone 4

kf = 0; %Milestone 3 -->kappa fwd,do not require rev cause we assume they are equal

beta_r = 0; % Real Beta value Milestone 2 Real component adds a "Twist"
beta_i = 0; % Imaginary Beta Value Milestone 2 Imj component adds gain to the waveform

kappa = kappa0 * ones(size(z));%Milestone 3 --> initalizing kappa
kappa(z<L * kappaStart) = 0;% Milestone 3 --> establishing grating boundaries
kappa(z>L * kappaStop) = 0;

Ef = zeros(size(z)); % Electric field forward
Er = zeros(size(z));% Electric field Reflected

Pf = zeros(size(z)); %Milestone4
Pr = zeros(size(z)); %Milestone4

Efp = Ef; %Milestone 4
Erp = Er; %Milestone 4
Pfp = Pf; %Milestone 4
Prp = Pr; %Milestone 4

Ef1 = @SourceFct; % Generates source for electric field
ErN = @SourceFct;

%Time Loop
t = 0; % initalizing time of simulation
time(1) = t;

InputL(1) = Ef1(t, InputParasL);
InputR(1) = ErN(t, InputParasR);

OutputR(1) = Ef(Nz);
OutputL(1) = Er(1);

Ef(1) = InputL(1);
Er(Nz) = InputR(1);

% Setup For Simulation Figures
figure('Name','Fields')
subplot(3,1,1)
plot(z*1000, real(Ef), 'r');
hold off
xlabel('z(\mum')
ylabel('E_f')
subplot(3,1,2)
plot(z*1000, real(Er), 'b');
xlabel('z(\mum)')
ylabel('E_f')
hold off
subplot (3,1,3)
plot(time * 1e2, real(InputL), 'r'); hold on
plot(time * 1e2, real(OutputR), 'r--'); 
plot(time * 1e2, real(InputR), 'b'); hold on
plot(time * 1e2, real(OutputR), 'b--'); 
xlabel('time(ps)')
ylabel('E')


hold off

%Main Time Stepping Loop iterativle increases time, soures and outputs
for i = 2:Nt
    t = dt *(i-1);
    time(i) = t;

    InputL(i) = Ef1(t, InputParasL);
    InputR(i) = ErN(t, 0);

    Ef(1) = InputL(i) + RL*Er(1);
    Er(Nz) = InputR(i) + RR*Ef(Nz);

    beta = ones(size(z)) * (beta_r + 1i * beta_i); %Milestone 2 Beta Equation
    exp_det = exp(-1i * dz * beta); %Milestone 2 Exponent

    Ef(2:Nz) = fsync * exp_det(1:Nz-1) .* Ef(1:Nz-1) + 1i*dz*kappa(1:Nz-1) .* Er(1:Nz-1); % Edited Milestone 2 & 3
    Er(1:Nz-1) = fsync * exp_det(2:Nz).* Er(2:Nz) + 1i*dz*kappa(1:Nz-1) .* Ef(2:Nz) ; % Edited Milestone 2 & 3
        
    Pf(1) = 0;
    Pf(Nz) = 0;
    Pr(1) = 0;
    Pr(Nz) = 0;
    Cw0 = -LGamma + 1i*Lw0;

    Tf = LGamma * Ef(1:Nz-2) + Cw0 * Pfp(2:Nz-1) + LGamma * Efp(1:Nz-2);
    Pf(2:Nz-1) = (Pfp(2:Nz-1) + 0.5 * dt * Tf) ./ (1 - 0.5 * dt * Cw0);
    Tr = LGamma * Er(3:Nz) + Cw0 * Prp(2:Nz-1) + LGamma * Erp(3:Nz);
    Pr(2:Nz-1) = (Prp(2:Nz-1) + 0.5 * dt * Tr) ./ (1 - 0.5 * dt * Cw0);

    Ef(2:Nz-1) = Ef(2:Nz-1) - LGain*(Ef(2:Nz-1) - Pf(2:Nz-1));
    Er(2:Nz-1) = Er(2:Nz-1) - LGain*(Er(2:Nz-1) - Pr(2:Nz-1));

    Efp = Ef;
    Erp = Er;
    Pfp = Pr;
    Prp = Pr;

  

    OutputR(i) = Ef(Nz)*(1-RR);
    OutputL(i) =  Er(1)*(1-RL);

    if mod(i, plotN) == 0
        subplot(3,1,1)
        plot(z*10000, real(Ef), 'r'); hold on 
        plot(z*10000, imag(Ef), 'r--'); hold off
        xlim(XL*1e4)
        ylim(YL)
        xlabel('z(\mum)')
        ylabel('E_f')
        legend('\Re', '\Im')
        hold off

        subplot(3,1,2)
        plot(z*10000, real(Er), 'b'); hold on 
        plot(z*10000, imag(Er), 'b--'); hold off
        xlim(XL*1e4)
        ylim(YL)
        xlabel('z(\mum)')
        ylabel('E_r')
        legend('\Re', '\Im')
        hold off

        subplot(3,1,3)
        plot(time*1e12, real(InputL), 'r'); hold on
        plot(time*1e12, real(OutputR), 'g');
        plot(time*1e12, real(InputR), 'b');
        plot(time*1e12, real(OutputL), 'w');

        xlim([0, Nt*dt*1e12])
        ylim(YL)

        xlabel('time(ps)')
        ylabel('0')
        legend('Left Input', 'Right Output', 'Right Input', 'Left Output', ...
            'Location','east')
        hold off
        pause(0.01)
    end
end

fftOutput = fftshift(fft(OutputR));
fftOutput2 = fftshift(fft(OutputL));
fftInput = fftshift(fft(InputL));
omega =  fftshift(wspace(time));



figure('Name', 'Frequency Domain Analysis MAG')

plot(omega, 20*log10(abs(fftOutput)), 'r', omega, 20*log10(abs(fftInput)), 'b')
xlabel('Frequency (rad/s)')
ylabel('Magnitude')
title('Mag Spectrum, R = Out, B = In')
grid on

%subplot(3,1,2)
%plot(omega, abs(fftInput), 'b')
%xlabel('Frequency (rad/s)')
%ylabel('Magnitude')
%title('Mag Spectrum IN')
%grid on

figure('Name', 'Frequency Domain Analysis Phase')

plot(omega, unwrap(angle(fftOutput)), 'r', omega, unwrap(angle(fftInput)), 'b' )
xlabel('Frequency (rad/s)')
ylabel('Phase(Rads)')
title('Phase Spectrum, R = Out, B = In')
grid on

%subplot(3,1,2)
%plot(omega, unwrap(angle(fftInput)), 'b')
%xlabel('Frequency (rad/s)')
%ylabel('Phase(Rads)')
%title('Phase Spectrum IN')
%grid on

%subplot(3,1,3)
%plot(omega, real(fftOutput), 'g', omega, imag(fftOutput), 'm')
%xlabel('Frequency (rad/s)')
%ylabel('Amplitude')
%title('Real and Imaj parts FFTOUTPUT')
%legend('Real', 'Imaj')
%grid on

%figure('Name', 'Kappa vs Z')
%plot(z, kappa, 'p')
%xlabel('z')
%ylabel('Kappa')
%grid on