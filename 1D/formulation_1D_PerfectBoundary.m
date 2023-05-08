%% Housekeeping
clc
close all;
clear all;

load('colormap_BWR.txt','CMAP');

%% Constants
c_0 = 299792458;
eps_0 = 8.854187817e-12;
mu_0 = 1.256637061e-6;
eta_0 = sqrt(mu_0/eps_0);

% Simulation Parameters

% Compute Default Grid Resolution
lambda = [5e-4 1e-3 1e-2 1e-1]; % Array of wavelengths of interest, in meter
n_max = 3.927; % Maximum refractive index of material in use
N_lambda = 20; % Resolution of one wavelength in number of cells, typical value: 10

d_min = 1e-5; % Size of smallest feature of the structure of interest, in meter
N_d = 6; % Resolution of smallest feature in number of cells, typical value: 4

dz1 = min(lambda)/n_max/N_lambda;
dz2 = d_min/N_d;
dz = min(dz1,dz2);

% Snap Grid to Criticcal Dimensions
dc = 8e-3; % Critical Dimension
N = ceil(dc/dz);
dz = dc/N;
Nz = N;

% Compute Time Step
n_bc = 1; % Refractive index at grid boundaries
dt = n_bc*dz/2/c_0;

% Compute Gaussian Source Parameters
f_c = c_0/min(lambda); % Upper limit of frequency for analysis
tau = 0.5/f_c; % Pulse-Width of Gaussian Source
t_0 = 6*tau; % Delay of Gaussian Source to avoid abrupt change of fields

% Compute Required Simulation Time
t_prop = n_max*Nz*dz/c_0; % Time it takes for wave to propagate across the grid once
n_bounce = 5; % Number of bounces of the wave in the simulation
T_sim = 12*tau + n_bounce*t_prop; % Allow for the entire pulse to propagate, along with allowance for n_bounce number of bounces for the wave

N_Steps = ceil(T_sim/dt);

% Initialize Materials to Free Space
eps_r = ones(1,Nz);
mu_r = ones(1,Nz);

% Insert Material
n_material_0 = 1.00293;
n_material_1 = 3.927;
n_material_2 = 1.33;
n_material_3 = 2.41;
eps_r(:) =  n_material_0^2;
eps_r(ceil(0.3*Nz):ceil(0.4*Nz)) = n_material_3^2;
eps_r(ceil(0.4*Nz):ceil(0.5*Nz)) = n_material_2^2;
eps_r(ceil(0.5*Nz):ceil(0.6*Nz)) = n_material_1^2;
%eps_r(ceil(0.6*Nz):ceil(0.7*Nz)) = n_material_3^2;
%eps_r(ceil(0.7*Nz):ceil(0.8*Nz)) = n_material_2^2;
%eps_r(ceil(0.8*Nz):ceil(0.9*Nz)) = n_material_1^2;

% Compute Source
nz_src = 3; % z-index of Source, typically set as the first cell of Total Field region, ie z=3
n_src = n_material_0^2; % Refractive index at Source
t = (0:N_Steps-1)*dt;
del_t = 0.5*(n_src*dz/c_0 + dt);
A = -sqrt(eps_r(nz_src)/mu_r(nz_src));

Ey_src = exp(-((t-t_0)/tau).^2);
Hx_src = A*exp(-((t-t_0+del_t)/tau).^2);

% Source Parameters for calculating Fourier Transforms
fMax_FFT = 1.0*f_c;
N_FFT = 2048;
f_FFT = linspace(0,fMax_FFT,N_FFT);

% Initialize Fourier Transforms
kernel_FFT = exp(-1i*2*pi*f_FFT*dt);
reflectance_FFT = zeros(1,N_FFT);
transmittance_FFT = zeros(1,N_FFT);
source_FFT = zeros(1,N_FFT);


% Initialize Fields to zero
Ey = zeros(1,Nz);
Ey1 = Ey2 = 0; % z-High boundary fields
Hx = zeros(1,Nz);
Hx1 = Hx2 = 0; % z-Low boundary fields



% Compute Update Coefficients
m_Ey = (c_0*dt)./eps_r;
m_Hx = (c_0*dt)./mu_r;


%% Main FDTD Loop
disp("Starting Simulation:");
disp(["Total steps: ", num2str(N_Steps)]);

z = (0:Nz-1)*dz;
figure(1)
h = plot((0:Nz-1)*dz*1e12/c_0,Ey,'-b', 'LineWidth',2);
hold on
h = plot((0:Nz-1)*dz*1e12/c_0,Hx,'-r', 'LineWidth',2);
hold off
grid on
plotYLim = [-2 2];

hParent = get(h,'parent');
set(hParent, 'ylim', plotYLim);
set(hParent, 'xlim', [0 Nz-1]*dz*1e12/c_0);
set(hParent, 'FontSize', 14, 'LineWidth', 2);


for T = 1:N_Steps
  % Inject Soft Source (Not needed with TFSF Source)
  %Ey(nz_src) = Ey(nz_src) + Ey_src(T);
  %Hx(nz_src-1) = Hx(nz_src-1) + Hx_src(T);

  % Update H from E(Perfect Boundary Conditions)
  Hx2 = Hx1;
  Hx1 = Hx(1);

  for nz = 1:Nz-1
    Hx(nz) = Hx(nz) + m_Hx(nz)*( Ey(nz+1) - Ey(nz) )/dz;
  endfor
  Hx(Nz) = Hx(Nz) + m_Hx(Nz)*( Ey2 - Ey(Nz) )/dz; % Perfectly Matched Load at z-High Boundary

  % Handle H-field source
  Hx(nz_src-1) = Hx(nz_src-1) - m_Hx(nz_src-1)*Ey_src(T)/dz;


  % Update E from H (Perfect Boundary Conditions)
  Ey2 = Ey1;
  Ey1 = Ey(Nz);

  Ey(1) = Ey(1) + m_Ey(1)*( Hx(1) - Hx2 )/dz; % Perfectly Matched Load at z-Low Boundary
  for nz = 2:Nz
    Ey(nz) = Ey(nz) + m_Ey(nz)*( Hx(nz) - Hx(nz-1) )/dz;
  endfor

  % Handle E-field source
  Ey(nz_src) = Ey(nz_src) - m_Ey(nz_src)*Hx_src(T)/dz;

  % Update Fourier Transforms
  kernel_FFT_T = kernel_FFT.^T;
  reflectance_FFT = reflectance_FFT + kernel_FFT_T*Ey(1);
  transmittance_FFT = transmittance_FFT + kernel_FFT_T*Ey(Nz);
  source_FFT = source_FFT + kernel_FFT_T*Ey_src(T);

  % Visualize
  % Update plot every T_plotUpdate steps
  T_plotUpdate = 50;
  if mod(T,T_plotUpdate) == 0 || T == N_Steps
    clc
    disp("Starting Simulation:");
    disp(["Total steps: ", num2str(N_Steps)]);
    disp(["Progress: ", num2str(T/N_Steps*100), "% (t = ", num2str(T*dt*1e12), " pico-second)"]);

    figure(1)
    h = plot((0:Nz-1)*dz*1e12/c_0,Ey,'-b', 'LineWidth',2);
    hold on
    h = plot((0:Nz-1)*dz*1e12/c_0,Hx,'-r', 'LineWidth',2);
    hold off
    grid on
    xlabel("z (in light-picoseconds)")
    title(["Progress: ", num2str(T/N_Steps*100), "% (t = ", num2str(T*dt*1e12), " pico-second)"]);
    hParent = get(h,'parent');
    set(hParent, 'ylim', plotYLim);
    set(hParent, 'xlim', [0 Nz-1]*dz*1e12/c_0);
    set(hParent, 'FontSize', 14, 'LineWidth', 2);
    refresh

  endif



endfor

% Finalize Fourier Transform
reflectance_FFT = reflectance_FFT*dt;
transmittance_FFT = transmittance_FFT*dt;
source_FFT = source_FFT*dt;

% Normalize Fourier Transforms wrt Source Spectrum
reflectance_FFT = abs(reflectance_FFT./source_FFT).^2;
transmittance_FFT = abs(transmittance_FFT./source_FFT).^2


%% Visualize Steady State Reflectance and Transmittance
xlim_fft = [0 fMax_FFT*1e-12];
ylim_fft = [-30 3]
ylim_fft_2 = [-0.5 0.5];

% Set Tick Markings
xm = [min(xlim_fft)*1e-12:0.25:fMax_FFT*1e-12];
xt = strtrim(cellstr(num2str(xm','%3.2f'))');

ym = [min(ylim_fft):3:max(ylim_fft)];
yt = strtrim(cellstr(num2str(ym','%2.1f'))');


fig_FFT = figure(2);

subplot(311)
h = plot(f_FFT*1e-12,10*log10(reflectance_FFT),'-b', 'LineWidth',2);
grid on
xlim(xlim_fft)
ylim(ylim_fft)
xlabel('Frequency (THz)')
ylabel('Reflectance R(f) [dB]')
set(get(h,'Parent'), 'FontSize', 14, 'LineWidth', 2, 'XTick',xm, 'XTickLabel',xt, 'YTick',ym, 'YTickLabel',yt);

subplot(312)
h = plot(f_FFT*1e-12,10*log10(transmittance_FFT),'-r', 'LineWidth',2);
grid on
xlim(xlim_fft)
ylim(ylim_fft)
xlabel('Frequency (THz)')
ylabel('Transmittance T(f) [dB]')
set(get(h,'Parent'), 'FontSize', 14, 'LineWidth', 2, 'XTick',xm, 'XTickLabel',xt, 'YTick',ym, 'YTickLabel',yt);

subplot(313)
h = plot(f_FFT*1e-12,10*log10(reflectance_FFT+transmittance_FFT),'-k', 'LineWidth',2);
grid on
xlim(xlim_fft)
ylim(ylim_fft_2)
xlabel('Frequency (THz)')
ylabel('R(f)+T(f) [dB]')
set(get(h,'Parent'), 'FontSize', 14, 'LineWidth', 2, 'XTick',xm, 'XTickLabel',xt);


