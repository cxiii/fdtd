%% Housekeeping
clc
close all;
clear all;

load('colormap_BWR.txt','CMAP');

%% Constants
c_0 = 299792458;
eps_0 = 8.854187817e-12;
mu_0 = 1.256637061e-6;
eta_0 = sqrt(mu_0/eps_0)

% Simulation Parameters
dt = 1e-3;
Nz = 100;
N_Steps = 1000;

% Initialize Fields to zero
Ey = zeros(1,Nz);
Hx = zeros(1,Nz);

% Initialize Materials to Free Space
eps_r = ones(1,Nz);
mu_r = ones(1,Nz);

% Compute Update Coefficients
m_Ey = (c_0*dt)./eps_r;
m_Hx = (c_0*dt)./mu_r;


%% Main FDTD Loop

for T = 1:N_Steps;
  % Update H from E(Dirichlet Boundary Conditions)
  for nz = 1:Nz-1
    Hx(nz) = Hx(nz) + m_Hx(nz)*( Ey(nz+1) - Ey(nz) )/dz;
  endfor
  Hx(Nz) = Hx(Nz) + m_Hx(Nz)*( 0 - Ey(nz) )/dz; % Ey(Nz) = 0 [Perfect Electric Conductor at z-High Boundary]

  % Update E from H (Dirichlet Boundary Conditions)
  Ey(1) = Ey(1) + m_Ey(1)*( Hx(1) - 0 )/dz; % Hx(1) = 0 [Perfect Magnetic Conductor at z-Low Boundary]
  for nz = 2:Nz
    Ey(nz) = Ey(nz) + m_Ey(nz)*( Hx(nz) - Hx(nz-1) )/dz;
  endfor

endfor



