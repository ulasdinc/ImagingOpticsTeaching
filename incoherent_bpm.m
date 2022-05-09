%% Split-step Fourier beam propagation method (SSF-BPM)

% For nonlinear propagation in linear or chi-3 media, according
% to the nonlinear Schrodinger equation (NLSE).

% To play with the code, you normally only need to edit the values or lines
% labeled with the following tag: <---### EDIT HERE ###

%% SETUP THE ENVIRONMENT AND DEFINE THE COMPUTATION WINDOW

clear all;	% clear all variables in current environment
close all;	% close all figures
clc

% Computation domain discretization
Nz       = 200;         % Number of steps in the z direction <---### EDIT HERE ###
Nx       = 512;         % x-direction size of computational grid
Ny       = Nx;          % x-direction size of computational grid. The computation domain is square
Navg     = 200;


% Physical dimension of the computation space. Physical values are denoted with an underscore. The corresponding normalized value are written without underscore.
% We use SI units
Lx_ = 200e-6;       % width of the computation window [m]  <---### EDIT HERE ###
Ly_ = Lx_;          % height of the computation window [m] <---### EDIT HERE ###
Lz_ = 800e-6;       % propagation distance [m]             <---### EDIT HERE ###

n0_ = 1;        % linear refractive index of background    <---### EDIT HERE ###

lambda0_ = 400e-9;      % free space wavelength [m]        <---### EDIT HERE ###
delta_   = 1.0;			% normalization parameter (see documentation on SSF)
n2_      = 2.4e-19;	    % nonlinear coefficient [m2/W]     <---### EDIT HERE ###
V = zeros(Ny, Nx);      % Index potential. This correspond to the refractive index difference with respect to background
% for a homogeneous medium, V = 0. <---### EDIT HERE ###

%% SETUP THE SSF-BPM VARIABLES
% Normally, you shouldn't need to edit this section at all

% Physical constants
mu0 = 4.0e-7 * pi;       % free space magnetic permeability [Vs/Am]
c0  = 2.99792458e+8;     % free space light speed [m/s]

epsilon0 = 1.0 / (mu0 * c0^2);      % free space permittivity [As/Vm]
eta0     = sqrt(mu0 / epsilon0);    % free space impedance [ohm]

% Derived parameters
n2_el   = n0_ * n2_ / (2 * eta0);       % nonlinear refractive index [m2/V2]
k0_     = 2 * pi / lambda0_;			% free space wavenumber [m-1]
k_      = n0_ * k0_;                    % medium wavenumber [m-1]
lambda_ = lambda0_ / n0_;               % medium wavelength [m]

% Normalization coefficients
% The equation can be normalized to a dimensionless form
% spatial normalization factor in the x-y plane
spatial_transverse_scale = 1/(k0_ * sqrt(2 * n0_ *  delta_));
% spatial normalization factor in the z direction
spatial_longitudinal_scale = 1/(delta_ * k0_);

scale_ratio = spatial_longitudinal_scale/spatial_transverse_scale; % = sqrt(2*n0_/delta_)
% normalization factor for the electric field
field_scale = sqrt(delta_ / n2_el);

% ************* Normalized parameters *************
Lx = Lx_ / spatial_transverse_scale;             % normalized model width
Ly = Ly_ / spatial_transverse_scale;             % normalized model height
Lz = Lz_ / spatial_longitudinal_scale;           % normalized propagation distance
k  = 2*pi * spatial_transverse_scale / lambda_;  % normalized light k-vector

% ************ Numeric model parameters ***********
dx_ = Lx_/Nx;                       % normalized discretization step in x
dx  = Lx/Nx;                        % discretization step in x
x_  = dx_ * (-Nx/2:Nx/2-1);         % x dimension vector
x   = dx * (-Nx/2:Nx/2-1);          % normalized x dimension vector
dkx = 2*pi/Lx;                      % discretization in the spatial spectral domain along the y direction
kx  = dkx * [0:Nx/2-1, -Nx/2:-1];	% spatial frequencies vector in the x direction
% Note that we swap the left and right part of this vector because the FFT algorithm puts the
% low frequencies in the edges of the transformed vector and the high frequencies in the middle.
% There is a fftshift function that perform this swap in MATLAB. Here we prefer to redefine k so
% that we don't need to call fftshift everytime we use fft.

% We do the same in the y and z direction
dy  = Ly/Ny;
dy_ = Ly_/Ny;
y   = dy * (-Ny/2:Ny/2-1)';
y_  = dy_ * (-Ny/2:Ny/2-1)';
dky = 2*pi/Ly;
ky  = dky * [0:Ny/2-1 -Ny/2:-1]';

dz  = Lz/Nz;
dz_ = Lz_/Nz;
z   = dz * (1:Nz);
z_  = dz_ * (1:Nz);

% Here we create the spatial computation grid (physical and normalized)
[X_, Y_]    = meshgrid(x_, y_);
[Xz_, Z_]   = meshgrid(x_, z_);
[X, Y]      = meshgrid(x, y);
[Xz, Z]     = meshgrid(x, z);

% The same for the spatial frequencies domain
[Kx, Ky]    = meshgrid(kx, ky);

K2 = Kx.^2 + Ky.^2; % Here we define some variable so that we don't need to compute them again and again

%% DEFINITION OF INPUT FIELD AND LENSES
% Here, we define the complex input field either by generating here or
% loading from another file

I_ = 1e+10;                 % Peak intensity of the input beam [W/m2] <---### EDIT HERE ###
A_ = sqrt(2 * eta0 * I_);	% Peak field amplitude of the input beam [V/m]
A = A_ / field_scale;		% normalized peak field amplitude of the input beam


circular_aperture  = zeros(Nx,Ny);

aperture_diameter  = 10e-6;                        
circular_aperture(X_.^2+Y_.^2<=(aperture_diameter/2)^2)=1;

%figure,imagesc(circular_aperture)

F_lens=Lz_/4;
lens=exp(-1i*pi/lambda0_/F_lens*(X_.^2+Y_.^2));

grating_period=12e-6;
amplitude_grating=0.5+0.5*cos(2*pi*X_/grating_period);

square_aperture = zeros(Nx,Ny);
aperture_width  = 150e-6;                   % <---### EDIT HERE ###
square_aperture(end/2-round(aperture_width/Lx_/2*Nx)+1:end/2+round(aperture_width/Lx_/2*Nx),end/2-round(aperture_width/Ly_/2*Ny)+1:end/2+round(aperture_width/Ly_/2*Ny))=1; 
% figure, imagesc(square_aperture)

u=amplitude_grating.*square_aperture;

% We declare this variable to monitor propagation (in the y-z plane)
profile=zeros(Nx,Nz);

% Variables that allow us to monitor the progress of the algorithm
progress_start = 0;
progress_diff = 100;
% Sometimes, the computation takes a long time (especially when you have a large discretized window)
% and you want to know how long it will take (i.e. if mankind will still exist when the computation ends).
%% PROPAGATION ROUTINE (CORE OF THE STORY)

output_each=zeros(Nx,Ny,Navg);
output=zeros(Nx,Ny);
input_each=zeros(Nx,Ny,Navg);

u_in = u;
% Now we are ready to propagate
tic         % to monitor how much time the propagation takes to simulate
count=0;
scrsz = get(0,'ScreenSize'); % get screen size of your computer
figure('Position',scrsz/1.5,'MenuBar','none','ToolBar','none','resize','off') 
for index_avg = 1:Navg
 
    count=count+1;
    
    Noise                         = imresize(randn(20,20),[Nx,Ny]);%wgn(Ny, Nx, 0);
    incident                      = (Noise) .* u_in;
    input_each(:,:,index_avg)     = abs(incident).^2;
    u0        = incident;%ifftshift(incident);
    
    % Paraxial code
    u0 = ifft2(fft2(u0) .* exp(-1i * K2 * 1 * dz * Nz/4));	
    u0 = u0.*lens;
    u0 = ifft2(fft2(u0) .* exp(-1i * K2 * 1 * dz * Nz/4));		
    u0 = u0.*circular_aperture;
    u0 = ifft2(fft2(u0) .* exp(-1i * K2 * 1 * dz * Nz/4));
    u0 = u0.*lens;
    u0 = ifft2(fft2(u0) .* exp(-1i * K2 * 1 * dz * Nz/4));
    
    
    % Monitor the progress
    % For MATLAB:
    fprintf('\nProgress: %f %% \n\n', progress_start + progress_diff * index_avg / Navg);
    % For Octave:
    % printf('\nProgress: %f %% \n', progress_start + progress_diff * index_z / Nz)
    % fflush(stdout);
    output_each(:,:,index_avg)=u0;
    output      = output + abs(u0).^2;
    
    

    %figure(2)
    subplot(1, 2, 2), imagesc(1e+6*x_,1e+6*y_,output),axis image, colormap gray
    title('Intensity of the output building up')
    xlabel('x coordinate [\mum]')
    ylabel('y coordinate [\mum]')
    subplot(1, 2, 1), imagesc(1e+6*x_,1e+6*y_,abs(u0).^2),axis image, colormap gray
    xlabel('x coordinate [\mum]')
    ylabel('y coordinate [\mum]')    
    title(['Intensity of current iteration:',num2str(index_avg)])
    drawnow;

end
toc


%% DISPLAY RESULTS

figure, imagesc(1e+6*x_,1e+6*y_,(abs(u).^2))
set(colorbar,'FontName','Times New Roman','FontSize',20,'LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',1)
xlabel('x coordinate [\mum]')
ylabel('y coordinate [\mum]')
title('Intensity of the input')
axis image
% print(gcf,'input.png','-dpng','-r300');    % save image 


figure, imagesc(1e+6*x_,1e+6*y_,(output))
set(colorbar,'FontName','Times New Roman','FontSize',20,'LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',1)
xlabel('x coordinate [\mum]')
ylabel('y coordinate [\mum]')
title('Intensity of the output')
axis image
% print(gcf,'output.png','-dpng','-r300');    % save image 