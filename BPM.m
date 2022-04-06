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
Nz       = 100;         % Number of steps in the z direction <---### EDIT HERE ###
Nx       = 512;         % x-direction size of computational grid
Ny       = Nx;          % x-direction size of computational grid. The computation domain is square


% Physical dimension of the computation space. Physical values are denoted with an underscore. The corresponding normalized value are written without underscore.
% We use SI units
Lx_ = 2e-3;     % width of the computation window [m]  <---### EDIT HERE ###
Ly_ = Lx_;		% height of the computation window [m] <---### EDIT HERE ###
Lz_ = 40e-3;    % propagation distance [m]             <---### EDIT HERE ###

n0_ = 1;        % linear refractive index of background    <---### EDIT HERE ###

lambda0_ = 532e-9;      % free space wavelength [m]        <---### EDIT HERE ###
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

%% LINEARITY, PARAXIAL APROXIMATION, BOUNDARY CONDITIONS

% Below, factor defining whether we want a medium that is linear (nonlinearity = 0), focusing (n2 > 0, nonlinearity = 1)
% or defocusing (n2 < 0, nonlinearity = -1).
nonlinearity = 0;   % Either -1 or 0 or 1 <---### EDIT HERE ###

% Variable that allows to switch between the paraxial (nonparaxial = 0) and nonparaxial (nonparaxial = 1) algorithm.
nonparaxial  = 0;   % <---### EDIT HERE ###

% For the side boundaries, if our field gets there, it will appear from
% other side due to periodic boundary condition from FFT. To prevent that,
% we can either make the computation window bigger if we have room in terms
% of time and computation power OR we can apply absorbing boundaries by
% lowering the intensity towords boundaries. However we have to do it
% smoothly, otherwise numerical reflections will appear and results would
% be messed-up. To do that, we will multiply the field at every step with a
% rectangular super gaussian function that is mostly 1 but goes to 0 
% towards boundaries in a smooth manner. 
% IF YOU WANT TO HAVE absorbing boundaries, set absorbing_boundary=1 below, 
% otherwise set it to 0.
absorbing_boundary = 1;
super_gaussian=exp(-((X_ / (0.85*Lx_/(2*sqrt(log(2)))) ).^20 + (Y_ / (0.85*Ly_/(2*sqrt(log(2)))) ).^20));
% figure,imagesc(super_gaussian)
%% DEFINITION OF INPUT FIELD
% Here, we define the complex input field either by generating here or
% loading from another file

I_ = 1e+10;                 % Peak intensity of the input beam [W/m2] <---### EDIT HERE ###
A_ = sqrt(2 * eta0 * I_);	% Peak field amplitude of the input beam [V/m]
A = A_ / field_scale;		% normalized peak field amplitude of the input beam

% Generate input field. Below, an example with Gaussian beam as the input
beam_fwhm_ = 0.05e-3;                        % <---### EDIT HERE ###
beam_scale_ = beam_fwhm_ / (2*sqrt(log(2)));
gaussian_beam = 1 .* exp(-(X_/beam_scale_).^2 - (Y_/beam_scale_).^2);
% figure, imagesc(abs(gaussian_beam)) ,colorbar

% Generate angle of the input beam. We do it by multiplying the input field 
% with a blazed grating
theta_x = 0  ;      % [degree] angle in x direction <---### EDIT HERE ###
theta_y = 0  ;      % [degree] angle in y direction <---### EDIT HERE ###
blazed_grating = exp(1i*2*pi/lambda_*(sind(theta_x)*X_+sind(theta_y)*Y_));

% Generate aperture for the input
rectangular_aperture = zeros(Nx,Ny);
width_x  = 200e-6;        % [m] <---### EDIT HERE ###
width_y  = 200e-6;        % [m]<---### EDIT HERE ###
x_shift = 0e-6;           % [m]<---### EDIT HERE ###
y_shift = 0e-6;           % [m]<---### EDIT HERE ###
rectangular_aperture(end/2+round(y_shift/Ly_*Ny)-round(width_y/Ly_/2*Ny)+1:end/2+round(y_shift/Ly_*Ny)+round(width_y/Ly_/2*Ny),end/2+round(x_shift/Lx_*Nx)-round(width_x/Lx_/2*Nx)+1:end/2+round(x_shift/Lx_*Nx)+round(width_x/Lx_/2*Nx))=1; 
% figure, imagesc(rectangular_aperture)


circular_aperture = zeros(Nx,Ny);
aperture_diameter  = 0;                % <---### EDIT HERE ###
circular_aperture(X_.^2+Y_.^2<=(aperture_diameter/2)^2)=1;
% figure, imagesc(circular_aperture)

% u is our input field. If you want to load a file, do not forget to assign it as u
% load('einsteins5x5');          % <---### EDIT & UNCOMMENT HERE ### if you want to load a input file
u = rectangular_aperture;        % <---### EDIT HERE ###

% Variables that allow us to monitor the progress of the algorithm.
progress_start = 0;
progress_diff = 100;
% Sometimes, the computation takes a long time (especially when you have a large discretized window)
% and you want to know how long it will take (i.e. if mankind will still exist when the computation ends).

% We declare this variable to monitor propagation (in the y-z plane)
fields = zeros(Nx,Ny,Nz);

% Assign the first step input value to be your input field
u0 = u;

%% PROPAGATION ROUTINE (CORE OF THE STORY)

% Now we are ready to propagate
tic         % to monitor how much time the propagation takes to simulate
count=0;
for index_z = 1:Nz
    count=count+1;
    if nonparaxial == 0
        % Paraxial code
        u1 = ifft2(fft2(u0) .* exp(-1i * K2 * 0.5 * dz));				% First linear half step
        u2 = u1 .* exp(1i * dz * (nonlinearity .* abs(u1).^2 + V));     % Nonlinear step
        u3 = ifft2(fft2(u2) .* exp(-1i * K2 * 0.5 * dz));				% Second linear step
    else
        % Nonparaxial code
        u1 = ifft2(fft2(u0) .* exp(-1i * K2 * 0.5 * dz * scale_ratio./ (k + sqrt(k^2 - K2))));
        u2 = u1 .* exp(1i * dz * (nonlinearity .* abs(u1).^2 + V));
        u3 = ifft2(fft2(u2) .* exp(-1i * K2 * 0.5 * dz * scale_ratio ./ (k + sqrt(k^2 - K2))));
    end
    
    % Monitor the progress
    % For MATLAB:
    fprintf('\nProgress: %f %% \n\n', progress_start + progress_diff * index_z / Nz);
    % For Octave:
    % printf('\nProgress: %f %% \n', progress_start + progress_diff * index_z / Nz)
    % fflush(stdout);
    
    % Get ready for the next step
    if absorbing_boundary==1
        u0 = u3.*super_gaussian;
    else
        u0 = u3;
    end
    
    fields(:,:,count)=u3;
end
toc
% Let's store the result in a variable with a more explicit name
output = u3;

%% DISPLAY RESULTS

% This monitors the physical input power of your beam
input_physical_power = n0_/(2 * eta0) * sum(sum(abs(u*field_scale).^2)) *dx_*dy_;

% Now you can plot your data (we convert the distances into millimeters for clearer display):

% Input object:
figure,
imagesc(1000*x_,1000*y_,abs(u))
set(colorbar,'FontName','Times New Roman','FontSize',12,'LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1)
xlabel('x coordinate [mm]')
ylabel('y coordinate [mm]')
title('Abs of the input field')
axis image
% Output field:
figure,imagesc(1000*x_,1000*y_,(abs(output)))
set(colorbar,'FontName','Times New Roman','FontSize',12,'LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',12,'LineWidth',1)
xlabel('x coordinate [mm]')
ylabel('y coordinate [mm]')
title('Abs of the output field')
axis image
%print(gcf,'output.png','-dpng','-r300');   % save image with 300 dpi resolution 

% Propagation profile in the middle of the transverse window:
profile=squeeze(fields(:,end/2,:));
figure,imagesc(1000*z_,1000*x_,abs(profile))
set(colorbar,'FontName','Times New Roman','FontSize',20,'LineWidth',2)
set(gca,'FontName','Times New Roman','FontSize',20,'LineWidth',2)
xlabel('z coordinate [mm]')
ylabel('y coordinate [mm]')
title('Propagation profile')
%print(gcf,'propagation.png','-dpng','-r300');    % save image with 300 dpi resolution 


%% MAKE YOUR OWN PROPAGATION MOVIE
% WAIT UNTIL the figure closes itself and you see that 'video is ready' line appears in the command window 
% to save the video
% vidObj = VideoWriter('propagation','MPEG-4');
% vidObj.FrameRate = 10;
% vidObj.Quality = 100;
% open(vidObj);
% 
% scrsz = get(0,'ScreenSize'); % get screen size of your computer
% figure('Position',scrsz,'MenuBar','none','ToolBar','none','resize','off') % fullscreen
% 
% for k=1:Nz
%     A=zeros(length(x_),length(z_));
%     A(1:length(x_)/4,k)=max(max(abs(profile)));
%     
%     subplot(2,1,1),imagesc(1000*x_' , 1000*y_',interp2(abs(fields(:,:,k)))), colorbar, caxis([min(min(abs(profile))) max(max(abs(profile)))])
%     xlabel('x coordinate [mm]')
%     ylabel('y coordinate [mm]')
%     title('Abs of the field in XY plane')
%     axis image
% 
%     subplot(2,1,2),imagesc(1000*z_' , 1000*x_', (abs(profile)+A)); colorbar, caxis([min(min(abs(profile))) max(max(abs(profile)))])
%     xlabel('z coordinate [mm]')
%     ylabel('y coordinate [mm]')
%     title('Abs of the field in YZ plane')
%   
%     pause(0.05)
%     frame = getframe(gcf);
%     writeVideo(vidObj,frame);  
% end
% close
% close(vidObj);
% fprintf('\nvideo is ready\n');