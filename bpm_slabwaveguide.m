%% Split-step Fourier beam propagation method (SSF-BPM)

% For nonlinear propagation in chi-3 media, according
% to the nonlinear Schr??dinger equation (NLSE).
% We use SI units

% To play with the code, you nromally only need to edit the values or lines
% labeled with the following tag: <---### EDIT HERE ###

%% SETUP THE ENVIRONMENT AND DEFINE THE COMPUTATION WINDOW

%clear all;	% clear all variables in current environment
close all;	% close all figures

% Computation domain discretization
% Note that for optimal use of fast Fourier transform (FFT),
% we take the size of the computation window to be a power of two.

exponent = 11;          % <---### EDIT HERE ### 
Nz = 500;               % Number of steps in the z direction <---### EDIT HERE ###
Nx = 2^exponent;        % x-direction size of computational grid

% Physical dimension of the computation space. Physical values are denoted with an underscore. The corresponding normalized value are written without underscore. 

Lx_ = 30e-6;		% width of the computation window [m]  <---### EDIT HERE ###
Lz_ = 100e-6;	% propagation distance [m]             <---### EDIT HERE ###

lambda0_ = 800e-9;	% free space wavelength [m] <---### EDIT HERE ###
delta_ = 1.0 ;			% normalization parameter (see documentation on SSF)
n2_ = 2.5e-21;			% nonlinear coefficient [m2/W] (2.4e-19m2/W for acetone) <---### EDIT HERE ###
n0_ = 1.5    ;          % n_cladding or background medium

%% SETUP THE SSF-BPM VARIABLES
% Physical constants

mu0 = 4.0e-7 * pi;	% free space magnetic permeability [Vs/Am]
c0 = 2.99792458e+8;	% free space light speed [m/s]

epsilon0 = 1.0 / (mu0 * c0^2);	% free space permittivity [As/Vm]
eta0 = sqrt(mu0 / epsilon0);		% free space impedance [ohm]

% Derived parameters

n2_el = n0_ * n2_ / (2 * eta0);	% nonlinear refractive index [m2/V2]
k0_ = 2 * pi / lambda0_;			% free space wavenumber [m-1]
k_ = n0_ * k0_;						% medium wavenumber [m-1]
lambda_ = lambda0_ / n0_;			% medium wavelength [m]

% Normalization coefficients
% The equatsion can be normalized to a dimensionless form

% spatial normalization factor in the x-y plane
spatial_transverse_scale = 1/(k0_ * sqrt(2 * n0_ *  delta_));
% spatial normalization factor in the z direction
spatial_longitudinal_scale = 1/(delta_ * k0_);

scale_ratio = spatial_longitudinal_scale/spatial_transverse_scale; % = sqrt(2*n0_/delta_)
% normalization factor for the electric field
field_scale = sqrt(delta_ / n2_el);

% ************* Normalized parameters *************

Lx = Lx_ / spatial_transverse_scale;	% normalized model width
Lz = Lz_ / spatial_longitudinal_scale;	% normalized propagation distance
k = 2*pi * spatial_transverse_scale / lambda_; % normalized light k-vector

% ************ Numeric model parameters ***********

dx_ = Lx_/Nx;		% normalized discretization step in x
dx = Lx/Nx;			% discretization step in x
x_ = dx_ * (-Nx/2:Nx/2-1);	% x dimension vector
x = dx * (-Nx/2:Nx/2-1);	% normalized x dimension vector
dkx = 2*pi/Lx;		% discretization in the spatial spectral domain along the y direction
kx = dkx * [0:Nx/2-1, -Nx/2:-1];	% spatial frequencies vector in the x direction
% Note that we swap the left and right part of this vector because the FFT algorithm puts the
% low frequencies in the edges of the transformed vector and the high frequencies in the middle.
% There is a fftshift function that perform this swap in MATLAB. Here we prefer to redefine k so
% that we don't need to call fftshift everytime we use fft.

dz = Lz/Nz;
dz_ = Lz_/Nz;
z = dz * (1:Nz);
z_ = dz_ * (1:Nz);

% Here we create the spatial computation grid (physical and normalized)
[Xz_, Z_] = meshgrid(x_, z_);
% The same for the spatial frequencies domain
K2 = kx.^2; % Here we define some variable so that we don't need to compute them again and again

%% WAVEGUIDE
% Define waveguide

Length_wg=1*Lz_;
d = 1e-6;

n_core = 1.53;
n_cladding = n0_;

NA= sqrt(n_core^2 - n_cladding^2);
v = NA*(pi*d/lambda0_) ;     % V mumber of the waveguide

V = (n_core - n_cladding)* (double(abs(x_) < d/2));         

x_clad_1 = [];
x_clad_2 = [];
x_core = [];

for ind = 1:Nx/2
    if V(ind) < abs(n_core - n_cladding)
        x_clad_1 = [x_clad_1, x_(ind)];
    end
end

for ind = Nx/2:Nx
    if V(ind) < abs(n_core - n_cladding)
        x_clad_2 = [x_clad_2, x_(ind)];
    end
end

for ind = 1:Nx
    if V(ind) > 0
        x_core = [x_core, x_(ind)];
    end
end

%% CALCULATION OF THE PROPAGATION CONSTANTS AND THE WAVEGUIDE MODES

%plot the circle 
f=@(h,q) h^2+q^2-v^2;
b1=ezplot(f,[0 10 0 10]);
set(b1,'color','red')
axis square;

g=@(h,q) h*tan(h)-q;
hold on;
b2=ezplot(g,[0 4*pi 0 10]);
set(b2,'color','blue')

m=@(h,q) -q- h*cot(h);
hold on;
b3=ezplot(m,[0 4*pi 0 10]);
title('Graphic Solution of transverse electric field')
xlabel('u');
ylabel('v');
set(b3,'color','green')
legend('v number','TE','TM')

%%%% EVEN MODES
Number = 3;         % Number of iterations 
Se = zeros(Number,1);

syms x 
fe = abs(x)*tan(abs(x))-sqrt(v^2-x^2);
for n = 1:Number 
  Se(n,1) = vpasolve(fe,x, [0 5],'random',true);
  fd = unique(fe,'sorted');
end
Su_e = unique(Se);      % Unique solutions of the equation fe

up_e = Su_e;            % Parameters u and v
vp_e = sqrt(v^2-up_e.^2);

qp_e = 2*vp_e(1)/d;        % Wavenumbers q and h 
hp_e = 2*up_e(1)/d;

beta_ = sqrt((n_core * k0_)^2 - hp_e^2);

B_e = 1;
D_e = cos(hp_e*d/2)/exp(-qp_e*d/2)*B_e;
C_e = D_e;

mode_0e = [D_e * exp(qp_e * x_clad_1), B_e * cos(hp_e .* x_core), C_e * exp(-qp_e * x_clad_2)];
beam_e = mode_0e;

%%%% ODD MODES

if v>pi/2
Number = 3;
So = zeros(Number,1);
syms x 
fo = abs(x)*cot(abs(x))+sqrt(v^2-x^2);
for n = 1:Number 
  So(n,1) = vpasolve(fo,x, [0 5],'random',true);
  fdo = unique(fo,'sorted');
end
Su_o = unique(So);      % Unique solutions of the equation fo

up_o = Su_o;
vp_o = sqrt(v^2-up_o.^2);


qp_o = 2*vp_o/d;
hp_o = 2*up_o/d;

beta_o = sqrt((n_core * k0_)^2 - hp_o^2);

B_1 = 1;
D_1 = sin(hp_o*d/2)/exp(-qp_o*d/2)*B_1;
C_1 = D_1;

mode_1 = [-D_1 * exp(qp_o * x_clad_1), B_1 * sin(hp_o .* x_core), C_1 * exp(-qp_o * x_clad_2)];
beam_o = mode_1;


figure(),plot(x_,(beam_o)); hold on
         plot(x_,(beam_e))
         
end
%figure, plot(x_,(beam_e))

%% DEFINITION OF INPUT FIELD

% Factor defining whether we want a medium that is linear (nonlinearity = 0), focusing (n2 > 0, nonlinearity = 1)
% or defocusing (n2 < 0, nonlinearity = -1).
nonlinearity = 0;           % <---### EDIT HERE ###   
% Variable that allows to switch between the paraxial (nonparaxial = 0) and nonparaxial (nonparaxial = 1) algorithm.
nonparaxial = 0;            % <---### EDIT HERE ###
% Variable for absorbin boundaries (if 1, it exists. It does not exist otherwise)
absorbing_boundary=1;

I_ = 1e+16;                 % Peak intensity of the input beam [W/m2] <---### EDIT HERE ###
A_ = sqrt(2 * eta0 * I_);	% Peak field amplitude of the input beam [V/m]
A = A_ / field_scale;		% normalized peak field amplitude of the input beam

beam_fwhm_ = 1e-6;%1.37e-6;
beam_scale_ = beam_fwhm_ / (2*sqrt(log(2)));
gaussian_beam = exp(-(x_/beam_scale_).^2);
%figure,plot(x_,gaussian_beam)

rect_input=zeros(1,Nx);
rect_input(x_>=-d/2)=1;
rect_input(x_>=d/2)=0;
%figure,plot(x_,rect_input)

window = exp(-(x_/(0.5*Lx_)).^20);
%figure,plot(x_,window)

u = beam_e; 
%u = gaussian_beam;
%u=rect_input;


% Variables that allow you to monitor the progress of the algorithm. 
progress_start = 0;	% This means you start at 0%
progress_diff = 100;		% This means that the following procedure will constitute 100% of the execution time
% Note that these values may change if you perform several propagation in the same code.
% Sometimes, the computation takes a long time (especially when you have a large discretized window)
% and you want to know how long it will take (i.e. if mankind will still exist when the computation ends).

% You declare this variable if you want to record the propagation profile (in the y-z plane)

% Assign the first step input value to be your input field
u0 = u;
V_ = V;
%% PROPAGATION ROUTINE (CORE OF THE STORY)
profile  = zeros(Nx,Nz);
profile_v= zeros(Nx,Nz);
% Now we are ready to propagate
for index_z = 1:Nz

	u1 = ifft(fft(u0) .* exp(-1i * K2 * dz * scale_ratio./ (k + sqrt(k^2 - K2))));
 	if index_z*dz_ <= Length_wg       
		u3 = u1 .* exp(1i * dz * (nonlinearity .* abs(u1).^2 + V_));
		profile_v(:,index_z) = V_.';          % [profile_v V.'];
	else
		profile_v(:,index_z) = zeros(Nx, 1); % [profile_v zeros(Nx, 1)];
		u3 = u1;
	end
	
	% Monitor the progress
	% For MATLAB:
	fprintf('\nProgress: %f %% \n\n', progress_start + progress_diff * index_z / Nz);
	% For Octave:
	% printf('\nProgress: %f %% \n', progress_start + progress_diff * index_z / Nz)
	% fflush(stdout);
	
    % Get ready for the next step
    if absorbing_boundary==1
        u0 = u3.*window;
    else
        u0 = u3;
    end
	
	% Record the profile (note that you can record it along a different plane)
	profile(:,index_z) = u3.';%[profile u3.'];
end

% Let's store the result in a variable with a more explicit name
output = u3;

%% DISPLAY RESULTS

figure,imagesc(1e+6*z_,1e+6*x_,profile_v)
set(colorbar,'FontName','Times New Roman','FontSize',16,'LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',1)
xlabel('z coordinate [\mum]')
ylabel('x coordinate [\mum]')
title('slab waveguide')
print(gcf,'input.png','-dpng','-r300');    % save image 

figure, plot(1e+6*x_,abs(u));
set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',1)
xlabel('x coordinate [\mum]')
ylabel('abs of the input')
title('input field')

figure, imagesc(1e+6*z_,1e+6*x_,abs(profile));
set(colorbar,'FontName','Times New Roman','FontSize',16,'LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',1)
xlabel('z coordinate [\mum]')
ylabel('x coordinate [\mum]')
title('abs of the propagation profile')

figure, imagesc(1e+6*z_,1e+6*x_,angle(profile));
set(colorbar,'FontName','Times New Roman','FontSize',16,'LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',1)
xlabel('z coordinate [\mum]')
ylabel('x coordinate [\mum]')
title('Phase of the propagation profile')

figure, plot(1e+6*x_,abs(output));
set(gca,'FontName','Times New Roman','FontSize',16,'LineWidth',1)
xlabel('x coordinate [\mum]')
ylabel('abs of the output')
title('output field')

input_physical_power = n0_/(2 * eta0) * sum(abs(u .* field_scale).^2) *dx_;
output_physical_power = n0_/(2 * eta0) * sum(abs(output .* field_scale).^2) *dx_;