% In this script we initialize the values of the experimental set-up and
% the details concerning the sample:



%% sample details 

%NW details
corewidth= .2; %edge-to-edge distance %NW 22
%corewidth = 2*cosd(50)* .087; %take estimate from SEM image
NW_length = 2;
phi=0; %misalignment of NW from vertical

% real space
a_latparam = 4.24e-4; %for InGaAs, 14.3% Ga fraction
c_latparam = 6.93e-4; 
mplane_spacing = sind(60)*a_latparam;
aplane_spacing = a_latparam/2;



% reciprocal space 
q_mplane = 2*pi/(mplane_spacing*1e4);
q_aplane = 2*pi/aplane_spacing *1e-4;
q_cplane = 2*pi/c_latparam * 1e-4;
%th_mplane110 = asind(lam/(2*mplane_spacing));
%th_mplane220 = asind(2*lam/(2*mplane_spacing));
%th_aplane112 = asind(lam/(2*a_latparam/2));


%% ZP details
zpdiam = 240;
outerzone = .04;
bsdiam = 50;
binaryprobe_flag =0; %use binary probe defined at FWHM of probe

%cutrad = .12;
meshdata = 0;
cutrad = .06;
edgepad = 1.25;
mncntrate = 1;

%% Phase retrieval parameters:

% Iteration parameters:
Niter_rho = 2000;
Niter_pos = 1;
Niter_theta = 1;
freq_pos = 1;
%freq_rho = 10;
%freq_store = 10;

% Beta adaptative step parameters:

tau_backtrack_rho = 0.1;
beta_ini_rho = 1;
counter_max_rho = 10;

tau_backtrack_theta = 0.1;
beta_ini_theta = 1e-5;
counter_max_theta = 15;



