% This script contains the value of the experimental parameterse
% when the [1010] reflection was measured.

% experimental setup details
pixsize = 55; %microns, for Merlin
lam = etolambda(10400)*1e-4;

Npix = 516;
detdist = 0.33e6; % in micrometers
d2_bragg = detdist * lam /(Npix*pixsize);
depth = 100; 
defocus = 0;

% expected values of the motors for m-plane
th = -9.3;
del = -18.6; %in plane
gam = 0; %out of plane

% rocking curve scans
thscanvals =  [-9.42:.02:-9.13];
delta_thscanvals = thscanvals-th;



edgepad = 1.3; %for support during phase ret.
