% This script contains the value of the experimental parameterse
% when the [2110] reflection was measured.

% experimental setup details
pixsize = 55; %microns, for Merlin
lam = etolambda(10400)*1e-4;
 
Npix = 200;
detdist = 0.529e6; % in micrometers
d2_bragg = detdist * lam /(Npix*pixsize);
depth = 60;%60;
defocus = 0;

% expected values of the motors for m-plane
th = 73.3;
del = -32.6; %in plane
gam = 0; %out of plane


% rocking curve scans
thscanvals =  [72.8:1.0/64:73.8-1.0/64];;
delta_thscanvals = thscanvals-th;
    

 
  
   
    

    

    mncntrate = 0.1;
    %usesimI = 1;