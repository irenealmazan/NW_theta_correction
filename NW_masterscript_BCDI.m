global d2_bragg X Y Z ki_o kf_o

warning off;


addpath(genpath('/Users/ialmazn/Box Sync/Nanowire_ptychography/NSLS II/NSLS II March 2017/Analysis_end_of_beamtime'));
addpath(genpath('/Users/ialmazn/Documents/MATLAB/ptycho/m_scripts/'))
addpath(genpath('./calc_functions'));
addpath(genpath('./display_functions'));

%% Flags:
noiseflag = 1; 
if(noiseflag) display('ADDING NOISE');end
addNWstrain = 0; 
if(addNWstrain) display('ADDING STRAIN FIELD');end
addNWsf = 0; 
if(addNWsf) display('ADDING STACKING FAULTS');end
plotdqshift = 0; 
if(plotdqshift) display('PLOTTING DQ ROCKING CURVE'); end
usesimI = 1; 
if(usesimI) display('USING SIMULATED DATA');end
flagContinue = 0;
if flagContinue == 0
    display('STARTING A NEW PHASE RETRIEVAL OPERATION');

    %% initialize parameters relative to experimental setup and sample
    NW_initial_and_experimental_variables;

    %% Scattering condition:
    scatgeo = 2110; %for strain image
    %scatgeo = 1010; %for stacking-fault image
    switch scatgeo   
        case 1010 %SF
            NW_scatgeo_1010;
        case 2110 %strain
            NW_scatgeo_2110;
    end

    %% Create sample
    NW_diff_vectors_BCDI; % does both the vectors ki and kf and creates the object
    NW_make_InGaAs_nocoreshell_BCDI;
    NW_plot_diff_vectors_sample_BCDI;
    probe = ones(size(X));

    %% Calculate diffraction patterns
    NW_calc_dp_BCDI;
    NW_add_dp_noise;

    else
        display('CONTINUING A PHASE RETRIEVAL OPERATION');   
end

%% Phase retrieval algorithm
NW_ph_retrieval_BCDI;
