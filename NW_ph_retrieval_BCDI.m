% NW_ph_retrieval_Irene_v3 tests the combination of PIE with the steepest descent method


% Iteration parameters:
Niter_rho = 2000;
Niter_pos = 1;
Niter_theta = 1;
freq_pos = 100;
freq_rho = 10;
freq_store = 10;


display(['set # rho iterations to ' num2str(Niter_rho) ' temp, position frequency to ' num2str(Niter_pos) 'per rho iteration'])

midsl = round(depth/2); % index of the section of the object represented in fig. 5
printind = round( [10:10:100]*(numel(data_exp)/100));

% support:
support = abs(NW);

% Is this the continuation of a phase retrieval process or does this start
% from scracth?
if flagContinue == 0
    cnt_ntheta = 1;
    cnt_store = 1;
    
    % initial list of angles:
    angles_list = zeros(numel(data_exp),1);
    for ii = 1:numel(data_exp)
        angles_list(ii) = data_exp(ii).dth_iter;
    end
    
    % initial guess and initial error:
    %rho_ini = rand(Npix,Npix,depth).* exp(i*2*pi*rand(Npix,Npix,depth));
    rho_ini = NW;
    [scale_fact,err_scale_fact] = Phretrieval_functions.ini_guess_scalefactor(probe, rho_ini,angles_list,data_exp,[350]*mncntrate/mn,ki_o,kf_o,X,Y,Z);
    rho = rho_ini.*scale_fact .* support;
    
    fprintf('initial  error: %4.4d \n',err_scale_fact);
    errlist = [min(err_scale_fact)];

end
 
show_rho_theta_update(5,errlist,rho,midsl,angles_list,delta_thscanvals'+dth_disp,1)

%% Iterative engine:

for nrho = 1:Niter_rho

    tic;
   err=0;
   fprintf('PIE iter %i: ', nrho);

    %RHO ITERATIONS

    if(1)

        [rho,beta_rho] = Phretrieval_functions.rho_update(probe, rho,angles_list,support, data_exp,depth,errlist(end),ki_o,kf_o,X,Y,Z);
        
        [err] = DiffractionPatterns.calc_error_multiangle(probe, rho, data_exp,angles_list,ki_o,kf_o,X,Y,Z);
        fprintf('\n     error: %4.4d \n', err);
        errlist = [errlist err];

        show_rho_theta_update(5,errlist,rho,midsl,angles_list,delta_thscanvals'+dth_disp,0)


        % store the current reconstruction:
        rho_store(nrho).rho_square = rho(:,:,midsl);
        rho_store(nrho).rho_hex = squeeze(rho(100,:,:));
        rho_store(nrho).beta_rho = beta_rho;
    end

     % THETA ANNEALING
    %%{
    tic;
    if mod(nrho,freq_pos) == 0

        [angles_list,dq_shift,grad_final_theta,beta_theta] = Phretrieval_functions.theta_update(probe, rho,data_exp,Niter_theta,th,errlist(end));

        % store the updated theta list
        for ii = 1:numel(data_exp)%index_to_distort%
           data_exp(ii).dqshift(:) = dq_shift(ii,:); 
           data_exp(ii).dth_iter = angles_list(ii);
           data_exp(ii).theta_iter(cnt_ntheta).grad_final = grad_final_theta(ii);
           data_exp(ii).theta_iter(cnt_ntheta).beta = beta_theta;
           data_exp(ii).theta_iter(cnt_ntheta).dth_new_iter = angles_list(ii);
           data_exp(ii).theta_iter(cnt_ntheta).dqshift(:) = dq_shift(ii,:);
        end

        [err] = DiffractionPatterns.calc_error_multiangle(probe, rho, data_exp,angles_list,ki_o,kf_o,X,Y,Z);
        fprintf('     error: %4.4d \n', err);
        errlist = [errlist err];

        % plot
        show_rho_theta_update(5,errlist,rho,midsl,angle_list,delta_thscanvals'+dth_disp,0)


        cnt_ntheta = cnt_ntheta + 1;
    end

    %}

end

   



