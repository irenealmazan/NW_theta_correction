% Phase retrieval algorithm which combines the retrieval of the object and
% the angles simultaneously. The variation of the object and the angles is
% based on the gradient of the error metric and how far we correct in the
% given direction is determined by a linear search algorithm 

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
    [scale_fact,err_scale_fact] = Phretrieval_functions.ini_guess_scalefactor(probe, rho_ini,angles_list,data_exp,250*mncntrate/mn,ki_o,kf_o,X,Y,Z);
    rho = rho_ini.*scale_fact .* support;
    
    % initial value of the gradient in rho and theta, assumed to be zero
    norm_grad_rho = zeros(Niter_rho,1);
    beta_rho = zeros(Niter_rho,1);
    norm_grad_theta = zeros(round(Niter_rho/freq_pos),1);
    beta_theta = zeros(round(Niter_rho/freq_pos),1);
    
    fprintf('initial  error: %4.4d \n',err_scale_fact);
    errlist = [min(err_scale_fact)];

end
 
DisplayResults.show_rho_theta_update(5,errlist,rho,midsl,angles_list,delta_thscanvals'+dth_disp,norm_grad_rho(1),beta_rho(1),norm_grad_theta(1:cnt_ntheta),beta_theta(1:cnt_ntheta),'Ini');

%% Iterative engine:

for nrho = 1:Niter_rho

    tic;
   err=0;
   fprintf('PIE iter %i: ', nrho);

    %RHO ITERATIONS

    if(1)

        [rho,beta_rho(nrho),norm_grad_rho(nrho)] = Phretrieval_functions.rho_update(probe, rho,angles_list,support, data_exp,depth,errlist(end),tau_backtrack_rho,beta_ini_rho,counter_max_rho,ki_o,kf_o,X,Y,Z);
        
        [err] = DiffractionPatterns.calc_error_multiangle(probe, rho, data_exp,angles_list,ki_o,kf_o,X,Y,Z);
        fprintf('\n     error: %4.4d \n', err);
        errlist = [errlist err];

        DisplayResults.show_rho_theta_update(5,errlist,rho,midsl,angles_list,delta_thscanvals'+dth_disp,norm_grad_rho(1:nrho),beta_rho(1:nrho),norm_grad_theta(1:cnt_ntheta),beta_theta(1:cnt_ntheta),'rho');


        % store the current reconstruction:
        rho_store(nrho).rho_square = rho(:,:,midsl);
        rho_store(nrho).rho_hex = squeeze(rho(100,:,:));
        rho_store(nrho).beta_rho = beta_rho(nrho);
    end

     % THETA ANNEALING
    %%{
    tic;
    if mod(nrho,freq_pos) == 0

        [angles_list,dq_shift,grad_final_theta,norm_grad_theta(cnt_ntheta),beta_theta(cnt_ntheta)] = Phretrieval_functions.theta_update(probe, rho,angles_list,data_exp,Niter_theta,errlist(end),tau_backtrack_theta,beta_ini_theta,counter_max_theta,ki_o,kf_o,X,Y,Z);

        % store the updated theta list
        for ii = 1:numel(data_exp)%index_to_distort%
           data_exp(ii).dqshift(:) = dq_shift(ii,:); 
           data_exp(ii).dth_iter = angles_list(ii);
           data_exp(ii).theta_iter(cnt_ntheta).grad_final = grad_final_theta(ii);
           data_exp(ii).theta_iter(cnt_ntheta).beta = beta_theta(cnt_ntheta);
           data_exp(ii).theta_iter(cnt_ntheta).dth_new_iter = angles_list(ii);
           data_exp(ii).theta_iter(cnt_ntheta).dqshift(:) = dq_shift(ii,:);
        end

        [err] = DiffractionPatterns.calc_error_multiangle(probe, rho, data_exp,angles_list,ki_o,kf_o,X,Y,Z);
        fprintf('     error: %4.4d \n', err);
        errlist = [errlist err];

        % plot
        DisplayResults.show_rho_theta_update(5,errlist,rho,midsl,angles_list,delta_thscanvals'+dth_disp,norm_grad_rho(1:nrho),beta_rho(1:nrho),norm_grad_theta(1:cnt_ntheta),beta_theta(1:cnt_ntheta),'theta');


        cnt_ntheta = cnt_ntheta + 1;
    end

    %}

end

   



