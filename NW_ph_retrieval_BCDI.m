% NW_ph_retrieval_Irene_v3 tests the combination of PIE with the steepest descent method


% Iteration parameters:
Niter_rho = 2000;
Niter_pos = 1;
Niter_theta = 1;
freq_pos = 1;
freq_rho = 10;
freq_store = 10;


display(['set # rho iterations to ' num2str(Niter_rho) ' temp, position iterations to ' num2str(Niter_pos) ' and cyles ' num2str(Ncycles)])

midsl = round(depth/2); % index of the section of the object represented in fig. 5
printind = round( [10:10:100]*(numel(data_exp)/100));

% support:
support = abs(NW);

% Is this the continuation of a phase retrieval process or does this start
% from scracth?
if flagContinue == 0
    cnt_ntheta = 1;
    cnt_store = 1;
    
    % initial guess and initial error:
    rho_ini = rand(Npix,Npix,depth).* exp(i*2*pi*rand(Npix,Npix,depth));
    [scale_fact,err_scale_fact] = Phretrieval_functions.ini_guess_scalefactor(probe, rho_ini, data_exp,[1]);
    scale_fact = 1;%1.5e-4; % scale factor for random ini
    rho = rho_ini.*scale_fact .* support;
    
    fprintf('initial  error: %4.4d \n',err_scale_fact);
    errlist = [errlist err_scale_fact];
end
 
figure(5); clf; setfigsize(gcf, 1000,500); pause(.1);

%% Iterative engine:

for nrho = 1:Niter_rho

    tic;
   err=0;
   fprintf('PIE iter %i: ', nrho);

    %RHO ITERATIONS

    if(1)

        [rho] = Phretrieval_functions.rho_update(probe, rho, data_exp,depth);

        [err] = DiffractionPaterns.calc_error_multiangle(probe, rho, data_exp);
        fprintf('\n     error: %4.4d \n', err);
        errlist = [errlist err];

        % plot
        subplot(131); imagecomp(rho(:,:,midsl)); colorbar; axis image; %zoom(1.5);
        subplot(132); imagecomp(squeeze(rho(100,:,:))); colorbar; axis image; %zoom(1.5);
        subplot(133); plot(log10(errlist),'LineWidth',3.0);
        drawnow;

        % store the current reconstruction:
        rho_store(nrho).rho_square = rho(:,:,midsl);
        rho_store(nrho).rho_hex = squeeze(rho(100,:,:));
        rho_store(nrho).beta_rho = beta_rho;
    end

     % THETA ANNEALING
    %%{
    tic;
    if mod(nrho,freq_pos) == 0

        [dth_new,dq_shift,grad_final_theta,beta_theta] = Phretrieval_functions.theta_correction(probe, rho,data_exp,Niter_theta,index_to_distort,del/2,errlist(end));

        % store the shift
        for ii = 1:numel(data_exp)%index_to_distort%
           data_exp(ii).dqshift(:) = dq_shift(ii,:); 
           data_exp(ii).dth_iter = dth_new(ii);
           data_exp(ii).theta_iter(cnt_ntheta).grad_final = grad_final_theta(ii);
           data_exp(ii).theta_iter(cnt_ntheta).beta = beta_theta;
           data_exp(ii).theta_iter(cnt_ntheta).dth_new_iter = dth_new(ii);
           data_exp(ii).theta_iter(cnt_ntheta).dqshift(:) = dq_shift(ii,:);
        end

        [err] = DiffractionPaterns.calc_error_multiangle(probe, rho, data_exp);
        fprintf('     error: %4.4d \n', err);
        errlist = [errlist err];

        % plot
        figure(5);subplot(132); plot(log10(errlist));drawnow;

        cnt_ntheta = cnt_ntheta + 1;
    end

    %}

end

   



