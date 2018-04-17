
clear data_sim data_exp

% initialize
dth_disp = zeros(numel(delta_thscanvals),1);
rock_curve = zeros(numel(delta_thscanvals),1);
dq_shift_nominal = zeros(numel(delta_thscanvals),3);
dq_shift_real = zeros(numel(delta_thscanvals),3);
mxI = zeros(size(delta_thscanvals));
im_sum_sim = zeros(Npix);

% distortion:
index_to_distort = [34:numel(delta_thscanvals)];%[1:numel(fly2Danglist)];%[28:1:45];%[28 31 35 36];%[36];%[28:1:33];%randi(numel(fly2Danglist),1,number_angles_distort);
dth_disp(index_to_distort) = [0.002];%[0.008];%[0.017 -0.008 -0.013 0.005];%[0.017];% -0.008 -0.013 0.005];%[0.002];%[-.003:.001:.003];

 [dq_shift_nominal] = DiffractionPatterns.calc_dqshift_for_given_th(delta_thscanvals,ki_o,kf_o,kf_o-ki_o);
     
 [dq_shift_real] = DiffractionPatterns.calc_dqshift_for_given_th(delta_thscanvals + dth_disp',ki_o,kf_o,kf_o-ki_o);
    
 [ simI,rock_curve,Proj_vol,FT_Proj_vol] = DiffractionPatterns.calc_dp(dq_shift_real,probe,NW,X,Y,Z);

figure(27);clf;
for ii = 1:numel(delta_thscanvals)
    
     mxI(ii) = max(max(squeeze(simI(ii,:,:))));
     
     im_sum_sim = im_sum_sim + squeeze(simI(ii,:,:));    
    
     % store angles, dqshifts and diffraction patterns in structure
    data_exp(ii).dth_real = delta_thscanvals(ii)+dth_disp(ii);
    data_exp(ii).dth_nominal = delta_thscanvals(ii);
    data_exp(ii).dth_iter = delta_thscanvals(ii);
    data_exp(ii).dth_disp = dth_disp(ii);
    data_exp(ii).dshift_nominal = dq_shift_nominal;
    data_exp(ii).dqshift_real = dq_shift_real(ii,:);
    data_exp(ii).dqshift = dq_shift_real(ii,:); % initial value of dq
    data_exp(ii).simI = squeeze(simI(ii,:,:));
    data_exp(ii).rock = rock_curve(ii);
    
    % plot:
    if(mod(ii,1)==0) 
        display(['simulating dp, ' num2str(ii) ' of ' num2str(numel(data_exp))]); 
        subplot(121); 
        imagecomp(squeeze(Proj_vol(ii,:,:))); 
        axis image; 
        
        subplot(122); 
        imagesc(squeeze(simI(ii,:,:)));
        axis image;
        
        drawnow;
    end

   
end

middpind = round(numel(data_exp)/2);

