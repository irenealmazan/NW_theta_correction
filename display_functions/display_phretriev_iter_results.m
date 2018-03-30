% this script runs through data_exp structure and displays the different
% values of the gradient, the second derivative, the angle step at
% different iterations values

addpath(genpath('/Users/ialmazn/Documents/MATLAB/ptycho/m_scripts/'))
addpath(genpath('/Users/ialmazn/Box Sync/Nanowire_ptychography/NSLS II/NSLS II March 2017/Analysis_end_of_beamtime'));



fig1 = figure;


index_to_show = [34];


rows_num = numel(index_to_show);
cols_num = 5; %things to plot
odd_index = [11:10:numel(errlist)]; % at this index, a theta correction is performed

counter = 1;


for ii = index_to_show
    
   
    
    for jj = 1:numel(data_exp(ii).theta_iter)
       grad_final(ii,jj) = data_exp(ii).theta_iter(jj).grad_final;
       alpha(ii,jj) = data_exp(ii).theta_iter(jj).beta;
       theta_step(ii,jj) = data_exp(ii).theta_iter(jj).dth_new_iter;
      
      
    end
    
    
     step_theta_altern = diff(theta_step(ii,:));
     
    figure(fig1);
    subplot(rows_num,cols_num,(counter-1)*cols_num+1);
    plot(squeeze(grad_final(ii,:)));
    xlabel('iterations');
    title(['gradient, theta index ' num2str(ii)]);
    
   
    
    subplot(rows_num,cols_num,(counter-1)*cols_num+2);    
    plot( step_theta_altern,'*r');
    xlabel('iterations');
    title('theta step');
    
    subplot(rows_num,cols_num,(counter-1)*cols_num+3);
    plot(squeeze(theta_step(ii,:)),'-xb');
    hold on;
    plot(ones(numel(data_exp(ii).theta_iter),1).*(data_exp(ii).dth+data_exp(ii).dth_delta),'*r');
    xlabel('iterations');
    title('theta new');
    
    subplot(rows_num,cols_num,(counter-1)*cols_num+4);
    plot(squeeze(alpha(ii,:)));
    xlabel('iterations');
    title('beta parameter (adaptative step coef)');
    
     subplot(rows_num,cols_num,(counter-1)*cols_num+5);
    plot(errlist(odd_index));
    xlabel('iterations');
    title('error');
    
   
    
    counter = counter + 1;
end

%return;

fig2 = figure;
figurepath = '/Users/ialmazn/Box Sync/ptycho_pos_grad/pictures/results_02262018/day/';

for jj = [1 10 20 30 40 50 60 numel(rho_store)]%numel(rho_store)%[1:numel(rho_store)]
    subplot(141)
   imagecomp(rho_store(jj).rho_square);
   axis image;
   title(['iteration ' num2str(jj)]);
   
   subplot(142);
   imagecomp(rho_store(jj).rho_hex);
   axis image;
   
   beta_rho(jj) = rho_store(jj).beta_rho;
   
   subplot(143);
   plot(real(beta_rho(1:jj)),'ob');
   title(['real beta rho ']);
  
   subplot(144);
   plot(imag(beta_rho(1:jj)),'xr');
   title(['imaginary beta rho ']);
  
   figurename = ['rho_iter' num2str(jj)];
   
   print([figurepath figurename],'-dpdf','-bestfit');
   %pause(.5);
end

return ;



