function show_rho_theta_update(figure_num,errlist,rho,midsl,anglelist,trueangles,flagINI)
% this function creates the figure where we plot the improvements of the
% rho and the angles during the phase retrieval iterations. flagINI
% indicates if the figure needs to be created (initial state) or only
% updated

    figure(figure_num);
    if flagINI == 1
        clf;
        setfigsize(gcf, 1000,500);
    end
    % plot
    subplot(141); imagecomp(rho(:,:,midsl)); colorbar; axis image; %zoom(1.5);
    subplot(142); imagecomp(squeeze(rho(100,:,:))); colorbar; axis image; %zoom(1.5);
    subplot(143); plot(log10(errlist),'LineWidth',3.0);
    subplot(144); plot(anglelist,'ob');
    hold on; plot(trueangles,'*r');

    drawnow;

end