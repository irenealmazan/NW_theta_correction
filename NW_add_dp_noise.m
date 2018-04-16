% In this script we add Poison noise to the simulated data and we set the
% signal to background ratio with mncntrate; we also overwrite the
% experimental data with the simulated data.

mn = mean(data_exp(middpind).simI(:));

if (noiseflag)
    for ii=1:numel(data_exp)
        
        if(mod(ii,50)==0) display(['adding noise to sim dp ' num2str(ii) ' of ' num2str(numel(data_exp))]); end
        data_exp(ii).noiseI = poisrnd( data_exp(ii).simI /mn * mncntrate);
        %data_exp(ii).simI = poisrnd( data_exp(ii).simI);
        %data_exp(ii).simI = data_exp(ii).simI;
    end
end
%%

%%
%if use simulated data, overwrite the experimental I with the noisy simI
%and calculates the noisy rocking curve

rock_curve_noise = zeros(numel(data_exp),1);

if(usesimI)
    display(['overwriting experimental I with sim I']);
    for ii=1:numel(data_exp)
        if (noiseflag)
            data_exp(ii).I = data_exp(ii).noiseI;            
        else
            data_exp(ii).I = data_exp(ii).simI;            
        end
        rock_curve_noise(ii) = sum(sum(data_exp(ii).I));
        
    end
end

figure(7); 
clf; setfigsize(gcf, 800,400);
colormap jetvar;
hax1 = axes('position', [0 0 .5 1]);
hax2 = axes('position', [0.5 0 .5 1]);
ca = [0 2];
for ii=1:numel(data_exp)
    %imagesc(hax1, log10(data_exp(ii).I));
    imagesc(hax1, (data_exp(ii).I));
    %caxis(hax1, ca); 
    set(hax1, 'xtick', Npix/2+1, 'ytick', Npix/2+1); grid(hax1, 'on');
    drawnow;
    %imagesc(hax2, log10(data_exp(ii).simI)); axis image;
    imagesc(hax2, (data_exp(ii).simI)); axis image;
    %caxis(hax2, ca); 
    set(hax2, 'xtick', Npix/2+1, 'ytick', Npix/2+1); grid(hax2, 'on');
    %caxis(ca);
    display(['dth: ' num2str(data_exp(ii).dth_nominal)]);
    pause(.1);
end
