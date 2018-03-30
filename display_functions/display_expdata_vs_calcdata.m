% this scripts compares for a  given angle of the rocking curve, the
% experimental and the calculated diffraction pattern


theta_indx = 34;

data = data_exp(theta_indx);

% calculated diffraction pattern

 Qterm_display = exp(i* data.dqshift(1) * X) .* ...
            exp(i* data.dqshift(2) * Y) .* ...
            exp(i* data.dqshift(3) * Z);
        
   

temp = sum( rho_guess.*probe.*Qterm_display, 3);
Psij_display = fftshift(fftn(fftshift(temp)));

Psij_mod_display = Psij_display.*conj(Psij_display);


figure;

subplot(1,2,1);
imagesc(data.I);
axis image;
title(['experiental at ' num2str(data.dth) ' and index ' num2str(theta_indx) ]);

subplot(1,2,2);
imagesc(Psij_mod_display);
axis image;
title([' Calculated ' ]);


