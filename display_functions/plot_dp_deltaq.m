% this script allows to plot the diffraction pattern at slightly different
% theta values

Qterm_center =  exp(1i* dq_shift(1) * X) .* ...
        exp(1i* dq_shift(2) * Y) .* ...
        exp(1i* dq_shift(3) * Z);
    
err_center = calc_error_theta_singlepos(rho,probe,data,Qterm_center);
  
    
Qterm_grid_1 = exp(1i* dq_shift_grid(1,1) * X) .* ...
exp(1i* dq_shift_grid(1,2) * Y) .* ...
exp(1i* dq_shift_grid(1,3) * Z);


Psij_1 = probe.*rho.*Qterm_grid_1;
Psij_1 = sum(Psij_1,3);                     %Radon oper
Psij_1 = fftshift(fftn(fftshift(Psij_1)));
Psij_1_conj = conj(Psij_1);

Psij_1_mod = Psij_1.*Psij_1_conj;

Qterm_grid_2 = exp(1i* dq_shift_grid(2,1) * X) .* ...
exp(1i* dq_shift_grid(2,2) * Y) .* ...
exp(1i* dq_shift_grid(2,3) * Z);

Psij_2 = probe.*rho.*Qterm_grid_2;
Psij_2 = sum(Psij_2,3);                     %Radon oper
Psij_2 = fftshift(fftn(fftshift(Psij_2)));
Psij_2_conj = conj(Psij_2);

Psij_2_mod = Psij_2.*Psij_2_conj;

figure;
subplot(121);
imagesc(log(Psij_1_mod));
axis image;

subplot(122);
imagesc(log(Psij_2_mod));
axis image;

% calcualte the error:

 th_fine_grid = [-0.05:1e-3 :0.05];
 err_grid = zeros(size(th_fine_grid,1),1);  
 
   for jj = 1:numel(th_fine_grid)
        dth_grid(jj) = dth_nominal + dth_delta +  th_fine_grid(jj);
        
        Ry = [cosd(-dth_grid(jj)) 0 sind(-dth_grid(jj));
            0 1 0;
            -sind(-dth_grid(jj)) 0 cosd(-dth_grid(jj))];
        
        ki = (Ry * ki_o.').';
        kf = (Ry * kf_o.').';
                
        dq_shift_grid(jj,:) = kf - ki - qbragg;
        
        Qterm_grid_1 = exp(1i* dq_shift_grid(jj,1) * X) .* ...
            exp(1i* dq_shift_grid(jj,2) * Y) .* ...
            exp(1i* dq_shift_grid(jj,3) * Z);
        
        err_grid(jj) = calc_error_theta_singlepos(rho,probe,data,Qterm_grid_1);
   end

   figure;
   subplot(121);plot(dth_grid,err_grid);
   hold on;      
   cstx = err_center-grad_final_theta*(dth_nominal + dth_delta );
   plot(dth_grid,cstx + grad_final_theta* dth_grid);
   plot(dth_grid(find(th_fine_grid==0)),err_center,'or');
   plot(dth_grid(find(th_fine_grid==0)),error_theta_integ,'xg');

 
  



