
function [error_theta_integ_up,error_theta_integ_down] = show_Qterm_effect(data,probe,rho,dth_up,dth_down)

global  ki_o kf_o X Y Z

    % theta + delta_theta
%     
    
%     
     Ry = [cosd(-dth_up) 0 sind(-dth_up);
        0 1 0;
        -sind(-dth_up) 0 cosd(-dth_up)];
    
    ki = (Ry * ki_o.').';
    kf = (Ry * kf_o.').';

    qbragg = kf_o - ki_o;
    
    dqtest_up = kf - ki - qbragg;


     Qterm_up = exp(1i* dqtest_up(1) * X) .* ...
            exp(1i* dqtest_up(2) * Y) .* ...
            exp(1i* dqtest_up(3) * Z);
        
     error_theta_integ_up = calc_error_theta_singlepos(rho,probe,data,Qterm_up);

% Qterm_up=  exp(i* dq_shift(1) * X) .* ...
%             exp(i* dq_shift(2) * Y) .* ...
%             exp(i* dq_shift(3) * Z);

% Qterm_up=  exp(i* dqshift(1) * X) .* ...
%             exp(i* dqshift(2) * Y) .* ...
%             exp(i* dqshift(3) * Z);
        
        
    %Psij = probe.*rho.*Qterm_up;
%     Psij = probe.*rho.*Qterm_up;
%     Psij = sum(Psij,3);                     %Radon oper
%     Psij = fftshift(fftn(fftshift(Psij)));
%     Psij_conj = conj(Psij);
%     
%     Psij_mod_up = Psij.*Psij_conj;
%     mn_up = mean(Psij_mod_up(:));
%     
    %%%%%%%%
   
    
     Ry = [cosd(-dth_down) 0 sind(-dth_down);
        0 1 0;
        -sind(-dth_down) 0 cosd(-dth_down)];
    
    ki = (Ry * ki_o.').';
    kf = (Ry * kf_o.').';

    qbragg = kf_o - ki_o;
    
    dqtest_down = kf - ki - qbragg;


     Qterm_down = exp(1i* dqtest_down(1) * X) .* ...
            exp(1i* dqtest_down(2) * Y) .* ...
            exp(1i* dqtest_down(3) * Z);
        
    error_theta_integ_down = calc_error_theta_singlepos(rho,probe,data,Qterm_down);
%         

% Qterm_down=  exp(i* data.dqshift_delta(1) * X) .* ...
%             exp(i* data.dqshift_delta(2) * Y) .* ...
%             exp(i* data.dqshift_delta(3) * Z);

% Qterm_down=  exp(i* data.dqshift_delta(1) * X) .* ...
%             exp(i* data.dqshift_delta(2) * Y) .* ...
%             exp(i* data.dqshift_delta(3) * Z);
%         
%     Psij = probe.*rho.*Qterm_down;
%     Psij = sum(Psij,3);                     %Radon oper
%     Psij = fftshift(fftn(fftshift(Psij)));
%     Psij_conj = conj(Psij);
%     
%     Psij_mod_down = Psij.*Psij_conj;
%     mn_down = mean(Psij_mod_down(:));
%     
%     figure;
%     subplot(131);
%     imagesc((Psij_mod_down./mn_down));
%     axis image;
%     title('true dth');
% 
%     subplot(132);
%     imagesc((Psij_mod_up./mn_up));
%     axis image;
%     title('wrong dth');
%     
%     subplot(133);
%     imagesc(((Psij_mod_down./mn_down) - (Psij_mod_up./mn_up)));
%     axis image;
%     title('difference dth');
    
end
