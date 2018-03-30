function display_calc_dqshift_deriv2(qbragg,dth_nominal,dqshift_deriv2_calc,dq_shift_deriv0,thBragg)

    global ki_o kf_o

    %      if dth_nominal < 0
%         th_fine_grid = [dth_nominal+dth_nominal/2:abs(dth_nominal)/10:dth_nominal-dth_nominal/2];%
%      else
% 
%         th_fine_grid = [dth_nominal-dth_nominal/2:abs(dth_nominal)/10:dth_nominal+dth_nominal/2];%
%      end

      th_fine_grid  = [-0.1:1e-2:0.1];%[-180:10:180];
     
     dq_shift_deriv = zeros(numel(th_fine_grid),3);
     
     h1 = figure;
     hold on;
     for jj = 1:numel(th_fine_grid)

        dth_grid(jj) = dth_nominal  +  th_fine_grid(jj);

        [dq_shift_deriv(jj,:)] = calc_dq_deriv_analytical(qbragg,dth_grid(jj),thBragg);
        
        %%% in the case that do you want to check dq_shift_deriv is well calculated
        %%{
        Ry = [cosd(-dth_grid(jj)) 0 sind(-dth_grid(jj));
            0 1 0;
            -sind(-dth_grid(jj)) 0 cosd(-dth_grid(jj))];
        
        ki = (Ry * ki_o.').';
        kf = (Ry * kf_o.').';
        
        dq_shift = kf - ki - qbragg;
        
        
        display_calc_dqshift_deriv(qbragg,dth_grid(jj),dq_shift_deriv(jj,:),dq_shift,thBragg,h1)
        %}
     end

    figure;
    
    subplot(1,2,1);
    plot(dth_grid,dq_shift_deriv(:,1));
    hold on;
    cstx = dq_shift_deriv0(1)-dqshift_deriv2_calc(1)*(dth_nominal);%dth_nominal + dth_delta-dq_shift_x_deriv*(dth_nominal + dth_delta);
    plot(dth_grid,cstx +dqshift_deriv2_calc(1).*dth_grid);
   title(['Graph of the second derivative of delta q (x component is left and z component is right) with respect to theta at a nominal angle of ' num2str(dth_nominal)])

    
    
    subplot(1,2,2);
    plot(dth_grid,dq_shift_deriv(:,3));
    hold on;
    cstz = dq_shift_deriv0(3)-dqshift_deriv2_calc(3)*(dth_nominal);%dth_nominal + dth_delta-dq_shift_x_deriv*(dth_nominal + dth_delta);
    plot(dth_grid,cstz +dqshift_deriv2_calc(3).*dth_grid);




end