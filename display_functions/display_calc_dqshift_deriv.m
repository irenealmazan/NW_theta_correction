function display_calc_dqshift_deriv(qbragg,dth_nominal,dqshift_deriv_calc,dq_shift0,thBragg,h1)

    % h1 is the handle for the figure
    
    global  ki_o kf_o

    %      if dth_nominal < 0
    %         th_fine_grid = [dth_nominal+dth_nominal/2:abs(dth_nominal)/10:dth_nominal-dth_nominal/2];%
    %      else
    % 
    %         th_fine_grid = [dth_nominal-dth_nominal/2:abs(dth_nominal)/10:dth_nominal+dth_nominal/2];%
    %      end

     th_fine_grid  = [-1:1e-2:1];%[-180:10:180];

     dq_shift = zeros(numel(th_fine_grid),3);

     for jj = 1:numel(th_fine_grid)

        dth(jj) = dth_nominal  +  th_fine_grid(jj);

         Ry = [cosd(dth(jj)) 0 -sind(dth(jj));
            0 1 0;
            sind(dth(jj)) 0 cosd(dth(jj))];

        ki = (Ry * ki_o.').';
        kf = (Ry * kf_o.').';                                             

        dq_shift(jj,:) = kf - ki - qbragg;

        [dq_shift_analytic(jj,:)] = calc_dq_analytical(qbragg,dth(jj),thBragg);
     end

     figure(h1);
     
   subplot(1,2,1);
   hold on;
    plot( dth,dq_shift(:,1),'b');
    hold on;
    plot( dth,dq_shift_analytic(:,1),'.b');
    cstx = dq_shift0(1)-dqshift_deriv_calc(1)*(dth_nominal);%dth_nominal + dth_delta-dq_shift_x_deriv*(dth_nominal + dth_delta);
    plot( dth,cstx +dqshift_deriv_calc(1).*dth,'r');
    plot(dth_nominal,dq_shift0(1),'.r');
    title(['Graph of the derivative of delta q (x component is left and z component is right) with respect to theta at a nominal angle of ' num2str(dth_nominal)])
    legend('dq\_{shift} calculated with the matrix ','dq\_{shift} calculated analytically ','derivative of dq\_{shift} calculated analytically')
    
    subplot(1,2,2);
    hold on;
    plot( dth,dq_shift(:,3),'b');
    hold on;
    plot( dth,dq_shift_analytic(:,3),'.b');
    plot(dth_nominal,dq_shift0(3),'.r');
    cstz = dq_shift0(3)-dqshift_deriv_calc(3)*(dth_nominal);%dth_nominal + dth_delta-dq_shift_x_deriv*(dth_nominal + dth_delta);
    plot( dth,cstz +dqshift_deriv_calc(3).*dth,'r');




end