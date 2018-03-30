classdef GeneralGradient
    % This library contains all the functions which allow us to rotate the
    % detector and the sample
    properties(Constant)
    end
    
    
    methods(Static)
        
        
        
        function [gradtot] = calc_grad_multiangle(probe, rho, data,X,Y,Z)
            
            gradtot = zeros(size(rho));
            
            for ii=1:numel(data)
                
                
                [~,~,Psij,~,Qterm] = DiffractionPatterns.calc_single_dp(data(ii).dqshift,probe,rho,X,Y,Z);
                
                Psig = flipud(sqrt(data(ii).I)).*exp(i*angle(Psij));%
                grad = Psij - Psig; %%%% NEW % Psij - Psig;
                grad = fftshift(ifftn(fftshift(grad)));
                grad = repmat(grad, [1 1 size(probe,3)]);
                grad = grad .* conj(Qterm);
                grad = grad .* conj(probe);
                grad = 2*grad; %add factor of 2 to make sure that gradient is equivalent
                
                gradtot = gradtot + grad;
            end
        end
    
        function [grad_final_theta,error_theta_integ] = calc_grad_theta(probe, rho, data, dth_nominal, dth_delta,thBragg,X,Y,Z,ki,kf)
            % this function calculates the gradient of the error metric with
            % respect to theta. dth_delta = 0 if we are using the function to
            % refine the angular position.    
            
            qbragg = kf - ki;
            
            %%% calculate the derivative of dq_shift with respect to theta:
            [dq_shift_deriv] = GeneralGradient.calc_dq_deriv_analytical(qbragg,dth_nominal,thBragg);
          
            % calculate the current momentum transfer:
            [dq_shift] = DiffractionPatterns.calc_dqshift_for_given_th(-(dth_nominal + dth_delta),ki,kf,qbragg);
           
            %h2 = figure;
            %display_calc_dqshift_deriv(qbragg,dth_nominal,dq_shift_deriv,dq_shift,thBragg,h2)
      
            [Psij_mod,~,~,FT_Psij,Qterm] = DiffractionPatterns.calc_single_dp(dq_shift,probe,rho,X,Y,Z);
                      
            %%% derivative of Qterm with respect to dq_shift
            deriv_Qterm_theta = (dq_shift_deriv(1).*X + dq_shift_deriv(3).*Z).*Qterm;                     
            
            % derivative of the diffraction pattern with respect to dq_shift
            Psij_deriv_theta = probe.*rho.*deriv_Qterm_theta;
            
            Psij_deriv_theta = sum(Psij_deriv_theta,3);
            
            Psij_deriv_theta = fftshift(fftn(fftshift(Psij_deriv_theta)));
                        
            Im_Psij_deriv_theta = imag(conj(FT_Psij).*Psij_deriv_theta);
            
            % calculate the gradient:
          
            grad_theta = (1 -sqrt(flipud(data.I)./(Psij_mod))).* Im_Psij_deriv_theta; %%%% NEW
            
            grad_final_theta = -2*sum(sum(grad_theta))/(numel(grad_theta));
                       
            %%% check that you calculate the gradient:
            error_theta_integ = DiffractionPatterns.calc_error_multiangle(probe, rho, data,data.dth_iter,ki,kf,X,Y,Z);
            %display_calc_grad_theta(probe, rho, data, dth_nominal,grad_final_theta,error_theta_integ)
            %display_calc_grad2_theta(probe, rho, data, dth_nominal,grad2_final_theta,grad_final_theta)
            
        end
   
        function [dq_shift_deriv] = calc_dq_deriv_analytical(qbragg,dth_nominal,thBragg)
            % This function calculates the first derivative of the vecotr dq
            % (linking the different positions of the detector at different theta
            % angles (differents points in the rocking curve).
            
            
            qbragg_mod = sqrt(qbragg*qbragg');

            dq_shift_mod_analytic = 2*qbragg_mod*sind(abs(dth_nominal/2));%- 2*sind(dth(jj)/2);%
            
            if dth_nominal< 0
                
                dq_shift_mod_deriv = -qbragg_mod*(pi/180)*cosd(abs(dth_nominal/2));
               
                % see summary slide 50:
                dq_shift_x_analytic_unit = -sind(thBragg+dth_nominal/2);%sind(thBragg-abs(dth(jj))/2);I don't understand why do I need a minus sign there!
                dq_shift_z_analytic_unit = cosd(thBragg+dth_nominal/2);%cosd(thBragg-abs(dth(jj))/2);%2*cosd(thBragg-dth(jj)/2)*sind(-dth(jj)/2)*qbragg_mod;%
                
                dq_shift_x_analytic_unit_deriv = -(pi/180)*0.5*cosd(thBragg+dth_nominal/2);
                dq_shift_z_analytic_unit_deriv = -(pi/180)*0.5*sind(thBragg+dth_nominal/2);
            else
                
                
                dq_shift_mod_deriv = qbragg_mod*(pi/180)*cosd(abs(dth_nominal/2));
                
                
                % see summary slide 50:
                
                dq_shift_x_analytic_unit = -sind(abs(thBragg)-abs(dth_nominal)/2);
                dq_shift_z_analytic_unit = -cosd(abs(thBragg)-abs(dth_nominal)/2);%2*cosd(thBragg-dth(jj)/2)*sind(-dth(jj)/2)*qbragg_mod;%
                
                dq_shift_x_analytic_unit_deriv = +(pi/180)*0.5*cosd(abs(thBragg)-abs(dth_nominal)/2);
                dq_shift_z_analytic_unit_deriv = -(pi/180)*0.5*sind(abs(thBragg)-abs(dth_nominal)/2);
                
            end
            
            
            dq_shift_x_deriv = dq_shift_mod_deriv*dq_shift_x_analytic_unit+dq_shift_mod_analytic*dq_shift_x_analytic_unit_deriv;
            
            dq_shift_z_deriv =  dq_shift_mod_deriv*dq_shift_z_analytic_unit+dq_shift_mod_analytic*dq_shift_z_analytic_unit_deriv;
            
            dq_shift_deriv = [dq_shift_x_deriv 0 dq_shift_z_deriv];
            
        end
        
        function [alpha_iter,err_plusalpha,alpha_track] = calc_beta_adaptative_step(probe, rho, data,gradtot,err_0,direction,flag,ki,kf,X,Y,Z)
            % this function calculates the adaptative step length for the
            % correction of the rho. We folllow Nocedal "Backtracking Line Search
            % The imputs are the following:
            %   probe: information about the beam (important in ptychography)
            %   rho: the object in its current state (before the update)
            %   data: structure containing the angles
            %   gradtot: vector with all the gradients for each angle
            
            
            
            % Armijo's condition: if err(alpha_i) > err(alpha=0) +
            % c1*deriv_err_with_theta *alpha_i
            
            % calculate the value of alpha_ini that we choose as the value for
            % which the linear approximation of the error metric becomes positive
            c1 = 1e-3; % see Nocedal
            counter = 1; % counter to track the evolution of the error, alpha and the approximation of the error with the iterations
            counter_max = 5;
            tau_backtrack = 0.01;
            err_linear_aprox(counter) = - err_0;
            
            slope_alpha_0 = c1*real(direction(:)'*gradtot(:));
            alpha_ini = 1;%-err_0/(slope_alpha_0);
            
            % rho/theta at initial alapha
            theta_alpha = zeros(numel(data),1);
            switch flag
                
                case 'rho'
                    rho_alpha = rho + alpha_ini * direction;
                    for ii=1:numel(data)
                        theta_alpha(ii) = data(ii).dth_iter;
                    end
                case 'theta'
                    rho_alpha = rho;
                    for ii=1:numel(data)
                        theta_alpha(ii) = data(ii).dth_iter +alpha_ini * gradtot(ii);
                    end
            end
            
            % estimate the value of the error at rho + alpha_ini*gradtot_rho
            err_plusalpha(counter) =  DiffractionPatterns.calc_error_multiangle(probe, rho_alpha, data,theta_alpha,ki,kf,X,Y,Z);
            
            Delta_err(counter) = err_plusalpha(counter) - err_0;
            
            alpha_iter = alpha_ini;
            
            alpha_track(counter) = alpha_iter;
            
            while( Delta_err(counter) >= err_linear_aprox(counter))
                
                % update the counter
                counter = counter + 1;
                
                % update alpha
                alpha_iter = alpha_iter*tau_backtrack;
                
                % move rho
                switch flag
                    
                    case 'rho'
                        rho_alpha = rho + alpha_iter * direction;
                        
                    case 'theta'
                        for ii=1:numel(data)
                            theta_alpha(ii) = data(ii).dth_iter +alpha_iter*gradtot(ii);
                        end
                end
                
                
                % estimate the value of the error at rho + alpha_ini*gradtot_rho
                err_plusalpha(counter) =  DiffractionPatterns.calc_error_multiangle(probe, rho_alpha, data,theta_alpha,ki,kf,X,Y,Z);
                
                
                %recalculate the difference between error metric at alpha_iter and
                %error metric at alpha = 0
                Delta_err(counter) = err_plusalpha(counter) - err_0;
                
                % calculate the linear approximation of the error metric
                err_linear_aprox(counter) =  alpha_iter*slope_alpha_0;
                
                display(['err(' num2str(alpha_iter) ' ) - err(0) = ' num2str(Delta_err(counter)) ' and linear approx. ' num2str(err_linear_aprox(counter))])
                
                alpha_track(counter) = alpha_iter;
                
                if counter > counter_max
                    display('beta not found');
                    break;
                    
                end
                
            end
            
            
            
            
        end
    end
        
   
        
        
end