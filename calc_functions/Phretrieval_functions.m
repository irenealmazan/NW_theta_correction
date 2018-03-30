classdef Phretrieval_functions
    % This library contains all the functions which allow us to rotate the
    % detector and the sample
    properties(Constant)
    end
    
    
    methods(Static)
        
        function [support] = make_support(corners,facet_spacing,edgepad)
            % This function calculates the dqshift connecting the diffraction
            % patterns at th = thBragg and at th = thBragg + dth
            
            
            support_edge = facet_spacing * edgepad;
            
            %make parallel planes parallel to each pair of facets
            v1 = [corners(2,:) - corners(8,:)]; v1 = v1/norm(v1);
            v2 = [corners(1,:) - corners(8,:)]; v2 = v2/norm(v2);
            v3 = cross(v1,v2); v3 = v3/norm(v3);
            T1 = v3(1)*X + v3(2)*Y + v3(3)*Z;
            T1 = (T1>-support_edge/2 & T1<support_edge/2); %two parallel lines
            
            v1 = [corners(2,:) - corners(8,:)]; v1 = v1/norm(v1);
            v2 = [corners(3,:) - corners(8,:)]; v2 = v2/norm(v2);
            v3 = cross(v1,v2); v3 = v3/norm(v3);
            T2 = v3(1)*X + v3(2)*Y + v3(3)*Z;
            T2 = (T2>-support_edge/2 & T2<support_edge/2); %two parallel lines
            
            v1 = [corners(4,:) - corners(9,:)]; v1 = v1/norm(v1);
            v2 = [corners(3,:) - corners(9,:)]; v2 = v2/norm(v2);
            v3 = cross(v1,v2); v3 = v3/norm(v3);
            T3 = v3(1)*X + v3(2)*Y + v3(3)*Z;
            T3 = (T3>-support_edge/2 & T3<support_edge/2); %two parallel lines
            
            v3 = [0 1 0];
            T4 = v3(1)*X + v3(2)*Y + v3(3)*Z;
            T4 = (T4>-support_edge/2 & T4<support_edge/2); %two parallel lines
            
            
            support = T1&T2&T3&T4;
            
            %support = T3;
            %support = abs(NW);
            
            
        end
        
        function [scale_fact,errtot] = ini_guess_scalefactor(probe_BCDI, rho_ini, data,scale_fact_guess)
            % This function estimates the scaling factor by which the initial guess
            % has to be calculated in order to establish good initial conditions
            % for phase retrieval.
            
            errtot = zeros(numel(scale_fact_guess),1);
            
            for jj=1:numel(scale_fact_guess)
                
                rho_ini_scale = rho_ini.*scale_fact_guess(jj);
                
                
                [errtot(jj)] = calc_error_multiangle_Irene(probe_BCDI, rho_ini_scale, data);
                
            end
            
            scale_fact = scale_fact_guess(find(min(errtot)==errtot));
            
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
            
            % rho at initial alpha
            rho_alpha = rho + alpha_ini * direction;
            
            % estimate the value of the error at rho + alpha_ini*gradtot_rho
            err_plusalpha(counter) =  DiffractionPaterns.calc_error_multiangle(probe, rho_alpha, data,ki,kf,X,Y,Z);
            
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
                            data(ii).dth_iter = data(ii).dth_new +alpha_iter*gradtot(ii);
                        end
                end
                
                
                % estimate the value of the error at rho + alpha_ini*gradtot_rho
                err_plusalpha(counter) =  DiffractionPaterns.calc_error_multiangle(probe, rho_alpha, data,ki,kf,X,Y,Z);;
                
                
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
        
        function [rho_new] = rho_update(probe, rho, data_exp,depth)
            % this functions updates rho
            
            % gradient calculation
            [gPIEiter] = GeneralGradient.calc_grad_multiangle(probe, rho, data_exp,X,Y,Z);
            
            % direction of the scaled steepest descent:
            D = 1/(max(max(max(abs(probe).^2))));
            
            %       direction_rho = -conj(gPIEiter);
            direction_rho = - (D/depth)*gPIEiter;
            
            % calculate the adaptative step length
            [beta_rho] = Phretrieval_functions.calc_beta_adaptative_rhostep(probe, rho, data_exp,gPIEiter,errlist(end),direction_rho,'rho',ki,kf,X,Y,Z);
            
            % update the object:
            rho = rho + beta_rho * direction_rho;
            
            % multiply by the support:
            rho_new = rho .*support;
            
            
            
        end
        
        function [dth_new,dq_shift,grad_final_theta,beta] = theta_correction(probe, rho,data_exp,Niter_theta,index_to_distort,dthBragg,error_0)
            %%% this function calculates the gradient of the error metric with
            %%% respect to the position of the angles analytically, and
            %%% the correction theta  step
            
            global ki_o kf_o
            
            qbragg = kf_o - ki_o;
                     
            % initialize varibles:
            
            dth_new = zeros(numel(data_exp),1);
            for ii=1:numel(data_exp)
                dth_new(ii) = data_exp(ii).dth_iter;
            end
            
            dq_shift = zeros(numel(data_exp),3);
            
            grad_final_theta = zeros(numel(data_exp),Niter_theta);
            
            for ntheta = 1:Niter_theta
                                
                orderrandom = randperm(numel(data_exp));%randperm(numel(index_to_distort));
                
                for ii = orderrandom%index_to_distort(orderrandom)%index_to_distort%1:numel(data_exp)%

                    [grad_final_theta(ii,ntheta)] = GeneralGradient.calc_grad_theta(probe, rho, data_exp(ii), dth_new(ii),0,dthBragg);
                    
                    display(['dth_new = ' num2str(dth_new(ii)) 'gradient = ' num2str(grad_final_theta(ii,ntheta))]);
                    
                end
                
                % define the direction of the descent:
                direction = -grad_final_theta;
                
                % calculate an adaptative step size:
                [beta] = Phretrieval_functions.calc_beta_adaptative_thetastep(bmtemp, rho, data_exp,grad_final_theta,error_0,direction,'theta');
                
                % corrected theta :
                dth_new = dth_new + beta* direction;
                
                % corrected dqshift:               
                [dq_shift] = calc_dqshift_for_given_th(dth_new,ki_o,kf_o,qbragg);                                 
                
                
            end
            
            
            
            
            %%% TEST OF  GRADIENT
            %{
            for jj= dthsearchind(2:end-1)
                grad_manual = test_grad_theta_manually(jj,thscan,fitQerr,data_exp(ii),bmtemp,rho,grad_final_theta(ii).grad(jj));

            end
        
            %}
            
            
        end
        
    end
        
   
        
        
end
    

