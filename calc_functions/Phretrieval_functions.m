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
        
        function [scale_fact,errtot] = ini_guess_scalefactor(probe_BCDI, rho_ini,angle_list_ini, data,scale_fact_guess,ki,kf,X,Y,Z)
            % This function estimates the scaling factor by which the initial guess
            % has to be calculated in order to establish good initial conditions
            % for phase retrieval.
            
            errtot = zeros(numel(scale_fact_guess),1);
                      
            for jj=1:numel(scale_fact_guess)
                rho_ini_scale = rho_ini.*scale_fact_guess(jj);                
                [errtot(jj)] = DiffractionPatterns.calc_error_multiangle(probe_BCDI, rho_ini_scale, data,angle_list_ini,ki,kf,X,Y,Z);                
            end
            
            scale_fact = scale_fact_guess(find(min(errtot)==errtot));
            
        end
                
        function [rho_new,beta_rho] = rho_update(probe, rho,angles_list,support, data_exp,depth,err_0,ki,kf,X,Y,Z)
            % this functions updates rho
            
            % gradient calculation
            [gPIEiter] = GeneralGradient.calc_grad_multiangle(probe, rho,angles_list,data_exp,ki,kf,X,Y,Z);
            
            % direction of the scaled steepest descent:
            D = 1/(max(max(max(abs(probe).^2))));
            
            %       direction_rho = -conj(gPIEiter);
            direction_rho = - (D/depth)*(gPIEiter.*support);
            
            % calculate the adaptative step length
            [beta_rho] = GeneralGradient.calc_beta_adaptative_step(probe, rho,angles_list,data_exp,gPIEiter,err_0,direction_rho,'rho',ki,kf,X,Y,Z);
            
            % update the object:
            rho_new = rho + beta_rho * direction_rho;
                                 
        end
        
        function [dth_new,dq_shift,grad_final_theta,beta] = theta_update(probe, rho,angles_list,data_exp,Niter_theta,dthBragg,error_0,ki,kf,X,Y,Z)
            %%% this function calculates the gradient of the error metric with
            %%% respect to the position of the angles analytically, and
            %%% the correction theta  step
                        
            qbragg = kf - ki;
                     
            % initialize varibles:
            dth_new = angles_list;
            
            dq_shift = zeros(numel(data_exp),3);
            
            grad_final_theta = zeros(Niter_theta,numel(data_exp));
            
            for ntheta = 1:Niter_theta
                                
                [grad_final_theta(ntheta,:)] = GeneralGradient.calc_grad_theta(probe, rho, data_exp, dth_new,ki,kf,X,Y,Z);
                
                for ii = 1:numel(data_exp)%
                    display(['dth_new = ' num2str(dth_new(ii)) 'gradient = ' num2str(grad_final_theta(ntheta,ii))]);
                end
                
                % define the direction of the descent:
                direction = -squeeze(grad_final_theta(ntheta,:)');
                
                % calculate an adaptative step size:
                [beta] = GeneralGradient.calc_beta_adaptative_step(probe, rho, dth_new,data_exp,grad_final_theta,error_0,direction,'theta',ki,kf,X,Y,Z);

                % corrected theta :
                dth_new = dth_new + beta* direction;
                
                % corrected dqshift:
                [dq_shift] = DiffractionPatterns.calc_dqshift_for_given_th(dth_new,ki,kf,qbragg);
                
                
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
    

