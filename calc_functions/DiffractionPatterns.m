classdef DiffractionPatterns
    % This library contains all the functions which allow us to rotate the
    % detector and the sample
    properties(Constant)
    end
    
    
    methods(Static)
        
        function [dqshift,ki,kf] = calc_dqshift_for_given_th(dth,ki_ini,kf_ini,qbragg)
            % This function calculates the dqshift connecting the diffraction
            % patterns at th = thBragg and at th = thBragg + dth
            
            ki = zeros(numel(dth),3); 
            kf = zeros(numel(dth),3); 
            dqshift = zeros(numel(dth),3);
            
            for ii = 1:numel(dth)
                
                [Ry,Ry_deriv] = RotationMatrix.rock_curve(dth(ii));
                
                ki(ii,:) = (Ry * ki_ini.').';
                kf(ii,:) = (Ry * kf_ini.').';
                
                dqshift(ii,:) = (kf(ii,:)-ki(ii,:))-qbragg;
            end
            
        end
        
        function [dq_shift_deriv] = calc_dq_deriv_analytical(dth,ki,kf,qbragg)
            % This function calculates the first derivative of the vector dq
            % (linking the different positions of the detector at different theta
            % angles - differents points in the rocking curve)
            
            dq_shift_deriv = zeros(numel(dth),3);
            
            for jj = 1:numel(dth)
                [~,Ry_rock_deriv] = RotationMatrix.rock_curve(dth(jj));
                dq_shift_deriv(jj,:) = (Ry_rock_deriv * kf' - Ry_rock_deriv * ki');
                %dq_shift_deriv(jj,:) = (Ry_rock_deriv * qbragg');
            end
            
         end
        
       function [dq_shift_deriv] = calc_dq_deriv_manually(dth,ki_ini,kf_ini,qbragg)
            % This function calculates the first derivative of the vector dq
            % (linking the different positions of the detector at different theta
            % angles - differents points in the rocking curve)
            
            dq_shift_deriv = zeros(numel(dth),3);
            
            dth_infinitesimal = 1e-4;
            
            for jj = 1:numel(dth)
                dq_shift_plus = DiffractionPatterns.calc_dqshift_for_given_th(dth(jj)+ dth_infinitesimal,ki_ini,kf_ini,qbragg);
                dq_shift_minus = DiffractionPatterns.calc_dqshift_for_given_th(dth(jj)- dth_infinitesimal,ki_ini,kf_ini,qbragg);

                dq_shift_deriv(jj,:) = (dq_shift_plus - dq_shift_minus)/(2*dth_infinitesimal);
            end
            
        end
         
        function [Psi_mod,rock_curve,Psij,FT_Psij,Qterm] = calc_dp(dq_shift,probe,rho,X,Y,Z)
            % This function calculates a diffraction pattern corresponding
            % to one single point of the rocking curve and one single dq
            % shift
            %
            Psij = zeros(size(dq_shift,1),size(rho,1),size(rho,2));
            FT_Psij = zeros(size(dq_shift,1),size(rho,1),size(rho,2));
            Psi_mod = zeros(size(dq_shift,1),size(rho,1),size(rho,2));
            rock_curve = zeros(size(dq_shift,1),1);
            
            for jj = 1:size(dq_shift,1)
                % calculate the Q term which depends on the theta angle distorted!!! :
                Qterm(jj).Qterm = exp(1i* dq_shift(jj,1) * X) .* ...
                    exp(1i* dq_shift(jj,2) * Y) .* ...
                    exp(1i* dq_shift(jj,3) * Z);
                
                % R(Q*P_j*rho): projected volume
                Psij(jj,:,:) = sum( rho.*probe.*Qterm(jj).Qterm, 3);
                
                % 2D_FT: diffracted wave
                FT_Psij(jj,:,:) = fftshift(fftn(fftshift( squeeze(Psij(jj,:,:)))));
                
                % modulus square of the diffracted wave
                Psi_mod(jj,:,:) =  squeeze(FT_Psij(jj,:,:)).* conj( squeeze(FT_Psij(jj,:,:)));
                
                % integrated intensity:
                rock_curve(jj) = sum(sum(squeeze(Psi_mod(jj,:,:)),1));
            end
            
        end
                
        function [errtot] = calc_error_multiangle(probe, rho, data,angle_list,ki,kf,X,Y,Z)
            % This function calculates the diffracted intensities
            % differences for the set of data in the the structure "data"
            
            qbragg = kf - ki;
            
            errtot=0;
            
            
            [dq_shift] = DiffractionPatterns.calc_dqshift_for_given_th(angle_list,ki,kf,qbragg);
            
            [Psij_mod] = DiffractionPatterns.calc_dp( dq_shift,probe,rho,X,Y,Z);
            
            for ii=1:numel(data)
                err = sqrt(data(ii).I) - sqrt(squeeze(Psij_mod(ii,:,:)));
                err = sum(err(:).^2)/numel(err);
                errtot = errtot + err;
            end
            
            
            
        end
        
        
    end
end
    

