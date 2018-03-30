classdef DiffractionPatterns
    % This library contains all the functions which allow us to rotate the
    % detector and the sample
    properties(Constant)
    end
    
    
    methods(Static)
        
        function [dqshift,ki,kf] = calc_dqshift_for_given_th(dth,ki_ini,kf_ini,qbragg)
            % This function calculates the dqshift connecting the diffraction
            % patterns at th = thBragg and at th = thBragg + dth
            
            for ii = 1:numel(dth)
                
                Ry = [cosd(-dth(ii)) 0 sind(-dth(ii));
                    0 1 0;
                    -sind(-dth(ii)) 0 cosd(-dth(ii))];
                
                ki(ii,:) = (Ry * ki_ini.').';
                kf(ii,:) = (Ry * kf_ini.').';
                
                dqshift(ii,:) = (kf(ii,:)-ki(ii,:))-qbragg;
            end
            
        end
        
        
        function [Psi_mod,rock_curve,Psij,FT_Psij,Qterm] = calc_single_dp(dq_shift,probe,rho,X,Y,Z)
            % This function calculates a diffraction pattern corresponding
            % to one single point of the rocking curve and one single dq
            % shift
            %
            
            % calculate the Q term which depends on the theta angle distorted!!! :
            Qterm = exp(1i* dq_shift(1) * X) .* ...
                exp(1i* dq_shift(2) * Y) .* ...
                exp(1i* dq_shift(3) * Z);
            
            % R(Q*P_j*rho): projected volume
            Psij = sum( rho.*probe.*Qterm, 3);
            
            % 2D_FT: diffracted wave
            FT_Psij = fftshift(fftn(fftshift(Psij)));            
           
            % modulus square of the diffracted wave
            Psi_mod =  FT_Psij.* conj( FT_Psij);
            
            % integrated intensity:
            rock_curve = sum(sum(Psi_mod,1));
            
        end
        
        
        
        function [errtot] = calc_error_multiangle(probe, rho, data,ki,kf,X,Y,Z)
           % This function calculates the diffracted intensities
           % differences for the set of data in the the structure "data"
           
            qbragg = kf - ki;
           
            errtot=0;
            
            dq_shift = zeros(numel(data),3);
            
            for ii=1:numel(data)
                
                [dq_shift(ii,:)] = DiffractionPaterns.calc_dqshift_for_given_th(data(ii).dth_iter,ki,kf,qbragg);
                
                [Psij_mod] = DiffractionPaterns.calc_single_dp( dq_shift(ii,:),probe,rho,X,Y,Z);
                                
                err = sqrt(data(ii).I) - sqrt(Psij_mod);
                err = sum(err(:).^2)/numel(err);
                errtot = errtot + err;
            end
            
            
            
        end
        
        
    end
end
    

