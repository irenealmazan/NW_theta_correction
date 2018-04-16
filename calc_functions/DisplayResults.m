classdef DisplayResults
    % This library contains all the functions which allow us to rotate the
    % detector and the sample
    properties(Constant)
    end
    
    
    methods(Static)
        
        function [theta_iter] = read_angles_iterations(data,ini_angle_val)
            
            % read values for the different angles and the iterations
            for jj = 1:numel(data)
                for kk = 1:numel(data(jj).theta_iter)
                    theta_iter(jj).dth_new_iter(kk) = data(jj).theta_iter(kk).dth_new_iter;
                    theta_iter(jj).beta(kk) = data(jj).theta_iter(kk).beta;
                    theta_iter(jj).grad_final(kk) = data(jj).theta_iter(kk).grad_final;
                    theta_iter(jj).dqshift(kk,:) = data(jj).theta_iter(kk).dqshift;
                    theta_iter(jj).delta_dth(kk) = data(jj).theta_iter(kk).dth_new_iter-ini_angle_val(jj);
                end
            end
            
            
        end
               
        function [h] = display_angles_iterations(theta_iter,data,angles_to_plot)
            
            fig_num = [1:numel(angles_to_plot)];
            
            for ff = 1:numel(fig_num)
                figure(fig_num(ff))
                clf;
                
                for jj = 1:numel(theta_iter)
                    
                    if jj == angles_to_plot(ff)%mod(jj,5) == 0
                        subplot(3,2,1);
                        hold on;
                        plot(theta_iter(jj).dth_new_iter,'-.');
                        plot(data(jj).dth_real*ones(numel(data(jj).theta_iter),1),'r');
                        xlabel('iterations');ylabel('theta');
                        title(['theta vs iterations @ theta = ' num2str(data(jj).dth_nominal)]);
                        
                        
                        subplot(3,2,2);
                        hold on;
                        plot(theta_iter(jj).beta)
                        xlabel('iterations');ylabel('beta');
                        title('beta vs iterations');
                        
                        subplot(3,2,3);
                        hold on;
                        plot(theta_iter(jj).grad_final)
                        xlabel('iterations');ylabel('gradient');
                        title('gradient vs iterations');
                        
                        subplot(3,2,4);
                        hold on;
                        plot(theta_iter(jj).dqshift(:,1),'.');
                        plot(data(jj).dqshift_real(1)*ones(numel(data(jj).theta_iter),1),'r');
%                         plot(theta_iter(jj).dqshift(:,2),'.');
%                         plot(data(jj).dqshift_real(2)*ones(numel(data(jj).theta_iter),1),'r');
%                         plot(theta_iter(jj).dqshift(:,3),'.');
%                         plot(data(jj).dqshift_real(3)*ones(numel(data(jj).theta_iter),1),'r');
                        xlabel('iterations');ylabel('dq components');
                        title('dqshift vs iterations');
                        legend('qx','real qx');%,'qy','real qy','qz','real qz');
                        
                         subplot(3,2,5);
                        hold on;
                        plot(theta_iter(jj).delta_dth)
                        xlabel('iterations');ylabel('delta_dth');
                        title('delta_dth vs iterations');
                        
                    end
                    
                end
              
                
            end
            
            
        end
        
        function [] = show_rho_theta_update(figure_num,errlist,rho,midsl,anglelist,trueangles,norm_grad_rho,beta_rho,norm_grad_theta,beta_theta,flag)
                       % this function creates the figure where we plot the improvements of the
            % rho and the angles during the phase retrieval iterations. flagINI
            % indicates if the figure needs to be created (initial state) or only
            % updated
            
            figure(figure_num);
            if strcmp(flag,'Ini')
                clf;
                setfigsize(gcf, 1000,500);
            end
            % plot
            subplot(231); imagecomp(rho(:,:,midsl)); colorbar; axis image;title('Longitudinal section'); %zoom(1.5);
            subplot(232); imagecomp(squeeze(rho(100,:,:))); colorbar; axis image;title('Transversal section') %zoom(1.5);
            subplot(233); plot(log10(errlist),'LineWidth',3.0);title('Error metric');xlabel('Iterations')
            subplot(234); plot(anglelist,'ob');title('Angles');xlabel('angle index');
            hold on; plot(trueangles,'*r');
            
            if strcmp(flag,'rho')
                subplot(235); plot(norm_grad_rho,'ob');title('norm of grad\_rho and beta');
                yyaxis right;
                plot(beta_rho,'*k');
                ax = gca;
                set(ax.YAxis(1),'Color','b');
                yyaxis left;
                ylabel('norm of gradien');
                yyaxis right;
                ylabel('beta');
                xlabel('iterations');
            elseif strcmp(flag,'theta')    
                subplot(236); plot(norm_grad_theta,'ob');title('norm of grad\d_theta and beta');
                yyaxis right;
                plot(beta_theta,'*k');
                set(ax.YAxis(1),'Color','b');
                yyaxis left;
                ylabel('norm of gradien');
                yyaxis right;
                ylabel('beta');
                xlabel('iterations');
            end
            
            drawnow;
            
        end
        
        function [h,lgd] = show_rocking_curve(angles,rock_curve,flag,fig_num,color,legend_str)
           % this function displays the rocking curve with respect to a specified angular grid
           
           if strcmp(flag,'new')
               figure(fig_num);
               clf;
               h = bar(angles,rock_curve,'FaceColor',color,'BarWidth',0.2,'EdgeColor',color);
               lgd = legend(legend_str);
           elseif strcmp(flag,'hold')
               figure(fig_num);
               ax = gca;
               hold on;
               h = bar(angles,rock_curve,'FaceColor',color,'BarWidth',0.2,'EdgeColor',color);
               ax.Legend.String{size(ax.Legend.String,2)} = legend_str;
               legend(ax.Legend.String);
           end
           
        end
        
    end
end