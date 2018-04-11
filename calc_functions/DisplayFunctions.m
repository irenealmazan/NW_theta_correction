classdef DisplayFunctions
    % This library contains all the functions which allow us to rotate the
    % detector and the sample
    properties(Constant)
    end
    
    
    methods(Static)
        
        function h = display_diff_geom(rho,ki,kf,qbragg,fig_num,X,Y,Z)
            
            if ishandle(fig_num)
                handleFLAG = get(fig_num,'type');
                switch handleFLAG
                    case 'figure'
                        figure(fig_num);
                        
                    case 'axes'
                        subplot(fig_num);
                end
            else
                figure(fig_num)
            end
            
            clf;
            
            hold on;
            
            % the sample 
            h=di(angle(rho), -.5, 'y', X,Y,Z); alpha(h,.5);
            
            % the wavevectors
            quiver3(0,0,0, ki(1), ki(2), ki(3), 'r');
            quiver3(0,0,0, kf(1), kf(2), kf(3), 'k');
            quiver3(0,0,0, qbragg(1), qbragg(2), qbragg(3), 'b');
            
            % the detector 
            [Xd Yd] = meshgrid([-.1 .1]);
            surf(Xd,Yd,ones(size(Xd)));
            view(-2,53);
            
            xlabel('x');ylabel('y'); zlabel('z');
            axis image
            colormap default
            colorbar
            title('units of displacement (microns)')
                      
        end
        
         function h = display_strain(strptsmid,phvals,fig_num)
            
             
             if ishandle(fig_num)
                 handleFLAG = get(fig_num,'type');
                 switch handleFLAF
                     case 'figure'
                         figure(fig_num);
                         
                     case 'axes'
                         subplot(fig_num);
                 end
             else
                 figure(fig_num)
             end
             
             
            clf;
            
            hold on;
            
            % the strain
            h = scatter3(strptsmid(:,1), strptsmid(:,2), strptsmid(:,3), [], phvals(:));

                      
        end
        
    end
end