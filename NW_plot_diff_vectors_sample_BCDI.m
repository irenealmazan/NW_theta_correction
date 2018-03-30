% In this script we plot the incoming, outgoing wavevectors, the
% corresponding momentum transfer and the shifts in momentum transfer
% associated to the rocking curve measurement

% plot the frame:
%%{
[X_square,Y_square,Z_square] = meshgrid([-Npix/2 Npix/2-1]*d2_bragg,[-Npix/2 Npix/2-1]*d2_bragg,[-depth/2 depth/2-1]*d2_bragg);
X_square_toplot = X_square(:);
Y_square_toplot = Y_square(:);
Z_square_toplot = Z_square(:);

% display the geometry:



figure(3);
clf; 
hold on;
quiver3(0,0,0, ki(1), ki(2), ki(3), 'r');
quiver3(0,0,0, kf(1), kf(2), kf(3), 'k');
quiver3(0,0,0, qbragg(1), qbragg(2), qbragg(3), 'b');

if plotdqshif
    for ii=1:numel(delta_thscanvals)   
        quiver3(0,0,0, ki_list(ii,1), ki_list(ii,2), ki_list(ii,3), 'r--');
        quiver3(0,0,0, dqlist(ii,1), dqlist(ii,2), dqlist(ii,3), 'r--');            
    end   
end

[Xd Yd] = meshgrid([-.1 .1]);
surf(Xd,Yd,ones(size(Xd)));
view(-2,53);
scatter3(X_square_toplot,Y_square_toplot,Z_square_toplot)
h=di(NW, -.5, 'y', X,Y,Z); alpha(h,.5);
axis image
colormap hot
hold on;
xlabel('x');ylabel('y'); zlabel('z');
