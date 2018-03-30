% this script plots the distance dq^2 = Q^2_true - Q^2_achieved vs 
% vs the angular positions at different iterations


figure;plot(fly2Danglist-thBragg,mod_ini_square,'xb');
hold on
plot(fly2Danglist-thBragg,mod_fin_square_inter(1,:),'.m')
plot(fly2Danglist-thBragg,mod_fin_square_inter(2,:),'*g')
plot(fly2Danglist-thBragg,mod_fin_square_inter(10,:),'ok')
plot(fly2Danglist-thBragg,mod_fin_square_inter(18,:),'*y')
plot(fly2Danglist-thBragg,mod_fin_square_inter(30,:),'og')
plot(fly2Danglist-thBragg,mod_fin_square_inter(40,:),'og')
plot(fly2Danglist-thBragg,mod_fin_square_inter(47,:),'og')
plot(fly2Danglist-thBragg,mod_fin_square_inter(46,:),'og')

legend('initial position','1 angle annealing','2 angle annealing','10 angle annealings','18 annealings','30 annealings','40 annealings','46 annealings (end)')

xlabel('rocking angle');ylabel('\Delta q^2 in Angstroms^{-2}');

savefig('results/dqdistane_vs_angularpos_nois')

%print -dpdf '/Users/ialmazn/Box Sync/ptycho_pos_grad/pictures/dq_vs_angularposition_approach_45angles_noise' '-bestfit'
