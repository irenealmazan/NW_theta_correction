% this script plots the distance dq^2 = Q^2_true - Q^2_achived vs 
% vs the number of angle annealign


figure;
hold on; plot(mod_fin_square_inter(:,index_to_distort(1)),'LineWidth',3.0);
hold on; plot(mod_fin_square_inter(:,index_to_distort(10)),'LineWidth',3.0);
hold on; plot(mod_fin_square_inter(:,index_to_distort(20)),'LineWidth',3.0);
hold on; plot(mod_fin_square_inter(:,index_to_distort(30)),'LineWidth',3.0);
hold on; plot(mod_fin_square_inter(:,index_to_distort(40)),'LineWidth',3.0);
hold on; plot(mod_fin_square_inter(:,index_to_distort(46)),'LineWidth',3.0);

xlabel('angle iterations');ylabel('dq^2 in Angstroms^{-2}');


legend(['angle' num2str(fly2Danglist(index_to_distort(1))-thBragg)] ,...
    ['angle' num2str(fly2Danglist(index_to_distort(10))-thBragg)],...
    ['angle' num2str(fly2Danglist(index_to_distort(20))-thBragg)],...
    ['angle' num2str(fly2Danglist(index_to_distort(30))-thBragg)],...
    ['angle' num2str(fly2Danglist(index_to_distort(40))-thBragg)],...
    ['angle' num2str(fly2Danglist(index_to_distort(46))-thBragg)]);

%savefig('results/dqdistance_vs_annealingsiter_cstangle_noise');

%print -dpdf '/Users/ialmazn/Box Sync/ptycho_pos_grad/pictures/dqdistance_vs_annealingsiter_cstangle_noise' '-bestfit'