% this scripts compare the evolution of the theta_angle over the angle
% annealings for different distorted points in the rocking curve

figure; 
subplot(321);
jj = 1;
plot(data_exp(index_to_distort(jj)).dth_new_iter,'LineWidth',3.0);
hold on;
plot(repmat(data_exp(index_to_distort(jj)).dth + data_exp(index_to_distort(jj)).dth_delta,numel(data_exp(index_to_distort(jj)).dth_new_iter)),'*r');
title(['angle' num2str(fly2Danglist(index_to_distort(jj))-thBragg)]);

subplot(322);
jj = 10;
plot(data_exp(index_to_distort(jj)).dth_new_iter,'LineWidth',3.0);
hold on;
plot(repmat(data_exp(index_to_distort(jj)).dth + data_exp(index_to_distort(jj)).dth_delta,numel(data_exp(index_to_distort(jj)).dth_new_iter)),'*r');
title(['angle' num2str(fly2Danglist(index_to_distort(jj))-thBragg)]);



subplot(323);
jj = 20;
plot(data_exp(index_to_distort(jj)).dth_new_iter,'LineWidth',3.0);
hold on;
plot(repmat(data_exp(index_to_distort(jj)).dth + data_exp(index_to_distort(jj)).dth_delta,numel(data_exp(index_to_distort(jj)).dth_new_iter)),'*r');
title(['angle' num2str(fly2Danglist(index_to_distort(jj))-thBragg)]);

subplot(324);
jj = 30;
plot(data_exp(index_to_distort(jj)).dth_new_iter,'LineWidth',3.0);
hold on;
plot(repmat(data_exp(index_to_distort(jj)).dth + data_exp(index_to_distort(jj)).dth_delta,numel(data_exp(index_to_distort(jj)).dth_new_iter)),'*r');
title(['angle' num2str(fly2Danglist(index_to_distort(jj))-thBragg)]);

subplot(325);
jj = 40;
plot(data_exp(index_to_distort(jj)).dth_new_iter,'LineWidth',3.0);
hold on;
plot(repmat(data_exp(index_to_distort(jj)).dth + data_exp(index_to_distort(jj)).dth_delta,numel(data_exp(index_to_distort(jj)).dth_new_iter)),'*r');
title(['angle' num2str(fly2Danglist(index_to_distort(jj))-thBragg)]);

subplot(326);
jj = 46;
plot(data_exp(index_to_distort(jj)).dth_new_iter,'LineWidth',3.0);
hold on;
plot(repmat(data_exp(index_to_distort(jj)).dth + data_exp(index_to_distort(jj)).dth_delta,numel(data_exp(index_to_distort(jj)).dth_new_iter)),'*r');
title(['angle' num2str(fly2Danglist(index_to_distort(jj))-thBragg)]);