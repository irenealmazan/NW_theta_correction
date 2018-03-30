% this scripts displays the evolution of the angular position as a function
% of the iterations of the phase retrieval algorithm


for jj = index_to_distort
    
   dth_new_iter = [data_exp(jj).dth data_exp(jj).dth_new_iter];
   dth_steps = diff(dth_new_iter);
   dth_true = ones(size(dth_new_iter,2),1).*(data_exp(jj).dth + data_exp(jj).dth_delta);
   
   figure;
   subplot(121);
   hold on;
   plot(dth_new_iter,'-xb');
   plot(dth_true,'*r');
   title('angle vs iterations')
   
   subplot(122);
   plot(dth_steps,'-ob');
   title('step\_angles vs iterations');
    
end