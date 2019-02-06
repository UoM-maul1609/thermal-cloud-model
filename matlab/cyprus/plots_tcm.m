nc=netcdf('/tmp/tcm_output_mixing_theta_q_cons.nc');
figure('position',[1         503        1204         163]);
k=1;
for i=1:10:600
   subplot(131);
   pcolor(nc{'x'}(:),nc{'z'}(:)+400,nc{'q'}(i,:,:,4)'./1e6);shading flat;caxis([0 400]);
   h=colorbar;
   ylabel(h,'N_d (cm^{-3})');
   xlabel('x (m)');
   ylabel('z (m)');
   h=streamslice(nc{'x'}(:),nc{'z'}(:)+400,nc{'u'}(i,:,:)',nc{'w'}(i,:,:)');
   set(h,'color','r');
   
   subplot(132);
   pcolor(nc{'x'}(:),nc{'z'}(:)+400,nc{'q'}(i,:,:,2)'.*1e3);shading flat;caxis([0 1.2]);
   h=colorbar;
   ylabel(h,'q_l (g m^{-3})');
   xlabel('x (m)');
   ylabel('z (m)');
   h=streamslice(nc{'x'}(:),nc{'z'}(:)+400,nc{'u'}(i,:,:)',nc{'w'}(i,:,:)');
   set(h,'color','r');
   text(0.1,0.9,['time=',num2str(nc{'time'}(i)),' s'],'units','normalized','color',[1 1 1]);
   
   subplot(133);
   pcolor(nc{'x'}(:),nc{'z'}(:)+400,nc{'precip'}(i,:,:,1)');shading flat;%caxis([0 1.2]);
   h=colorbar;
   ylabel(h,'P (mm hr^{-1})');
   xlabel('x (m)');
   ylabel('z (m)');
   h=streamslice(nc{'x'}(:),nc{'z'}(:)+400,nc{'u'}(i,:,:)',nc{'w'}(i,:,:)');%caxis([0 1]);
   set(h,'color','r');
   print('-dpng', ['/tmp/mix_400_theta_q',num2str(k,'%03d'),'.png'])
   k=k+1;
end

close(nc);
