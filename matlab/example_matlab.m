
printflag=false;

nc1=netcdf('/tmp/output1.nc')

for i=10:300
    subplot(121)
    contourf(nc1{'x'}(:),nc1{'z'}(:),nc1{'precip'}(i,:,:,1)',[0:5:30]);shading flat
    h=streamslice(nc1{'x'}(:),nc1{'z'}(:),nc1{'u'}(1,:,:)',nc1{'w'}(1,:,:)');
    set(h,'color',[1 1 1]);
    caxis([0 30]);
    colorbar
    title(['time: ',num2str(nc1{'time'}(i))])
    
    
    subplot(122)
    contourf(nc1{'x'}(:),nc1{'z'}(:),nc1{'q'}(i,:,:,7)',[0:10:50].*1000);shading flat
    caxis([0 50000])
    colorbar
    pos=get(gca,'position');
    h2=axes('position',pos,'visible','off');
    [c,h]=contour(nc1{'x'}(:),nc1{'z'}(:),nc1{'q'}(i,:,:,2)',[0:.2:2]./1000,'linecolor',[1 1 1])
    set(h2,'color','none');
    caxis([0 2e-3]);

    
    if ~printflag
        pause(0.05)
    end
    
    if printflag
        eval(['print -dpng /tmp/hm_picture', num2str(i,'%03d'),'.png']);
    end
    delete([h2]);
end

close(nc1);
close(nc2);
