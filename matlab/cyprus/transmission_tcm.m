m1=nc{'q'}(i,:,:,2)'./nc{'q'}(i,:,:,4)';
d1=(m1.*6./(pi.*1000)).^(1./3);
beta1=nc{'q'}(i,:,:,4)'.*d1.^2.*pi./2;
beta1(find(isnan(beta1(:))))=0;       
albedo=sqrt(3).*0.15*cumsum(nanmean(beta1(end:-1:1,1:260)'.*10))./(2.+sqrt(3)*0.15.*cumsum(nanmean(beta1(end:-1:1,1:260)'.*10)));

trans=1-albedo;
z=nc{'z'}(:);
z1=z(end:-1:1)+400;
plot(trans,z1,'k','linewidth',3)