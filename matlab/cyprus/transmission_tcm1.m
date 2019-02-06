iqc=2;
inc=4;
iqc=12;
inc=11;
z=nc{'z'}(:);
x=nc{'x'}(:);
z1=z(end:-1:1)+400;
% data1=zeros([600,length(z),length(x)]);
data1=zeros([600,length(z)]);
for i=1:600
    m1=nc{'q'}(i,:,:,iqc)'./nc{'q'}(i,:,:,inc)';
    d1=(m1.*6./(pi.*1000)).^(1./3);
    beta1=nc{'p'}(i,:,:)'./nc{'t'}(i,:,:)'./287.*nc{'q'}(i,:,:,inc)'.*d1.^2.*pi./2;
    beta1(find(isnan(beta1(:))))=0;       
%     albedo=sqrt(3).*0.15*cumsum((beta1(end:-1:1,1:end).*10))./(2.+sqrt(3)*0.15.*cumsum((beta1(end:-1:1,1:end).*10)));
    albedo=sqrt(3).*0.15*cumsum(prctile(beta1(end:-1:1,1:end)'.*10,50))./(2.+sqrt(3)*0.15.*cumsum(prctile(beta1(end:-1:1,1:end)'.*10,50)));

    trans=1-albedo;
%     data1(i,:,:)=trans;
    data1(i,:)=trans;
end

fillx=[data1(1,:), data1(1,end), data1(600,end:-1:1)];
filly=[z1(1:end);z1(end);z1(end:-1:1)];
   
h=fill(fillx,filly,'c')
set(h,'facealpha',0.5)