function [sice,qv,ql,qi,ni,t]=analyse_lw_vs_w(fn)
Ra=287;
Rv=461;
cp=1005;
nc=netcdf(fn);

qv=nc{'q'}(:,:,:,1);
ql=nc{'q'}(:,:,:,2);
qi=nc{'q'}(:,:,:,6);
ni=nc{'q'}(:,:,:,7);
th=nc{'theta'}(:,:,:);
p=nc{'p'}(:,:,:);
t=th.*(p./100000.).^(Ra./cp);

eps1=Ra./Rv;
sice=qv./(eps1.*svp(t,'buck2','ice')./p);


close(nc);
