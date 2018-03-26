function [x,z,xn,zn,phi,u,w,wij,uij]=convection_analytical
% quiver(x(1:end-1,1:end-1),z(1:end-1,1:end-1),real(u(:,1:end-1)),real(w(1:
% end-1,:)));
%pcolor(x(1:end-1,1:end-1),z(1:end-1,1:end-1),real(u(:,1:end-1)));shading
%flat
%pcolor(x(1:end-1,1:end-1),z(1:end-1,1:end-1),real(w(1:end-1,:)));shading
%flat

resx=50;
resz=50;

xs=0:resx:10000;
zs=0:resz:7000;


xns=xs+resx./2;
zns=zs+resz./2;

[x,z]=meshgrid(xs,zs);
[xn,zn]=meshgrid(xns,zns);
zc=xns(1);
% some definitions:
Md=29e-3;
Mv=18e-3;
beta=Md./Mv-1.;
ttr=273.15;
alpha=1./ttr;

% some constants from figure 1:
k=0.001875577703636; % changes the width
dsm_by_dz_z_eq_zc=-1.6e-6;
b=1e-6; % range from 0 to 5e-6, default 1e-6
del_gamma_mac=0.5e-4;
del_c_s=0;
del_c_T=.5;  % changing alters height
epsilon=3e-7; % changing this alters the height too
cp=1005;
L=2.5e6;
K=epsilon.*cp./L;
grav=9.81;


% equation 32
N_bar_mac=sqrt(grav.*(alpha.*del_gamma_mac+beta.*(dsm_by_dz_z_eq_zc+b)));


% equation 33
z_bar=(alpha.*del_c_T+beta.*del_c_s)./ ...
    (alpha.*del_gamma_mac+beta.*(dsm_by_dz_z_eq_zc+b));

k1=0.5.*(alpha.*epsilon-beta.*K)./ ...
    (alpha.*del_gamma_mac+beta.*(dsm_by_dz_z_eq_zc+b));

% equation 39 of ZZRA
Z=k.*N_bar_mac.* ...
    sqrt(2.*(z-zc).*(z_bar+(z-zc)./2-k1./3.*(z-zc).^2));

% equation 41 of ZZRA
X=1./k.^2.*cos(k.*x);

% equation 42 of ZZRA
phi=(Z).*X;

% u=-(phi(2:end,:)-phi(1:end-1,:))./resz;
% w=(phi(:,2:end)-phi(:,1:end-1))./resx;
small=1e-30;

u=-N_bar_mac./k.*(z_bar+(z-zc)-k1.*(z-zc).^2) ./ ...
    sqrt(small+2.*(z-zc).*(z_bar+(z-zc)./2-k1./3.*(z-zc).^2)).*sin(k.*xn);
% u(1,:)=0;

w=(-N_bar_mac.* ...
    sqrt(2.*(zn-zc).*(z_bar+(zn-zc)./2-k1./3.*(zn-zc).^2)).*-cos(k.*x));

wij=real(-N_bar_mac.* ...
    sqrt(2.*(z-zc).*(z_bar+(z-zc)./2-k1./3.*(z-zc).^2)).*-cos(k.*x));

uij=zeros(size(wij));
[r,c]=size(wij);
for i=2:c
    uij(2:end,i)=uij(2:end,i-1)-(wij(2:end,i-1)-wij(1:end-1,i-1)).*resx./resz;
end


