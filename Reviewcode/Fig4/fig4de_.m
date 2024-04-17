clear; clc; close all;

nx=300; ny=nx;  %numbers of x and y grids
np=6; % number of species
Lx=30; Ly=30;  % length of x and y domains
max_it=100000000; % maximum of time
x=linspace(0,Lx,nx+1)'; x=x(1:end-1);
y=linspace(0,Ly,ny+1); y=y(1:end-1);
h=x(2)-x(1); % interval of x grid
dt=10*h^2; % interval of time
t=dt*(1:1:max_it);
scale=0.02;%scaler
%kinetics
k = [10 8.3 1.7]*100*scale/9; % kf: binding constant, kb: unbinding constant kcat: catalytic constant
km=(k(2)+k(3))/k(1); % Michaelis-Menten (MM) constant

% initial u
u=cell(1,np); % E D S S_P ES DS_P
iu{1}=5*cos(2*x).*ones(1,ny)+100; %E
iu{2}=100*ones(nx,ny);%D
iu{3}=40*ones(nx,1).*cos(2*y)+100; %S
iu{4}=0*iu{1}; iu{5}=0*iu{1}; iu{6}=0*iu{1};

st=mean(iu{3},"all");
%%%
fig1=figure(1);
TitleStr={'D','S/D','E/D'};
subplot(1,3,1);
surf(x,y,iu{2}'); shading interp;
title(TitleStr{1}); set(gca,'FontSize',15); alpha(0.9);
xlabel('x');ylabel('y');view(0,90);
subplot(1,3,2);
surf(x,y,(iu{3}./iu{2})'); shading interp;
title(TitleStr{2}); set(gca,'FontSize',15); alpha(0.9);
xlabel('x');ylabel('y');view(0,90);
subplot(1,3,3);
surf(x,y,(iu{1}./iu{2})'); shading interp;
title(TitleStr{3}); set(gca,'FontSize',15); alpha(0.9);
xlabel('x');ylabel('y');view(0,90);
colormap(flip(magma_white));
set(fig1,'position',[100 100 900 300])


% diffusion & reaction functions
lap=-(2*pi*x).^2-(2*pi*y).^2;
nrea1=@(ET,kcat,km,s)(kcat*ET.*s./(km+s));%sQ
nrear1=@(ET,kcat,km,s_hat)( 0.5*kcat*(ET+km+s_hat-sqrt((ET+km+s_hat).^2-4*ET.*s_hat)) );%tQ

%%% pde solution
D=0.2*scale*ones(np,1);D(1)=0;D(2)=0;
fu=iu; tu=iu; su=iu;
for it=1:max_it

    %%% sQSSA
    for ip=1:np
        su{ip}=real(idct2( dct2(su{ip})./(1-dt*D(ip)*lap) ));
    end
    et1=su{1}+su{5};
    dt1=su{2}+su{6};
    r1=nrea1(et1,k(3),km,su{3});
    r2=nrea1(dt1,k(3),km,su{4});
    su{3}=su{3}+dt*(-r1+r2);
    su{4}=su{4}+dt*(r1-r2);

    %%% tQSSA
    for ip=1:np
        tu{ip}=real(idct2( dct2(tu{ip})./(1-dt*D(ip)*lap) ));
    end
    et2=tu{1}+tu{5};
    dt2=tu{2}+tu{6};
    s_hat=tu{3}+tu{5};
    sp_hat=tu{4}+tu{6};
    r3=nrear1(et2,k(3),km,s_hat);
    r4=nrear1(dt2,k(3),km,sp_hat);
    tu{3}=tu{3}+dt*(-r3+r4);
    tu{4}=tu{4}+dt*(r3-r4);

    %%% Full Model
    for ip=1:np
        fu{ip}=real(idct2( dct2(fu{ip})./(1-dt*D(ip)*lap) ));
    end
    r1=-k(1)*fu{1}.*fu{3}+(k(2)+k(3))*fu{5};
    r2=-k(1)*fu{2}.*fu{4}+(k(2)+k(3))*fu{6};
    r3=-k(1)*fu{1}.*fu{3}+k(2)*fu{5}+k(3)*fu{6};
    r4=-k(1)*fu{2}.*fu{4}+k(2)*fu{6}+k(3)*fu{5};
    r5=-r1;
    r6=-r2;
    fu{1}=fu{1}+dt*r1;
    fu{2}=fu{2}+dt*r2;
    fu{3}=fu{3}+dt*r3;
    fu{4}=fu{4}+dt*r4;
    fu{5}=fu{5}+dt*r5;
    fu{6}=fu{6}+dt*r6;
end

fig3=figure(3);
subplot(1,3,1);
surf(x,y,(su{4}./st)'); shading interp;
title('sQSSAp'); set(gca,'FontSize',15); %alpha(0.9);
set(gca,'XTick',[0 0.5 1]);set(gca,'YTick',[0 0.5 1]);
xlabel('x');ylabel('y');
axis image;view(0,90);caxis([0,1]);colorbar;
subplot(1,3,2);
surf(x,y,(tu{4}./st)'); shading interp;
title('tQSSAp'); set(gca,'FontSize',15); %alpha(0.9);
set(gca,'XTick',[0 0.5 1]);set(gca,'YTick',[0 0.5 1]);
axis image;view(0,90);caxis([0,1]);colorbar;colormap(deep);
xlabel('x');ylabel('y');
subplot(1,3,3);
surf(x,y,((fu{4}+fu{6})./st)'); shading interp;
title('Full'); set(gca,'FontSize',15); %alpha(0.9);
set(gca,'XTick',[0 0.5 1]);set(gca,'YTick',[0 0.5 1]);
axis image;view(0,90);caxis([0,1]);colorbar;colormap(deep);
xlabel('x');ylabel('y');


