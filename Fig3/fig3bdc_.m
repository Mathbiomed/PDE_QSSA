clear; clc; close all;

np=4;% number of species
nx=300; % number of x grids
L=30; % length of domain
x=linspace(0,L,nx);
h=x(2)-x(1); % interval of x grid
dt=0.1*h^2; % interval of time
max_it=11921; % maximum of time
t=dt*(1:1:max_it);

%kinetics
k(1)=0.67; % kf: binding constant
k(2)=0.53; % kb: unbinding constant
k(3)=0.13; % kcat: catalytic constant

km=(k(2)+k(3))/k(1); % Michaelis-Menten (MM) constant


% initial condition
iu=zeros(np,nx); %E S P C
et=5; st=et*8-km; data=zeros(1,nx);
data= et/mean(normpdf(x, 15, 0.2))*normpdf(x, 15, 0.2);
iu(1,:)=data;
iu(2,:)=st*ones(1,nx);

% diffusion & reaction functions
lap=-((0:nx-1)*2*pi/L).^2;
rea1=@(e,s,c)(-k(1)*e.*s+k(2)*c);
rea2=@(c)(k(3)*c);
nrea1=@(et1,s,p,c)(k(3)*et1.*s./(km+s));
nrear1=@(et2,st2,p)(.5*k(3)*(et2+km+st2-p-sqrt((et2+km+st2-p).^2-4*et2.*(st2-p))));

%%% pde solution
D=0.2*ones(np,1); D(1)=0;D(4)=0;
% return;
fu=iu;tu=iu;su=iu;
for it=1:max_it
    %%% sQSSA
    su=real(idct( dct(su')./(1-dt*D'.*lap') ))'; % dct for diffusion
    et1=sum(su([1,np],:)); % total enzyme concentration
    r1=nrea1(et1,su(2,:),su(3,:),su(np,:));
    su(2:3,:)=su(2:3,:)+dt*[-r1;r1];
    su(4,:)=et1.*su(2,:)./(su(2,:)+km); % C calculation
    su(1,:)=et1-su(4,:); % E=ET-C
    sup(:,it)=su(3,:); % P concentration save

    %%% tQSSA
    tu=real(idct( dct(tu')./(1-dt*D'.*lap') ))'; % dct for diffusion
    et2=sum(tu([1,np],:)); % total enzyme concentration
    st2=sum(tu([2,3,np],:)); % total substrate concentration
    r2=nrear1(et2,st2,tu(3,:));

    tu(2:3,:)=tu(2:3,:)+dt*[-r2;r2];

    tu(4,:)=nrear1(et2,st2,tu(3,:))/k(3); % C calculation
    tu(1,:)=et2-tu(4,:); % E=ET-C
    tu(2,:)=st2-tu(4,:)-tu(3,:); %S=ST-C-P
    tup(:,it)=tu(3,:); % P concentration save

    %%%Full Model
    fu=real(idct( dct(fu')./(1-dt*D'.*lap') ))'; % dct for diffusion
    r1=rea1(fu(1,:),fu(2,:),fu(np,:));
    r2=rea2(fu(np,:));
    fu(1:np,:)=fu(1:np,:)+dt*[ r1+r2;r1; r2; -r1-r2];
    fup(:,it)=fu(3,:); % P concentration save
end
save fig3bcd
%%
figure(1);clf;
plot(x,iu(1,:),':k','LineWidth',3);hold on;plot(x,iu(2,:)+km,'-k','LineWidth',3)
xlabel('x (\mum)');ylabel('Concentration (\muM)');
legend('$E$','$S+K_M$','Interpreter','latex');legend boxoff;
axis([0 30 0 400]);set(gca,'YTick',[0 200 400]);set(gca,'XTick',[0 10 20 30]);
set(gca,'fontsize',15);set(gca,'TickDir','out')
%%
figure(2);
plot(t,mean(fup),'-', 'color',"#C3C3C3",'LineWidth',3);hold on;
plot(t,mean(tup),'--', 'color',"#EA3680",'LineWidth',3);
plot(t,mean(sup),'--', 'color',"#22B14C",'LineWidth',3);
axis([0 12 0 10])
legend('Full','tQSSAp','sQSSAp','location','southeast');legend boxoff;
xlabel('time (s)');
ylabel('$\bar{P}$($\mu$M)','Interpreter','latex');
set(gca,'fontsize',15);
set(gca,'XTick',[0 12]);set(gca,'YTick',[0 5 10]);
set(gca,'TickDir','out')
%%
figure(3);
load('magma.mat');colormap(flip(magma_white));
subplot(1,3,1)
mesh(x,t,sup');shading interp;view(0,90);
xlabel('x (\mum)');ylabel('time (s)');title('sQSSAp');
axis([0 30 0 max_it*dt]);
set(gca,'XTick',[0 15 30]);set(gca,'YTick',[0 12]);
colorbar;caxis([0 30]);

subplot(1,3,2)
mesh(x,t,fup');shading interp;view(0,90);
xlabel('x (\mum)');ylabel('time (s)');title('Full');
set(gca,'FontSize',15);axis([0 30 0 max_it*dt]);
set(gca,'XTick',[0 15 30]);set(gca,'YTick',[0 12]);
colorbar;caxis([0 30]);

subplot(1,3,3)
mesh(x,t,tup');shading interp;view(0,90);
xlabel('x (\mum)');ylabel('time (s)');title('tQSSAp');
set(gca,'FontSize',15);axis([0 30 0 max_it*dt]);
set(gca,'XTick',[0 15 30]);set(gca,'YTick',[0 12]);
colorbar;caxis([0 30]);

