clear; clc; close all;

np=4;
nx=101; L=10; x=linspace(0,L,nx); h=x(2)-x(1);
dt=.0001*h; max_it=500; t=dt*(1:1:max_it);

%kinetics
k(1)=1.7*10^(4);
k(2)=3*10^(4);
k(3)=1.6*10^(3);

km=(k(2)+k(3))/k(1);
% initial u
iu=zeros(np,nx);
iu(1,:)=1.85*(tanh(20*(x-9.5))+1)+0.00001;
iu(2,:)=1.85*(tanh(20*(x-9.5))+1)+0.00001;

figure(1)
semilogy(x,iu(1,:),'k--','LineWidth',3);hold on;
semilogy(x,iu(2,:)+km,'k-','LineWidth',3)
xlabel('x(\mum)');ylabel('Concentration (mM)');
axis([0 10 10^(-5.5) 10])
% legend('E','S','S+K_M','location','best');legend boxoff;
set(gca,'fontsize',22)
set(gca,'TickDir','out')
set(gca,'XTick',[0 5 10]);set(gca,'YTick',[1e-5 1e-3 1e-1 10])
box off

% diffusion & reaction functions
lap=-((0:nx-1)*2*pi/L).^2;
rea1=@(e,s,c)(-k(1)*e.*s+k(2)*c);
rea2=@(c)(k(3)*c);
nrea1=@(et1,s,p,c)(k(3)*et1.*s./(km+s));
nrear1=@(et2,st2,p)(.5*k(3)*(et2+km+st2-p-sqrt((et2+km+st2-p).^2-4*et2.*(st2-p))));
for im=1:10
    multi=10^(im-1);

    D=1e-2/9*multi*ones(np,1);
    % return;
    fu=iu;tu=iu;su=iu;
    for it=1:max_it
        %%% sQ
        su=real(idct( dct(su')./(1-dt*D'.*lap') ))';
        et1=sum(su([1,np],:));
        r1=nrea1(et1,su(2,:),su(3,:),su(np,:));
        su(2:3,:)=su(2:3,:)+dt*[-r1;r1];
        su(4,:)=et1.*su(2,:)./(su(2,:)+km);
        su(1,:)=et1-su(4,:);
        sup(:,it)=su(3,:);

        %%% TQ
        tu=real(idct( dct(tu')./(1-dt*D'.*lap') ))';
        et2=sum(tu([1,np],:));
        st2=sum(tu([2,3,np],:));
        r2=nrear1(et2,st2,tu(3,:));

        tu(2:3,:)=tu(2:3,:)+dt*[-r2;r2];

        tu(4,:)=nrear1(et2,st2,tu(3,:))/k(3);
        tu(1,:)=et2-tu(4,:);
        tu(2,:)=st2-tu(4,:)-tu(3,:);
        tup(:,it)=tu(3,:);


        %%%Full Model
        fu=real(idct( dct(fu')./(1-dt*D'.*lap') ))';
        r1=rea1(fu(1,:),fu(2,:),fu(np,:));
        r2=rea2(fu(np,:));
        fu(1:np,:)=fu(1:np,:)+dt*[ r1+r2;r1; r2; -r1-r2];
        fup(:,it)=fu(3,:);

    end
    %%
    figure(2);
    maxy=0.21;
    subplot(2,5,im);
    plot(t,mean(fup),'-', 'color',"#C3C3C3",'linewidth',3);hold on;
    plot(t,mean(sup), '--', 'color',"#22B14C",'linewidth',3);hold on;
    plot(t,mean(tup),'--', 'color',"#EA3680",'linewidth',3);hold on;
    % legend( 'Full', 'sQSSAp', 'tQSSAp', 'location','best');legend boxoff;
    axis([0 max_it*dt 0 maxy])
    xlim([0 max_it*dt])
%     xlabel('Time (s)');ylabel('$\bar{P}$ (mM)','Interpreter','latex');
    set(gca,'FontSize',15)
    set(gca,'TickDir','out')
    % title('D=',mean(D))
    set(gca,'XTick',[0 0.01])
    set(gca,'YTick',[0 0.1 0.2])
%     fontname(gcf,"Arial")
    box off;

end