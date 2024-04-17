clear; clc; close all;

np=4;% number of species
nx=300; % number of x grids
L=30; % length of domain
x=linspace(0,L,nx);
h=x(2)-x(1); % interval of x grid
dt=0.1*h^2; % interval of time
max_it=11920; % maximum of time
t=dt*(1:1:max_it);

%kinetics
k(1)=0.67; % kf: binding constant
k(2)=0.53; % kb: unbinding constant
k(3)=0.13; % kcat: catalytic constant

km=(k(2)+k(3))/k(1); % Michaelis-Menten (MM) constant


% initial condition
iu=zeros(np,nx); %E S P C
et=5; st=et*8-km*10;
iu(2,:)=st*ones(1,nx);


l=0;
% standard deviation = sigma
sigmalist=[25 12.5 10 8 7 6 5.5 5 4.5 4 3.5 3 2.5 2 1.5 1 0.9 0.8 0.7 0.6 0.5 0.3 0.1];

for is=1:length(sigmalist)
    sigma=sigmalist(is);
    data= et/mean(normpdf(x,15,sigma))*normpdf(x, 15, sigma);
    iu(1,:)=data; % iu(1,:)=et*ones(1,nx);
    iu(2,:)=st*ones(1,nx);

    % diffusion & reaction functions
    lap=-((0:nx-1)*2*pi/L).^2;
    rea1=@(e,s,c)(-k(1)*e.*s+k(2)*c);
    rea2=@(c)(k(3)*c);
    nrea1=@(et1,s,p,c)(k(3)*et1.*s./(km+s));
    nrear1=@(et2,st2,p)(.5*k(3)*(et2+km+st2-p-sqrt((et2+km+st2-p).^2-4*et2.*(st2-p))));


    %%% pde solution
    D=0.2*ones(np,1); D(1)=0;D(4)=0; %diffusion coefficients
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
    l=l+1;
    mfp(l,:)=mean(fup);msp(l,:)=mean(sup);mtp(l,:)=mean(tup);

end
%%
inittime=max_it;t(inittime)
for il=1:length(sigmalist)
    initialvelocity(il,:)=[(mfp(il,inittime)-mfp(il,1))/(t(inittime)-t(1)),...
        (mtp(il,inittime)-mtp(il,1))/(t(inittime)-t(1)),...
        (msp(il,inittime)-msp(il,1))/(t(inittime)-t(1))];
end

%%

figure(4);clf;
plot(10*log10(25./sigmalist),initialvelocity(:,1),'-', 'color',"#C3C3C3",'linewidth',3);hold on;
plot(10*log10(25./sigmalist),initialvelocity(:,2),'--', 'color',"#EA3680",'linewidth',3);
plot(10*log10(25./sigmalist),initialvelocity(:,3),'--', 'color',"#22B14C",'linewidth',3);
xlabel('homogeneity');ylabel('initial velocity (\muM^2/{s})');xlim([0 25]);
legend('Full','tQSSA_p','sQSSA_p','location','southwest');legend boxoff;
set(gca,'XTick',[0 12.5 25]);
set(gca,'YTick',[0.25 0.45 0.65]);
set(gca,'FontSize',15);
axis([0 25 0.25 0.65]);
box off;
set(gca, 'TickDir','out');
