%% estimate errors in sheet area due to dense tubular matrix
% treat as a packing of circular holes separated by tubule diameter

rh = 150; % radius of holes
rt = 50; % tubule radius
pack = 0.9069; % packing efficiency

% surface area of inner half for each donut
S1 = 2*pi*(rh+rt-2*rt/pi)*pi*rt
Sflat = pi*(rh+rt^2)*2

% true area relative to apparent area
Ascl = pack*(pi*rt)*(rh+rt-2*rt/pi)/(rh+rt)^2 + (1-pack)

%%
rtlist = linspace(20,60);
rhlist = 150;
Ascl = pack*(pi*rtlist).*(rhlist+rtlist-2*rtlist/pi)./(rhlist+rtlist).^2 + (1-pack);

rhlist = 200-rtlist;
Ascl2 = pack*(pi*rtlist).*(rhlist+rtlist-2*rtlist/pi)./(rhlist+rtlist).^2 + (1-pack);
% plot R_estimate / R_true as the tube radius changes
%plot(rtlist,1./Ascl)
% plot R_estimate as the tube radius changes
plot(rtlist,rtlist./Ascl, 'b', rtlist,rtlist./Ascl2,'r', 'LineWidth',2)
set(gca,'FontSize',14,'defaultTextInterpreter','latex','TickLabelInterpreter','latex')
xlabel('actual tube radius')
ylabel('our estimated radius')
legend('$r_{hole} = 150nm$', '$r_{hole}=200nm-r_{tube}$','Interpreter','latex')

%% given the estimated measurements, what are the true ones?
Asclfunc = @(rh,rt) pack*(pi*rt).*(rh+rt-2*rt/pi)./(rh+rt).^2 + (1-pack);

Rest = linspace(40,60);
clear rtrue1 rtrue2
for rc = 1:length(Rest) % cycle over estimated measurements
    rtrue1(rc) = fzero(@(rt) rt/Asclfunc(150,rt) - Rest(rc), Rest(rc));
    rtrue2(rc) = fzero(@(rt) rt/Asclfunc(200-rt,rt) - Rest(rc), Rest(rc));
end

plot(Rest,rtrue1,'b',Rest,rtrue2,'r','LineWidth',2)
set(gca,'FontSize',14,'defaultTextInterpreter','latex','TickLabelInterpreter','latex')
xlabel('our estimated radius')
ylabel('true tube radius')
legend('$r_{hole} = 150nm$', '$r_{hole}=200nm-r_{tube}$','Interpreter','latex')

%% relative measurements
Rmid = 54;
rtrue1mid = interp1(Rest,rtrue1,Rmid);
rtrue2mid = interp1(Rest,rtrue2,Rmid);
plot(Rest/Rmid,rtrue1/rtrue1mid,'b',Rest/Rmid,rtrue2/rtrue2mid,'r','LineWidth',2)
set(gca,'FontSize',14,'defaultTextInterpreter','latex','TickLabelInterpreter','latex')
xlabel('our estimated radius, relative')
ylabel('true tube radius, relative')
legend('$r_{hole} = 150nm$', '$r_{hole}=200nm-r_{tube}$','Interpreter','latex')

% ------
%% now consider the perforated sheet case
rh = 0.13/2; % inner radius of holes
rho = 10; % hole density (per um^2)
rs = 0.05/2; % half sheet thickness

Ascl = 1-rho*(pi*(rh+rs)^2 - pi^2*rs*(rh+rs-2*rs/pi))


%% hexagonal tubule packing
rh = 190;
rt = 10;

alpha = rh/rt;
a = (rh+2*rt)/sqrt(3)
Ascl = (12*pi*(alpha-1)/sqrt(3) + (2/3*(4*sqrt(3)-2*pi) + 4*pi*(pi-1)))/(4*sqrt(3)*(alpha+1)^2)

rt./Ascl
%%

%% hexagonal tubule packing for different true tube radius
rtlist = linspace(30,60);
rh = 150;

alpha = rh/50;
Ascl0 =  (12*pi*(alpha-1)/sqrt(3) + (2/3*6*(4*sqrt(3)-2*pi) + 4*pi*(pi-1)))/(4*sqrt(3)*(alpha+1)^2);

for rc = 1:length(rtlist)
    rt = rtlist(rc);
    alpha = rh./rt;
    Ascl(rc) = (12*pi*(alpha-1)/sqrt(3) + (2/3*6*(4*sqrt(3)-2*pi) + 4*pi*(pi-1)))/(4*sqrt(3)*(alpha+1)^2)
    
    alpha = (200-rt)./rt;
    Ascl2(rc) = (12*pi*(alpha-1)/sqrt(3) + (2/3*6*(4*sqrt(3)-2*pi) + 4*pi*(pi-1)))/(4*sqrt(3)*(alpha+1)^2)            
end

plot(rtlist,rtlist./Ascl0,'k',rtlist,rtlist./Ascl, 'b', rtlist,rtlist./Ascl2,'r', 'LineWidth',2)
set(gca,'FontSize',14,'defaultTextInterpreter','latex','TickLabelInterpreter','latex')
xlabel('actual tube radius')
ylabel('our estimated radius')
legend('$r_t$ constant', '$r_h = 150nm$', '$r_h=200nm-r_t$','Interpreter','latex')

title('hexagonal matrix of tubules')

%% relative values
Rmid = 54;
rest1mid = interp1(rtlist,rtlist./Ascl,Rmid);
rest2mid = interp1(rtlist,rtlist./Ascl2,Rmid);
plot(rtlist/Rmid,rtlist./Ascl/rest1mid,'b',rtlist/Rmid,rtlist./Ascl2/rest2mid,'r','LineWidth',2)
set(gca,'FontSize',14,'defaultTextInterpreter','latex','TickLabelInterpreter','latex')
xlabel('our estimated radius, relative')
ylabel('true tube radius, relative')
legend('$r_{hole} = 150nm$', '$r_{hole}=200nm-r_{tube}$','Interpreter','latex')


%% given the estimated measurements, what are the true ones?
Asclfunc = @(alpha) (12*pi*(alpha-1)/sqrt(3) + (2/3*(4*sqrt(3)-2*pi) + 4*pi*(pi-1)))./(4*sqrt(3)*(alpha+1).^2)

Rest = linspace(40,60);
clear rtrue1 rtrue2
for rc = 1:length(Rest) % cycle over estimated measurements
    rtrue1(rc) = fzero(@(rt) rt/Asclfunc(150./rt) - Rest(rc), Rest(rc));
    rtrue2(rc) = fzero(@(rt) rt/Asclfunc((200-rt)./rt) - Rest(rc), Rest(rc));
end

plot(Rest,rtrue1,'b',Rest,rtrue2,'r','LineWidth',2)
set(gca,'FontSize',14,'defaultTextInterpreter','latex','TickLabelInterpreter','latex')
xlabel('our estimated radius')
ylabel('true tube radius')
legend('$r_{hole} = 150nm$', '$r_{hole}=200nm-r_{tube}$','Interpreter','latex')


%% Compare volumes for different models

rtlist = linspace(30,60);
rh = 150;
t = 4; %membrane thickness
% expected sheet height 
h = 50;
for rc = 1:length(rtlist)
    rt = rtlist(rc);
    
    % for dense tubules
    rh = 200-rt;
    ell = (2*rh-rt)/sqrt(3);
    a = 2*(rh+rt)/sqrt(3);    
    %Vtruematrix(rc) = (3*pi*(rt-t)^2*ell + 4*sqrt(3)*(rt-t)^3)/(2*sqrt(3)*(rh+rt)^2*h);
     
    aj = (4*sqrt(3)-2*pi)*rt^2;
    vj = 2*pi*(rt+(rt-t)-4*(rt-t)/3/pi)*pi/2*(rt-t)^2; % inner part of torus
    Vtruematrix(rc) = (3*pi*(rt-t)^2*ell + aj*(h-2*t)*6/3+vj)/(2*sqrt(3)*(rh+rt)^2*h);
    
    %Atruematrix(rc) = ((2*pi*rt)*3*ell + 2*sqrt(3)*rt^2*2)/(3*sqrt(3)*a^2);
    alpha = rh/rt;
    Atruematrix(rc) = (12*pi*(alpha-1)/sqrt(3) + (2/3*6*(4*sqrt(3)-2*pi) + 4*pi*(pi-1)))/(4*sqrt(3)*(alpha+1)^2)       
    
    ab(rc) = Atruematrix(rc)/Vtruematrix(rc);
    
    Rest(rc) = rtlist(rc)/Atruematrix(rc);
end

%plot(rtlist,Vtruematrix, rtlist,Atruematrix)
plot(rtlist,rtlist/h,rtlist,ab/h.*rtlist,'LineWidth',2)
set(gca,'FontSize',14,'defaultTextInterpreter','latex','TickLabelInterpreter','latex')
xlabel('radius (nm)')
ylabel('$(I_\mathrm{m,s}/I_\mathrm{m,t})/(I_\mathrm{l,s}/I_\mathrm{l,t})$')
legend('sheets','tubular matrices','Interpreter','latex')

%%
plot(rtlist,Rest)
