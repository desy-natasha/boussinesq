clc
clear all
close all

% input parameters
L=200;     % channel length (m)
T=15;      % observation time(s)

dx=0.65;
x=0:dx:L;
Nx=length(x);

Ld=100; % dam location (m)
hu=1;   % initial water depth in  reservoir (m)
hd=0;   % initial water depth in  channel (m)

g=9.81;
cu=sqrt(g*hu);
cd=sqrt(g*hd);

% analytical solution for h(x)
for i=1:Nx
    if x(i)<=(Ld-T*cu)
        h(i)=hu;
    elseif x(i)>=(Ld-T*cu) && x(i)<=(Ld+2*T*cu)
        h(i)=4/(9*g)*(cu-((x(i)-Ld)/(2*T)))^2;
    elseif x(i)>=(Ld+2*T*cu)
        h(i)=hd;
    end
end

plot(x,h,'k','LineWidth',1);
axis([0 200 0 1.1])
xlabel('Distance (m)','fontsize',12)
ylabel('Flow depth (m)','fontsize',12)
legend('Analytical Solution','fontsize',12)
