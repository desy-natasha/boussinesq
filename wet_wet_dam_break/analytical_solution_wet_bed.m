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
h1=1;   % initial water depth in  reservoir (m)
h0=0.5; % initial water depth in  channel (m)

g=9.81;
c1=sqrt(g*h1);
c0=sqrt(g*h0);

% calculate e
syms e
eqn= e==2*c1+(g*h0/(4*e)*(1+sqrt(1+8*e^2/(g*h0))))-(2*g*h0*(sqrt(1+8*e^2/(g*h0))-1))^(0.5);
e_array = double(solve(eqn,e));
for i=1:size(e_array,1);
    if e_array(i)>c0 && e_array(i)<c1
        edot=e_array(i);
    end
end

% calculate u2 dan h2
h2 = h0*0.5*(sqrt(1+(8*edot^2/(g*h0)))-1);
u2 = edot - (g*h0/(4*edot)*(1+sqrt(1+(8*edot^2/(g*h0)))));

% analytical solution for h(x)
for i=1:Nx
    if x(i)<=(Ld-T*c1)
        h(i)=h1;
    elseif x(i)>(Ld-T*c1) && x(i)<=(Ld+T*(u2-sqrt(g*h2)))
        h(i)=4/(9*g)*(c1-((x(i)-Ld)/(2*T)))^2;
    elseif x(i)>(Ld+T*(u2-sqrt(g*h2)))&& x(i)<(Ld+T*edot)
        h(i)=h2;
    elseif x(i)>=(Ld+T*edot)
        h(i)=h0;
    end
end

plot(x,h,'k','LineWidth',1)
axis([0 200 0.4 1.1])
xlabel('Distance (m)','fontsize',12)
ylabel('Flow depth (m)','fontsize',12)
legend('Analytical Solution','fontsize',12)
