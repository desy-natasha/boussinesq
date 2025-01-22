clc
clear all
close all

% input parameters
L=200;     % channel length (m)
Tmax=15;   % t_max(s)

dx=0.65;   % step size(m) 
x=0:dx:L;
Nx=length(x);

Ld=100; % dam location (m)
hu=1;   % initial water depth in  reservoir (m)
hd=0.5; % initial water depth in  channel (m)
g=9.81; % gravity acceleration (m/s^2)
n=0;    % manning coefficient
So=0;   % bed slope
K=0.5;

% boussinesq terms (adjust the value to apply the term in the equation)
% 1: boussinesq term is on
% 0: boussinesq term is off
B1=1; 
B2=1;
B3=1;
Cn=0.6; % courant number 

% initial condition (t=0)
t=0;
for i=1:Nx
    u(i)=0;
    if i<=floor(Nx/2)
        h(i) = hu;
    else
        h(i) = hd;
    end
end
dt = Cn*dx/max(abs(u)+sqrt(g*h)); % stability criterion

t=t+dt;

while t<=Tmax
    %% PREDICTOR
    % boundary condition
    hp(1:2)=hu;
    up(1:2)=0;

    hp(Nx-1:Nx)=hd;
    up(Nx-1:Nx)=0;
  
    for i=3:Nx-1
        b2(i)= -B2*h(i)^3*u(i)/3*((u(i+1)-2*u(i)+u(i-1))/dx^2);
        b3(i)= B3*h(i)^3/3*((u(i+1)-u(i-1))/(2*dx))^2;
        T(i)=u(i)^2*h(i)+(g*h(i)^2/2)+b2(i)+b3(i);
    end
    
    for i=3:Nx-2
        Sf = u(i)^2*n^2/h(i)^(4/3);
        hp(i)= h(i)+ dt/dx*(u(i)*h(i)-u(i+1)*h(i+1));
        up(i)= (1/hp(i))*(u(i)*h(i)+ dt/dx*(T(i)-T(i+1))+dt*g*h(i)*(So-Sf));
    end

    %% CORRECTOR
    hc=hp;uc=up;
    
    for i=2:Nx-3
        b2(i)= -B2*hp(i)^3*up(i)/3*((up(i+1)-2*up(i)+up(i-1))/dx^2);
        b3(i)= B3*hp(i)^3/3*((up(i+1)-up(i-1))/(2*dx))^2;
        T(i) = up(i)^2*hp(i)+(g*hp(i)^2/2)+b2(i)+b3(i);
    end

    for i=3:Nx-2
        Sf = up(i)^2*n^2/hp(i)^(4/3);
        hc(i)= hp(i)+dt/dx*(up(i-1)*hp(i-1)-up(i)*hp(i));
        uc(i)= (1/hc(i))*(up(i)*hp(i)+ dt/dx*(T(i-1)-T(i))+dt*g*hp(i)*(So-Sf));
    end
    
    %% INTERMEDIATE
    hi=hc;ui=uc;
        
    for i=3:Nx-2
        hi(i)=(h(i)+hc(i))/2;
        ui(i)=(u(i)+uc(i))/2;
    end
    
    %% FINAL
    hf=hi;uf=ui;
    
    for i=2:Nx-1
        b1(i)=-B1*hi(i)^3/(6*dx*dt)*(ui(i+1)-ui(i-1)-u(i+1)+u(i-1));
    end
    
    for i=3:Nx-2
        uf(i)=ui(i)-dt/hi(i)*(b1(i+1)-b1(i-1))/(2*dx);
    end
    
    % boundary condition
    hf(1:2)=hu;
    uf(1:2)=0;

    hf(Nx-1:Nx)=hd;
    uf(Nx-1:Nx)=0;
    
    % ARTIFICIAL VISCOSITY (see Jameson, et al. 1981)
    v=zeros(1,Nx);
    ha=hf;ua=uf;
    for i=2:Nx-1
        v(i) = abs(hf(i+1)-2*hf(i)+hf(i-1))/(abs(hf(i+1))+2*abs(hf(i))+abs(hf(i-1)));
    end
    
    for i=3:Nx-2
        vp = K*max(v(i+1), v(i));
        vm = K*max(v(i), v(i-1));
        ha(i) = hf(i) + vp*(hf(i+1)-hf(i)) - vm*(hf(i)-hf(i-1));
        ua(i) = uf(i) + (vp*(uf(i+1)-uf(i)) - vm*(uf(i)-uf(i-1)));
    end
    
    h=ha;
    u=ua;
    
    dt = Cn*dx/max(abs(u)+sqrt(g*h)); 
    t=t+dt;
end

plot(x,h,'b:','LineWidth',1)
axis([0 200 0.4 1.1])
xlabel('Distance (m)','fontsize',12)
ylabel('Flow depth (m)','fontsize',12)
legend('Modified Mohapatra-Chaudhry','fontsize',12)