clc
clear all
close all

% input parameters
L=200;     % channel length (m)
Tmax=15;   % t_max(s)

dx=0.65;   % step size(m) 
x=0:dx:L;
Nx=length(x);

eps=10^(-242);
hu=1;   % initial water depth in  reservoir (m)
hd=eps; % initial water depth in  channel (m)
g=9.81; % gravity acceleration (m/s^2)
n=0;    % manning coefficient
So=0;   % bed slope
K = 0.5;

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
dt = 0.7*Cn*dx/max(abs(u)+sqrt(g*h)); % stability criterion

uh=u.*h;
t=t+dt;

while t<=Tmax
    %% PREDICTOR
    % boundary condition
    hp(1:3)=hu;
    uhp(1:3)=uh(1:3);

    hp(Nx-2:Nx)=hd;
    uhp(Nx-2:Nx)=uh(Nx-2:Nx);
  
    for i=4:Nx-1
        b2(i)= -B2*h(i)^3*u(i)/3*((u(i+1)-2*u(i)+u(i-1))/dx^2);
        b3(i)= B3*h(i)^3/3*((u(i+1)-u(i-1))/(2*dx))^2;
        T(i)=u(i)^2*h(i)+(g*h(i)^2/2)+b2(i)+b3(i);
    end
    
    for i=4:Nx-3
        Sf = u(i)^2*n^2/h(i)^(4/3);
        hp(i)= h(i)+ dt/dx*(uh(i)-uh(i+1));
        if hp(i)>eps
            uhp(i)=uh(i)+ dt/dx*(T(i)-T(i+1))+dt*g*h(i)*(So-Sf);
        else
            uhp(i)=0;
        end
    end
    
    for j=1:Nx
        if hp(j)>eps
            up(j) = uhp(j)/hp(j);
        else
            up(j)=0;
        end
    end
    % boundary condition
    up(1:3)=0;
    up(Nx-2:Nx)=0;

    %% CORRECTOR
    hc=hp;uhc=uhp;
    
    for i=2:Nx-3
        b2(i)= -B2*hp(i)^3*up(i)/3*((up(i+1)-2*up(i)+up(i-1))/dx^2);
        b3(i)= B3*hp(i)^3/3*((up(i+1)-up(i-1))/(2*dx))^2;
        T(i) = up(i)^2*hp(i)+(g*hp(i)^2/2)+b2(i)+b3(i);
    end

    for i=4:Nx-3
        Sf = up(i)^2*n^2/hp(i)^(4/3);
        hc(i)= hp(i)+dt/dx*(up(i-1)*hp(i-1)-up(i)*hp(i));
        if hc(i)>eps
            uhc(i)=uhp(i)+ dt/dx*(T(i-1)-T(i))+dt*g*hp(i)*(So-Sf);
        else
            uhc(i)=0;
        end
    end
    
    for j=1:Nx
        if hc(j)>eps
            uc(j) = uhc(j)/hc(j);
        else
            uc(j)=0;
        end
    end
    % boundary condition
    uc(1:3)=0;
    uc(Nx-2:Nx)=0;

    %% INTERMEDIATE
    hi=hc;uhi=uhc;
        
    for i=4:Nx-3
        hi(i)=(h(i)+hc(i))/2;
        if hi(i)>eps
            uhi(i)=(uh(i)+uhc(i))/2;
        else
            uhi(i)=0;
        end
    end
    
    for i=1:Nx
        if hi(i)>eps
            ui(i)=uhi(i)/hi(i);
        else
            ui(i)=0;
        end
    end
    
    % boundary condition
    ui(1:3)=0;
    ui(Nx-2:Nx)=0;
    
    %% FINAL
    hf=hi;uhf=uhi;
    
    for i=3:Nx-2
        b1(i)=-B1*hi(i)^3/(6*dx*dt)*(ui(i+1)-ui(i-1)-u(i+1)+u(i-1));
    end
    
    for i=4:Nx-3
        if hf(i)>eps
            uhf(i)=uhi(i)-dt*(b1(i+1)-b1(i-1))/(2*dx);
        else
            uhf(i)=0;
        end
    end
    
    for i=1:Nx
        if hf(i)>eps
            uf(i)=uhf(i)/hf(i);
        else
            uf(i)=0;
        end
    end
    
    % boundary condition
    uf(1:3)=0;
    uf(Nx-2:Nx)=0;
    
    % ARTIFICIAL VISCOSITY (see Jameson, et al. 1981)
    ha=hf;uha=uhf;
    v=zeros(1,Nx);
    for i=2:Nx-1
        v(i) = abs(hf(i+1)-2*hf(i)+hf(i-1))/(abs(hf(i+1))+2*abs(hf(i))+abs(hf(i-1)));
    end
    
    for i=2:Nx-1
        vp = K*max(v(i+1), v(i));
        vm = K*max(v(i), v(i-1));
        ha(i) = hf(i) + vp*(hf(i+1)-hf(i)) - vm*(hf(i)-hf(i-1));
        if hf(i)>eps
            uha(i) = uhf(i) + vp*(uhf(i+1)-uhf(i)) - vm*(uhf(i)-uhf(i-1));
        else
            uha(i)=0;
        end
    end
    
    h=ha;
    uh=uha;
    
    for j=1:Nx
        if h(j)>eps
            u(j) = uh(j)/h(j);
        else
            u(j)=0;
        end
    end

    t=t+dt;
end

plot(x,h,'b:','LineWidth',1)
axis([0 L 0 1.1])
xlabel('Distance (m)','fontsize',12)
ylabel('Flow depth (m)','fontsize',12)
legend('Modified Mohapatra-Chaudhry','fontsize',12)