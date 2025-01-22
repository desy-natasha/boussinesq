clc
close all
clear all

% boussinesq terms
b1 = 1;  
b2 = 1;  
b3 = 1;  

% input parameters
L      = 77892; % meter
Tfinal = 30.7*60*60; % second 

g = 9.81;      
K = 0.5;  
Cn= 0.6;      
S0=0;

dx = 50;
x  = 0:dx:L;        % discretized space
Nx = length(x);     % number of space points

dt = 0.5;
Nt=floor(Tfinal/dt)+1;

% initial condition
d0=2.2;
for i=1:Nx
    u(i) = 0;
    e(i) = 0;
    h(i) = e(i)+d0;
end

t=0;
uh=u.*h;

for n=1:Nt
    % boundary condition 
    e(1:3)=2.5*cos((n*dt)*((2*pi)/44712)+3.37703);
    h(1:3)=e(1:3)+d0;
   
    %% PREDICTOR
    for i=4:Nx-3
        hp(i)=h(i)+ dt/(6*dx) *(uh(i+2) -8*uh(i+1) +7*uh(i));
    end
    
    for i=2:Nx-1
        Tp(i) = u(i)^2*h(i) + 0.5*g*h(i)^2 ...
            -b2*(1/3*u(i)*(h(i)^3) * (u(i+1)-2*u(i)+u(i-1))/dx^2) ...
            +b3*(1/3*h(i)^3 * ((u(i+1)-u(i-1))/(2*dx))^2);
    end
    
    for i=4:Nx-3
        uhp(i)=uh(i)+ dt/(6*dx) *(Tp(i+2) -8*Tp(i+1) +7*Tp(i))-dt*g*h(i)*S0;
    end
    
    % update left boundary
    hp(1:3)=e(1:3)+d0;
    uhp(1:3)=uh(1:3);
    % update right boundary
    hp(Nx-2:Nx)=e(Nx-2:Nx)+d0;
    uhp(Nx-2:Nx)=uh(Nx-2:Nx);
    
    for j=1:Nx
        if hp(j)>d0
            up(j) = uhp(j)/hp(j);
        else
            up(j)=0;
        end
    end
    
    ep = hp - d0;
    
    %% CORRECTOR
    hc=hp; uhc=uhp;
    for i=4:Nx-3
        hc(i)=hp(i)+ dt/(6*dx) *(-uhp(i-2) +8*uhp(i-1) -7*uhp(i));
    end
    
    for i=2:Nx-1
        Tc(i) = up(i)^2*hp(i) + 0.5*g*hp(i)^2 ...
            -b2*(1/3*up(i)*(hp(i)^3) * (up(i+1)-2*up(i)+up(i-1))/dx^2) ...
            +b3*(1/3*hp(i)^3 * ((up(i+1)-up(i-1))/(2*dx))^2);
    end
    
    for i=4:Nx-3
        if hc(i)>d0
           uhc(i)=uhp(i)+ dt/(6*dx) *(-Tc(i-2) +8*Tc(i-1) -7*Tc(i))-dt*g*hp(i)*S0;
        else
           uhc(i)=0;
        end
    end
    
    % update left boundary
    hc(1:3)=e(1:3)+d0;
    uhc(1:3)=uh(1:3);
    % update right boundary
    hc(Nx-2:Nx)=e(Nx-2:Nx)+d0;
    uhc(Nx-2:Nx)=uh(Nx-2:Nx);
    
    for j=1:Nx
        if hc(j)>d0
            uc(j) = uhc(j)/hc(j);
        else
            up(j)=0;
        end
    end
    
    ec = hc - d0;
    
    %% INTERMEDIATE
    hi= 0.5*(h + hc);
    for i=4:Nx-3
        if hi(i)>d0
            uhi(i)= 0.5*(uh(i) + uhc(i));
        else
            uhi(i)=0;
        end
    end
    
    % update left boundary
    hi(1:3)=e(1:3)+d0;
    uhi(1:3)=uh(1:3);
    % update right boundary
    hi(Nx-2:Nx)=e(Nx-2:Nx)+d0;%10^(-9);
    uhi(Nx-2:Nx)=uh(Nx-2:Nx);
    
    for j=1:Nx
        if hi(j)>d0
            ui(j) = uhi(j)/hi(j);
        else
            ui(j)=0;
        end
    end
    
    ei = hi - d0;
    
    %% FINAL
    hf = hi;   
    
    for i=2:Nx-1
        B1(i) = b1*-1*hi(i)^3/(6*dx*dt)*(ui(i+1)-ui(i-1)-u(i+1)+u(i-1)); 
    end
    % update left boundary for B1
    B1(1:3)=0;
    % update right boundary for B1
    B1(Nx-2:Nx)=0;
    
    for i=4:Nx-3
        if hf(i)>d0
            uhf(i)= uhi(i) - dt/(2*dx)*(B1(i+1)-B1(i-1));
        else
            uhf(i)=0;
        end
    end
    
    % update left boundary
    hf(1:3)=e(1:3)+d0;
    uhf(1:3)=uh(1:3);
    % update right boundary
    hf(Nx-2:Nx)=e(Nx-2:Nx)+d0;
    uhf(Nx-2:Nx)=uh(Nx-2:Nx);
    
    for j=1:Nx
        if hf(j)>d0
            uf(j) = uhf(j)/hf(j);
        else
            uf(j)=0;
        end
    end
    
    ef = hf - d0;
    
    % ARTIFICIAL VISCOSITY (see Jameson, et al. 1981)
    v=zeros(1,Nx);
    for i=2:Nx-1
        v(i) = abs(hf(i+1)-2*hf(i)+hf(i-1))/(abs(hf(i+1))+2*abs(hf(i))+abs(hf(i-1)));
    end
    
    for i=2:Nx-1
        vp = K*max(v(i+1), v(i));
        vm = K*max(v(i), v(i-1));
        ha(i) = hf(i) + vp*(hf(i+1)-hf(i)) - vm*(hf(i)-hf(i-1));
        if hf(i)>d0
            uha(i) = uhf(i) + vp*(uhf(i+1)-uhf(i)) - vm*(uhf(i)-uhf(i-1));
        else
            uha(i)=0;
        end
    end
    
    % update left boundary
    ha(1:3)=e(1:3)+d0;
    uha(1:3)=uh(1:3);
    % update right boundary
    ha(Nx-2:Nx)=e(Nx-2:Nx)+d0;
    uha(Nx-2:Nx)=uh(Nx-2:Nx);
    
    h = ha;
    uh= uha;
    
    for j=1:Nx
        if h(j)>d0
            u(j) = uh(j)/h(j);
        else
            u(j)=0;
        end
    end
    
    e = h-d0;
    
    for k=1:Nx
        HH(k,n) = h(k);
    end

    t = t+dt;
end

t1 = 0:dt:Tfinal;

% Measurement data (Bayu, 2019)
t2=xlsread('data_sungai_kampar.xlsx','A1:A101');
E2=xlsread('data_sungai_kampar.xlsx','B1:B101');

plot((t1/3600)-5.7,HH(round(60000/dx),:),'r-','LineWidth',1)
hold on
plot(t2,E2,'k-','LineWidth',1)
axis([0 25 0 5])
xlabel('Hours','fontsize',12)
ylabel('Elevation (m)','fontsize',12)
legend('Mohapatra-Chaudhry','Measurement Result','fontsize',12)
