clc
close all
clear all

% boussinesq terms
b1 = 1;
b2 = 1;  
b3 = 1; 

% input parameters
L = 250;      % space domain
time = [0;3.19;6.39;9.58;12.77;15.96];
S = 1;
Tfinal = time(S);
g = 9.81;      
K = 0.5;        

dx = 0.65;
x  = 0:dx:L;        
Nx = length(x);     

dt = 0.01;
Nt = floor(Tfinal/dt)+1; 

d0 = 1;
S0 = 0;

a = 5;
uo= 0.1;

% initial condition
for n=1:Nx
    u(n)=0.5*uo*sqrt(g*d0)*(1-tanh(((n-1)*dx-15)/a));
    e(n)=(u(n)*sqrt(g*d0)+0.25*u(n)^2)/g;
    h(n)=e(n)+d0;
end

t  = 0;
uh = u.*h;

for n=1:Nt
    %% PREDICTOR
    
    % left boundary
    hp(1:3)=e(1:3)+d0;
    uhp(1:3)=uh(1:3);
    % right boundary
    hp(Nx-2:Nx)=e(Nx-2:Nx)+d0;
    uhp(Nx-2:Nx)=uh(Nx-2:Nx);
    
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
        if hp(j)>0
            up(j) = uhp(j)/hp(j);
        else
            up(j)=0;
        end
    end
    
    ep = hp - d0;
    
    %% CORRECTOR
    hc=hp;uhc=uhp;
    for i=4:Nx-3
        hc(i)=hp(i)+ dt/(6*dx) *(-uhp(i-2) +8*uhp(i-1) -7*uhp(i));
    end
    
    for i=2:Nx-1
        Tc(i) = up(i)^2*hp(i) + 0.5*g*hp(i)^2 ...
            -b2*(1/3*up(i)*(hp(i)^3) * (up(i+1)-2*up(i)+up(i-1))/dx^2) ...
            +b3*(1/3*hp(i)^3 * ((up(i+1)-up(i-1))/(2*dx))^2);
    end
    
    for i=4:Nx-3    
       uhc(i)=uhp(i)+ dt/(6*dx) *(-Tc(i-2) +8*Tc(i-1) -7*Tc(i))-dt*g*hp(i)*S0;
    end
    
    % update left boundary
    hc(1:3)=e(1:3)+d0;
    uhc(1:3)=uh(1:3);
    % update right boundary
    hc(Nx-2:Nx)=e(Nx-2:Nx)+d0;
    uhc(Nx-2:Nx)=uh(Nx-2:Nx);
    
    for j=1:Nx
        if hc(j)>0
            uc(j) = uhc(j)/hc(j);
        else
            up(j)=0;
        end
    end
    
    ec = hc - d0;
    
    %% INTERMEDIATE
    hi= 0.5*(h + hc);
    for i=4:Nx-3
        uhi(i)= 0.5*(uh(i) + uhc(i));
    end
    
    % update left boundary
    hi(1:3)=e(1:3)+d0;
    uhi(1:3)=uh(1:3);
    % update right boundary
    hi(Nx-2:Nx)=e(Nx-2:Nx)+d0;%10^(-9);
    uhi(Nx-2:Nx)=uh(Nx-2:Nx);
    
    for j=1:Nx
        if hi(j)>0
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
        uhf(i)= uhi(i) - dt/(2*dx)*(B1(i+1)-B1(i-1));
    end
    
    % update left boundary
    hf(1:3)=e(1:3)+d0;
    uhf(1:3)=uh(1:3);
    % update right boundary
    hf(Nx-2:Nx)=e(Nx-2:Nx)+d0;%10^(-9);
    uhf(Nx-2:Nx)=uh(Nx-2:Nx);
    
    for j=1:Nx
        if hf(j)>0
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
        uha(i) = uhf(i) + vp*(uhf(i+1)-uhf(i)) - vm*(uhf(i)-uhf(i-1));
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
        if h(j)>0
            u(j) = uh(j)/h(j);
        else
            u(j)=0;
        end
    end
    
    e = h-d0;

    t = t+dt;
end

% Experimental data (Frazao, 2002)
x2=xlsread('Data_UndularBore_Frazao_New.xlsx',S,'A1:A695');
E2=xlsread('Data_UndularBore_Frazao_New.xlsx',S,'B1:B695');

figure(1)
plot(x,h,'r--','LineWidth',1.5);
hold on
str=['t= ',num2str(time(S)),' s'];
annotation('textbox',[0.15,0.12,0.1,0.1],'String',str);
plot(x2,E2(:),'k-','LineWidth',1.5);
grid on
legend('Mohapatra-Chaudhry','MUSCL4-Frazao');
xlabel('x (m)');
ylabel('h (m)');
axis([(S-1)*10 (S+2)*10 1 1.15]);
