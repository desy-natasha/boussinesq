clc
close all
clear all 

L  = 250;
time = [0;3.19;6.39;9.58;12.77;15.96];
S  = 1;
T  = time(S);
g  = 9.81;

dx = 0.1;
x  = 0:dx:L;
Nx = length(x);

dt = 0.01;
Nt = floor(T/dt)+1;

d  = 1;
a  = 5;
uo = 0.1;

% input parameters
for n=1:Nx
    u(n,1)=0.5*uo*sqrt(g*d)*(1-tanh(((n-1)*dx-15)/a));
    e(n,1)=(u(n,1)*sqrt(g*d)+0.25*u(n,1)^2)/g;
    h(n,1)=e(n,1)+d;
end
h(Nx+1,1)=h(Nx,1);
u(Nx+1,1)=0;

% x1/2 x1 x3/2 x2 ..... xj-1/2 xj xj+1/2 ..... xNx xNx+1/2
% u1/2 eta1 u3/2 eta2 .... uj-1/2 etaj uj+1/2 ... etaNx uNx+1/2
% grid matlab:
% u1 eta1 u2 eta2... uj etaj uj+1 .. etaNx uNx+1
% size from u[Nx+1,Nt] and eta[Nx,Nt]

for j=2:Nt+1
    q(1)=u(1,j-1)*h(1,j-1);
    for i=2:Nx
        if (u(i,j-1)>0)
            h1(i,j-1)=h(i-1,j-1);
        else
            h1(i,j-1)=h(i,j-1);
        end
        q(i) = u(i,j-1)*h1(i,j-1);
    end
    q(Nx+1)=u(Nx+1,j-1)*h1(Nx,j-1);
    
    for i=1:Nx
        e(i,j)=e(i,j-1)+dt/dx*(q(i)-q(i+1));
        h(i,j)=e(i,j)+d;
    end
    
    u(1,j)=u(1,j-1);
    for i=2:Nx
        F(i,j-1)=u(i,j-1)-(h(i,j)^2*(u(i+1,j-1)-2*u(i,j-1)+u(i-1,j-1))/(3*dx^2));
        B2(i,j-1)=-1*(u(i,j-1)*(h(i,j)^3)*(u(i+1,j-1)-2*u(i,j-1)+u(i-1,j-1))/(3*dx^2));
        B3(i,j-1)=1*(h(i,j)^3 *((u(i+1,j-1)-u(i-1,j-1))/(2*dx))^2)/3;
    end

    for i=2:Nx
        q_bar(i)=(q(i)+q(i-1))/2;
        h_bar(i)=(h(i,j)+h(i-1,j))/2;
    end
    
    for i=2:Nx-1
        if h(i,j)>0||h(i-1,j)>0
            % calculate uu_x
            if q_bar(i)>=0
                uux(i,j-1)=q_bar(i)*(u(i,j-1)-u(i-1,j-1))/(h_bar(i)*dx);
            else
                uux(i,j-1)=q_bar(i+1)*(u(i+1,j-1)-u(i,j-1))/(h_bar(i)*dx);
            end
        else
            uux(i,j-1)=0;
        end
    end
    
    % calculate F_n+1
    for i=2:Nx-1
        if h_bar(i)>0
            F(i,j)=F(i,j-1)-dt*((g*(e(i,j)-e(i-1,j))/dx)+uux(i,j-1)+...
                ((B2(i,j-1)+B3(i,j-1)-B2(i-1,j-1)-B3(i-1,j-1))/(dx*h_bar(i))));
        else
            F(i,j)=0;
        end
    end
    F(2,j)=F(2,j)+h(1,j)^2/(3*dx^2)*u(1,j-1);
    F(Nx,j)=F(Nx,j)+h(Nx,j)^2/(3*dx^2)*u(Nx+1,j-1); 
    
    for i=1:Nx-1
        a(i)= -h(i,j)^2/(3*dx^2);      %A(i,i-1)
        b(i)= 1+(2*h(i,j)^2)/(3*dx^2); %A(i,i)
        c(i)= a(i);                    %A(i,i+i)
    end
    
    n=Nx-1;
    Fn = zeros(1,n);
    for i=1:n
        Fn(i)=F(i+1,j); % reshape F into 1 x Nx-1
    end
    
    for k=1:Nx-2
        if b(k)>0
            p=a(k+1)/b(k);
        else
            p=0;
        end
        b(k+1)=b(k+1)-p*c(k);
        Fn(k+1)=Fn(k+1)-p*Fn(k);
        a(k+1)=0;
    end
    
    ubaru(n)=Fn(n)/b(n);
    for k=Nx-2:-1:1
        ubaru(k)=(Fn(k)-c(k)*ubaru(k+1))/b(k);
    end
    
    % update u value
    u(1,j)=u(1,j-1);
    u(2:Nx,j)=ubaru(1:Nx-1);
    u(Nx+1,j)=u(Nx+1,j-1);
   
    for i=1:Nx
        H(i) = h(i,j);
        E(i) = e(i,j);
    end

end

% Experimental data (Frazao, 2002)
x2=xlsread('Data_UndularBore_Frazao_New.xlsx',S,'A1:A695');
E2=xlsread('Data_UndularBore_Frazao_New.xlsx',S,'B1:B695');

figure(1)
plot(x,H,'m-.','LineWidth',1.5);
hold on
str=['t= ',num2str(time(S)),' s'];
annotation('textbox',[0.15,0.12,0.1,0.1],'String',str);
plot(x2,E2(:),'k-','LineWidth',1.5);
grid on
legend('Finite Volume Scheme','MUSCL4-Frazao');
xlabel('x (m)');
ylabel('h (m)');
axis([(S-1)*10 (S+2)*10 1 1.15]);
pause(0.001);