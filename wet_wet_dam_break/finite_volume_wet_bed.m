clc
close all
clear all 

% input parameters
L = 200;       % space domain
Tfinal = 15;   % final time
g = 9.81;      % gravity
d = 0;

hu= 1;
hd= 0.5;

dx= 0.65;
x = 0:dx:L;       % discretized space
Nx= length(x);    % number of space points

% x1/2 x1 x3/2 x2 ..... xj-1/2 xj xj+1/2 ..... xNx xNx+1/2
% u1/2 eta1 u3/2 eta2 .... uj-1/2 etaj uj+1/2 ... etaNx uNx+1/2
% grid matlab:
% u1 eta1 u2 eta2... uj etaj uj+1 .. etaNx uNx+1
% size for u[Nx+1,Nt] and eta[Nx,Nt]

% initial condition
for i =1:Nx 
    u(i,1)=0;
    
    if x(i)<=100
        e(i,1) = hu;
    else
        e(i,1) = hd;
    end
    
    h(i,1)=e(i,1)+d;
end

Cn= 0.6;
dt= dx/sqrt(g*max(h));
t = 0:dt:Tfinal;  % discretized time
Nt= length(t);    % number of space points

% boundary condition
u(Nx+1,1)=0;
u(1,1:Nt)=0;
u(Nx+1,1:Nt)=0;

eps=0.15;

for j=2:Nt
    
    q(1)=u(1,j-1)*h(1,j-1);
    for i=2:Nx
        if u(i,j-1)>=0
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
    
    % filter (e_t = eps*e_xx)
    for i=2:Nx-1
        e(i,j)=e(i,j)+dt*(eps*(e(i+1,j)-2*e(i,j)+e(i-1,j)))/(dx^2);
    end
    
    for i=2:Nx
        % if boussinesq term is on
        F(i,j-1)=u(i,j-1)-(h(i,j)^2*(u(i+1,j-1)-2*u(i,j-1)+u(i-1,j-1))/(3*dx^2));
        B2(i,j-1)=-1*(u(i,j-1)*(h(i,j)^3)*(u(i+1,j-1)-2*u(i,j-1)+u(i-1,j-1))/(3*dx^2));
        B3(i,j-1)=1*(h(i,j)^3 *((u(i+1,j-1)-u(i-1,j-1))/(2*dx))^2)/3;
        
%         % if boussinesq term is off
%         F(i,j-1)=u(i,j-1); %b1 off
%         B2(i,j-1)=0;
%         B3(i,j-1)=0;
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
        Fn(i)=F(i+1,j);% reshape F into 1 x Nx-1
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
   
    % filter (u_t = eps*u_xx)
    for i=2:Nx
        u(i,j)=u(i,j)+dt*(eps*(u(i+1,j)-2*u(i,j)+u(i-1,j)))/(dx^2);
    end
    
    for i=1:Nx
        H(i) = h(i,j);
    end

end

plot(x,H,'k-','LineWidth',1)
axis([0 200 0.4 1.1])
xlabel('Distance (m)','fontsize',12)
ylabel('Flow depth (m)','fontsize',12)
legend('Finite Volume Model','fontsize',12)
