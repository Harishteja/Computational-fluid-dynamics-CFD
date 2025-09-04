close all
clear all
clc
L=1;
H=0.5;
error=10;
errormax=0.001;
Nx=80;
Ny=80;
S.F=input('enter strech_factor -0, for normal mesh or enter any number for refined mesh:');
if S.F ~= 0
     L1=0.8;
    L2=0.2;
    Nx=96;
    Ny=96;
    Nx1=64;
    Nx2=32;

    for i=1:Nx1
        dx(i)=L1/Nx1;
    end
    for i=Nx1+1:Nx
        dx(i)=L2/Nx2;
    end
    for i=1:Ny
        dy(i)=H/Ny;
    end
else
    
    for i =1:Nx
        dx(i)=L/Nx;
    end
    for i=1:Ny
        dy(i)=H/Ny;
    end
end

GLx(1)=0;
for i=2:Nx+1
    GLx(i)=GLx(i-1)+dx(i-1);
end 

GLy(1)=0;            
for i=2:Ny+1
    GLy(i)=GLy(i-1)+dy(i-1);
end  

nodex(1)=0;
nodex(Nx+2)=GLx(Nx+1);
for i=2:Nx+1
    nodex(i)=(GLx(i-1)+GLx(i))/2; 
end  

nodey(1)=0;
nodey(Ny+2)=GLy(Ny+1);
for i=2:Ny+1
    nodey(i)=(GLy(i-1)+GLy(i))/2; 
end

for i=2:Nx+1
    A(i) = dx(i-1)*dy(i-1)
end
for i=1:Nx+1
    k(i)=(16*(2-(GLy(i)/H)));
end


for i=1:Nx+2
    k2(i)=(16*(2-(nodey(i)/H)));
end


for i=2:Nx+1
    DXe(i)=nodex(i+1)-nodex(i);
    DXw(i)=nodex(i)-nodex(i-1);
end

for i=2:Ny+1
    DYn(i)=nodey(i+1)-nodey(i);   
    DYs(i)=nodey(i)-nodey(i-1);
end

T = ones(Nx+2,Ny+2);



T(Nx+2,:)=15;
for i=1:Nx+2
    T(i,Ny+2)= 5*(nodey(i)/H) + 15*sin(pi*(H-nodey(i))/H);
end
condition=input('give 1 for boundary condition to be as in given problem  or give any other number for changing Neumann to Dirichlet:');
if condition==1
    T(1,:)=10;
    T(:,1)=T(:,2);
else
    T(1,:)=T(2,:)
    T(:,1)=10
end
itr = 0;
while error>errormax
    T_old=T;
    residual=0;
    for i=2:Nx+1
        for j=2:Ny+1
            Ae(i)=k(i)*A(i)/DXe(i);
            Aw(i)=k(i)*A(i)/DXw(i);
            An(i)=k(i)*A(i)/DYn(i);
            As(i)=k(i)*A(i)/DYs(i);
            Ap(i)=Ae(i)+ Aw(i)+An(i)+As(i)+(1.5*dx(i-1)*dy(i-1))/T_old(i,j);
            
            if condition~=1
                T(1,:)=T(2,:);
            else
                T(:,1)=T(:,2);
            end
            T(i,j)=(Ae(i)*T(i,j+1)+Aw(i)*T(i,j-1)+An(i)*T(i-1,j)+As(i)*T(i+1,j))/Ap(i);
        end
    end
      
    disp(T)

    temp=flip(T);
    error=0;
   for i=2:Nx+1
      for j=2:Ny+1
        error = error+abs(T_old(i,j)-T(i,j));
      end
   end 
   itr=itr+1;
   err3(itr)=error;
end
figure('Name','MESH GRID')
[x,y]=meshgrid(GLx,GLy);
plot(x,GLy)
hold on
plot(GLx,y')
set(gca,'XTick',0:0.1:1)
[x,y]=meshgrid(nodex,nodey);
scatter(x(:),y(:),'.');
figure('Name','MESH INDEPENDENCE')
load nodex8.mat;
load T_8.mat;
plot(nodex8,T_8);
hold on
load nodex16.mat;
load T_N16.mat;
plot(nodex16,T_N16);
axis([0 1 0 15])
legend('80*80','160 *160')


gradientx=zeros(Nx+2);	
gradienty=zeros(Ny+2);

for i=2:Nx+1
    for j=2:Ny+1    
        gradientx(i,j)=(temp(i,j+1)-temp(i,j-1))/(nodex(i+1)-nodex(i-1));
        gradienty(i,j)=(temp(i+1,j)-temp(i-1,j))/(nodey(i+1)-nodey(i));
    end
end
for j=1:Nx+1
        i=1;
        gradientx(i,j)=(temp(i,+1)-temp(i,j))/(nodex(j+1)-nodex(j));
end
for j=1:Nx+2
        i=1;
        gradienty(i,j)=(temp(i+1,j)-temp(i,j))/(nodey(i+1)-nodey(i));
end
for i=1:Nx+2
        j=Nx+2;
        gradientx(i,j)=(temp(i,j)-temp(i,j-1))/(nodex(j)-nodex(j-1));
end
for j=1:Nx+2
        i=Nx+2;
        gradienty(i,j)=(temp(i,j)-temp(i-1,j))/(nodey(i)-nodey(i-1));
end
for i=2:Ny+1
        j=1;
        gradienty(i,j)=(temp(i+1,j)-temp(i-1,j))/(nodey(i+1)-nodey(i));
end
for i=2:Ny+1
        j=Nx+2;
        gradienty(i,j)=(temp(i+1,j)-temp(i-1,j))/(nodey(i+1)-nodey(i));
end
fluxx=zeros(Nx+2);
fluxy=zeros(Ny+2);
for i=1:Nx+2
    for j=1:Ny+2    
        fluxx(i,j)=(-k2(i))*gradientx(i,j);
        fluxy(i,j)=(-k2(i))*gradienty(i,j);
    end
end

figure('Name','HEAT FLUX VECTOR')
quiver(nodex,nodey,fluxx,fluxy,20);
xlim([0 1]);
ylim([0 0.5]);
figure('Name','TEMPERATURE CONTOUR')
contourf(nodex,nodey,temp)

iteration3 = linspace(1,itr,itr);
figure('Name','ERROR vs ITERATIONS')
plot(iteration3,err3,'o')
hold on
plot(iteration5,err5,'-')
legend('errorE-3','errorE-5')