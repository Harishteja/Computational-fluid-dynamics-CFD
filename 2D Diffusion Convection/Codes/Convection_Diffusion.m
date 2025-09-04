
close all
clear
clc
% read xc
load xc.dat
nim1=length(xc);
% nim1 = ni-1 = number of grid lines. Number of cell nodes = ni
% For a 10x10 mesh, ni=no of nodes=12,nim1=no of grid lines=11
ni=nim1+1;
%
% read yc
load yc.dat
njm1=length(yc);
nj=njm1+1;
%
% read u
load u.dat
u2d=reshape(u,ni,nj);
% read v
load v.dat
v2d=reshape(v,ni,nj);
%
% xc and yc are the coordinates of the grid
%
%       o---------------o  yc(j)
%       |               |
%       |               |
%       |       P       |
%       |               |
%       |               |
%       o---------------o  yc(j-1)
%
%     xc(i-1)          xc(i)          
%
% The cell-node (P) has the coordinates xp(i), yp(j)
%
% compute the x-coordinates of the cell centres
for i=2:nim1
   xp(i)=0.5*(xc(i)+xc(i-1));
end
xp(1)=xc(1);
xp(ni)=xc(nim1);
%
% take the transpose of x
xp=xp';

% compute the y-coordinates of the cell centres
for j=2:njm1
   yp(j)=0.5*(yc(j)+yc(j-1));
end
yp(1)=yc(1);
yp(nj)=yc(njm1);
%
% take the transpose of y
yp=yp';

% plot the velocity field
% length of vectors = vec
vec= 5;
figure
quiver(xp,yp,u2d,v2d,vec)
axis('equal');
xlabel('x'); ylabel('y')
%print vectxy.ps -deps

%To plot the mesh
figure
[x,y]= meshgrid(xc,yc);
plot(x,yc);
hold on;
plot(xc,y');
xlabel('x'); ylabel('y')

%To compute the delta_xe:dist b/w node 
for i=1:nim1
d_x(i)= xp(i+1)-xp(i);
end

for i=1:njm1
d_y(i)= yp(i+1)-yp(i);
end
%To compute the cell distance
for i=1:nim1-1
delx(i)= xc(i+1)-xc(i);
end

for i=1:njm1-1
dely(i)= yc(i+1)-yc(i);
end
%To compute Peclet number, Pe
rho=1;gamma=1/50;
for i=2:nim1
    Pe=rho*u2d(i)*d_x(i)/gamma;
end
for i=1:nim1-1
fx(i)= 0.5*delx(i)/d_x(i);
end
for i=1:njm1-1
fy(i)= 0.5*dely(i)/d_y(i);
end
for i=1:nim1-1
fxe(i)= 0.5*delx(i)/d_x(i+1);
fxw(i)= 0.5*delx(i)/d_x(i);
end
for j=1:njm1-1
fxs(j)= 0.5*dely(j)/d_y(j);
fxn(j)= 0.5*dely(j)/d_y(j+1);
end

% u3d=ones(nim1-1,njm1-1);
u3d(1:nim1-1,1)=u2d(2:ni-1,1);
u3d(1:nim1-1,njm1)=u2d(2:ni-1,nj);
for i=1:nim1-1
    for j=2:njm1-1                                                                                                                                                                                                                                                                                                                                                                                                       
        u3d(i,j) = fx(j)*u2d(i+1,j+1)+(1-fx(j))*u2d(i+1,j);
    end
end
% v3d=ones(nim1-1,njm1-1);
v3d(1,1:njm1-1)=v2d(1,2:nj-1);
v3d(nim1,1:njm1-1)=v2d(ni,2:nj-1);
for i=2:nim1-1
    for j=1:njm1-1
        v3d(i,j) = fy(j)*v2d(i+1,j+1)+(1-fy(j))*v2d(i,j+1);
    end
end
H=yc(nim1)-yc(1);
h=0.068*H;
% finding no of nodes in the boundary for boundary condition
L=0;
num=0;
for i=1:nim1-1
    L=L+dely(i);
    if L<h 
        num=num+1;
    end
end
% Versteeg hybrid scheme
for i=1:nim1-1
    for j=1:njm1-1
        Fe(i,j)=rho*u3d(i,j+1);
        Fw(i,j)=rho*u3d(i,j);
        Dw(i,j)=gamma/d_x(j);  
        De(i,j)=gamma/d_x(j+1);
        
    end
end


for j=1:nim1-1
    for i=1:njm1-1
         Fs(i,j)=rho*v3d(i,j);
         Fn(i,j)=rho*v3d(i+1,j);
         Dn(i,j)=gamma/d_y(i+1);
         Ds(i,j)=gamma/d_y(i);
    end
end


for i=1:nim1-1
    for j=1:njm1-1
  aw1(i,j)=max(Fw(i,j)*dely(i),(Dw(i,j)+(Fw(i,j)*fxw(j)))*dely(i));
  a_w_B(i,j)=max(aw1(i,j),0);
  ae1(i,j)=max(-1*Fe(i,j)*dely(i),(De(i,j)-(Fe(i,j)*fxe(j)))*dely(i));
  a_e(i,j)=max(ae1(i,j),0);
    end
end


for i=(ni-(num+1)):(ni-2)
    fx_w_B(i+1)=0.5;
    fy_n_B(i+1)=0.5*(yc(i)-yc(i+1))/(yp(i)-yp(i+1));
    fy_s_B(i+1)=0.5*(yc(i)-yc(i+1))/(yp(i+1)-yp(i+2)); 
end


for i=(ni-(num+1)):(ni-2)
    D_w_B(i+1)=gamma/(xp(nj)-xp(nj-1));
    F_w_B(i+1)=rho*u2d(i,nj-1);
    D_n_B(i+1)=gamma/(yp(i)-yp(i+1));
    D_s_B(i+1)=gamma/(yp(i+1)-yp(i+2));
    F_n_B(i+1)=rho*v3d(i,njm1-1);
    F_s_B(i+1)=rho*v3d(i+1,njm1-1);
end

for i=(ni-(num+1)):(ni-2)
    a_w1_B(i+1)=max(F_w_B(i+1)*(yc(i)-yc(i+1)),0);
    a_w_B(i+1)=max(a_w1_B(i+1),(D_w_B(i+1)+F_w_B(i+1)*fx_w_B(i+1))*(yc(i)-yc(i+1)));

    a_n1_B(i+1)=max(-1*F_n_B(i+1)*d_x(j),0);
a_n_B(i+1)=max(a_n1_B(i+1),(D_n_B(i+1)-F_n_B(i+1)*fy_n_B(i+1))*(xp(nj)-xp(nj-1)));

    a_s1_B(i+1)=max(F_s_B(i+1)*d_x(j),0);
    a_s_B(i+1)=max(a_s1_B(i+1),(D_s_B(i+1)+F_s_B(i+1)*fy_s_B(i+1))*(xp(nj)-xp(nj-1)));
end

for i=(ni-(num+1)):(ni-2) 
    a_p_B(i+1)=a_w_B(i+1)+a_n_B(i+1)+a_s_B(i+1);
end


for i=1:nim1-1
    for j=1:njm1-1
        as1(i,j)=max(Fs(i,j)*delx(i),((Ds(i,j)+Fs(i,j)*fxs(j)))*delx(i));
        a_s(i,j)=max(as1(i,j),0);
        an1(i,j)=max(-1*Fn(i,j)*delx(i),((Dn(i,j)-Fn(i,j)*fxn(j)))*delx(i));
        a_n(i,j)=max(an1(i,j),0);
    end
end
for i=1:nim1-1
    for j=1:njm1-1
        a_p(i,j)=a_w_B(i,j)+a_e(i,j)+a_s(i,j)+a_n(i,j);
    end
end

T=ones(ni,nj);
D=ones(nim1,njm1);
T(ni-num:ni,1)=20;
T(num+2:ni,nj)=50;
itr=0;
error=1;
errormax=0.01;
 while (error>errormax)
    itr=itr+1;
    T_old=T;
    D_old=D;
    resid=0;
    res=0;
    residual(itr)=0;

for i=2:nim1
     for j=2:njm1
   
          D(i,j)=a_n(i-1,j-1)*T(i+1,j)+a_s(i-1,j-1)*T(i-1,j);
          B(i,j)=a_e(i-1,j-1);
          A(i,j)=a_p(i-1,j-1);
          C(i,j)=a_w_B(i-1,j-1);
          if j==2
             P(i,j)=B(i,j)/A(i,j);
             Q(i,j)=(D_old(i,j)+C(i,j)*T(i,j-1))/A(i,j);
          else
             P(i,j)=B(i,j)/(A(i,j)-C(i,j)*P(i,j-1));
             Q(i,j)=(C(i,j)*Q(i,j-1)+D_old(i,j))/(A(i,j)-C(i,j)*P(i,j-1));
          end

     end
     
end
for i=nim1:-1:2
       for j=njm1:-1:2
          T(i,j)=P(i,j)*T(i,j+1)+Q(i,j);
          T(1:ni-(num+1),1)=T_old(1:ni-(num+1),2);
          % CHANGE OF BOUNDARY CONDITION ON THE LEFT SIDE 
          %T(1:ni-(num+1),1) =50;

                   for G=(ni-(num+1)):(ni-2)
                        T(G+1,nj)=((a_w_B(G+1)*T(G,nj-1)) + (a_n_B(G+1)*T(G-1,nj)) + (a_s_B(G+1)*T(G+1,nj))) / (a_p_B(G+1));
                    end
          
          T(1,1:nj)=T_old(2,1:nj);
          T(ni,1:nj)=T_old(ni-1,1:nj);
          %T(ni,1:nj) = 50;

       end
end
for i=2:nim1-1
     for j=2:njm1-1
       resid=resid+abs((a_e(i-1,j-1)*T(i,j+1))+(a_w_B(i-1,j-1)*T(i,j-1))+(a_n(i-1,j-1)*T(i+1,j))+(a_s(i-1,j-1)*T(i-1,j))-(a_p(i-1,j-1)*T(i,j)));
     end
end
R=nj;
for i=1:6
    T(i,nj)=T(R,nj);
    R=R-1;
end
T(ni-num:ni,nj)=50;
T(1,1:nj)=T_old(2,1:nj);
residual(itr)=resid;
temp=T(1:num+1,nj);
del_T=abs(20-mean(temp));
f=rho*1*h*del_T;
error=residual(itr)/f
 end
TEMP_F = flip(T);

 for i = 2:nj-1
     D_T_IN(i) = (u2d(i,1)*(T(i,2)-T(i,1))) -(gamma* (T(i,2)-T(i,1))/(xp(2)-xp(1)));
 end

%  for i = 1:num+1
%      D_T_OT(i) = (u2d(i,nj)*(T(i,nj)-T(i,nj-1))) - (gamma*T(i,nj)-T(i,nj-1)/(xp(27)-xp(26)));
%  end
   for i = 2:ni-1
     D_T_B(i) = (-gamma*(T(i,nj)-T(i,nj-1))/(xp(nj)-xp(nj-1)))+(u2d(i,nj)*(T(i,nj)-T(i,nj-1)));
   end
   for i=2:ni-1
        flux(i) = D_T_IN(i)-D_T_B(i);
   end

 figure
contourf(xp,yp,T)
figure
contour(xp,yp,T,10000)
