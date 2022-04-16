% MATLAB program for Linear MPC: Reducing online computation (method 2)
clear all;
close all
% System parameters and simulation parameters
A=[0.9 0.2;-0.4 0.8];
B=[0.1;0.01];
NT=50;N=5;NC=2;n=2;m=1; 
Q=eye(n); QN=Q; R=eye(m);
Fx=[1 0;0 1;-1 0;0 -1];gx=[10;10;10;10];
Fu=[1;-1];gu=[1;1];


x0=[10;5];
x=zeros(n,NT+1); x(:,1)=x0;
Xk=zeros(n*(N+1),1); Xk(1:n,1)=x0;
u=zeros(m,NT);
Uk=zeros(m*NC,1);
zk=Uk;

% constructing AX,BU,QX,RU
for i=1:N+1
      AX((i-1)*n+1:i*n,:)=A^(i-1);
  for j=1:N
      if i>j
          BU((i-1)*n+1:i*n,(j-1)*m+1:j*m)=A^(i-j-1)*B;
      else
          BU((i-1)*n+1:i*n,(j-1)*m+1:j*m)=zeros(n,m);
      end    
  end
end
QX=Q;RU=R;
FX=Fx;gX=gx;FU=Fu;gU=gu;
for i=1:N-1
  QX=blkdiag(QX,Q); RU=blkdiag(RU,R);
  FX=blkdiag(FX,Fx);gX=[gX;gx];
  FU=blkdiag(FU,Fu);gU=[gU;gu];
end
QX=blkdiag(QX,QN);
FX=blkdiag(FX,Fx);
gX=[gX;gx];
H=BU(:,1:NC)'*QX*BU(:,1:NC)+RU(1:NC,1:NC);


% simulating system with MPC
for k=1:NT
   xk=x(:,k);
   qk=2*xk'*AX'*QX*BU(:,1:NC);rk=xk'*AX'*QX*AX*xk;
   fun = @(z)z'*H*z+qk*z+rk;
   F=[FX*BU(:,1:NC);FU(:,1:NC)];g=[gX-FX*AX*xk;gU];Feq=[];geq=[];
   lb=[];ub=[];
   z=fmincon(fun,zk,F,g,Feq,geq,lb,ub);
   u(:,k)=z(1:m,1);
   x(:,k+1)=A*x(:,k)+B*u(:,k);
   zk=z;
end     

% plotting response
figure(1)
time = (0:NT);
subplot(2,1,1)
plot(time,x(1,:),'r.-','LineWidth',.7) 
hold on
plot(time,x(2,:),'k.-','LineWidth',.7) 
legend('$x_1$','$x_2$','Interpreter','latex');
%axis ([0 50-10 10])
xlabel('$k$','Interpreter','latex');ylabel('$\textbf{x}_{k}$','Interpreter','latex');
grid on
ax = gca;
set(gca,'xtick',[0:5:50])
set(gca,'ytick',[-10:5:10])
ax.GridAlpha = 1
ax.GridLineStyle = ':'
subplot(2,1,2)
stairs(time(1:end-1),u,'r.-','LineWidth',.7)
%axis([0 50 -10 0])
xlabel('$k$','Interpreter','latex');ylabel('${u}_{k}$','Interpreter','latex');
grid on
ax = gca;
set(gca,'xtick',[0:5:50])
set(gca,'ytick',[-1:.5:1])
ax.GridAlpha = 1
ax.GridLineStyle = ':'
print -dsvg lmpc4


