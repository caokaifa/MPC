close all;
clear;
L=76;
yaw=pi/9;

x=[];
y=[];
xL=[];
yL=[];
xR=[];
yR=[];
d=2;
N = 640;
kappa = (yaw*(2*pi))/N;

ds = L / N;

for i=1:N
   yaw=yaw+ kappa ;
    co = cos(yaw);
    si = sin(yaw);
    f = [co;si];
    n = [-si;co];
    if i==1
        x=[x ds*f(1)];
        y=[y ds*f(2)];
    else
        x=[x x(end)+ds*f(1)];
        y=[y y(end)+ds*f(2)];
    end
    xL=[xL  x(end)+d*n(1)];
    yL=[yL  y(end)+d*n(2)];
    xR=[xR  x(end)-d*n(1)];
    yR=[yR  y(end)-d*n(2)];
end

figure
plot(x,y,'b')
hold on
plot(xL,yL,'.r')
hold on
plot(xR,yR,'.y')

% cx=[]
% cy=[]
% for i=1:20
%     cx=[cx 2*cos(sqrt(i))+2*sin(sqrt(i))*sqrt(i)];
%     cy=[cy 2*sin(sqrt(i))-2*cos(sqrt(i))*sqrt(i)];
% end
% figure
% plot(cx,cy)
%     



