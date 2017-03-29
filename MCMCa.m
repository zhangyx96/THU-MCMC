close all;clear;clc;
tic;
N=1;
p=zeros(1,N);
for m=1:N
cur_x=[5;10];
mu=[5;10];
sigma=[1,1;1,4];
cur_f=fx(cur_x,mu,sigma);
n=200000;  %模拟次数
save_x=[];
 r1=unifrnd(-1,1,1,n);
 r2=unifrnd(-1,1,1,n);
 r3=unifrnd(0,1,1,n);
 r=[7.5*r1;15*r2;r3];
for i=1:n
    new_x=cur_x+r(1:2,i);
    new_f=fx(new_x,mu,sigma);
    if new_f/cur_f>r(3,i)
        cur_f=new_f;
        cur_x=new_x;
        save_x=[save_x cur_x];
    end
end
x1=save_x(1,500:end);
x2=save_x(2,500:end);
z=corrcoef(x1,x2);
p(m)=z(2,1);
end
e=100*sum(abs(p-0.5)/0.5)/N;  %平均相对误差
lg=length(x1);
x_theory=mvnrnd(mu,sigma,lg);
subplot(1,2,1),
plot(x_theory(:,1),x_theory(:,2),'.'),
set(gca,'XLim',[0,10],'YLim',[0,20]),
xlabel('x'),ylabel('y'),
title('理论计算结果','FontSize',14);
subplot(1,2,2),
plot(x1(1:end),x2(1:end),'.'),
set(gca,'XLim',[0,10],'YLim',[0,20]),
xlabel('x'),ylabel('y'),
title('Metropolis-Hasting算法','FontSize',14);
toc;