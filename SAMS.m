function Z=SAMS(W1,W2,W3,W4,b1,b2,b3,b4,a1,a2,a3,a4)
v_num=[length(b1),length(b2),length(b3),length(b4)];
h_num=[length(a1),length(a2),length(a3),length(a4)];
W=zeros(max(v_num),max(h_num),length(v_num));
W(:,1:h_num(1),1)=W1 ;
W(:,1:h_num(2),2)=W2;
W(:,1:h_num(3),3)=W3;
W(:,1:h_num(4),4)=W4;
b=zeros(4,max(v_num));
b(1,:)=b1;
b(2,:)=b2;
b(3,:)=b3;
b(4,:)=b4;
a=zeros(4,h_num(4));
a(1,1:h_num(1))=a1 ;
a(2,1:h_num(2))=a2;
a(3,1:h_num(3))=a3;
a(4,1:h_num(4))=a4;
G=[0,1,0,0;0.5,0,0.5,0;0,0.5,0,0.5;0,0,1,0];
pai=[0.25,0.25,0.25,0.25];
Z=zeros(1,4);
m=4;
t0=randi([2*m,100*m]);
L=1;
v=zeros(1,v_num(1));
N=100000;
% start the interate procedure
for t=1:N
    % generate j from Gamma(L,:)
    j=find(cumsum(G(L,:))>rand,1,'first');
    % calculate the accept probability of j
    W_j=W(1:v_num(j),1:h_num(j),j);
    b_j=b(j,:);
    a_j=a(j,1:h_num(j));
    log_q_j=b_j*v'+sum(log(1+exp(v*W_j+a_j)));
    log_p_j = log(1/m) - Z(j) + log_q_j;
    
    W_L=W(1:v_num(L),1:h_num(L),L);
    b_L=b(L,:);
    a_L=a(L,1:h_num(L));
    log_q_L=b_L*v'+sum(log(1+exp(v*W_L+a_L)));
    log_p_L = log(1/m) - Z(L) + log_q_L;
    
    temp = min([1,G(j,L)/G(L,j)*exp(log_p_j - log_p_L)]);
    t_rand = rand;
    L = (t_rand<temp)*j+(t_rand>=temp)*L;
    % generate v_(t+1) using Gibbs sample
    W_L=W(1:v_num(L),1:h_num(L),L);
    b_L=b(L,:);
    a_L=a(L,1:h_num(L));
    
    prob_h_v=sigm(v*W_L+a_L);
    h=prob_h_v>rand(size(prob_h_v));
    prob_v_h=sigm(h*W_L'+b_L);
    v=prob_v_h>rand(size(prob_v_h));
    
    % update Z
    delta=((1:m)==L);
    gamma_t=t0/max([t0,t]);
    Z=Z+gamma_t*(delta-pai);
    Z=Z-Z(1);
end
end