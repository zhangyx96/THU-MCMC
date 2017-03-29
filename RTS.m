function logZZ=RTS(W,b_B,a_B,beta,batchdata)
a_A=zeros(1,length(a_B));
b_A=baserate(batchdata);
[v_num,h_num]=size(W);
N=100;
K=length(beta);
logZA=sum(log(1+exp(b_A)))+sum(log(1+exp(a_A)));
logZ=zeros(1,K);  %Zk初始化为0

%算法执行200次，每次得到更精确的ZK初始值带入
for n=1:300
    c=zeros(1,K);
    v=zeros(1,v_num);
    %每次从beta中随机抽取一个betak
    beta_next=beta(randi([1,K]));
    for i=1:N
        prob_hB_v=1./(1+exp(-beta_next*(v*W+a_B)));
        hB=prob_hB_v>rand(size(prob_hB_v));
        prob_v_h=1./(1+exp(-(1-beta_next)*(+b_A)-beta_next*(hB*W'+b_B)));
        v=prob_v_h>rand(size(prob_v_h));
        log_fk=zeros(1,K);
        v_bA=v*b_A';v_bB=v*b_B';Wh=v*W+a_B;
        for k=1:K
            log_fk(k)=(1-beta(k))*v_bA+sum(log(1+exp((1-beta(k))*a_A)))+beta(k)*v_bB+sum(log(1+exp(beta(k)*Wh)));
            log_fk(k)=log_fk(k)-logZ(k)-(1-beta(k))*logZA;
        end
%         log_fk=(1-beta)'*v_bA+log(prod(1+exp((1-beta)'*a_A),2))+log(prod(1+exp(beta'*Wh),2));
%         log_fk=(log_fk'-logZ-(1-beta)*logZA);
        prob_beta_v=exp(log_fk-logsum(log_fk(:)));
        beta_next=beta(randi([1,K]));
        c=c+prob_beta_v/N;
    end
    logZ=[logZ(1),logZ(2:K)+log(c(2:K)/c(1))];
end
logZZ=logZ(end);
