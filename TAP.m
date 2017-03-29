function logZZ=TAP(W,a,b)
[v_num,h_num]=size(W);
mh=0.5*ones(1,h_num);
mv=0.5*ones(1,v_num);
N=15000;
for i=1:N
    mh=sigm(b+mv*W-(mh-1/2).*((mv-mv.^2)*(W.^2)));
    mv=sigm(a+mh*W'-(mv-1/2).*((mh-mh.^2)*(W.^2)'));
end
mh = mh + 1e-6;
mv = mv + 1e-6;
s=mh*log(mh)'+(1-mh)*log(1-mh)'+(mv*log(mv)'+(1-mv)*log(1-mv)');
logZZ=abs(-s+a*mv'+b*mh'+mv*W*mh'+(mv-mv.^2)*((W.^2)/2)*(mh-mh.^2)');
end
