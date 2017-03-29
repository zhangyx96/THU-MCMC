function vb_base=baserate(a)
[numcases,numdims,numbatches]=size(a);
count_int=zeros(numdims,1);
for batch=1:numbatches
    xx=sum(a(:,:,batch));
    count_int=count_int+xx';
end
lp=5;
p_int=(count_int+lp*numbatches)/(numcases*numbatches+lp*numbatches);
log_base_rate = log(p_int) - log(1-p_int);
vb_base=log_base_rate';

