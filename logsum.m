function ls = logsum(xx,dim)

if(length(xx(:))==1) ls=xx; return; end

xdims=size(xx);
if(nargin<2) 
  dim=find(xdims>1);
end
alpha = max(xx,[],dim)-log(realmax)/2;
repdims=ones(size(xdims));
repdims(dim)=xdims(dim);
ls = alpha+log(sum(exp(xx-repmat(alpha,repdims)),dim));

