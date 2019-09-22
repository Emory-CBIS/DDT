%datum--correlation matrices
%design--design matrix
%    -
%model--adjust (model=1), unadjusted(model=0)
function [Dinter_new, Sign_new]=adjusted_pvalue_lin_table(datum,design,model)
%Input
%  datum:   struct object containing the correlation matrices
%  design:  design matrix of covariates; first column must be group
%           membership (1-condition, 0 control); DO NOT include a column of
%           1s
%  model:   indicator of whether to adjust for covariates in design matrix
%           (model=1) or ignore covariates in design matrix (model=0)
%Output
%  Dinter_new: estimated matrix of p-values
%  Sign_new:   Matrix of the sign of the group coefficient at each edge
N=size(datum(1).data,1);
Dinter_new=zeros(N,N);
Sign_new=zeros(N,N);
if model==1
    Xlin=design;
else
    Xlin=design(:,1);
end
for j=1:(N-1)
        for k=(j+1):N
            latt=zeros(size(datum,2),1);
           for h=1:size(datum,2)
              latt(h)=datum(h).data(j,k);
           end
           fit=fitlm(Xlin, latt);
           Dinter_new(j,k)=fit.Coefficients.pValue(2);
           Sign_new(j,k)=sign(fit.Coefficients.Estimate(2));
        end       
end