function [ Null, Check, thresh] = HQS_fun_nonpar( C,M )
%HQS function
%   Returns null "Corr" matrices which respect the distribution of the
%   psuedo correlation matrix C.  
%   C--observed PSD matrix
%   M--number of null networks desired
%   CC--transform elements of C:  log(cij/(1-cij)) in (-Inf,Inf)
CC=log(C./(1-C));

%for computational ease...
locspos=find(CC>=12);
locsneg=find(CC<=-12);
neg=-12+(-13+12)*rand(1); pos=12+(13-12)*rand(1);
CC(locspos)=pos; %random number between 9 and 10...
CC(locsneg)=neg; %random number between -10 and -9...

Null=struct('null',[]);
Check=struct('null',[]);
quant=zeros(M,1);
N=size(C,2);
e=(sum(CC(:))-trace(CC))/(N*(N-1));
  D1=CC;
  D1(eye(N)~=0)=nan;
v=nanvar(D1(:)); %variance of off diagonal elements
ebar=trace(CC)/N;
m=max(2,floor((ebar^2-e^2)/v));
Mu=sqrt(e/m);
sigsq=-Mu^2+sqrt(Mu^4+(v/m));


for j=1:M
    X=Mu+sqrt(sigsq)*randn(N,m);
    XX=X*X';
    Check(j).null=XX;
    Null(j).null=exp(XX)./(1+exp(XX));
    quant(j)=prctile(XX(find(tril(ones(size(XX)),-1))),97.5);
end

thresh=exp(max(real(quant)))/(1+exp(max(real(quant))));

end


