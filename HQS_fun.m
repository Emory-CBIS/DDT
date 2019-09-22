function [ Null, Check, thresh] = HQS_fun( C,M )
%HQS function
%   Returns null "Corr" matrices which respect the distribution of the
%   psuedo correlation matrix C.  f_ij sampled from truncated normal dist
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
    %quant(j)=prctile(XX(find(triu(ones(size(XX)),1) )),95);
end
pd = makedist('Normal');
ll=icdf(pd,(1-2*.05/(N*(N-1))));
%threshun=e+sqrt(v)*ll; %threshold on (-Inf,Inf) scale
threshun=e+sqrt(v)*1.65;
%thresh=exp(ll)/(1+exp(ll));
thresh=exp(ll)/(1+exp(ll));
%Power analysis approach
%sig_tilde=(exp(2*e)*sigsq)/((1+exp(e))^4);
%thresh=.9 + sqrt(sig_tilde/M)*(norminv((1-(2*.05/(N*(N-1)))),0,1)-sqrt(2)*erfinv(-.6));
end