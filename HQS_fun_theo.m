function [ Null, Check, thresh] = HQS_fun_theo( C,M )
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
end
df=m;%degrees of freedom
ncp=abs(m*(4*(Mu)^2)/(2*sigsq));%non-centrality parameter; abs() bc it was producing negative values-->NAN values for pd
mcon=2*sigsq/4; %constant 
pd=mcon*ncx2rnd(df,ncp,1000000,1)-mcon*chi2rnd(df,1000000,1);
%vec=zeros(1000,1);
%for i=1:1000
%    vec(i)=quantile(pd(:,i),.975);
%end
%ll=mean(vec);
%thresh=exp(ll)/(1+exp(ll));
%display(num2str(thresh))
ll=quantile(pd,.975);
thresh=exp(ll)/(1+exp(ll));
end