function [ pvals_final,estvec,mean_null,sd_null ] = nodediff_piesttest( R_undernull,Dtilde_obs, bb)
%nodediffest: nodal connections that differ across the two populations
N=size(Dtilde_obs,2); %number of nodes
M=size(R_undernull,2);%number of null networks
    for m = 1:M
        LL=R_undernull(m).null;
        LL(find(eye(N)==1))=0;%set diagonals to 0 so they don't affect the counts
        MM(:,m)=sum(LL>=(bb),2);
    end
mean_null=mean(MM,2);
sd_null=std(MM,0,2);


estvec=(1/(M*(N-1)))*sum(MM,2); %vector of estimated p_star_i, for node i under the null


%
% BINOMIAL TEST PROCEDURE
%
Dtilde_obs(find(eye(N)==1))=0;

poss=0:(N-1);
und_null=repmat(nan,N,(N+1));%N nodes by (N+1)-->(0,..,N-1) possible values that can be taken, last column will be the observed probability

for i=1:N;
und_null(i,1:N)=binopdf(poss,(N-1),estvec(i));
end

obs_degree=sum(Dtilde_obs>=bb,2);
obs_prob=binopdf(obs_degree,repmat((N-1),N,1),estvec); %observed probability under the null

und_null(:,(N+1))=obs_prob;
pvals_final=repmat(nan,N,1);
for i=1:N;
    interest=und_null(i,1:N)<=und_null(i,(N+1));
    pvals_final(i)=sum(und_null(i,interest));
end


end

