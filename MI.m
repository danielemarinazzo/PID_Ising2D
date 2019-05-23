%        configurazione spin
%
%              s4-------s6
%               |        |
%               |        |
%       s2------s1------s3
%               |
%               |
%              s5
% input:   i        drivers 1
%          j        drivers 2
%          p(2^ns)  spin correlations  ns number spins
% output:  mi       Mutual Informazio (i,j)
%          mib1     Mutual Informazio (i,1)
%          mib2     Mutual Informazio (j,1)
%          p1ij     spin corrlations (1 i j)
function [mi23 mib1 mib2 p1ij p1i p1j]=MI(i,j,p)
mi23=0;mib1=0;mib2=0;
np=length(p);
ns=log2(np);
p=reshape(p,2,2,2,2,2,[]);
ind=setdiff(1:ns,[1 i j]);
p1ij=sum(p,ind(1));
for j=2:length(ind);
    p1ij=sum(p1ij,ind(j));
end
p1ij=squeeze(p1ij);
pij=squeeze(sum(p1ij,1));
p1i=squeeze(sum(p1ij,3));
p1j=squeeze(sum(p1ij,2));
pi=sum(pij,2);
for s1=1:2
    for si=1:2
        mib1=mib1+p1i(s1,si)*log(p1i(s1,si)/(pi(s1)*pi(si)));
        mib2=mib2+p1j(s1,si)*log(p1j(s1,si)/(pi(s1)*pi(si)));
        for sj=1:2
            mi23=mi23+p1ij(s1,si,sj)*log(p1ij(s1,si,sj)/(pi(s1)*pij(si,sj)));
        end
    end
end

