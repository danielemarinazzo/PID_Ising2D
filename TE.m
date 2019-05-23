function[data p1ijt p1t p1 p1ij num1 num2 den1 den2]=TE(i,j,p,beta,N);
%        configurazione spin
%
%              s4-------s6
%               |        |
%               |        |
%       s2------s1------s3
%               |
%               |
%              s5
%lo stato futuro dello spin 1 ? t
%
%calcoliamoci ora la causalit? (ij)->1
%serve        p(s(1),s(i),s(j),s(t))* p(s(1)
%        log(--------------------------------)
%             p(s(1),s(t))*p(s(1),s(i),s(t))
%calcoliamo la bivariata
%serve   <log(p(s1 s2 s6) p(s1)/p(s1s6)p(s2 s1))>
p1(1:2)=0;
p1t(1:2,1:2)=0;
p1i(1:2,1:2)=0;
p1j(1:2,1:2)=0;
p1ij(1:2,1:2,1:2)=0;
p1it(1:2,1:2,1:2)=0;
p1jt(1:2,1:2,1:2)=0;
p1ijt(1:2,1:2,1:2,1:2)=0;
p1234(1:2,1:2,1:2,1:2)=0;
p12346(1:2,1:2,1:2,1:2,1:2)=0;
H=0;Hi=0;Hj=0;Hg=0;H3=0;
np=length(p);
ns=log2(np);
for k=1:np
    s=anticodice(k,ns)+1;
    h=2*(s(2)+s(3)+s(4)+s(5))-12;
    p1(s(1))=p1(s(1))+p(k);
    p1i(s(1),s(i))=p1i(s(1),s(i))+p(k);
    p1j(s(1),s(j))=p1j(s(1),s(j))+p(k);
    p1ij(s(1),s(i),s(j))=p1ij(s(1),s(i),s(j))+p(k);

    pflip=1/(N*(1+exp(2*beta*(2*s(1)-3)*h)));
    p1t(s(1),s(1))=p1t(s(1),s(1))+p(k)*(1-pflip);
    p1t(s(1),3-s(1))=p1t(s(1),3-s(1))+p(k)*pflip;
    
    p1it(s(1),s(i),s(1))=p1it(s(1),s(i),s(1))+p(k)*(1-pflip);
    p1it(s(1),s(i),3-s(1))=p1it(s(1),s(i),3-s(1))+p(k)*pflip;
    p1jt(s(1),s(j),s(1))=p1jt(s(1),s(j),s(1))+p(k)*(1-pflip);
    p1jt(s(1),s(j),3-s(1))=p1jt(s(1),s(j),3-s(1))+p(k)*pflip;
    p1ijt(s(1),s(i),s(j),s(1))=p1ijt(s(1),s(i),s(j),s(1))+p(k)*(1-pflip);
    p1ijt(s(1),s(i),s(j),3-s(1))=p1ijt(s(1),s(i),s(j),3-s(1))+p(k)*pflip;
    p1234(s(1),s(2),s(3),s(4))=p1234(s(1),s(2),s(3),s(4))+p(k);
    p12346(s(1),s(2),s(3),s(4),s(1))=p12346(s(1),s(2),s(3),s(4),s(1))+p(k)*(1-pflip);
    p12346(s(1),s(2),s(3),s(4),3-s(1))=p12346(s(1),s(2),s(3),s(4),3-s(1))+p(k)*pflip;
end
for k=1:np
    
    s=anticodice(k,ns)+1;
    h=2*(s(2)+s(3)+s(4)+s(5))-12;
    pflip=1/(N*(1+exp(2*beta*(2*s(1)-3)*h)));
    Hg=Hg+p(k)*(1-pflip)*log((1-pflip)*p1(s(1))/p1t(s(1),s(1)))+...
          p(k)*pflip*log(pflip*p1(s(1))/p1t(s(1),3-s(1)));

    H3=H3+p(k)*(1-pflip)*log(p12346(s(1),s(2),s(3),s(4),s(1))*p1(s(1))/(p1t(s(1),s(1))*p1234(s(1),s(2),s(3),s(4))));
    H3=H3+p(k)*(pflip)*log(p12346(s(1),s(2),s(3),s(4),3-s(1))*p1(s(1))/(p1t(s(1),3-s(1))*p1234(s(1),s(2),s(3),s(4))));
    Hi=Hi+p(k)*(1-pflip)*log(p1it(s(1),s(i),s(1))*p1(s(1))/(p1t(s(1),s(1))*p1i(s(1),s(i))));
    Hi=  Hi+p(k)*(pflip)*log(p1it(s(1),s(i),3-s(1))*p1(s(1))/(p1t(s(1),3-s(1))*p1i(s(1),s(i))));
    Hj=Hj+p(k)*(1-pflip)*log(p1jt(s(1),s(j),s(1))*p1(s(1))/(p1t(s(1),s(1))*p1j(s(1),s(j))));
    Hj=Hj+p(k)*(pflip)*log(p1jt(s(1),s(j),3-s(1))*p1(s(1))/(p1t(s(1),3-s(1))*p1j(s(1),s(j))));
    num1=p1ijt(s(1),s(i),s(j),s(1))*p1(s(1));
    num2=p1ijt(s(1),s(i),s(j),3-s(1))*p1(s(1));
    den1=p1t(s(1),s(1))*p1ij(s(1),s(i),s(j));
    den2=p1t(s(1),3-s(1))*p1ij(s(1),s(i),s(j));
    H=H+p(k)*(1-pflip)*log(num1/den1)+p(k)*pflip*log(num2/den2);
end
data=N*[H Hi Hj Hg H3];
