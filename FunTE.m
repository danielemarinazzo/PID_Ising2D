function f=FunTE(v,Z,p1)
p=p1+Z*v;
f=0;
if isempty(find(p<0,1))
    p0(1:8)=0;
    for k=1:16
        k3=mod(k-1,8)+1;
        p0(k3)=p0(k3)+p(k);
        f=f-p(k)*log(p(k));
    end
    for k=1:8
        f=f+p0(k)*log(p0(k));
    end
end