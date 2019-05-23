function s=anticodice(i,l)
s=zeros(l,1);
j=i-1;
for h=1:l
    c=2^(l-h);
    f=mod(j,c);
    s(l-h+1,1)=(j-f)/c;
    j=f;
end