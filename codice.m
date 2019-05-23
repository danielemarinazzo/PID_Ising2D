function k=codice(s)
k=1;
for i=1:length(s)
    k=k+2^(i-1)*s(i);
end
