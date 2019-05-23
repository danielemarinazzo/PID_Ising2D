function [S R U1 U2 data Fm v]=TE_red(i,j,p,beta,N)
Z=[[0, 0, 0, 1];[ 0, 0, 1, 0];[0, 0, 0,-1];[0, 0,-1, 0];[0, 0, 0,-1]; [ 0, 0,-1, 0];...
   [0, 0, 0, 1];[ 0, 0, 1, 0];[0, 1, 0, 0];[1, 0, 0, 0];[0,-1, 0, 0]; [-1, 0, 0, 0];...
   [0,-1, 0, 0];[-1, 0, 0, 0];[0, 1, 0, 0];[1, 0, 0, 0]];
[data p1ijt p1t p1]=TE(i,j,p,beta,N);
p11=reshape(p1ijt,16,1);
v_start=zeros(4,1);
fun=@(x)-FunTE(x,Z,p11);
warning('off');
options=optimset('Display','off','TolFun',1.e-12,'TolX',1.e-12,'Algorithm','interior-point');
[v fmax,exitflag,output]=fminunc(fun,v_start,options);
warning('on');
Fm=-fmax;
p_red=p11+Z*v;
p1ij_red=sum(reshape(p_red,8,2),2);
red=0;
for k=1:16
    s=anticodice(k,4)+1;
    k3=codice(s(1:3)-1);
    red=red+p_red(k)*log(p_red(k)*p1(s(1))/(p1t(s(1),s(4))*p1ij_red(k3)));
end
red=red*N;
S=data(1)-red;
R=data(2)+data(3)+S-data(1);
U1=data(2)-R;
U2=data(3)-R;

