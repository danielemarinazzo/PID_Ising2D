function [S R U1 U2 mi1 mi2 mi v]=MI_red(i,j,p)
[mi mi1 mi2 p123]=MI(i,j,p);
pp=reshape(p123,8,1);
v_start=zeros(2,1);
warning('off');
fun=@(x)FunMI(x,pp);
options=optimset('Display','off','TolFun',1.e-4,'TolX',1.e-8);
[v red]=fminunc(fun,v_start,options);
warning('on');
S=mi-red;
R=mi1+mi2-red;
U1=mi1-R;
U2=mi2-R;
