% The information decomposition of the mutual information between a target
% spin and two neighbors is evaluated according to the approach described
% in Bertschinger et al., Entropy 16(4), pp. 2161-2183, 2014.
%As an example, the file G6spin_256_0.mat cointans the probability distribution of six
%spins, arranged as in the figure below, in the 2D Ising model on a 256x256 lattice with Glauber dynamics, for several values of the coupling beta.

%Spin 1 is the target; the two driving spins can be chosen among the
%remaining 5. In the output figure, the synergy S, the redundancy R, and
%unique information terms U_1 and U_2 are depicted versus beta.


clear;clc;load('G6spin_256_0.mat');
betac=0.440687;
ind=find(beta<betac);
x=beta(ind);
Pm=squeeze(mean(P,2));
warning('off');
drs=[2 3];%the two driving spins, the be chosen among 2,3,4,5,6
target=1;
for ib=1:length(ind)
    [Sm(ib) Rm(ib) U1m(ib) U2m(ib)]=MI_red(drs(1),drs(2),Pm(:,ib));
end
warning('on');
s1=[num2str(drs(1)) '->1'];
s2=[num2str(drs(2)) '->1'];
s=[num2str(drs(1)) num2str(drs(2)) '->1' ];
ymn=min(0,min([Rm Sm U1m U2m]));
ymx=max([Rm Sm U1m U2m]);
dy=0.03*(ymx-ymn);
figure(drs(1));plot(x, Sm,'*',x,Rm,'s',x,U1m,'>',x,U2m,'<');
line([betac betac],[ymn ymx],'LineStyle',':','color','k', 'LineWidth', 2);
ylim([ymn ymx+dy]);
xlabel('\beta','FontSize',20,'FontWeight','bold');
ylabel('Mutual Information','FontSize',16,'FontWeight','bold');
legend({'S','R','U_1' ,'U_2'},'Location','NorthWest' ,'FontSize',18,'FontWeight','bold');
