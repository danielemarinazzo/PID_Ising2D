% The information decomposition of the transfer entropy from two neighboring spins  to a target
% spin is evaluated according to the approach described
% in Bertschinger et al., Entropy 16(4), pp. 2161-2183, 2014.
%As an example, the file G6spin_256_0.mat cointans the probability distribution of six
%spins, arranged as in the figure below, in the 2D Ising model on a 256x256 lattice with Glauber dynamics, for several values of the coupling beta.

%Spin 1 is the target; the two driving spins can be chosen among the
%remaining 5. In the output figure, the synergy S, the redundancy R, and
%unique information terms U_1 and U_2 are depicted versus beta.




%        spin configuration
%               4 6
%             2 1 3
%               5
%
clear;clc;load('G6spin_256_0.mat');
betac=0.440687;
ind=find(beta<betac);
Pm=squeeze(mean(P,2));
x=beta(ind);
N=256;
drs=[2 3]; %the two driving spins, the be chosen among 2,3,4,5,6
for ib=1:length(ind)
   [SS(ib) RR(ib) UU1(ib) UU2(ib)]=TE_red(drs(1),drs(2),Pm(:,ib),beta(ib),N);
end
ymn=min(0,min([RR SS UU1 UU2]));
ymx=max([RR SS UU1 UU2]);
dy=0.03*(ymx-ymn);
figure(drs(1));plot(x, SS,'*',x,RR,'s',x,UU1,'>',x,UU2,'<');
line([betac betac],[ymn ymx],'LineStyle',':','color','k', 'LineWidth', 2);
ylim([ymn ymx+dy]);
xlabel('\beta','FontSize',20,'FontWeight','bold');
ylabel('Transfer Entropy','FontSize',16,'FontWeight','bold');
legend({'S','R','U_1' ,'U_2'},'Location','NorthWest' ,'FontSize',18,'FontWeight','bold');