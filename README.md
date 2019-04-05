# PID_Ising2D

The information decomposition of the mutual information, as well as those of the transfer entropy, between a target
spin and two neighbors is evaluated according to the approach described in Bertschinger et al., Entropy 16(4), pp. 2161-2183, 2014.
As an example, the file G6spin_256_0.mat cointans the probability distribution of six spins, arranged as in the figure below, in the 2D Ising model on a 256x256 lattice with Glauber dynamics, for several values of the coupling beta.

Spin 1 is the target; the two driving spins can be chosen among the remaining 5. In the output figure, the synergy S, the redundancy R, and unique information terms U_1 and U_2 are depicted versus beta.

To compute the PID of the mutual information, execute PID_MI; for the PID of the transfer entropy, execute PID_TE (drivers can be
changed modyfing the corresponding scripts).

References
Quantifying Unique Information https://www.mdpi.com/1099-4300/16/4/2161
Synergy as a warning sign of transitions: the case of the two-dimensional Ising model https://arxiv.org/abs/1901.05405
