# PID_Ising2D

The information decomposition of the mutual information, as well as those of the transfer entropy, between a target
spin and two neighbors is evaluated according to the approach described in Bertschinger et al., Entropy 16(4), pp. 2161-2183, 2014.

As an example, the file G6spin_256_0.mat contains the probability distribution of six spins, arranged as in the figure below, in the 2D Ising model on a 256x256 lattice with Glauber dynamics, for several values of the coupling beta.

![figure 1 of the paper](http://users.ugent.be/~dmarinaz/Ising_spins.png)

Spin 1 is the target; the two driving spins can be chosen among the remaining 5. In the output figure, the synergy S, the redundancy R, and unique information terms U_1 and U_2 are depicted versus beta.

To run the PID of the mutual information, execute PID_MI; for the PID of the transfer entropy, execute PID_TE (drivers can be
changed by modifying the scripts).

*References*

* Quantifying Unique Information https://www.mdpi.com/1099-4300/16/4/2161

* Synergy as a warning sign of transitions: the case of the two-dimensional Ising model, Phys. Rev. E 99, 040101(R) – Published 19 April 2019 - https://journals.aps.org/pre/abstract/10.1103/PhysRevE.99.040101. Also available on https://arxiv.org/abs/1901.05405
