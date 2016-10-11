# Choice of step size in the EMC program

## Introduction
In the manual of the EMC program the effective masses of GaAs and InP are calculated as a function of the step size. The step size determines the grid on which the second derivatives of energy with respect to the k-vector are performed. In Figure 1 we overlay the grid with a step size of 0.05 1/Bohr on the band structure of GaAs. It is clear that this step size is too large.

Figure 1: Band structure of GaAs along the ΓX direction (black lines). A parabola (red lines) is fitted through grid points (step size is 0.05 1/Bohr). Top of valence band is set to 0 eV.

This leads to sampling of the bands away from the parabolic region. The step sizes in the manual range from 0.005 to 0.1 1/Bohr. We also note that the calculations in the manual are performed on a conventional cell rather than the primitive one, which leads to "band folding". The Γ-point of the conduction band is sampled correctly, but the rest of the grid points actually lie on the band above the conduction band.

In Figure 2 we test the smallest step from the manual (0.005 1/Bohr). It turns out that this is the largest step size that leads to grids within the parabolic region. We use this step size to compare our results to a paper by the Kresse group [1]. Their calculations include spin-orbit coupling (SOC), which is necessary to correctly describe the split-off band in III-V semiconductors. Since our calculations also include SOC we obtain an excellent agreement (see Table below).

<table>
  <tr>
    <td colspan="2"></td><td colspan="2">GaAs</td><td colspan="2">InP</td>
  </tr>
  <tr>
    <td>Band</td><td>Description</td><td>Kresse</td><td>EMC</td><td>Kresse</td><td>EMC</td>
  </tr>
  <tr>
    <td>4</td><td>electron</td><td>0.030</td><td>0.029</td><td>0.054</td><td>0.054</td>
  </tr>
  <tr>
    <td>3</td><td>heavy hole</td><td>0.320</td><td>0.318</td><td>0.435</td><td><strong>0.371</strong></td>
  </tr>
  <tr>
    <td>2</td><td>light hole</td><td>0.036</td><td>0.035</td><td>0.073</td><td>0.074</td>
  </tr>
  <tr>
    <td>1</td><td>split-off hole</td><td>0.108</td><td>0.109</td><td>0.139</td><td>0.139</td>
  </tr>
</table>

The only problematic value is the heavy hole in InP. It turns out that the step is too small in this case. Since the band is rather flat the changes in the eigenvalues are too small and thus inaccurate. Recalculating with a larger step of 0.05 1/Bohr we obtain a value of 0.437 m0, close to the results of the Kresse group.

## Conclusions
Calculations on InP show that even in one material there does not exist a single reliable step value. It is thus advisable to overlay the grid on the band structure and check that one is within the parabolic region.

## Contributors
Karol Jarolimek (firstname.lastname@gmail.com)

## PDF version
<a href ="kjarolimek_step_benchmark.pdf" targe="_blank">kjarolimek_step_benchmark.pdf</a>

## References
1. Kim, Y.-S.; Marsman, M.; Kresse, G.; Tran, F.; Blaha, P. *Phys. Rev. B: Condens. Matter* **2010**, *82*, 205212.
