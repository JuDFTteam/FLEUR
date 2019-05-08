# CalculatingMCDSpectra

A possibility to calculate magnetic circular dichroism spectra was implemented similar to Wu et al. [J. Magn. Magn. Mater. 132 (1994) p. 103-127][2]. The input is provided by a file mcd\_inp that contains a lower and upper energy bound in hartree units. Transitions from all core-states within these bounds to the valence electrons are taken into account. Ff you set flags to generate DOS output, three files MCD\_SPEC.[+,-,0] are generated containing columnwise the MCD-spectra for all involved atoms for circular polarized light [+,-] and linear polarization [0].

 [2]: http://dx.doi.org/10.1016/0304-88539490305-0
