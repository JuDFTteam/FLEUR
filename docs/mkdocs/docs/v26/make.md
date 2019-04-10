Building FLEUR
=======================

The procedure is normally like this:

#  `imake`   (to generate the `Makefile`)
#  If needed, update the settings in the makefile
#  `make fleur.x inpgen.x`

When you later decide to change a preprocessor flag in the Imakefile or
Makefile, you should run

# ` imake ` to update the Makefile
# ` make rminv `  to remove preprocessor dependent object-files
# ` make fleur.x `

Then you have an appropriately changed fleur.x

