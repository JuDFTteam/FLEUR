Contributors guide
======================================

Everyone is very welcome to contribute to the enhancement of FLEUR.
Please use the [gitlab service] (https://iffgit.fz-juelich.de/fleur/fleur) to obtain the
latest development version of FLEUR.

##Coding rules for FLEUR:
In no particular order we try to collect items to consider when writing code for FLEUR

- Instead of 'stop' use calls to judft_error, judft_warn, judft_end

- Do not read and write any files. Files are neither replacements for common-blocks nor a
storage for status variables.
 Only exceptions: 
+ you create nice IO subroutines in the subdirectory io
+ you write to the typical FLEUR output files

