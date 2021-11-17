exclude: struct_input.f90
exclude: soc_or_ssdw.f90
exclude: test_FFTW.f90
src_dir: ../
output_dir: ./ford-doc
project_github: https://iffgit.fz-juelich.de/fleur/fleur
project_website: https://www.flapw.de
summary: The FLEUR code: All-electron full-potential augmented plane-wave method for DFT.
author: FLEUR development team 
author_description: Developers mostly from FZ-Juelich.
email: d.wortmann@fz-juelich.de
fpp_extensions: fpp
predocmark: >
media_dir: ./media
docmark_alt: #
predocmark_alt: <
display: public
         protected
         private
source: false
graph: true
search: true
macro: TEST
       LOGIC=.true.
extra_mods: json_module: http://jacobwilliams.github.io/json-fortran/
            futility: http://cmacmackin.github.io
license: by-nc
extra_filetypes: sh #

