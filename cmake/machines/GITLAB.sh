#NOSHOW
export FC=mpif90
FLEUR_INCLUDEDIR="/opt/include /usr/include"
#FLEUR_LIBRARIES="-L/opt/lib;-lxcf03;-lxc;-ldl;-L/usr/lib;-L/usr/lib/x86_64-linux-gnu;-lxml2;-lscalapack-openmpi;-lblacsF77init-openmpi;-lblacs-openmpi;-llapack;-lblas"
#FLEUR_LIBRARIES="-lwannier;-L/opt/lib;-ldl;-L/usr/lib;-L/usr/lib/x86_64-linux-gnu;-lxml2;-lscalapack-openmpi;-lblacsF77init-openmpi;-lblacs-openmpi"
FLEUR_LIBRARIES="-L/opt/lib;-ldl;-L/usr/lib;-L/usr/lib/x86_64-linux-gnu;-lscalapack-openmpi;-lfftw3"
#FLEUR_LIBRARIES="-L/opt/lib;-L/usr/lib"
