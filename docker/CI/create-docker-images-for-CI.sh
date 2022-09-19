#create the docker images for the CI and push them to the registry

for image in oneAPI NVHPC AOMP gfortran
do
   cd $image
   docker build -t iffregistry.fz-juelich.de/fleur/fleur:$image .
   docker push iffregistry.fz-juelich.de/fleur/fleur:$image
   cd ..
done
