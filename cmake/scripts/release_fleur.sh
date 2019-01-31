
echo " This script creates a tar file to be released on the website"
echo " Please make sure you have adjusted the release output in FLEUR"

git clone https://iffgit.fz-juelich.de/fleur/fleur.git 
rm -rf fleur/.git*

echo "Please enter name of release:"
read version

echo $version >fleur/version

tar czf fleur.release.tgz fleur

rm -rf fleur


