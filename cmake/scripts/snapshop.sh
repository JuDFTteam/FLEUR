#!/bin/bash
echo "Please make sure you are in a git directory on the develop branch"
echo "Also the pipeline status should be OK and your workdirectory should be clean"
echo "Press OK to continue"
read

TAG=`date "+%m.%y"`
TAG="snapshot-$TAG"

echo $TAG

git tag -a -m "Monthly snapshot" $TAG
git checkout stable
git merge --ff-only develop
git push origin $TAG
git push
