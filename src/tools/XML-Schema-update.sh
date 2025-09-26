#!/bin/bash

masci_tools_repo=https://github.com/JuDFTteam/masci-tools.git

in_schema=`grep "version=" FleurInputSchema.xsd |grep XMLSchema|cut -d= -f 3|sed -e s/[\"\>]//g`
out_schema=`grep "version=" FleurOutputSchema.xsd |grep XMLSchema|cut -d= -f 3|sed -e s/[\"\>]//g`

if [ -z "$in_schema" ] || [ -z "$out_schema" ]; then
  echo "Error: Unable to extract schema versions."
  exit 1
fi

echo "Input schema version: $in_schema"
echo "Output schema version: $out_schema"
#stop if not the same version
if [ "$in_schema" != "$out_schema" ]; then
  echo "Error: Input and output schema versions do not match."
  exit 1
fi  

#split version number into major and basic version
in_basic=`echo $in_schema | cut -d. -f 2`
in_major=`echo $in_schema | cut -d. -f 1`


#add current inputSchema.h.backup to inputSchema_old.h
sed -i .old  "s/FleurInputSchema_xsd/FleurInputSchema_${in_basic}_xsd/" inputSchema.h.backup
cat inputSchema.h.backup inputSchema_old.h > tmpSchema.h
mv tmpSchema.h inputSchema_old.h
#add current outputSchema.h.backup to outputSchema_old.h
sed -i .old  "s/FleurOutputSchema_xsd/FleurOutputSchema_${in_basic}_xsd/" outputSchema.h.backup
cat outputSchema.h.backup outputSchema_old.h > tmpSchema.h
mv tmpSchema.h outputSchema_old.h

#add number to schema version
in_new=`echo $in_major.$(($in_basic + 1))`
echo "New schema version: $in_new"

#update schema version in files
sed -i .old  "s/version=\"$in_schema\"/version=\"$in_new\"/g" FleurInputSchema.xsd
sed -i .old  "s/version=\"$out_schema\"/version=\"$in_new\"/g" FleurOutputSchema.xsd

#update version in xml file
sed -i .old  "s/enumeration value=\"$in_schema\"/enumeration value=\"$in_new\"/" FleurInputSchema.xsd
sed -i .old  "s/enumeration value=\"$out_schema\"/enumeration value=\"$in_new\"/" FleurOutputSchema.xsd

#create a new backup of the schema files
xxd -i FleurInputSchema.xsd > inputSchema.h.backup
xxd -i FleurOutputSchema.xsd > outputSchema.h.backup

#delete old backup files
rm -f *.old

#Now update the schema in masci-tools
python -m venv create masci_tools
source masci_tools/bin/activate
cd masci_tools
git clone $masci_tools_repo
pip install -e masci-tools
./bin/masci-tools fleur-schema add ../FleurInputSchema.xsd
./bin/masci-tools fleur-schema add ../FleurOutputSchema.xsd
./bin/masci-tools fleur-schema list

#update git repo 
cd masci-tools
git add masci_tools/io/parsers/fleur_schema/$in_new 
git commit -am "Update Fleur schemas to version $in_new"
#git push origin master  
git status
cd ../..
#deactivate and remove masci_tools
#deactivate
#rm -rf masci_tools

