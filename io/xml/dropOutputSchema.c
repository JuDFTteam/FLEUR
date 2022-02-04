#include <stdio.h>
#include "outputSchema.h"
#include "outputSchema_old.h"
#include <string.h>
/*
 * This method together with the variables defined in outputSchema.h
 * writes out the file FleurOutputSchema.xsd.
 *                                            GM'16
 */
int dropOutputSchema(char* version)
{
  char * xsd_txt;
  int xsd_len;
   if (strcmp(version,"0.35")==0){
     xsd_len=FleurOutputSchema_xsd_len;
     xsd_txt = FleurOutputSchema_xsd;
   }else if (strcmp(version,"0.34")==0){
     xsd_len=FleurOutputSchema0_34_xsd_len;
     xsd_txt = FleurOutputSchema0_34_xsd;
   }else{
     return 1;
   }
   char schemaString[xsd_len + 1];
   int i = 0;
   int errorCode = 0;
   FILE *file;
   for (i = 0 ; i < xsd_len ; ++i)
   {
      schemaString[i] = xsd_txt[i];
   }
   schemaString[xsd_len] = '\0';
   file = fopen("FleurOutputSchema.xsd", "w");
   errorCode = fprintf(file,"%s", schemaString);
   fclose(file);
   if(errorCode < 0) return 1;
   return 0;
}

/*
 * How to create the file outputSchema.h if the XML Schema has to be
 * changed:
 *
 * You have to write it by hand. ;)
 * ...But if you prefer an automatic generation just follow this recipe:
 *
 * 1. Generate the file FleurOutputSchema.xsd with the dropOutputSchema
 *    method in this file.
 * 2. Change the XML Schema file as desired.
 * 3. run: xxd -i FleurOutputSchema.xsd outputSchema.h
 *    (Note that the name of the file FleurOutputSchema.xsd determines
 *     the names of the variables)
 * 4. The new outputSchema.h is generated.
 */
