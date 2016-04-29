#include <stdio.h>
#include "inputSchema.h"

/*
 * This method together with the variables defined in inputSchema.h
 * writes out the file FleurInputSchema.xsd.
 *                                            GM'16
 */
int dropInputSchema()
{
   char schemaString[FleurInputSchema_xsd_len + 1];
   int i = 0;
   int errorCode = 0;
   FILE *file;
   for (i = 0 ; i < FleurInputSchema_xsd_len ; ++i)
   {
      schemaString[i] = FleurInputSchema_xsd[i];
   }
   schemaString[FleurInputSchema_xsd_len] = '\0';
   file = fopen("FleurInputSchema.xsd", "w");
   errorCode = fprintf(file,"%s", schemaString);
   fclose(file);
   if(errorCode < 0) return 1;
   return 0;
}

/*
 * How to create the file inputSchema.h if the XML Schema has to be
 * changed:
 *
 * You have to write it by hand. ;)
 * ...But if you prefer an automatic generation just follow this recipe:
 *
 * 1. Generate the file FleurInputSchema.xsd with the dropInputSchema 
 *    method in this file.
 * 2. Change the XML Schema file as desired.
 * 3. run: xxd -i FleurInputSchema.xsd inputSchema.h
 *    (Note that the name of the file FleurInputSchema.xsd determines
 *     the names of the variables)
 * 4. The new inputSchema.h is generated.
 */
