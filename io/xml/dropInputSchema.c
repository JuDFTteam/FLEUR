/*--------------------------------------------------------------------------------
 * Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
 * This file is part of FLEUR and available as free software under the conditions
 * of the MIT license as expressed in the LICENSE file in more detail.
 *--------------------------------------------------------------------------------
 */

#include <stdio.h>
#include "inputSchema.h"
#include "inputSchema_old.h"
#include <string.h>
/*
 * This method together with the variables defined in inputSchema.h
 * writes out the file FleurInputSchema.xsd.
 *                                            GM'16
 */
int dropInputSchema(char* version)
{
  char * xsd_txt;
  int xsd_len;
   if (strcmp(version,"0.35")==0){
     xsd_len=FleurInputSchema_xsd_len;
     xsd_txt = FleurInputSchema_xsd;
   }else if(strcmp(version,"0.34")==0){
     xsd_len=FleurInputSchema0_34_xsd_len;
     xsd_txt = FleurInputSchema0_34_xsd;
   }else if(strcmp(version,"0.33")==0){
     xsd_len=FleurInputSchema0_33_xsd_len;
     xsd_txt = FleurInputSchema0_33_xsd;
   }else if(strcmp(version,"0.32")==0){
     xsd_len=FleurInputSchema0_32_xsd_len;
     xsd_txt = FleurInputSchema0_32_xsd;
   }else if(strcmp(version,"0.31")==0){
     xsd_len=FleurInputSchema0_31_xsd_len;
     xsd_txt = FleurInputSchema0_31_xsd;
   }else if(strcmp(version,"0.27")==0){
     xsd_len=FleurInputSchema0_27_xsd_len;
     xsd_txt = FleurInputSchema0_27_xsd;
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
