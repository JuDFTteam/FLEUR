/*--------------------------------------------------------------------------------
 * Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
 * This file is part of FLEUR and available as free software under the conditions
 * of the MIT license as expressed in the LICENSE file in more detail.
 *--------------------------------------------------------------------------------
 */

#include <stdio.h>
#include "oxides_validation_econfig.h"

/*
 * This method together with the variables defined in oxides_validation_econfig.h
 * writes out the file oxides_validation.econfig.
 */
int dropOxidesValidationEConfig()
{
  int errorCode = 0;
  FILE *file;
  file = fopen("oxides_validation.econfig", "w");
  errorCode = fprintf(file,"%.*s",oxides_validation_econfig_len, oxides_validation_econfig);
  fclose(file);
  if(errorCode < 0) return 1;
  return 0;
}

/*
void main(){
  dropOxidesValidationEConfig();
}
*/
