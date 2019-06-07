/*--------------------------------------------------------------------------------
 * Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
 * This file is part of FLEUR and available as free software under the conditions
 * of the MIT license as expressed in the LICENSE file in more detail.
 *--------------------------------------------------------------------------------
 */

#include <stdio.h>
#include "default_econfig.h"

/*
 * This method together with the variables defined in default_econfig.h
 * writes out the file default.econfig.
 */
int dropDefaultEconfig()
{
  int errorCode = 0;
  FILE *file;
  file = fopen("default.econfig", "w");
  errorCode = fprintf(file,"%.*s",default_econfig_len, default_econfig);
  fclose(file);
  if(errorCode < 0) return 1;
  return 0;
}

/*
void main(){
  dropDefaultEconfig();
}
*/
