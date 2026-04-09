/*--------------------------------------------------------------------------------
 * Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
 * This file is part of FLEUR and available as free software under the conditions
 * of the MIT license as expressed in the LICENSE file in more detail.
 *--------------------------------------------------------------------------------
 */

#include <stdio.h>
#include "default2_econfig.h"

/*
 * This method together with the variables defined in default2_econfig.h
 * writes out the file default2.econfig.
 */
int dropDefault2EConfig()
{
  int errorCode = 0;
  FILE *file;
  file = fopen("default2.econfig", "w");
  errorCode = fprintf(file,"%.*s",default2_econfig_len, default2_econfig);
  fclose(file);
  if(errorCode < 0) return 1;
  return 0;
}

/*
void main(){
  dropDefault2EConfig();
}
*/
