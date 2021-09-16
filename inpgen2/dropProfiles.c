/*--------------------------------------------------------------------------------
 * Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
 * This file is part of FLEUR and available as free software under the conditions
 * of the MIT license as expressed in the LICENSE file in more detail.
 *--------------------------------------------------------------------------------
 */

#include <stdio.h>
#include "profileConfig.h"

/*
 * This method together with the variables defined in default_econfig.h
 * writes out the file default.econfig.
 */
int dropProfiles()
{
  int errorCode = 0;
  FILE *file;
  file = fopen("profile.config", "w");
  errorCode = fprintf(file,"%.*s",profile_config_len, profile_config);
  fclose(file);
  if(errorCode < 0) return 1;
  return 0;
}

/*
void main(){
  dropProfiles();
}
*/

