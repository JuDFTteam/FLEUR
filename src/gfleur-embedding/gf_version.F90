!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_version 
      use m_juDFT
          IMPLICIT NONE
!-----------------------------------------------                        
! DESC:Displays version information and checks basic preprocessor settin
!                 Daniel Wortmann, (05-08-12)                           
!-----------------------------------------------                        
      CONTAINS 
                                                                        
      !<-- S: gf_hello()                                                
      SUBROUTINE gf_hello() 
!-----------------------------------------------                        
! Say hello!                                                            
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_out 
      use m_juDFT
      IMPLICIT NONE 
      !<--Locals                                                        
                                   !version of FLEUR used               
      CHARACTER(LEN = 9)  :: ivers 
                                        !The Date of the CVS-version of 
      CHARACTER(LEN = 80) :: CVSVersion 
                                        !this file                      
      CHARACTER(LEN = 80) :: compiletime 
      CHARACTER(LEN = 80) :: compileplace 
      !>                                                                
                                                                        
      ivers="unkown" 
      compiletime="unkown" 
      compileplace=" " 
      CVSVersion="$Date: 2009/07/03 08:06:53 $" 
      !<--Make will hopefully replace some variables in here            
      !Do not edit the two comment lines                                
      !!make-replace                                                    
      ivers='fleur.25g' 
      compiletime='Fri Apr 16 17:05:40 CEST 2010' 
      compileplace='iff316 by wortmann' 
      !!make-replace                                                    
      !>


                                                                        
      !OPEN files and say Hello                                         
      CALL gf_out_create() 
                                                                        
      WRITE(6,*) "              This is GFLEUR" 
      WRITE(16,*)"              This is GFLEUR" 
      WRITE(6,*) ;WRITE(16,*) 
      WRITE(6,*) "Settings:" 
      WRITE(6,"(a,a50)") "CVS-Date   :",CVSVersion 
      WRITE(6,"(a,a50)") "Compiled at:",compiletime 
      WRITE(6,"(a,a50)") "         on:",compileplace 
      WRITE(6,"(a,a9)") "Based on   :",ivers 
#ifdef CPP_MPI                                                          
      WRITE(6,*) "          Parallel-version (CPP_MPI)" 
#endif                                                                  
#ifdef CPP_NOCO                                                         
      WRITE(6,*) "          Noco-version     (CPP_NOCO)" 
#ifdef CPP_SOC                                                          
      WRITE(6,*) "          SOC-version      (CPP_SOC)" 
#endif                                                                  
#else                                                                   
#ifdef CPP_SOC                                                          
       "Compile time error, No SOC without NOCO"
#endif                                                                  
#endif                                                                  
      !<--More checks for preprocessor switches                         
#ifndef CPP_APW                                                         
       "Compile time error, No GFLEUR without APW "
#endif                                                                  
#ifndef CPP_GF                                                          
       "Compile time error, No GLEUR without CPP_GF "
#endif                                                                  
#ifndef CPP_HDF                                                         
       "COMPILE time error, No GLEUR without HDF "
#endif                                                                  
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      END                                           
