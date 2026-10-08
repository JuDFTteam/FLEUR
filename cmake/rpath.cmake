#RPATH handling for the executables.
#FLEUR_LIBRARIES mostly contains plain linker flags (-L/-l), for which cmake adds
#no RPATH on its own. With 'configure.sh -rpath true' the RPATH of the build and the
#installed binaries contains the -L directories, the directories of shared libraries
#given by full path and the directories in which the shared libraries given as -lfoo
#are found (searched in the -L directories and in LIBRARY_PATH, as done by the linker).
#By default the cmake defaults are left untouched.

if (DEFINED CLI_FLEUR_USE_RPATH)
   set(FLEUR_USE_RPATH ${CLI_FLEUR_USE_RPATH})
else()
   set(FLEUR_USE_RPATH FALSE)
endif()

set(FLEUR_RPATH "")
if (FLEUR_USE_RPATH)
   #entries may also be separated by blanks if given via the environment
   string(REPLACE " " ";" rpath_items "${FLEUR_LIBRARIES}")
   set(rpath_implicit ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES} ${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}
                      ${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES} ${CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES})

   #search path of the linker for -lfoo: all -L directories first, then LIBRARY_PATH
   set(rpath_search "")
   foreach(item ${rpath_items})
      if (item MATCHES "^-L(.+)$")
         list(APPEND rpath_search "${CMAKE_MATCH_1}")
      endif()
   endforeach()
   if (DEFINED ENV{LIBRARY_PATH})
      string(REPLACE ":" ";" rpath_env "$ENV{LIBRARY_PATH}")
      foreach(dir ${rpath_env})
         get_filename_component(dir "${dir}" ABSOLUTE)
         list(APPEND rpath_search "${dir}")
         #the compilers report LIBRARY_PATH as implicit directories, these are kept
         list(REMOVE_ITEM rpath_implicit "${dir}")
      endforeach()
   endif()

   set(rpath_dirs "")
   foreach(item ${rpath_items})
      if (item MATCHES "^-L(.+)$")
         list(APPEND rpath_dirs "${CMAKE_MATCH_1}")
      elseif (item MATCHES "^-l(.+)$")
         set(rpath_lib "${CMAKE_MATCH_1}")
         if (rpath_lib MATCHES "^:(.+)$")
            #-l:filename links exactly this file
            set(rpath_names "${CMAKE_MATCH_1}")
         else()
            set(rpath_names "")
            foreach(suffix ${CMAKE_SHARED_LIBRARY_SUFFIX} ${CMAKE_FIND_LIBRARY_SUFFIXES})
               list(APPEND rpath_names "lib${rpath_lib}${suffix}")
            endforeach()
         endif()
         #the first directory containing the library is used; nothing to do for a static one
         foreach(dir ${rpath_search})
            set(rpath_found "")
            foreach(name ${rpath_names})
               if (NOT rpath_found AND EXISTS "${dir}/${name}")
                  set(rpath_found "${name}")
               endif()
            endforeach()
            if (rpath_found)
               if (NOT rpath_found MATCHES "\\.a$")
                  list(APPEND rpath_dirs "${dir}")
               endif()
               break()
            endif()
         endforeach()
      elseif (IS_ABSOLUTE "${item}" AND item MATCHES "\\.(so|dylib)(\\.[0-9.]+)?$")
         get_filename_component(rpath_dir "${item}" DIRECTORY)
         list(APPEND rpath_dirs "${rpath_dir}")
      endif()
   endforeach()

   foreach(rpath_dir ${rpath_dirs})
      if (IS_DIRECTORY "${rpath_dir}")
         get_filename_component(rpath_dir "${rpath_dir}" ABSOLUTE)
         list(FIND rpath_implicit "${rpath_dir}" rpath_index)
         if (rpath_index EQUAL -1)
            list(APPEND FLEUR_RPATH "${rpath_dir}")
         endif()
      endif()
   endforeach()
   if (FLEUR_RPATH)
      list(REMOVE_DUPLICATES FLEUR_RPATH)
   endif()

   set(CMAKE_BUILD_RPATH ${FLEUR_RPATH})
   set(CMAKE_INSTALL_RPATH ${FLEUR_RPATH})
   #also keep the directories of libraries cmake linked by full path
   set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()
