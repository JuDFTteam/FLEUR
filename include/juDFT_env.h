/*--------------------------------------------------------------------------------
 * Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
 * This file is part of FLEUR and available as free software under the conditions
 * of the MIT license as expressed in the LICENSE file in more detail.
 *--------------------------------------------------------------------------------
 */
       !The default is to call timestart only if requested
#define CPP_juDFT_timestart(name) call timestart(name,__FILE__,__LINE__)
#define CPP_juDFT_timestart_info(name)
#define CPP_juDFT_timestart_debug(name)
#define CPP_juDFT_timestop(name) call timestop(name)
#define CPP_juDFT_timestop_info(name)
#define CPP_juDFT_timestop_debug(name)

       !If CPP_INFO is defined we call timing more frequently
#if defined(CPP_INFO)|defined(CPP_DEBUG)
#define CPP_juDFT_timestart_info(name) call timestart(name,__FILE__,__LINE__)
#define CPP_juDFT_timestop_info(name) call timestop(name)
#endif
       !In debug mode even more often
#ifdef CPP_DEBUG
#define CPP_juDFT_timestart_debug(name) call timestart(name,__FILE__,__LINE__)
#define CPP_juDFT_timestop_debug(name) call timestop(name)
#endif

#define CPP_error(message) call juDFT_error(message,file=__FILE__,line=__LINE__)

       USE m_juDFT

