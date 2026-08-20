/*--------------------------------------------------------------------------------
 * Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
 * This file is part of FLEUR and available as free software under the conditions
 * of the MIT license as expressed in the LICENSE file in more detail.
 *--------------------------------------------------------------------------------*/

/* Portable process-memory queries used by m_judft_sysinfo.
 *
 *   peak_rss_bytes()    : peak (high-water-mark) resident set size, in bytes.
 *                         Uses getrusage(); works on Linux, macOS and BSD and
 *                         needs no polling - the kernel tracks the maximum.
 *   current_rss_bytes() : current resident set size, in bytes, or 0.0 if the
 *                         platform offers no cheap way to obtain it.
 *
 * Both return double (bytes) so the Fortran side needs no struct mirroring of
 * the platform-dependent `struct rusage`.
 */

#include <sys/resource.h>

#if defined(__linux__)
#include <stdio.h>
#include <unistd.h>
#elif defined(__APPLE__)
#include <mach/mach.h>
#endif

double peak_rss_bytes(void) {
    struct rusage ru;
    if (getrusage(RUSAGE_SELF, &ru) != 0) return 0.0;
#if defined(__APPLE__)
    return (double) ru.ru_maxrss;          /* macOS/BSD report bytes     */
#else
    return (double) ru.ru_maxrss * 1024.0; /* Linux reports kibibytes    */
#endif
}

double current_rss_bytes(void) {
#if defined(__linux__)
    long pages = 0, rss = 0;
    FILE *f = fopen("/proc/self/statm", "r");
    if (!f) return 0.0;
    /* field 1 = total program size (pages), field 2 = resident pages */
    if (fscanf(f, "%ld %ld", &pages, &rss) != 2) { fclose(f); return 0.0; }
    fclose(f);
    return (double) rss * (double) sysconf(_SC_PAGESIZE);
#elif defined(__APPLE__)
    mach_task_basic_info_data_t info;
    mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                  (task_info_t) &info, &count) != KERN_SUCCESS) return 0.0;
    return (double) info.resident_size;
#else
    return 0.0;
#endif
}
