#include <sched.h>
#ifndef __PGI
int sched_getcpu();

int findmycpu()
{
    int cpu = sched_getcpu();
    return cpu;
}
#endif
