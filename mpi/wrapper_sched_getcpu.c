#include <sched.h>

int sched_getcpu();

int findmycpu()
{
    int cpu = sched_getcpu();
    return cpu;
}
