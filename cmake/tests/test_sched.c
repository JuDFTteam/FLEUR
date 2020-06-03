#include <sched.h>
int sched_getcpu();

int main()
{
  int cpu = sched_getcpu();
  return (cpu);
}
