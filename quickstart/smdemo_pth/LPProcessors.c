#include <string.h>
#include <stdio.h>
int MPProcessors() {
/* Return the number of processors on the host computer
   on Linux platforms
local data                                              */
   int nproc = 0;
   char *cnerr = 0;
   char cline[82];
   FILE *unit;
   unit = fopen("/proc/cpuinfo","r");
/* Quit if file does not exist */
   if (!unit)
      return nproc;
/* Read next line */
   while ((cnerr = fgets(cline,81,unit))) {
      cline[9] = '\0';
      if (!strcmp(cline,"processor"))
         nproc += 1;
   }
   fclose(unit);
   return nproc;
}
