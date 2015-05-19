/* header file for LnxMP.c   */

void MP_Init(int *nproc, int *irc);
	       
void MP_Taskwait(int *taskid);

int MP_Sndsig(int *taskid);

int MP_Waitsig(int *taskid);

void MP_Killtask(int *taskid);

void MP_End();

void MP_Taskstart(int *taskid, void (*proc)(), int *nargs, ...);

void MP_Taskinit(int *taskid, void (*proc)(), int *nargs, ...);

void MP_Setstack(int stackval);

void MP_Taskbuild(int *taskid, void (*proc)(), int *nargs, ...);

void MP_Runtask(int *taskid);

void MP_Initialized(int *flag);

void prparms(int taskid);
