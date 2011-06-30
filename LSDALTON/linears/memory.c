#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>

void runps_(char* size, int* len)
{
 int  out_fd[2];
 char cmd[80];

 if(pipe(out_fd)==-1) {fprintf(stderr,"Error opening pipe!"); exit(1);}
 
 snprintf(cmd,80,"ps -o vsz= -o rss= -p %d >&%d",(int) getpid(),out_fd[1]);

 system(cmd);

 *len = (int) read(out_fd[0],size,(ssize_t) 30);

 close(out_fd[0]); close(out_fd[1]);

}

