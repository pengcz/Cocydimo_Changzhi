#include <silo.h>
#include <stdio.h>
int main(int argc,char**argv){puts("frame,time_s");for(int i=0;i<=51;i++){char p[4096];snprintf(p,sizeof p,"%s/LongElectrode_3D__%04d.silo",argv[1],i);DBfile*d=DBOpen(p,DB_UNKNOWN,DB_READ);if(!d)return 1;double t;if(!DBInqVarExists(d,"dtime"))return 2;DBReadVar(d,"dtime",&t);printf("%d,%.17g\n",i,t);DBClose(d);}}
