#include <silo.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
static double val(void *p,int type,int i){return type==DB_DOUBLE?((double*)p)[i]:((float*)p)[i];}
int main(int argc,char **argv){
 double cuts[]={1e-4,1e-3,1e-2};
 puts("frame,time_s,z_1e4_m,z_1e3_m,z_1e2_m,dz_m");
 for(int f=0;f<=51;f++){
 char path[4096];snprintf(path,sizeof path,"%s/LongElectrode_3D__%04d.silo",argv[1],f);
 DBfile *db=DBOpen(path,DB_UNKNOWN,DB_READ);if(!db)return 2;
 DBmultivar *mv=DBGetMultivar(db,"sigma_tot");if(!mv)return 3;
 double front[3]={NAN,NAN,NAN},dx=NAN,time=NAN;
 for(int b=0;b<mv->nvars;b++){
 DBquadvar *v=DBGetQuadvar(db,mv->varnames[b]);if(!v)return 4;
 int active=0;for(int i=0;i<v->nels;i++)if(val(v->vals[0],v->datatype,i)>=cuts[0]){active=1;break;}
 if(active){
 DBquadmesh *m=DBGetQuadmesh(db,v->meshname);if(!m)return 5;
 // Silo rectilinear coordinates; flatten using the variable's supplied strides.
 for(int k=m->min_index[2];k<m->max_index[2];k++)
 for(int j=m->min_index[1];j<m->max_index[1];j++)
 for(int i=m->min_index[0];i<m->max_index[0];i++){
 int idx=i*v->stride[0]+j*v->stride[1]+k*v->stride[2];
 double s=val(v->vals[0],v->datatype,idx);
 double z0=val(m->coords[2],m->datatype,k),z1=val(m->coords[2],m->datatype,k+1),z=.5*(z0+z1);
 for(int c=0;c<3;c++)if(s>=cuts[c]&&(!isfinite(front[c])||z>front[c])){front[c]=z;if(c==1)dx=z1-z0;}
 }
 DBFreeQuadmesh(m);
 }
 DBFreeQuadvar(v);
 }
 // Time is stored on the multimesh rather than block variables by Afivo.
 if(DBInqVarExists(db,"dtime")) DBReadVar(db,"dtime",&time);
 printf("%d,%.17g,%.17g,%.17g,%.17g,%.17g\n",f,time,front[0],front[1],front[2],dx);fflush(stdout);
 fprintf(stderr,"frame %d: z=%g blocks=%d time=%g\n",f,front[1],mv->nvars,time);
 DBFreeMultivar(mv);DBClose(db);
 }
}
