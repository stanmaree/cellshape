#include "minexcal.h"

unsigned char fixedCol[9][3] = {{0,0,0}, {255,255,255}, {255,0,0}, {0,255,0}, {0,0,255}, \
{255,255,0}, {255,0,255}, {100,100,100}, {0,255,255}};
int lastcolor=-1;
unsigned char **userCol=NULL;

extern int cellboundarycolor;

/********************************************functions*/

void AllocateUserCol(int color)
{
  int i,j;
  static int old=-1;

  if(color<0) {
    free(userCol[0]);
    free(userCol);
    old=-1;
    return;
  }
  
  if((userCol=realloc(userCol,(size_t)((color+1)*sizeof(unsigned char*))))== NULL) {
    fprintf(stderr,"AllocateUserCol: error in memory allocation\n");
    exit(EXIT_FAILURE);
  }
  if(old==-1)
    userCol[0]=NULL;
  if((userCol[0]=realloc(userCol[0],(size_t)((color+1)*4*sizeof(unsigned char))))== NULL) {
    fprintf(stderr,"AllocateUserCol: error in memory allocation\n");
    exit(EXIT_FAILURE);
  }
  for (i=1;i<=color;i++) {
    userCol[i] = userCol[i-1] + 4;
  }
  
  if(color>old)
    for(i=old+1;i<=color;i++)
      for(j=0;j<4;j++)
	userCol[i][j]=0;
  old=color;
  lastcolor=color;
}

int ColorTable(int s,int e, ...)
{
  va_list ap;        /* points to each unnamed arg in turn */
  int i,n,cnt,col;

  /* some checks, and now hope the user supplied enough arguments */
  if ( (s>e) || (s<0) ) {
    fprintf(stderr, "ColorTable: warning: start- or end-value out-of-range (%d, %d)\n",s,e);
    return -1;
  }

  if(e>lastcolor)
    AllocateUserCol(e);
  n = e - s + 1;
  va_start(ap,e);  /* make ap point to first unnamed arg */
  
  for(cnt=0;cnt<n;cnt++) {   /* loop over 'n' userdefined colors */
    col=va_arg(ap,int);      /* retrieve value from next unnamed arg */
    if (col<0 || col>8) {
      fprintf(stderr,"ColorTable: error: color-arg %d out-of-range (%d)\n",cnt+1,col);
      exit(EXIT_FAILURE);
    }
    userCol[s+cnt][0] = 1;
    for (i=0; i<3; i++) 
      userCol[s+cnt][i+1] = fixedCol[col][i];
  }
  va_end(ap); /* clean up */
  return s;
}
