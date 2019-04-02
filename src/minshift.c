#include "minexcal.h"

int (*boundaryfill)(TYPE **,int,int)=NULL;

extern int nrow;
extern int ncol;
extern int boundary;
extern TYPE boundaryvalue;
extern int (*newindex)(int,int);

/********************************************functions*/

int NewIndexWrap(int index, int bound)
{
  while ((index-1)<0)
    index+=bound;
  return(((index-1)%bound)+1);
}

int NewIndexEcho(int index, int bound)
{
  if (index<1)
    return 1;
  if (index>bound)
    return bound;
  return index;
}

int NewIndexFixed(int index, int bound)
{
  if (index<0)
    return 0;
  if (index>bound)
    return bound+1;
  return index;
}

TYPE **Boundaries(TYPE **a)
{
  int i,j;
  int nc,nr,nr1,nc1;
  TYPE bound;
  nr = nrow;
  nc = ncol;
  nr1 = nrow + 1;
  nc1 = ncol + 1;

  if(boundary==WRAP)
    newindex=NewIndexWrap;
  else if(boundary==ECHO)
    newindex=NewIndexEcho;
  else
    newindex=NewIndexFixed;

  if (boundary == WRAP) {
    for (i=1; i <= nr; i++) {
      a[i][0] = a[i][nc];
      a[i][nc1] = a[i][1];
    }
    for (j=0; j <= nc1; j++) {
      a[0][j] = a[nr][j];
      a[nr1][j] = a[1][j];
    }
  }
  else if (boundary == FIXED) {
    if(boundaryfill!=NULL) {
      for (i=1; i <= nr; i++) {
	a[i][0] = (*boundaryfill)(a,i,0);
	a[i][nc1] = (*boundaryfill)(a,i,nc1);
      }
      for (j=0; j <= nc1; j++) {
	a[0][j] = (*boundaryfill)(a,0,j);
	a[nr1][j] = (*boundaryfill)(a,nr1,j);
      }
    }
    else {
      bound = boundaryvalue;
      for (i=1; i <= nr; i++) {
	a[i][0] = bound;
	a[i][nc1] = bound;
      }
      for (j=0; j <= nc1; j++) {
	a[0][j] = bound;
	a[nr1][j] = bound;
      }
    }
  }
  else if (boundary == ECHO) {
    for (i=1; i <= nr; i++) {
      a[i][0] = a[i][1];
      a[i][nc1] = a[i][nc];
    }
    for (j=0; j <= nc1; j++) {
      a[0][j] = a[1][j];
      a[nr1][j] = a[nr][j];
    }
  }
  return a;
}
