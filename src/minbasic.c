#include "minexcal.h"

int nrow=100;
int ncol=100;
int boundary=0;//WRAP
TYPE boundaryvalue=0;
int cellboundarycolor=255; 

Position neighbours[9]={{1,0},{1,1},{-1,1},{-1,0},{-1,-1},{1,-1},{0,1},{0,-1},{0,0}};

extern int specialcelltype;
extern int targetarea;
extern int targetperimeter;
extern int perimeterconstraint;
extern int hexperimeterconstraint;
extern double dissipation;
extern double zerotemperature;
extern double temperature;
extern double dhcutoff;

/********************************************functions*/

TYPE **NewP()
{
  TYPE **a;
  a = (TYPE **)calloc((size_t)(nrow+2), sizeof(TYPE *));
  if (a == NULL) {
    fprintf(stderr,"NewP: error in memory allocation\n");
    exit(EXIT_FAILURE);
  }
  return a;
}

TYPE **New()
{
  TYPE **a;
  int i,j;
  a = NewP(); 
  a[0] = (TYPE *)calloc((size_t)((nrow+2)*(ncol+2)),sizeof(TYPE));
  if (a[0] == NULL) {
    fprintf(stderr,"New: error in memory allocation\n");
    exit(EXIT_FAILURE);
  }
  for (i=1,j=nrow+2; i < j; i++)
    a[i] = a[i-1] + ncol + 2;
  return a;
}

TYPE **NewPlane(int row,int col)
{
  TYPE **a;
  int i,j;

  a = (TYPE **)calloc((size_t)(row+2),sizeof(TYPE *));
  if (a == NULL) {
    fprintf(stderr,"NewPlane: error in memory allocation\n");
    exit(EXIT_FAILURE);
  }
  a[0] = (TYPE *)calloc((size_t)((row+2)*(col+2)),sizeof(TYPE));
  if (a[0] == NULL) {
    fprintf(stderr,"NewPlane: error in memory allocation\n");
    exit(EXIT_FAILURE);
  }
  for (i=1,j=row+2; i < j; i++)
    a[i] = a[i-1] + col + 2;
  
  return a;
}

TYPE **Fill(TYPE **a,TYPE c)
{
  PLANE(
	a[i][j] = c;
	);
  return a;
}

TYPE **Copy(TYPE **a,TYPE **b)
{
  if(a==b)
    return a;

   PLANE(
	 a[i][j] = b[i][j];
	 );
  return a;
}
