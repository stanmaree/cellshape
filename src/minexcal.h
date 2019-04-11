#ifndef _EXCAL_H
#define _EXCAL_H 

#ifdef __cplusplus
extern "C" {
#endif
#define VERSION 2.65
#ifdef _SMALL
  typedef short TYPE;
#else  
  typedef int TYPE;
#endif

/*********************************************general*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
/*****************************************end_general*/

/*********************************************include*/
#include <stdarg.h>
#include <math.h>
#include <float.h>
/*****************************************end_include*/
	
  typedef struct POSITION {
    int xx;
    int yy;
  } Position;

  typedef struct TREE {
    int cell;
    int count;
    struct TREE *next;
  } Tree;
  
  typedef struct SHAPE {
    int perimeter;
    int targetperimeter;
    Tree *neighbour;
    double sumx;
    double sumy;
    double sumxx;
    double sumxy;
    double sumyy;
    double meanx;
    double meany;
    double dx;
    double dy;
  } Shape;

  typedef struct CONTOUR {
    double x;
    double y;
  } Contour;

  typedef struct CELL {
    int maxcells;
    int maxtypes;
    int **J;
    int *area;
    int *targetarea;
    int *celltype;
    int chance0;
    int chance1;
    double *copyprob;
    Shape *shape;
    Contour **contour;
    int *contourlength;
    void *extras;
  } Cell;

#define PLANE(A) {int i,j,_ii,_jj;for(i=1,_ii=nrow+1;i<_ii;i++)for(j=1,_jj=ncol+1;j<_jj;j++){A};}
#define NEIGHBOURS(A) {extern Position neighbours[9];int neigh,_kk,x,y;	\
    for(neigh=0,_kk=9;neigh<_kk;neigh++){y=neighbours[neigh].yy;x=neighbours[neigh].xx;A};}
#define CELLS(C,A){int c,_c;for(c=0,_c=C.maxcells;c<_c;c++){A};}
  
/* #define RANDOM()	(rand()*(1.0/(RAND_MAX+1.0)))		*/
/* #define SEED(A)	(srand((unsigned int)A))	*/
  
#define RANDOM()	uniform()
#define SEED(A)		set_seed(A)

#define WRAP	   0
#define	FIXED	   1
#define	ECHO	   2

#define BLACK	   0

#define DEBUG 0

#define STRING   600

#define	max(a,b)	((a) > (b) ? (a) : (b))
#define	min(a,b)	((a) > (b) ? (b) : (a))
 
  /***********************************************basic*/
  TYPE **NewP();
  TYPE **New();
  TYPE **NewPlane(int,int);
  int PlaneFree(TYPE**);
  TYPE **Fill(TYPE**,TYPE);
  TYPE **Copy(TYPE**,TYPE**);

  /******************************************arithmetic*/
  TYPE **MultV(TYPE**, TYPE**, TYPE);

  /***********************************************color*/
  int ColorTable(int,int, ...);

  /***********************************************shift*/
  int NewIndexWrap(int,int);
  int NewIndexEcho(int,int);
  int NewIndexFixed(int,int);
  TYPE **Boundaries(TYPE**);
  
  /***************************************cellularpotts*/
  Cell *CNew(Cell*,int,int);
  void UpdateCFill(TYPE**,Cell*);
  void InitCellPosition(Cell*);
  void UpdateWideCellNeighbours(TYPE**,Cell*);
  void UpdateWideHexCellNeighbours(TYPE**,Cell*);
  void OneCellPosition(TYPE**,Cell*,int);
  void UpdateCellPosition(TYPE**,Cell*);

  /*********************************************locoefa*/
  void InitLOCOEFAContour(Cell*,int,int);
  void UpdateCellContour(TYPE**,TYPE**,Cell*);

  /*****************************************************/

#ifdef __cplusplus
}
#endif
#endif
