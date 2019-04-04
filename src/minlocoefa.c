#include "minexcal.h"

extern int nrow;
extern int ncol;
extern int (*newindex)(int,int);
extern int boundaryvalue;
extern int boundary;
extern int specialcelltype;

/********************************************functions*/

void InitLOCOEFAContour(Cell *cell,int cellnumber,int contourlength)
{
  if(cellnumber<0 && cellnumber>=cell->maxcells) {
    fprintf(stderr,"InitLOCOEFAContour: invalid cell number\n");
    exit(EXIT_FAILURE);
  }
  if(cell->contour[cellnumber]!=NULL) {
    //because of R issues
    //free(cell->contour[cellnumber]);
    cell->contour[cellnumber]=NULL;
  }
    
  if((cell->contour[cellnumber]=(Contour *)calloc((size_t)(contourlength+1),sizeof(Contour)))==NULL) {
    fprintf(stderr,"InitLOCOEFAContour: error in memory allocation\n");
    exit(EXIT_FAILURE);
  }
  cell->contourlength[cellnumber]=contourlength;
}

void UpdateCellContour(TYPE **state,TYPE **a,Cell *cell)
{
  //determine contour of all but state 0
  const int xval[8] = { 1, 1, 0,-1,-1,-1, 0, 1};
  const int yval[8] = { 0, 1, 1, 1, 0,-1,-1,-1};
  int removed,connections;
  int ii,jj,c,perimeter;
  int direction,connected;
  int next;
  int x,y;
  int tmpboundary,tmpboundaryvalue,tmpspecialcelltype;
  int (*tmpnewindex)(int,int);
  
  for(c=0;c<cell->maxcells;c++) {
    if(cell->contour[c]!=NULL) {
      //because of R issues
      //free(cell->contour[c]);
      cell->contour[c]=NULL;
      cell->contourlength[c]=0;
    }
  }

  tmpboundary=boundary;
  tmpboundaryvalue=boundaryvalue;
  tmpspecialcelltype=specialcelltype;
  tmpnewindex=newindex;
  
  boundary=FIXED;
  boundaryvalue=0;
  specialcelltype=0;
  newindex=NewIndexFixed;   

  Fill(state,0);
  Boundaries(a);
  Boundaries(state);
  Copy(state,a);
  MultV(state,state,-1);
  
  //we remove loose pixels but do not check for global connectivity, as SPM guarantees this
  do {
    removed=0;
    PLANE(
	  if(state[i][j]) {
	    connections=0;
	    NEIGHBOURS(
		       if(state[i+y][j+x]==state[i][j])
			 connections++;
		       );
	    if(connections<3) {
	      state[i][j]=0;
	      a[i][j]=0;
	      removed++;
	    }
	  }
	  );
    if(DEBUG)
      fprintf(stdout,"UpdateCellContour: %d weakly connected pixels removed\n",removed);
  } while(removed);

  UpdateCFill(a,cell);

  PLANE(
	if(state[i][j]<0) {
	  if(cell->contour[-state[i][j]]==NULL) {
	    c=state[i][j];
	    perimeter=1;
	    ii=i;
	    jj=j;
	    state[ii][jj]=perimeter;
	    direction=0;
	    connected=0;
	    do {
	      next=0;
	      do {
		x=xval[direction];
		y=yval[direction];
		if(state[newindex(ii+y,nrow)][newindex(jj+x,ncol)]==c || (state[newindex(ii+y,nrow)][newindex(jj+x,ncol)]>0 && a[newindex(ii+y,nrow)][newindex(jj+x,ncol)]==-c)) {
		  next=1;
		  ii=ii+y;
		  jj=jj+x;
		  if(state[ii][jj]<0) {
		    state[ii][jj]=++perimeter;
		  }
		  else {
		    connected=1;
		    InitLOCOEFAContour(cell,-c,perimeter);
		  }
		  direction=(direction+5)%8;
		}     
		else direction=(direction+1)%8;
	      } while(!next);
	    } while(!connected);
	  }
	  else {
	    state[i][j]=0;
	  }
	}
	);
  PLANE(
	if(state[i][j]) {
	  cell->contour[a[i][j]][state[i][j]-1].y=i;
	  cell->contour[a[i][j]][state[i][j]-1].x=j;
	  if(state[i][j]==1) {
	    cell->contour[a[i][j]][cell->contourlength[a[i][j]]].y=i;
	    cell->contour[a[i][j]][cell->contourlength[a[i][j]]].x=j;
	  }
	}
	);
  
  boundary=tmpboundary;
  boundaryvalue=tmpboundaryvalue;
  specialcelltype=tmpspecialcelltype;
  newindex=tmpnewindex;
}
