//libraries
#include "minexcal.h"

#define MYMAXDEPTH_ 10000
static int myoverflow=0;
static int mydepth;

typedef struct extras {
  int timepoint;
  int oldsigma;
  int minx;
  int maxx;
  int miny;
  int maxy;
  int width;
  int height;
} Extras;

#define EXTRAS ((Extras*)(originalcells.extras))
#define MAXOUTPUT 100000

static TYPE **image;
static TYPE **state;
static TYPE **outline;

static int totalcelltypes=2;
static Cell cells;
static Cell originalcells;
static Cell reconstructedcells;
static int totalcells=0;
static int maxcells=0;
static int ncells=0;
static int biggestCellOnly=1;
static int excludeEdgeCells=1;

//general
extern int (*newindex)(int,int);
extern TYPE boundaryvalue;
extern int nrow;
extern int ncol;
extern int specialcelltype;
extern int boundary;

/********************************************functions*/

void MyInitLOCOEFA(Cell *cell)
{
  //note that this MyInitLOCOEFA should never be called a second time for the same Cell struct
  if((cell->contour=(Contour **)calloc((size_t)cell->maxcells,sizeof(Contour*)))==NULL) {
    fprintf(stderr,"InitLOCOEFA: error in memory allocation\n");
    exit(EXIT_FAILURE);
  }
  
  if((cell->contourlength=(int *)calloc((size_t)cell->maxcells,sizeof(int)))==NULL) {
    fprintf(stderr,"InitLOCOEFA: error in memory allocation\n");
    exit(EXIT_FAILURE);
  }
}

void MyFindCell(TYPE **cellplane,TYPE **typeplane,TYPE **b,int i,int j,int c,int new)
{
  int ii,jj;

  if(boundary==WRAP)
    newindex=NewIndexWrap;
  else if (boundary==ECHO)
    newindex=NewIndexEcho;
  else
    newindex=NewIndexFixed;

  mydepth=(!new)*mydepth+1;
  if(mydepth<MYMAXDEPTH_) {
    b[i][j]=0;
    jj=newindex(j+1,ncol);
    if(b[i][jj] && typeplane[i][j]==typeplane[i][jj])
      MyFindCell(cellplane,typeplane,b,i,jj,c,0);
    ii=newindex(i+1,nrow);
    if(b[ii][j] && typeplane[i][j]==typeplane[ii][j])
      MyFindCell(cellplane,typeplane,b,ii,j,c,0);
    jj=newindex(j-1,ncol);
    if(b[i][jj] && typeplane[i][j]==typeplane[i][jj])
      MyFindCell(cellplane,typeplane,b,i,jj,c,0);
    ii=newindex(i-1,nrow);
    if(b[ii][j] && typeplane[i][j]==typeplane[ii][j])
      MyFindCell(cellplane,typeplane,b,ii,j,c,0);
    cellplane[i][j]=c;
  }
  else
    myoverflow=1;
}

TYPE **MyType2Cell(TYPE **cellplane,TYPE **typeplane,int destcell)
{
  TYPE **b;
  TYPE tmpboundaryvalue;
  int c;
 
  b=New();
  c=destcell;
  Fill(b,1);
  tmpboundaryvalue=boundaryvalue;
  boundaryvalue=0;
  Boundaries(b);
  Boundaries(typeplane);
  boundaryvalue=tmpboundaryvalue;
  PLANE(
        if(b[i][j]) {
          if(typeplane[i][j]<=specialcelltype) {
	cellplane[i][j]=typeplane[i][j];
	b[i][j]=0;
          }
          else {
	MyFindCell(cellplane,typeplane,b,i,j,c,1);
	if(myoverflow) {
	  do {
	    myoverflow=0;
	    PLANE(
	          if(cellplane[i][j]==c && !b[i][j]) {
		cellplane[i][j]=typeplane[i][j];
		MyFindCell(cellplane,typeplane,b,i,j,c,1);
	          }
	          );
	  } while(myoverflow);
	}
	c++;
	  }
	}
        );
  PlaneFree(b);
  return cellplane;
}

void RemoveProtrutions2D(TYPE **target, TYPE **tmpstate)
{
  int nswitched;
  int n,nself,delta,self;
  int round=0;

  int neighid[9];
  int neighsize[9];

  //fprintf(stdout,"RemoveProtrutions2D\n");
  nswitched=0;

  delta=1;
  while(delta) {
    delta=0;
    round++;
    PLANE(tmpstate[i][j]=target[i][j];);

    PLANE(
	  self=target[i][j];
	  nself=0;
	  NEIGHBOURS(
		     if(target[i+y][j+x]==self) {
		       nself++;
		       if(nself==5)
			 break;
		     }
		     );
	  if(nself==5)
	    continue;

	  for(n=1;n<=8;n++) {
	    neighid[n]=0;
	    neighsize[n]=0;
	  }

	  NEIGHBOURS(
		     if(!(i+y<1 || i+y>nrow || j+x<1 || j+x>ncol)) {
		       if(target[i+y][j+x]!=self) {
			 for(n=1;n<=8;n++) {
			   if(!neighid[n] || neighid[n]==target[i+y][j+x]) {
			     neighid[n]=target[i+y][j+x];
			     neighsize[n]++;
			     break;
			   }
			 }
		       }
		     }
		     );

	  for(n=1;n<=8;n++) {
	    if(neighsize[n]>=5) {
	      target[i][j]=neighid[n];
	      delta++;
	    }
	    else if(!neighid[n])
	      break;
	  }
	  );

    nswitched+=delta;
    //fprintf(stdout,"\tin round %d, %d pixels switched\n",round,delta);
  }
  //fprintf(stdout,"\tin total, %d pixels switched\n",nswitched);
  //fprintf(stdout,"\tcomplete\n\n");
}

int DetermineTotalNoOfCells(int *x,int *n)
{
  int c;
  TYPE **totaldifferentcells;
  int cellnumber=0;
  
  nrow=n[0];
  ncol=n[1];  

  image=New();
  state=New();
  outline=New();

  PLANE(
	//hack to deal with black on white images, this should be done differently
	image[nrow+1-i][j]=x[(j-1)*nrow+i-1]+1;
	);

  RemoveProtrutions2D(image,state);
  MyType2Cell(state,image,1);
  Copy(image,state);
  //determine ncells
  ncells=0;
  PLANE(ncells=max(ncells,image[i][j]););
  
  maxcells=ncells+1;
  //fprintf(stdout,"max cells in image(s): %d\n",maxcells);
  ncells=0;

  //what is highest number of cells?
  totaldifferentcells=NewPlane(1,maxcells-1);
  //register each cell that exists
  PLANE(
	if(image[i][j]>specialcelltype)
	  totaldifferentcells[1][image[i][j]]=1;
	);
  
  for(c=1;c<maxcells;c++)
    if(totaldifferentcells[1][c]) {
      cellnumber++;
    }
  
  //allocate cells
  CNew(&cells,maxcells,totalcelltypes);
  InitCellPosition(&cells);
  //becuase of R issue
  cells.shape->neighbour=NULL;
  MyInitLOCOEFA(&cells);

  return cellnumber;
}

void UpdateCells()
{
  //fprintf(stdout,"Initiating cells...\n");
  UpdateCellPosition(state,&cells);
  
  ncells=0;
  CELLS(cells,
	if(c>specialcelltype && cells.area[c])
	  ncells=c;
	);
  //fprintf(stdout,"\t...done, %d cells\n",ncells);
}

void SelectCells()
{
  int largestcell=0;
  int largestarea=0;
  //fprintf(stdout,"Selecting cells...\n");
  
  CELLS(cells,
	if(c<=specialcelltype) {
	  cells.contourlength[c]=0;
	}
	);

  if(excludeEdgeCells) {
    PLANE(
	  if(state[i][j]>specialcelltype && (i==1 || i==nrow || j==1 || j==ncol))
	    cells.contourlength[state[i][j]]=0;
	  );
  }
  
  if(biggestCellOnly) {
    CELLS(cells,
	  if(cells.contourlength[c]) {
	    if(cells.area[c]>largestarea) {
	      largestcell=c;
	      largestarea=cells.area[c];
	    }
	  }
	  );
    CELLS(cells,
	  if(c!=largestcell) {
	    cells.contourlength[c]=0;
	  }
	  );
  }
  
  ncells=0;
  CELLS(cells,
	if(cells.contourlength[c]) {
	  ncells++;
	}
	);
  //fprintf(stdout,"\t...done, %d cells selected in SelectCells.\n",ncells);
}

void SaveCellOutlines(int *out)
{
  int k,newsigma;
  
  CELLS(cells, 
	if(cells.contourlength[c]) {
	  //find last added cell with that number
	  newsigma=totalcells;
	  while(EXTRAS[newsigma].oldsigma!=c)
	    newsigma--;
	  if(originalcells.contourlength[newsigma]>=MAXOUTPUT)
	    originalcells.contourlength[newsigma]=MAXOUTPUT-1;
	  for(k=0;k<=originalcells.contourlength[newsigma];k++) {
	    out[k]=originalcells.contourlength[newsigma]+1;
	    out[k+originalcells.contourlength[newsigma]+1]=originalcells.contour[newsigma][k].x;
	    out[k+2*(originalcells.contourlength[newsigma]+1)]=originalcells.contour[newsigma][k].y;
	    out[k+3*(originalcells.contourlength[newsigma]+1)]=cells.area[c];
	  }
	}
	);
}

void AssignPerimeter(int celltype,int timepoint,int oldsigma)
{
  int i;

  totalcells++;
  originalcells.celltype[totalcells]=celltype;
  if(originalcells.celltype[totalcells]>=originalcells.maxtypes) {
    fprintf(stderr,"celltype larger or equal than totalcelltypes (%d>=%d). Exitting now!\n",originalcells.celltype[totalcells],originalcells.maxtypes);
    exit(EXIT_FAILURE);
  }

  EXTRAS[totalcells].timepoint=timepoint;
  EXTRAS[totalcells].oldsigma=oldsigma;
  originalcells.area[totalcells]=cells.area[oldsigma];

  InitLOCOEFAContour(&originalcells,totalcells,cells.contourlength[oldsigma]);
  
  for(i=0;i<=cells.contourlength[oldsigma];i++) {
    originalcells.contour[totalcells][i].x=cells.contour[oldsigma][i].x;
    originalcells.contour[totalcells][i].y=cells.contour[oldsigma][i].y;

    if(i==0) {//first datapoint
      EXTRAS[totalcells].minx=originalcells.contour[totalcells][i].x;
      EXTRAS[totalcells].maxx=originalcells.contour[totalcells][i].x;
      EXTRAS[totalcells].miny=originalcells.contour[totalcells][i].y;
      EXTRAS[totalcells].maxy=originalcells.contour[totalcells][i].y;
    }
    else {
      EXTRAS[totalcells].minx=min(EXTRAS[totalcells].minx,originalcells.contour[totalcells][i].x);
      EXTRAS[totalcells].maxx=max(EXTRAS[totalcells].maxx,originalcells.contour[totalcells][i].x);
      EXTRAS[totalcells].miny=min(EXTRAS[totalcells].miny,originalcells.contour[totalcells][i].y);
      EXTRAS[totalcells].maxy=max(EXTRAS[totalcells].maxy,originalcells.contour[totalcells][i].y);
    }
  }

  EXTRAS[totalcells].width=EXTRAS[totalcells].maxx-EXTRAS[totalcells].minx+1;
  EXTRAS[totalcells].height=EXTRAS[totalcells].maxy-EXTRAS[totalcells].miny+1;
}

void RunOutline2D(int *out)
{
  PLANE(
	state[i][j]=image[i][j];
	outline[i][j]=0;
	);

  UpdateCFill(state,&cells);
  UpdateCellContour(outline,state,&cells);
  /*
  CELLS(cells,
	if(cells.contourlength[c])
	  fprintf(stdout,"cell %d length %d\n",c,cells.contourlength[c]);
	);
  */
  
  UpdateCells();
  SelectCells();
  
  CELLS(cells, 
	if(cells.contourlength[c]) {
	  AssignPerimeter(1 /*celltype, currently always 1*/,1 /*timepoint*/,c /*sigma*/);
	}
	);
  
  if(ncells) {
    SaveCellOutlines(out);
  }
  else {
    printf("No cells selected\n\tMost likely either the matrix containing the image is empty\n\t(please check this matrix and how you assign an image to this matrix),\n\tor the shape to be analysed is touching the edge of the image\n\t(please check the image you are analysing).\n\tThis code will now kill your R session,\n\tin order not to cause trouble downstream due to corrupted data.\n");
    exit(EXIT_FAILURE);
  }
}   

void CellContour_(int *x,int *n,int *out)
{
  //now we need to know how many cells we will have, as well as setting of nooftimepoints
  //this will be overestimation of totalcells, since each picture will contain cells that are not valid
  //we therefore keep the more serious allocation to later on, limited to cells we are really interested in

  //so first time (after DetermineTotalNoOfCells) totalcells will be bigger than second time (after RunOutline2D)

  specialcelltype=-1;
  //protection step: under R, globally initialised variables are not guaranteed to be zero
  cells.shape=NULL;
  originalcells.shape=NULL;
  reconstructedcells.shape=NULL;

  totalcells=DetermineTotalNoOfCells(x,n);
  
  CNew(&originalcells,totalcells+1,totalcelltypes); //maxsigma+1 because we have to count 0 as well
  CNew(&reconstructedcells,totalcells+1,totalcelltypes); //maxsigma+1 because we have to count 0 as well
  if((originalcells.extras=(void *)calloc((size_t)originalcells.maxcells,sizeof(Extras)))==NULL) {
    fprintf(stderr,"ReadData: error in memory allocation\n");
    exit(EXIT_FAILURE);
  }

  MyInitLOCOEFA(&originalcells);
  
  totalcells=0;
  //now really read in
  RunOutline2D(out);
} 
