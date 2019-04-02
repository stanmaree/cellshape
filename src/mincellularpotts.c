#include "minexcal.h"

double dhcutoff=1e-9;
double temperature=6.;
double zerotemperature=0.5;
double dissipation=0.8;
int targetarea=30;
int targetperimeter=50;
int specialcelltype=0;
int *wideweight=NULL;
int *hexwideweight=NULL;
int perimeterconstraint=0;
int hexperimeterconstraint=0;

int maxneigh=8;
int hexmaxneigh;
Position *val=NULL;
Position *hexval=NULL;
int (*newindex)(int,int)=NULL;
int (*hamiltonianmask)(int,int,int)=NULL;
int (*hexhamiltonianmask)(int,int,int)=NULL;

extern int nrow;
extern int ncol;
extern int boundary;
extern unsigned char **userCol;
extern int lastcolor;
extern int cellboundarycolor;

extern void AllocateUserCol(int);

/********************************************functions*/

void TreeFreeMemory(Tree *ptree)
{
  if(ptree->next!=NULL)
    TreeFreeMemory(ptree->next);
  free(ptree);
}

void TreeAddition(Tree **pptree,int cell,int number)
{
  Tree *tmp;

  if(*pptree==NULL) {
    if((*pptree=(Tree *)malloc(sizeof(Tree)))==NULL) {
      fprintf(stderr,"TreeAddition: error in memory allocation\n");
      exit(EXIT_FAILURE);
    }
    (**pptree).cell=cell;
    (**pptree).count=number;
    (**pptree).next=NULL;
  }
  else if ((**pptree).cell < cell)
    TreeAddition(&(**pptree).next,cell,number);
  else if ((**pptree).cell == cell)
    (**pptree).count+=number;
  else {
    tmp=*pptree;
    if((*pptree=(Tree *)malloc(sizeof(Tree)))==NULL) {
      fprintf(stderr,"TreeAddition: error in memory allocation\n");
      exit(EXIT_FAILURE);
    }
    (**pptree).cell=cell;
    (**pptree).count=number;
    (**pptree).next=tmp;
  }
}

Cell *CNew(Cell *newcell, int maxcells,int maxtypes)
{
  int copyrange,i,j;

  newcell->maxcells=maxcells;
  newcell->maxtypes=maxtypes;

  if(boundary==WRAP)
    newindex=NewIndexWrap;
  else if (boundary==ECHO)
    newindex=NewIndexEcho;
  else
    newindex=NewIndexFixed;

  // copyprobability
  newcell->chance1  = (int)floor(-dissipation-DBL_EPSILON); //chance is 1,
							//even when temp==0
  if(temperature<DBL_EPSILON) { //temp==0
    if(fabs(dissipation+(double)(newcell->chance1+1))<DBL_EPSILON) 
      //we're at border, and need to use zerotemperature
      newcell->chance0 = newcell->chance1+2;
    else //no Delta H (which take int-values) is at border
      newcell->chance0 = newcell->chance1+1;
  }
  else
    newcell->chance0 = (int)ceil(-log(dhcutoff)*temperature-dissipation);

  copyrange = newcell->chance0-newcell->chance1-1;//are there non-0,
						  //non-1 values?
  if ((newcell->copyprob = (double*)calloc((size_t)copyrange,sizeof(double)))==NULL) {
    fprintf(stderr,"CNew: error in memory allocation\n");
    exit(EXIT_FAILURE);
  }
  
  for ( i = 0 ; i < copyrange ; i++ ) {  
    if(temperature<DBL_EPSILON)
      newcell->copyprob[i]=zerotemperature;//we must be at the border,
					   //otherwise copyrange would
					   //be zero
    else
      newcell->copyprob[i] = exp(-(((double)(i)+newcell->chance1+1+dissipation)/temperature));
  }

  // allocate
  if((newcell->J=(int **)calloc((size_t)maxtypes,sizeof(int *)))==NULL) {
    fprintf(stderr,"CNew: error in memory allocation\n");
    exit(EXIT_FAILURE);
  }
  
  if((newcell->J[0]=(int *)calloc((size_t)(maxtypes*maxtypes),sizeof(int)))==NULL) {
    fprintf(stderr,"CNew: error in memory allocation\n");
    exit(EXIT_FAILURE);
  }
  
  for(i=1;i<maxtypes;i++) 
    newcell->J[i]=newcell->J[i-1]+maxtypes;

  for(i=0;i<maxtypes;i++) 
    for(j=0;j<maxtypes;j++) 
      newcell->J[i][j]=0;
  
  if((newcell->area=(int *)calloc((size_t)maxcells,sizeof(int)))==NULL) {
    fprintf(stderr,"CNew: error in memory allocation\n");
    exit(EXIT_FAILURE);
  }
  
  if((newcell->targetarea=(int *)calloc((size_t)maxcells,sizeof(int)))==NULL) {
    fprintf(stderr,"CNew: error in memory allocation\n");
    exit(EXIT_FAILURE);
  }
  
  if((newcell->celltype=(int *)calloc((size_t)maxcells,sizeof(int)))==NULL) {
    fprintf(stderr,"CNew: error in memory allocation\n");
    exit(EXIT_FAILURE);
  }
  
  for(i=0;i<=specialcelltype;i++) {
    newcell->area[i]=0;
    newcell->targetarea[i]=0;
    newcell->celltype[i]=i;
  }
  
  if(lastcolor<cellboundarycolor)
    AllocateUserCol(cellboundarycolor);
  if(!userCol[cellboundarycolor][0])
    ColorTable(cellboundarycolor,cellboundarycolor,BLACK);
  if(perimeterconstraint || hexperimeterconstraint)
    InitCellPosition(newcell);
  return(newcell);
}

void UpdateCFill(TYPE** a,Cell* cell)
{
  int k;
  
  if(boundary==WRAP)
    newindex=NewIndexWrap;
  else if (boundary==ECHO)
    newindex=NewIndexEcho;
  else
    newindex=NewIndexFixed;

  for(k=0;k<cell->maxcells;k++)
    cell->area[k]=0;
  
  PLANE(
	if((int)a[i][j]>specialcelltype)
	  cell->area[a[i][j]]++;
	);
  if(perimeterconstraint)
    UpdateWideCellNeighbours(a,cell);
  else if(hexperimeterconstraint)
    UpdateWideHexCellNeighbours(a,cell);
  for(k=0;k<cell->maxcells;k++)
    if(cell->area[k]>0) { 
      if(cell->celltype[k]==0 && cell->maxtypes>1)
	cell->celltype[k]=1;
      if(cell->targetarea[k]==0)
	cell->targetarea[k]=targetarea;
      if((perimeterconstraint || hexperimeterconstraint) && cell->shape[k].targetperimeter==0)
	cell->shape[k].targetperimeter=targetperimeter;
    }
  Boundaries(a);
}

void InitCellPosition(Cell* cell)
{
  int i;
  double nr2,nc2;
  //static int warning=0;

  if(cell->shape==NULL) {
    nr2=(double)nrow/2.;
    nc2=(double)ncol/2.;
    
    if((cell->shape=(Shape *)calloc((size_t)cell->maxcells,sizeof(Shape)))==NULL) {
      fprintf(stderr,"InitCellPosition: error in memory allocation\n");
      exit(EXIT_FAILURE);
    }
    
    for(i=0;i<cell->maxcells;i++) {
      cell->shape[i].meanx=nc2;
      cell->shape[i].meany=nr2;
    }
  }
  /*
  else if(!warning) {
    fprintf(stderr,"InitCellPosition: warning: I assume that this is a repeated call to InitCellPosition.\n");
    fprintf(stderr,"InitCellPosition: warning: If not, the program will crash. In that case, change in your program\n");
    fprintf(stderr,"InitCellPosition: warning: `Cell cells;' into `static Cell cells;'.\n");
    warning=1;
  }
  */
}

void UpdateWideCellNeighbours(TYPE** a,Cell* cell)
{
  int i,j,k,x,y;
  int sigma2,sigmaneigh2;
  int nc,nr,nr1,nc1;

  if(val==NULL) {
    fprintf(stderr,"UpdateWideCellNeighbours: warning: no wide neighbourhood defined, command ignored\n");
    return;
  }
  
  if(boundary==WRAP)
    newindex=NewIndexWrap;
  else if (boundary==ECHO)
    newindex=NewIndexEcho;
  else
    newindex=NewIndexFixed;

  for(i=0;i<cell->maxcells;i++)
    if(cell->shape[i].neighbour!=NULL) {
      TreeFreeMemory(cell->shape[i].neighbour);
      cell->shape[i].neighbour=NULL;
      cell->shape[i].perimeter=0;
    }
  
  Boundaries(a);
  nr = nrow;
  nc = ncol;
  nr1 = nrow+1;
  nc1 = ncol+1;
  //--- for WRAP the boundary doesn't count: would mean counting twice
  //--- for ECHO the boundary doesn't count: the `fake' boundary does not
  //--- truly take part in the Hamiltonian
  //--- for FIXED the boundary has to be count: it is a truly separate entity
  //--- for FIXED, take care that the whole boundary is counted
  if(boundary==WRAP || boundary==ECHO) 
    for (i=1; i <= nr; i++)
      for (j=1; j <= nc; j++) {
	sigma2=(int)a[i][j];
	for(k=0;k<maxneigh;k++)
	  if(sigma2!=(sigmaneigh2=(int)a[newindex(i+val[k].yy,nr)][newindex(j+val[k].xx,nc)]))  {
	    if(hamiltonianmask!=NULL) {
	      if(wideweight[k]) {
		TreeAddition(&(cell->shape[sigma2].neighbour),sigmaneigh2,wideweight[k]);
		cell->shape[sigma2].perimeter+=wideweight[k];
	      }
	    }
	    else {
	      TreeAddition(&(cell->shape[sigma2].neighbour),sigmaneigh2,1);
	      cell->shape[sigma2].perimeter++;
	    }
	  }
      }
  
  else //--- if boundary==FIXED
    for (i=1; i <= nr; i++)
      for (j=1; j <= nc; j++) {
	sigma2=(int)a[i][j];
	for(k=0;k<maxneigh;k++)
	  if(sigma2!=(sigmaneigh2=(int)a[(y=newindex(i+val[k].yy,nr))][(x=newindex(j+val[k].xx,nc))])) {
	    if(hamiltonianmask!=NULL) {
	      if(wideweight[k]) {
		TreeAddition(&(cell->shape[sigma2].neighbour),sigmaneigh2,wideweight[k]);
		cell->shape[sigma2].perimeter+=wideweight[k];
		if(y==0 || y==nr1 || x==0 || x==nc1) {
		  TreeAddition(&(cell->shape[sigmaneigh2].neighbour),sigma2,wideweight[k]);
		  cell->shape[sigmaneigh2].perimeter+=wideweight[k];
		}
	      }
	    }
	    else {
	      TreeAddition(&(cell->shape[sigma2].neighbour),sigmaneigh2,1);
	      cell->shape[sigma2].perimeter++;
	      if(y==0 || y==nr1 || x==0 || x==nc1) {
		TreeAddition(&(cell->shape[sigmaneigh2].neighbour),sigma2,1);
		cell->shape[sigmaneigh2].perimeter++;
	      }
	    }
	  }
      }
}

void UpdateWideHexCellNeighbours(TYPE** a,Cell* cell)
{
  int i,j,k,x,y;
  int sigma2,sigmaneigh2;
  int nc,nr,nr1,nc1;
  int nneigh;

  if(hexval==NULL) {
    fprintf(stderr,"UpdateWideHexCellNeighbours: warning: no wide neighbourhood defined, command ignored\n");
    return;
  }

  if(boundary==WRAP)
    newindex=NewIndexWrap;
  else if (boundary==ECHO)
    newindex=NewIndexEcho;
  else
    newindex=NewIndexFixed;

  for(i=0;i<cell->maxcells;i++)
    if(cell->shape[i].neighbour!=NULL) {
      TreeFreeMemory(cell->shape[i].neighbour);
      cell->shape[i].neighbour=NULL;
      cell->shape[i].perimeter=0;
    }
  
  Boundaries(a);
  nr = nrow;
  nc = ncol;
  nr1 = nrow+1;
  nc1 = ncol+1;
  //--- for WRAP the boundary doesn't count: would mean counting twice
  //--- for ECHO the boundary doesn't count: the `fake' boundary does not
  //--- truly take part in the Hamiltonian
  //--- for FIXED the boundary has to be count: it is a truly separate entity
  //--- for FIXED, take care that the whole boundary is counted
  if(boundary==WRAP || boundary==ECHO) 
    for (i=1; i <= nr; i++) {
      nneigh=hexmaxneigh*(i%2);
      for (j=1; j <= nc; j++) {
	sigma2=(int)a[i][j];
	for(k=0;k<hexmaxneigh;k++)
	  if(sigma2!=(sigmaneigh2=(int)a[newindex(i+hexval[k].yy,nr)][newindex(j+hexval[k+nneigh].xx,nc)])) {
	    if(hexhamiltonianmask!=NULL) {
	      if(hexwideweight[k]) {
		TreeAddition(&(cell->shape[sigma2].neighbour),sigmaneigh2,hexwideweight[k]);
		cell->shape[sigma2].perimeter+=hexwideweight[k];
	      }
	    }
	    else {
	      TreeAddition(&(cell->shape[sigma2].neighbour),sigmaneigh2,1);
	      cell->shape[sigma2].perimeter++;
	    }
	  }
      }
    }
  
  else //--- if boundary==FIXED
    for (i=1; i <= nr; i++) {
      nneigh=hexmaxneigh*(i%2);
      for (j=1; j <= nc; j++) {
	sigma2=a[i][j];
	for(k=0;k<hexmaxneigh;k++)
	  if(sigma2!=(sigmaneigh2=(int)a[(y=newindex(i+hexval[k].yy,nr))][(x=newindex(j+hexval[k+nneigh].xx,nc))])) {
	    if(hexhamiltonianmask!=NULL) {
	      if(hexwideweight[k]) {
		TreeAddition(&(cell->shape[sigma2].neighbour),sigmaneigh2,hexwideweight[k]);
		cell->shape[sigma2].perimeter+=hexwideweight[k];
		if(y==0 || y==nr1 || x==0 || x==nc1) {
		  TreeAddition(&(cell->shape[sigmaneigh2].neighbour),sigma2,hexwideweight[k]);
		  cell->shape[sigmaneigh2].perimeter+=hexwideweight[k];
		}
	      }
	    }
	    else {
	      TreeAddition(&(cell->shape[sigma2].neighbour),sigmaneigh2,1);
	      cell->shape[sigma2].perimeter++;
	      if(y==0 || y==nr1 || x==0 || x==nc1) {
		TreeAddition(&(cell->shape[sigmaneigh2].neighbour),sigma2,1);
		cell->shape[sigmaneigh2].perimeter++;
	      }
	    }
	  }
      }
    }
}

void OneCellPosition(TYPE** a,Cell* cell,int cellnumber)
{
  int n,m,p,mmax;
  int lx,ly,hx,hy;
  int y,x,yy,xx;
  int count=0;
  int i,j,ii,jj;
  int nr,nc;

  if(cellnumber<=specialcelltype)
    return;
  if(!(count=cell->area[cellnumber]))
    return;

  nr=nrow;
  nc=ncol;
  i=(int)cell->shape[cellnumber].meany;
  j=(int)cell->shape[cellnumber].meanx;

  cell->shape[cellnumber].sumx=0.;
  cell->shape[cellnumber].sumy=0.;
  cell->shape[cellnumber].sumxx=0.;
  cell->shape[cellnumber].sumxy=0.;
  cell->shape[cellnumber].sumyy=0.;
  
  if(boundary==WRAP) {
    if(i==nrow/2 && j==ncol/2 && fabs(cell->shape[cellnumber].meany-(double)nrow/2.)<DBL_EPSILON && fabs(cell->shape[cellnumber].meanx-(double)ncol/2.)<DBL_EPSILON) {
      for(i=1,ii=nrow+1;i<ii;i++)
	for(j=1,jj=ncol+1;j<jj;j++)
	  if(a[i][j]==cellnumber)
	    goto furtherstart;
    }
  furtherstart:

    lx=-nc/2;
    hx=(nc-1)/2;
    ly=-nr/2;
    hy=(nr-1)/2;
  }
  else {
    lx=-j+1;
    hx=nc-j;
    ly=-i+1;
    hy=nr-i;
  }
  for(m=0,mmax=max(nr,nc);m<mmax;m++)
    for(n=max(-m,ly);n<=min(m,hy);n++) {
      for(p=max(-m,lx);p<=min(m,hx);p++) {
	if(abs(n)!=m && abs(p)!=m)
	  p=m-1;
	else {
	  y=NewIndexWrap((yy=i+n),nr);
	  x=NewIndexWrap((xx=j+p),nc);
	  if(a[y][x]==(TYPE)cellnumber) {
	    cell->shape[cellnumber].sumx+=(double)xx;
	    cell->shape[cellnumber].sumy+=(double)yy;
	    cell->shape[cellnumber].sumxx+=(double)(xx*xx);
	    cell->shape[cellnumber].sumxy+=(double)(xx*yy);
	    cell->shape[cellnumber].sumyy+=(double)(yy*yy);
	    count--;
	    if(count==0)
	      goto further;	
	  }
	}
      }
    }
 further:;
  if(count!=0)
    fprintf(stderr,"OneCellPosition: warning: %d pixels missing of cell %d\n",count,cellnumber);
 
  cell->shape[cellnumber].dx=cell->shape[cellnumber].sumx/(double)cell->area[cellnumber]-cell->shape[cellnumber].meanx;
  cell->shape[cellnumber].dy=cell->shape[cellnumber].sumy/(double)cell->area[cellnumber]-cell->shape[cellnumber].meany;
  cell->shape[cellnumber].meanx+=cell->shape[cellnumber].dx;
  cell->shape[cellnumber].meany+=cell->shape[cellnumber].dy;
}

void UpdateCellPosition(TYPE** a,Cell* cell)
{
  int k,kk;

  if(boundary==WRAP) {
    for(k=specialcelltype+1,kk=cell->maxcells;k<kk;k++)
      OneCellPosition(a,cell,k);
  }

  else {
    for(k=specialcelltype+1,kk=cell->maxcells;k<kk;k++)
      if(cell->area[k]) {
	cell->shape[k].sumx=0.;
	cell->shape[k].sumy=0.;
	cell->shape[k].sumxx=0.;
	cell->shape[k].sumxy=0.;
	cell->shape[k].sumyy=0.;
      }

    PLANE(
	  if((k=(int)a[i][j])>specialcelltype) {
	    cell->shape[k].sumx+=(double)j;
	    cell->shape[k].sumy+=(double)i;
	    cell->shape[k].sumxx+=(double)(j*j);
	    cell->shape[k].sumxy+=(double)(i*j);
	    cell->shape[k].sumyy+=(double)(i*i);
	  }
	  );
  
    for(k=specialcelltype+1,kk=cell->maxcells;k<kk;k++)
      if(cell->area[k]) {
	cell->shape[k].dx=cell->shape[k].sumx/(double)cell->area[k]-cell->shape[k].meanx;
	cell->shape[k].dy=cell->shape[k].sumy/(double)cell->area[k]-cell->shape[k].meany;
	cell->shape[k].meanx+=cell->shape[k].dx;
	cell->shape[k].meany+=cell->shape[k].dy;
      }
  }
}
