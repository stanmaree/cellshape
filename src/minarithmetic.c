#include "minexcal.h"

extern int nrow;
extern int ncol;

/********************************************functions*/

TYPE **MultV(TYPE **a,TYPE  **b,TYPE  c)
{
  PLANE(
        a[i][j] = b[i][j] * c;
        );
  return a;
}
