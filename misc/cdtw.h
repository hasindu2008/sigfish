#include <stdlib.h>


typedef struct Path
{
  int k;
  int *px;
  int *py;
} Path;


float std_dtw(float *x, float *y, int n, int m, float *cost, int squared);
int path(float *cost, int n, int m, int startx, int starty, Path *p);
void subsequence(float *x, float *y, int n, int m, float *cost);
int subsequence_path(float *cost, int n, int m, int starty, Path *p);
