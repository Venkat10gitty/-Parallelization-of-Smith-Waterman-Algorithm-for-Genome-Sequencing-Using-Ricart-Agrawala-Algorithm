#include <stdio.h> 
#include <stdlib.h>
#include <math.h> 
#include <omp.h> 
#include <time.h>
// DEFINING THE GLOBAL VARIABLES AND FIXED VARIABLES
#define RESET "\033[0m"
#define BOLDRED "\033[1m\033[31m" 
#define BOLDYELLOW "\033[1;33m"
#define GREEN   "\033[32m"
#define BLUE    "\033[34m"
#define CYAN    "\033[36m" 
#define BLACK   "\033[30m"
#define PATH -1
#define NONE 0
#define UP 1
#define LEFT 2
#define DIAGONAL 3
#define min(x, y) (((x) < (y)) ? (x) : (y)) 
#define max(a,b) ((a) > (b) ? a : b)
int c1=0,c2=0,c3=0,c4=0;
// FUNCTION DEFINITIONS
void similarityScore(long long int i, long long int j, int* H, int* P, long long int* maxPos);
int matchMissmatchScore(long long int i, long long int j); 
void backtrack(int* P, long long int maxPos);
void printMatrix(int* matrix);
void printPredecessorMatrix(int* matrix); 
void generate(void);
long long int nElement(long long int i);
void calcFirstDiagElement(long long int *i, long long int *si, long long int *sj);
// DEFINING THE MAIN VARIABLES
long long int m ; //Columns - Size of string a 
long long int n ; //Lines - Size of string b 
int matchScore = 5;
int missmatchScore = -3;
int gapScore = -4;
char *a, *b;
// MAIN FUNCTION FOR EXECUTION
int main(int argc, char* argv[]) 
{
  int thread_count = strtol(argv[1], NULL, 10); 
  m = strtoll(argv[2], NULL, 10);
  n = strtoll(argv[3], NULL, 10); 
  printf("\nMatrix[%lld][%lld]\n", n, m); 
  a = malloc(m * sizeof(char));
  b = malloc(n * sizeof(char)); 
  m++;
  n++;
  int *H;
  H = calloc(m * n, sizeof(int)); 
  int *P;
  P = calloc(m * n, sizeof(int)); 
  generate();
  long long int maxPos = 0;
  long long int i, j;
  double initialTime = omp_get_wtime(); 
  long long int si, sj, ai, aj;
  long long int nDiag = m + n - 3; 
  long long int nEle;
  #pragma omp parallel num_threads(thread_count) \
  default(none) shared(H, P, maxPos, nDiag) private(nEle, i, si, sj, ai, aj)
  {
    for (i = 1; i <= nDiag; ++i)
    {
      nEle = nElement(i); 
      calcFirstDiagElement(&i, &si, &sj); 
      #pragma omp for
      for (j = 1; j <= nEle; ++j)
      {
        ai = si - j + 1;
        aj = sj + j - 1;
        similarityScore(ai, aj, H, P, &maxPos);
      }
    }
  }
  backtrack(P, maxPos); 
  printf("\nSimilarity Matrix:\n"); 
  printMatrix(H); 
  printf(BOLDYELLOW);
  printf("\nInference Stage 1: \n");
  printf(RESET);
  printf("\nCommon bases count: A: %d, G: %d, T: %d, C: %d\n",c1,c2,c3,c4);
  int temp1=max(c1,c2);
  int temp2=max(c3,temp1);
  int temp3=max(c4,temp2);
  printf("Highest base count: %d\n",temp3);
  printf("Inference: Minimum number of common base pairs between the two sequences: %d\n\n",temp3);
  int ag=c1+c2;
  int tc=c3+c4;
  if(ag==tc)
  printf("Base pair equality satisfied for similar sequence\n");
  else
  {
    int dif=abs(ag-tc);
    float diff=((float)dif/temp3)*100;
    printf(CYAN);
    printf("Base pair equality deviates by %f percentage for the similar sequence\n",diff);
    printf(RESET);
  }
  printf("\nPredecessor Matrix:\n"); 
  printPredecessorMatrix(P);
  double finalTime = omp_get_wtime();
  printf("\nElapsed time: %f\n\n", finalTime - initialTime);
  free(H);
  free(P);
  free(a);
  free(b); 
  return 0;
}
long long int nElement(long long int i)
{ 
  if (i < m && i < n)
    return i;
  else if (i < max(m, n))
  { 
    long int min = min(m, n);
    return min - 1;
  }
  else
  {
    long int min = min(m, n);
    return 2 * min - i + abs(m - n) - 2;
  }
}
// FUNCTION TO CALCULATE THE ELEMENT FROM WHICH THE ALIGNMENT MUST START
void calcFirstDiagElement(long long int *i, long long int *si, long long int *sj)
{
  if (*i < n)
  {
    *si = *i;
    *sj = 1;
  }
  else
  {
    *si = n - 1;
    *sj = *i - n + 2;
  }
}
// FUNCTION TO CALCULATE THE SIMILARITY SCORES OF THE TWO SEQUENCES
void similarityScore(long long int i, long long int j, int* H, int* P, long long int* maxPos)
{
  int up, left, diag;
  long long int index = m * i + j; 
  up = H[index - m] + gapScore; 
  left = H[index - 1] + gapScore;
  diag = H[index - m - 1] + matchMissmatchScore(i, j); 
  int max = NONE;
  int pred = NONE;
  if (diag > max)
  { //same letter ↖
    max = diag;
    pred = DIAGONAL;
  }
  if (up > max)
  { //remove letter ↑ 
    max = up;
    pred = UP;
  }
  if (left > max)
  { //insert letter ← 
    max = left;
    pred = LEFT;
  }
  H[index] = max; 
  P[index] = pred;
  if (max > H[*maxPos])
  {
    #pragma omp critical
    *maxPos = index;
  }
}
// FUNCTION TO CALCULATE THE MISMATCH SCORES FOR THE TWO SEQUENCES
int matchMissmatchScore(long long int i, long long int j)
{ 
  if (a[j - 1] == b[i - 1])
  {
    if(i==j)
    {
      if(a[i-1]=='A')
        c1++;
        else if(a[i-1]=='G')
        c2++;
        else if(a[i-1]=='T')
        c3++;
        else 
        c4++;
    }
    return matchScore;
  }
  else
    return missmatchScore;
}
// THE MAIN BACKTRACKING FUNCTION FOR SEQUENCE ALIGNMENT
void backtrack(int* P, long long int maxPos)
{
  //hold maxPos value 
  long long int predPos; 
  do
  {
    if (P[maxPos] == DIAGONAL) 
    predPos = maxPos - m - 1;
    else if (P[maxPos] == UP) 
    predPos = maxPos - m;
    else if (P[maxPos] == LEFT) 
    predPos = maxPos - 1;
    P[maxPos] *= PATH;
    maxPos = predPos;
  }
  while (P[maxPos] != NONE);
}
// FUNCTION TO PRINT THE MATRIX
void printMatrix(int* matrix)
{ 
  long long int i, j; 
  printf("-\t-\t");
  for (j = 0; j < m-1; j++)
    printf("%c\t", a[j]);
  printf("\n-\t");
  for (i = 0; i < n; i++)
  { //Lines 
    for (j = 0; j < m; j++)
    {
      if (j==0 && i>0) printf("%c\t", b[i-1]); 
      printf("%d\t", matrix[m * i + j]);
    }
    printf("\n");
  }
}
// FUNCTION TO PRINT THE MAIN ALIGNMENT MATRIX
void printPredecessorMatrix(int* matrix)
{ 
  int flag1,flag2;
  long long int i, j, index;
  printf(" ");
  for (j = 0; j < m-1; j++)
    printf("%c ", a[j]);
  printf("\n ");
  for (i = 0; i < n; i++)
  { //Lines 
    for (j = 0; j < m; j++)
    {
      if (j==0 && i>0)
        printf("%c ", b[i-1]); 
      index = m * i + j;
      if (matrix[index] < 0)
      { 
        printf(BOLDRED),flag2++;
        if (matrix[index] == -UP) 
        printf("↑ "),flag1++;
        else if (matrix[index] == -LEFT) 
        printf("← "),flag1++;
        else if (matrix[index] == -DIAGONAL) 
        printf("↖ "),flag1++;
        else
        printf("- ");
        printf(RESET);
      }
      else
      {
        if (matrix[index] == UP) 
        printf("↑ "),flag1++;
        else if (matrix[index] == LEFT) 
        printf("← "),flag1++;
        else if (matrix[index] == DIAGONAL) 
        printf("↖ "),flag1++;
        else
        printf("- ");
      }
    }
    printf("\n");
  }
  float per = ((float)flag2/flag1)*100;
  printf(BOLDYELLOW);
  printf("\nInference Stage 2: \n");
  printf(RESET);
  printf("Possible number of similarities: %d\n",flag1);
  printf("Actual simialrity %d\n",flag2);
  printf("Percentage of similarity on a whole: %f\n\n",per);
  printf(CYAN);
  if(per>=93)
    printf("The two test organism might belong to the same species\n");
  else if(per>=80)
    printf("The two test organism might belong to the same genera\n");
  else
    printf("Proceed with characteristic comparison\n");
  printf(RESET);
}
// FUNCTION TO GENERATE THE SEQUENCES RANDOMLY
void generate()
{
  srand(0);
  long long int i;
  for (i = 0; i < m; i++)
  { 
    int aux = rand() % 4; 
    if (aux == 0)
    a[i] = 'A';
    else if (aux == 2)
    a[i] = 'C';
    else if (aux == 3)
    a[i] = 'G';
    else
    a[i] = 'T';
  }
  for (i = 0; i < n; i++)
  { 
    int aux = rand() % 4; 
    if (aux == 0)
    b[i] = 'A';
    else if (aux == 2)
    b[i] = 'C';
    else if (aux == 3)
    b[i] = 'G';
    else
    b[i] = 'T';
  }
}