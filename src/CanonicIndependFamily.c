#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

typedef struct {
  int n;
  int sValue;
  int grpsize;
  int **ListGroupElement;
  mpq_t **ListVectors;
  mpq_t **SoughtMat;
  mpq_t *TheV;
  int *eVect1;
  mpq_t tq1, tq2, tq3;
  int *CoordinateOfVector;
  int *CoordinateStatus;
} TheMeta;


void ReadData(FILE *f, TheMeta *TheMET)
{
  int n, sValue, grpsize, iElt, iVect;
  int eVal, iS, iDim, i, j;
  fscanf(f, "%d", &n);
  TheMET->n=n;
  fscanf(f, "%d", &sValue);
  TheMET->sValue=sValue;
  fscanf(f, "%d", &grpsize);
  TheMET->grpsize=grpsize;
  /* allocate the group */
  if ((TheMET->ListGroupElement = (int**)malloc(grpsize*sizeof(int*))) == 0)
    exit (EXIT_FAILURE);
  for (iElt=1; iElt<=grpsize; iElt++)
    if ((TheMET->ListGroupElement[iElt-1] = (int*)malloc(sValue*sizeof(int))) == 0)
      exit (EXIT_FAILURE);
  /* allocate the vectors */
  if ((TheMET->ListVectors = (mpq_t**)malloc(sValue*sizeof(mpq_t*))) == 0)
    exit (EXIT_FAILURE);
  for (iVect=1; iVect<=sValue; iVect++)
    if ((TheMET->ListVectors[iVect-1] = (mpq_t*)malloc(n*sizeof(mpq_t))) == 0)
      exit (EXIT_FAILURE);
  if ((TheMET->eVect1 = (int*)malloc(sValue*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((TheMET->TheV = (mpq_t*)malloc(n*sizeof(mpq_t))) == 0)
    exit (EXIT_FAILURE);
  for (i=1; i<=n; i++)
    mpq_init(TheMET->TheV[i-1]);
  mpq_init(TheMET->tq1);
  mpq_init(TheMET->tq2);
  mpq_init(TheMET->tq3);
  if ((TheMET->SoughtMat = (mpq_t**)malloc(n*sizeof(mpq_t*))) == 0)
    exit (EXIT_FAILURE);
  for (i=1; i<=n; i++)
    if ((TheMET->SoughtMat[i-1] = (mpq_t*)malloc(n*sizeof(mpq_t))) == 0)
      exit (EXIT_FAILURE);
  for (i=1; i<=n; i++)
    for (j=1; j<=n; j++)
      mpq_init(TheMET->SoughtMat[i-1][j-1]);
  if ((TheMET->CoordinateOfVector = (int*)malloc(n*sizeof(int))) == 0)
    exit (EXIT_FAILURE);
  if ((TheMET->CoordinateStatus = (int*)malloc(n*sizeof(int))) == 0)
    exit (EXIT_FAILURE);

  /* now reading */
  for (iElt=1; iElt<=grpsize; iElt++)
    {
      for (iS=1; iS<=sValue; iS++)
	{
	  fscanf(f, "%d", &eVal);
	  TheMET->ListGroupElement[iElt-1][iS-1]=eVal;
	}
    }
  for (iVect=1; iVect<=sValue; iVect++)
    {
      for (iDim=1; iDim<=n; iDim++)
	{
	  mpq_init(TheMET->ListVectors[iVect-1][iDim-1]);
	  mpq_inp_str(TheMET->ListVectors[iVect-1][iDim-1], f, 0);
	}
    }
}


void NormalizeMatrix(int *eVectINP, TheMeta *TheMET)
{
  int iVect, iRelCol, IsFinished, test, iCol;
  int TheRank, iS;
  TheRank=0;
  for (iCol=1; iCol<=TheMET->n; iCol++)
    TheMET->CoordinateStatus[iCol-1]=1;
  for (iS=1; iS<=TheMET->sValue; iS++)
    if (eVectINP[iS-1] == 1)
      {
	for (iCol=1; iCol<=TheMET->n; iCol++)
	  mpq_set(TheMET->TheV[iCol-1], TheMET->ListVectors[iS-1][iCol-1]);
	for (iVect=1; iVect<=TheRank; iVect++)
	  {
	    iRelCol=TheMET->CoordinateOfVector[iVect-1];
	    for (iCol=1; iCol<=TheMET->n; iCol++)
	      if (TheMET->CoordinateStatus[iCol-1] == 1)
		{
		  mpq_mul(TheMET->tq1, TheMET->TheV[iRelCol-1], TheMET->SoughtMat[iVect-1][iCol-1]);
		  mpq_sub(TheMET->TheV[iCol-1], TheMET->TheV[iCol-1], TheMET->tq1);
		}
	  }
	IsFinished=1;
	iRelCol=-47;
	for (iCol=1; iCol<=TheMET->n; iCol++)
	  if (TheMET->CoordinateStatus[iCol-1] == 1)
	    {
	      test=mpq_sgn(TheMET->TheV[iCol-1]);
	      if (test != 0)
		{
		  IsFinished=0;
		  iRelCol=iCol;
		  mpq_set(TheMET->tq3, TheMET->TheV[iRelCol-1]);
		}
	    }
	if (IsFinished == 1)
	  {
	    fprintf(stderr, "The vector configuration is not independent\n");
	    fprintf(stderr, "This is forbidden, please debug\n");
	    exit(1);
	  }
	for (iCol=1; iCol<=TheMET->n; iCol++)
	  if (TheMET->CoordinateStatus[iCol-1] == 0)
	    mpq_set_ui(TheMET->TheV[iCol-1], 0, 1);
	for (iCol=1; iCol<=TheMET->n; iCol++)
	  mpq_div(TheMET->TheV[iCol-1], TheMET->TheV[iCol-1], TheMET->tq3);
	for (iVect=1; iVect<=TheRank; iVect++)
	  {
	    mpq_set(TheMET->tq3, TheMET->SoughtMat[iVect-1][iRelCol-1]);
	    for (iCol=1; iCol<=TheMET->n; iCol++)
	      {
		mpq_mul(TheMET->tq1, TheMET->tq3, TheMET->TheV[iCol-1]);
		mpq_sub(TheMET->SoughtMat[iVect-1][iCol-1], TheMET->SoughtMat[iVect-1][iCol-1], TheMET->tq1);
	      }
	  }
	TheRank=TheRank+1;
	for (iCol=1; iCol<=TheMET->n; iCol++)
	  mpq_set(TheMET->SoughtMat[TheRank-1][iCol-1], TheMET->TheV[iCol-1]);
	TheMET->CoordinateOfVector[TheRank-1]=iRelCol;
	TheMET->CoordinateStatus[iRelCol-1]=0;
      }
}


int IsRankCorrect(TheMeta *TheMET, int iS, int siz)
{
  int iVect, iRelCol, IsFinished, test, iCol;
  for (iCol=1; iCol<=TheMET->n; iCol++)
    mpq_set(TheMET->TheV[iCol-1], TheMET->ListVectors[iS-1][iCol-1]);
  for (iVect=1; iVect<=siz; iVect++)
    {
      iRelCol=TheMET->CoordinateOfVector[iVect-1];
      for (iCol=1; iCol<=TheMET->n; iCol++)
	if (TheMET->CoordinateStatus[iCol-1] == 1)
	  {
	    mpq_mul(TheMET->tq1, TheMET->TheV[iRelCol-1], TheMET->SoughtMat[iVect-1][iCol-1]);
	    mpq_sub(TheMET->TheV[iCol-1], TheMET->TheV[iCol-1], TheMET->tq1);
	  }
    }
  IsFinished=1;
  iRelCol=-47;
  for (iCol=1; iCol<=TheMET->n; iCol++)
    if (TheMET->CoordinateStatus[iCol-1] == 1)
      {
	test=mpq_sgn(TheMET->TheV[iCol-1]);
	if (test != 0)
	  {
	    IsFinished=0;
	    iRelCol=iCol;
	    mpq_set(TheMET->tq3, TheMET->TheV[iRelCol-1]);
	  }
      }
  if (IsFinished == 1)
    return 0;
  else
    return 1;
}






void ActionVector(int *eVectImg, int *eVectINP, TheMeta *TheMET, int idxgrp)
{
  int iS, iSnew;
  for (iS=1; iS<=TheMET->sValue; iS++)
    {
      iSnew=TheMET->ListGroupElement[idxgrp-1][iS-1];
      eVectImg[iSnew-1]=eVectINP[iS-1];
    }
}


/* test if eVect1 is strictly smaller than eVect2 */
int IsStrictlySmaller(int *eVect1, int *eVect2, TheMeta *TheMET)
{
  int iDim;
  for (iDim=1; iDim<=TheMET->sValue; iDim++)
    {
      if (eVect1[iDim-1] == 1 && eVect2[iDim-1] == 0)
	return 1;
      if (eVect1[iDim-1] == 0 && eVect2[iDim-1] == 1)
	return 0;
    }
  /* they are equal */
  return 0;
}

int IsSmallestInOrbit(int *eVectINP, TheMeta *TheMET)
{
  int iElt, test;
  for (iElt=1; iElt<=TheMET->grpsize; iElt++)
    {
      ActionVector(TheMET->eVect1, eVectINP, TheMET, iElt);
      test=IsStrictlySmaller(TheMET->eVect1, eVectINP, TheMET);
      if (test == 1)
	return 0;
    }
  return 1;
}

  
void ReadOneVector(int *eVect, TheMeta *TheMET, FILE *f, int *maxval, int siz)
{
  int iDim, eVal;
  for (iDim=1; iDim<=TheMET->sValue; iDim++)
    eVect[iDim-1]=0;
  for (iDim=1; iDim<=siz; iDim++)
    {
      fscanf(f, "%d", &eVal);
      eVect[eVal-1]=1;
    }
  *maxval=eVal;
}

void PrintOneVector(int *eVect, TheMeta *TheMET, FILE *f)
{
  int iDim;
  for (iDim=1; iDim<=TheMET->sValue; iDim++)
    if (eVect[iDim-1] == 1)
      fprintf(f, " %d", iDim);
  fprintf(f, "\n");
}


void FullEnum(FILE *INPUTNB, FILE *INPUT, FILE *OUTPUTNB, FILE *OUTPUT, TheMeta *TheMET)
{
  int *eVectINP;
  int maxval;
  int nbVectorINPUT, nbVectorOUTPUT;
  int iVal, jVal, iVector;
  int siz;
  fscanf(INPUTNB, "%d", &siz);
  fscanf(INPUTNB, "%d", &nbVectorINPUT);
  nbVectorOUTPUT=0;
  if ((eVectINP = (int*)malloc(TheMET->sValue*sizeof(int))) == 0)
    exit (EXIT_FAILURE);  
  for (iVector=1; iVector<=nbVectorINPUT; iVector++)
    {
      ReadOneVector(eVectINP, TheMET, INPUT, &maxval, siz);
      NormalizeMatrix(eVectINP, TheMET);
      for (iVal=maxval+1; iVal<=TheMET->sValue; iVal++)
	{
	  for (jVal=maxval+1; jVal<=TheMET->sValue; jVal++)
	    eVectINP[jVal-1]=0;
	  eVectINP[iVal-1]=1;
	  if (IsSmallestInOrbit(eVectINP, TheMET) == 1)
	    if (IsRankCorrect(TheMET, iVal, siz) == 1)
	      {
		nbVectorOUTPUT=nbVectorOUTPUT+1;
		PrintOneVector(eVectINP, TheMET, OUTPUT);
	      }
	}
    }
  fprintf(OUTPUTNB, "%d\n", siz+1);
  fprintf(OUTPUTNB, "%d\n", nbVectorOUTPUT);
}


int main(int argc, char *argv[])
{

  FILE *PermanentDATA=NULL;
  FILE *INPUTNB=NULL;
  FILE *OUTPUTNB=NULL;
  FILE *INPUT=NULL;
  FILE *OUTPUT=NULL;
  TheMeta TheMET;
  if (argc != 6)
    {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "CanonicIndepFamily [PermanentDATA] [INPUTNB] [INPUT] [OUTPUTNB] [OUTPUT]\n");
      return -1;
    }
  PermanentDATA=fopen(argv[1], "r");
  if (PermanentDATA == NULL)
    {
      fprintf(stderr,"The file %s not found\n", argv[1]);
      return -1;
    }
  INPUTNB=fopen(argv[2], "r");
  if (INPUTNB == NULL)
    {
      fprintf(stderr,"The file %s not found\n", argv[2]);
      return -1;
    }
  INPUT=fopen(argv[3], "r");
  if (INPUT == NULL)
    {
      fprintf(stderr,"The file %s not found\n", argv[3]);
      return -1;
    }
  OUTPUTNB=fopen(argv[4], "w");
  if (OUTPUTNB == NULL)
    {
      fprintf(stderr,"The file %s could not be created\n", argv[4]);
      return -1;
    }
  OUTPUT=fopen(argv[5], "w");
  if (OUTPUT == NULL)
    {
      fprintf(stderr,"The file %s could not be created\n", argv[5]);
      return -1;
    }
  ReadData(PermanentDATA, &TheMET);
  FullEnum(INPUTNB, INPUT, OUTPUTNB, OUTPUT, &TheMET);
  fclose(INPUTNB);
  fclose(INPUT);
  fclose(OUTPUTNB);
  fclose(OUTPUT);
  return 1;
}
