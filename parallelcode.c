#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<omp.h>
#include<time.h>
FILE *ptr_file_1, *ptr_file_2;	

#define TRUE 1
#define FALSE 0
#define Match 3
#define MissMatch -3
#define Gap -2
#define GapExt 0
	
int inputI;
int StrLen1,StrLen2;
int intcheck = TRUE;

char holder, ch;
int filelen1 = 0;
int filelen2 = 0;
int i,j,k,l,m,n,lenA,lenB,compval;
char dash = '-';

char FASTA1[5000];			
char FASTA2[5000];			
int HiScore=0;			
int HiScorePos[2];		
int SWArray[5001][5001];	

char MaxA[5000];
char MaxB[5000];
char OptA[5000];
char OptB[5000];

int MaxAcounter = 1;	
int MaxBcounter = 1;	
int cont = TRUE;
int check;



int Align(int PosA,int PosB)
{
	
	int relmax = -1;		
	int relmaxpos[2];
	if(SWArray[PosA-1][PosB-1] == 0)
	{
		cont = FALSE;
	}
	while(cont == TRUE)
	{	
		
		for(i=PosA; i>0; --i)
		{
			if(relmax < SWArray[i-1][PosB-1])
			{
				relmax = SWArray[i-1][PosB-1];
				relmaxpos[0]=i-1;
				relmaxpos[1]=PosB-1;
			}
		}
	
		for(j=PosB; j>0; --j)
		{

			if(relmax < SWArray[PosA-1][j-1])
			{
				relmax = SWArray[PosA-1][j-1];
				relmaxpos[0]=PosA-1;
				relmaxpos[1]=j-1;
			}
		}

		
		if((relmaxpos[0] == PosA-1) && (relmaxpos[1] == PosB-1))
		{

			MaxA[MaxAcounter] = FASTA1[relmaxpos[0]-1];
			++MaxAcounter;
			MaxB[MaxBcounter] = FASTA2[relmaxpos[1]-1];
			++MaxBcounter;

		}
		else
		{
			if((relmaxpos[1] == PosB-1) && (relmaxpos[0] != PosA-1))
			{	
				
				for(i=PosA-1; i>relmaxpos[0]-1; --i)
				{	
						MaxA[MaxAcounter]= FASTA1[i-1];
						++MaxAcounter;
				}
				
				for(j=PosA-1; j>relmaxpos[0]; --j)
				{
					MaxB[MaxBcounter] = dash;
					++MaxBcounter;
				}

				MaxB[MaxBcounter] = FASTA2[relmaxpos[1]-1];	
				++MaxBcounter;
			}
			if((relmaxpos[0] == PosA-1) && (relmaxpos[1] != PosB-1))
			{	
			
				for(j=PosB-1; j>relmaxpos[1]-1; --j)
				{	
					MaxB[MaxBcounter] = FASTA2[j-1];
					++MaxBcounter;
				}
			
				for(i=PosB-1; i>relmaxpos[1]; --i)
				{		
					MaxA[MaxAcounter] = dash;
					++MaxAcounter;
				}
				MaxA[MaxAcounter] = FASTA1[relmaxpos[0]-1];
				++MaxAcounter;
			}
		}
	
		Align(relmaxpos[0],relmaxpos[1]);
	}
	return(cont);
}


int main()
{
    double s,e;
	s=clock();
    
    ptr_file_1 = fopen("s2.txt", "r");
    
    if(ptr_file_1 == NULL)
    {
        printf("Error opening 's2.txt'\n");
        system("PAUSE");
        exit(1);
    }
   
    ptr_file_2 = fopen("s1.txt", "r");
   
    if(ptr_file_2 == NULL)
    {
        printf("Error opening 's1.txt'\n");
        system("PAUSE");
        exit(1);
    }
 
    fgets(FASTA1, sizeof(FASTA1), ptr_file_1);
    fgets(FASTA2, sizeof(FASTA2), ptr_file_2);
    fclose(ptr_file_1);
    fclose(ptr_file_2);
	lenA = strlen(FASTA1);
	lenB = strlen(FASTA2);

	for(i=0;i<=lenA;++i)
	{
		SWArray[0][i]=0;
	}
	for(i=0;i<=lenB;++i)
	{
		SWArray[i][0]=0;
	}

	compval = 0;
//parallelization of for loop to minimizethe time if done serially 
	#pragma omp parallel for 
	for(i = 1; i <= lenA; ++i)
	{	
	    #pragma omp parallel for 
	    for(j = 1; j <= lenB; ++j)
	    {	
			if(FASTA1[i-1] == FASTA2[j-1])
			{				
				compval = (SWArray[i-1][j-1] + Match);	

			}
			if(compval < ((SWArray[i-1][j]) + Gap))
			{	    
				compval = ((SWArray[i-1][j]) + Gap);		
			}
			

			if(compval<(SWArray[i][j-1]+Gap))
			{	
			    compval=((SWArray[i][j-1])+Gap);   
			}
		

			if(compval < 0)
			{
				compval = 0;
			}

			

			if(FASTA1[i-1] != FASTA2[j-1])
			{					

				if(compval < (SWArray[i-1][j-1] + MissMatch))
				{	
					compval = SWArray[i-1][j-1] + MissMatch;
                		}

				if(compval < ((SWArray[i-1][j]) + Gap))
				{	   

					compval = ((SWArray[i-1][j]) + Gap);	
				}
	

				if(compval < ((SWArray[i][j-1]) + Gap))
				{	

					compval = ((SWArray[i][j-1]) + Gap);    
				}
			

				if(compval < 0)
				{
					compval = 0;
				}
			}
			SWArray[i][j] = compval;
			compval = 0;
		}
	}

	for(i=0; i<=lenA; ++i)
	{
	    
	    for(j=0; j<=lenB; ++j)
	    {
	    	if(SWArray[i][j] > HiScore)
	        {
			HiScore = SWArray[i][j];
			HiScorePos[0]=i;
			HiScorePos[1]=j;
		}
	    }
	}
	
	MaxA[0] = FASTA1[HiScorePos[0]-1];
	MaxB[0] = FASTA2[HiScorePos[1]-1];

	check = Align(HiScorePos[0],HiScorePos[1]);

	
	k=0;
	for(i = strlen(MaxA)-1; i > -1; --i)
	{
		OptA[k] = MaxA[i];
		++k;
	}

	k=0;
	for(j=strlen(MaxB)-1; j > -1; --j)
	{
		OptB[k] = MaxB[j];
		++k;
	}

	printf("%s\n%s	",OptA,OptB);
	e=clock();
	printf("\n\n%lf",(double)(e-s)/(double)CLOCKS_PER_SEC);
	return(0);
}
