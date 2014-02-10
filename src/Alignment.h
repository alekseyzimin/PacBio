#ifndef __ALIGNMENT_H
#define __ALIGNMENT_H
#include <iostream>
#include <string>
#include <string.h>
#include <stdint.h>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "time.h"
#include "BasicDataStructure.h"
#include "GraphConstruction.h"
using namespace std;


int score_M[5][5] = {
/*    N    A    C    G    T    */
	 -43,-44, -42, -42,  -42,  
	 -44, 91, -114, -31, -123,
	 -42,-123, 100, -125, -31, 
	 -42, -31, -125, 100, -114,
	 -42, -114, -31, -114,  91,
	
};


void str2int(char *s,int * out)
{
	int str_sz=strlen(s);
	for (int i=0;i<str_sz;++i)
	{
		char temp=s[i];
		switch(temp)
		{
		case 'A':
			out[i]=1;
			break;
		case 'C':
			out[i]=2;
			break;
		case 'G':
			out[i]=3;
			break;
		case 'T':
			out[i]=4;
			break;
		default:
			cout<<"error"<<endl;
			return;
		}
	}

}


int GlobalAlign(struct aln_t *aln_t,char *A,char *B,int match,int mismatch,int gap_cost,int band_width)
{
	int A_sz=strlen(A),B_sz=strlen(B);
	int Ai[500],Bi[500];
	memset(aln_t->score,0,sizeof(int)*(aln_t->MaxLen*aln_t->MaxLen));
	memset(aln_t->b,0,sizeof(uint8_t)*(aln_t->MaxLen*aln_t->MaxLen));
	
	str2int(A,Ai);
	str2int(B,Bi);
	for(int i=1;i<=A_sz;++i)
	{
		
		int idA=Ai[i-1];
		int upper_bound,lower_bound;
		if(B_sz>A_sz)
		{
			upper_bound=min(B_sz-A_sz+i+band_width,B_sz);
			lower_bound=max(1,i-band_width);
		}
		else
		{
			lower_bound=max(1,i-band_width-(A_sz-B_sz));
			upper_bound=min(i+band_width,B_sz);
		}
		for(int j=lower_bound;j<=upper_bound;++j)
		{
	
			int idB=Bi[j-1];
			int Sij=-1000000000,Sij_max=-100000000,from=0;
			if(j<upper_bound)
			{
			
				Sij=aln_t->score[i-1][j]+gap_cost;
				Sij_max=Sij;
				from=1;
				
			}
			
			if(j>lower_bound)
			{
				Sij=aln_t->score[i][j-1]+gap_cost;
				if(Sij>Sij_max)
				{
					Sij_max=Sij;
					from=2;
				}
			}
			if(idA==idB)
			{
				Sij=aln_t->score[i-1][j-1]+match;
			}
			else
			{
				Sij=aln_t->score[i-1][j-1]+mismatch;
			}
			if(Sij>Sij_max)
			{
				Sij_max=Sij;
				from=3;
			}
			
			aln_t->score[i][j]=Sij_max;
			aln_t->b[i][j]=from;
			
		}
		

	}
	
	/*
	ofstream o_sc("score.txt");
	ofstream o_path("path.txt");
	
	for(int i=1;i<=A_sz;++i)
	{
		for(int j=1;j<=B_sz;++j)
		{
			o_sc<<aln_t->score[i][j]<<" ";
		}
		o_sc<<endl;
	}
	for(int i=1;i<=A_sz;++i)
	{
		for(int j=1;j<=B_sz;++j)
		{
			int path=aln_t->b[i][j];
			o_path<<(path)<<" ";
		}
		o_path<<endl;
	}

	*/
	int max_score=-1000000;
	for(int i=B_sz;i>=max(i-band_width,0);--i)
	{
		if((aln_t->score[A_sz][i])>max_score)
		{max_score=aln_t->score[A_sz][i];}
	}
	return max_score;//aln_t->score[A_sz][B_sz];
	
	

}


void printAlign(struct aln_t *aln_t,char *A,char *B, string &A_aln,string &B_aln)
{
	A_aln.clear();
	B_aln.clear();
	int i,j,max_score=-100000000;
	int A_sz=strlen(A),B_sz=strlen(B);
	aln_t->deletions=0;
	aln_t->insertions=0;
	aln_t->substitutions=0;
	/*
	for(int k=1;k<=B_sz;++k)
	{
		if(aln_t->score[k][A_sz]>max_score)
		{i=k;}
	}
	*/
	i=A_sz;
	j=B_sz;

	max_score=-1000000;
	for(int k=B_sz;k>=0;--k)
	{
		if((aln_t->score[A_sz][k])>max_score)
		{max_score=aln_t->score[A_sz][k];j=k;}
	}


	while(1)
	{
	//	cout<<i<<" "<<j<<endl;//36//43

		if (i==0)
		{
			
			for(int k=j;k>=1;--k)
			{
				break;//A is used as reference, so we stop here
				A_aln.push_back('-');
				B_aln.push_back(B[k-1]);
			}
			
			break;
		}


		if (j==0)
		{
			//the head of B is missing here
			for(int k=i;k>=1;--k)
			{
				aln_t->insertions++;//A has more bases 
				B_aln.push_back('-');
				A_aln.push_back(A[k-1]);
			}
			break;
		}
		if(aln_t->b[i][j]==3||aln_t->b[i][j]==0)
		{

			A_aln.push_back(A[i-1]);
			B_aln.push_back(B[j-1]);
			if(A[i-1]!=B[j-1])
			{
				aln_t->substitutions++;
			}
			i--;j--;
			continue;
		}
		if(aln_t->b[i][j]==2)
		{
			if(i<A_sz)//trim B at the end to match the length of A
			{
				aln_t->deletions++;//A is missing a base			
				A_aln.push_back('-');
				B_aln.push_back(B[j-1]);
			}
			j--;
			continue;
		}
		if(aln_t->b[i][j]==1)
		{
			aln_t->insertions++;//A has one more base
			B_aln.push_back('-');
			A_aln.push_back(A[i-1]);
			i--;
			continue;
		}
	}
	reverse(A_aln.begin(),A_aln.end());
	reverse(B_aln.begin(),B_aln.end());
}





#endif

