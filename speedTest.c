#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include "kseq.h"
#include "fasthamming.hpp"

#define MAXDIFF 6
int fillSequence(char **a,int length);shoul
void fillRandomSequence(char **a,int length,int nSeqs);
void printSequences(char *sequences,int nSeqs,int nBases);


int main(int argc, char *argv[]){
	int nBases=8;
	int nSeqs=100000;
	int maxdist=MAXDIFF;
	clock_t start, diff;
 char *sequences;
 srand(1234);
	unsigned int  *encodedSequences;	
	long long  *encodedSequences_ll;
	if(argc >= 2){
		nBases=atoi(argv[1]);
	}
	if(argc>=3){
		nSeqs=atoi(argv[2]);
	}
	if(argc==4){
		maxdist=atoi(argv[3]);
	}
}	
int fillSequence(char **a,int length){
	int nSeqs=pow(5,length);
	char bases[]={'N','A','C','G','T'};
 char *myseq=malloc(nSeqs*length);
 if(!myseq){
		fprintf(stderr,"unable to allocate memory for sequences\n");
		exit(0);
	}
	//set up reverse odometer using %
	//multiply instead of divide
	int k=0;
	for	(int i=0;i<nSeqs;i++){
		int factor=1;
		for(int j=0;j<length;j++){
	  myseq[k++]=bases[(i/factor) % 5];
	  factor*=5;
		}
	}
	*a=myseq;
	return nSeqs;
}
void printSequences(char *sequences,int nSeqs,int nBases){
	for (int i=0;i<nSeqs*nBases;i++){
	 fputc(sequences[i],stdout);
	 if(i%nBases==nBases-1){
			fputc('\n',stdout);
		}	
	}	
}
void fillRandomSequence(char **a,int length,int nSeqs){
	char bases[]={'N','A','C','G','T'};
 char *myseq=malloc(nSeqs*length);
 if(!myseq){
		fprintf(stderr,"unable to allocate memory for sequences\n");
		exit(0);
	}
	int k=0;
	for	(int i=0;i<nSeqs*length;i++)
		myseq[i]=bases[rand()%5];
	*a=myseq;
}
