#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include "kseq.h"
#include "umitools.hpp"


void fillRandomSequence(char **a,int length,int nSeqs);
void printSequences(char *sequences,int nSeqs,int nBases);
void printBitsEndian(void *ptr, char nBytes, char endian);

int main(int argc, char *argv[]){
	int nBases=16;
	int nSeqs=1000;
	int maxdist=3;
	clock_t start;
 char *sequences;
 srand(1234);
	uint64_t *encodedSequences;	

	if(argc >= 2){
		nBases=atoi(argv[1]);
	}
	if(argc>=3){
		nSeqs=atoi(argv[2]);
	}
	if(argc==4){
		maxdist=atoi(argv[3]);
	}
	umipanel testPanel;
	
	
	exit();
	
	hamming<uint32_t> hammingEval=hamming<uint32_t>(nBases);
	const unsigned int chunkSize=hammingEval.chunkSize;
 fillRandomSequence(&sequences,nBases,nSeqs);
 nacgt<uint32_t> *encodedSeqs = new nacgt<uint32_t>[nSeqs];

 
	for(int i=0;i<nSeqs;i++){
  encodedSeqs[i]=nacgt<uint32_t>(sequences+i*nBases,nBases);
	}
	int lastDist;		
 fprintf(stderr,"checking correctness\n");
 start=clock();
 int dist;
 for(int i=0;i<nSeqs+1;i++){
		char *a=sequences+i*nBases;
		for(int j=i+1;j<nSeqs;j++){
			char *b=sequences+j*nBases;
			int d1,d2;
			d1=hammingEval.slowDistance(a,b,nBases,maxdist);
			d2=hammingEval.distance(encodedSeqs+i,encodedSeqs+j,maxdist);
			if(d1 !=d2){
				fprintf(stderr,"For index i=%d and j=%d  correct is %d hamming is= %d\n",i,j,d1,d2);

			}	
  }
	}
	start=clock();
	for(int i=0;i<nSeqs+1;i++){
		for(int j=i+1;j<nSeqs;j++){
		 dist=hammingEval.slowDistance(sequences+i*nBases,sequences+j*nBases,nBases,maxdist);
  }
	}
	fprintf (stderr,"slowHamming took %ld ms\n",(clock()-start)*1000/CLOCKS_PER_SEC);
	
	start=clock();
	for(int i=0;i<nSeqs+1;i++){
		for(int j=i+1;j<nSeqs;j++){
			dist=hammingEval.distance(encodedSeqs+i,encodedSeqs+j,maxdist);
  }
	}
	fprintf (stderr,"popCount took %ld ms %d\n",(clock()-start)*1000/CLOCKS_PER_SEC,dist);
	
	
	delete [] encodedSequences;
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
 char *myseq=(char*) malloc(nSeqs*length);
 if(!myseq){
		fprintf(stderr,"unable to allocate memory for sequences\n");
		exit(0);
	}
	for	(int i=0;i<nSeqs*length;i++)
		myseq[i]=bases[rand()%5];
	*a=myseq;
}
