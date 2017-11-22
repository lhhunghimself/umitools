#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include "kseq.h"
#include "umitools.hpp"

void fillRandomSequence(char **a,int length,int nSeqs);
void printSequences(char *sequences,int nSeqs,int nBases);
void printBitsEndian(void *ptr, char nBytes, char endian);
template <class T1,class T2>  void randomBarcodes(char **a,umipanel<T1,T2> *p,int nSeqs);
void obscureBarcodes(char **a,char *barcodes,int barcodeSize,int nSeqs,int nObscure);
void obscureBarcodesR(char **a,char *barcodes,int barcodeSize,int nSeqs,float obscureChance);
void testPanel(int nBases,int nSeqs,int maxdist);
//relative speeds
//for searching against panel hash: 32 ms fast: 352 ms slow: 1110 ms (1 million sequences - 1 N)
//~ 80:1 for fast vs slow for long sequences 

int main(int argc, char *argv[]){
	//usage speedTest <nBases> <nSeqs> <maxDist> <seed>
	int nBases=16;
	int nSeqs=1000;
	int maxdist=4*nBases;
	clock_t start;
 char *sequences;
 int seed=1234;
 
	if(argc >= 2){
		nBases=atoi(argv[1]);
	}
	if(argc>=3){
		nSeqs=atoi(argv[2]);
	}
	if(argc>=4){
		maxdist=atoi(argv[3]);
	}
	if(argc==5){
		seed=atoi(argv[4]);
	}
	srand(seed);
// testPanel(nBases,nSeqs,maxdist);

	
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
			char sa[64],sb[64];
			sa[nBases]=0;
			sb[nBases]=0;
			char *b=sequences+j*nBases;
			int d1,d2;
			d1=hammingEval.slowDistance(a,b,nBases,maxdist);
			d2=hammingEval.distance(encodedSeqs+i,encodedSeqs+j,maxdist);
			memcpy(sa,a,nBases);
			memcpy(sb,b,nBases);
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
	fprintf (stderr,"slowHamming took %ld ms lastDist %d\n",(clock()-start)*1000/CLOCKS_PER_SEC,dist);
	
	start=clock();
	for(int i=0;i<nSeqs+1;i++){
		for(int j=i+1;j<nSeqs;j++){
			dist=hammingEval.distance(encodedSeqs+i,encodedSeqs+j,maxdist);
  }
	}
	fprintf (stderr,"popCount took %ld ms lastDist %d\n",(clock()-start)*1000/CLOCKS_PER_SEC,dist);
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
template <class T1,class T2>  void randomBarcodes(char **a,umipanel<T1,T2> *p,int nSeqs){
 char *myseq=(char*) malloc(nSeqs*p->barcodeSize);
 if(!myseq){
		fprintf(stderr,"unable to allocate memory for sequences\n");
		exit(0);
	}
	for	(int i=0;i<nSeqs;i++){
		const char *newseq=p->sequences[rand()%p->nBarcodes].c_str();
		memcpy(myseq+i*p->barcodeSize,newseq,p->barcodeSize);
	}		
	*a=myseq;
}
void obscureBarcodesR(char **a,char *barcodes,int barcodeSize,int nSeqs,float obscureChance){
 char *myseq=(char*) malloc(nSeqs*barcodeSize);
 if(!myseq){
		fprintf(stderr,"unable to allocate memory for sequences\n");
		exit(0);
	}
	int ilimit=obscureChance*100000;
	memcpy(myseq,barcodes,nSeqs*barcodeSize);
	for	(int i=0;i<nSeqs*barcodeSize;i++){
		if(rand()%100000 < ilimit ) myseq[i]='N';
	}		
	*a=myseq;
}
void obscureBarcodes(char **a,char *barcodes,int barcodeSize,int nSeqs,int nObscure){
	unsigned char *masterDeck=new unsigned char[barcodeSize];
	unsigned char *deck=new unsigned char[barcodeSize];
	for(unsigned char i=0;i<barcodeSize;i++){
		masterDeck[i]=i;
	}	
 char *myseq=(char*) malloc(nSeqs*barcodeSize);
 if(!myseq){
		fprintf(stderr,"unable to allocate memory for sequences\n");
		exit(0);
	}
	memcpy(myseq,barcodes,nSeqs*barcodeSize);
	for	(int i=0;i<nSeqs;i++){
		memcpy(deck,masterDeck,barcodeSize);
		for(int j=0;j<nObscure;j++){
			int swap=rand()%(barcodeSize-j);
			myseq[i*barcodeSize+deck[swap]]='N';
			if(swap !=barcodeSize-1-j){
			 char temp=deck[barcodeSize-1-j];
			 deck[barcodeSize-1-j]=deck[swap];
			 deck[swap]=temp;
			}
		}
	}
	delete [] masterDeck;
	delete [] deck;		
	*a=myseq;
}


void testPanel(int nBases,int nSeqs,int maxdist){
	clock_t start;
	hamming<uint32_t> phammingEval=hamming<uint32_t>(nBases);
	umipanel <uint32_t,unsigned char> testPanel(string("/mnt/backup/barcodes_trugrade_96_set4.dat"),0,2);
 char *goodBarcodes,*badBarcodes;
	randomBarcodes<uint32_t,unsigned char>(&goodBarcodes,&testPanel,nSeqs);
 obscureBarcodes(&badBarcodes,goodBarcodes,nBases,nSeqs,maxdist);
 acgt<uint32_t> *encodedbcs = new acgt<uint32_t>[nSeqs];
 for(int i=0;i<nSeqs;i++){
  encodedbcs[i]=acgt<uint32_t>(badBarcodes+i*nBases,nBases);
	}

	char *dists=new char[nSeqs];
	int dist=0;
 for(int i=0;i<nSeqs;i++){
		int dist1=phammingEval.bestMatch(testPanel.sequences,badBarcodes+i*nBases,testPanel.nBarcodes);
  int dist2=phammingEval.bestMatch(&(testPanel.Aseqs[0]),encodedbcs+i,testPanel.nBarcodes);
  int dist3=testPanel.bestMatch(badBarcodes+i*nBases)-1;
  if(dist2 != dist3){
			fprintf(stderr,"mismatch %d %d %d\n",i,dist1,dist3);
			exit(1);
		}	
	}
	start=clock();
 for(int i=0;i<nSeqs;i++){
  dists[i]=phammingEval.bestMatch (&(testPanel.Aseqs[0]),encodedbcs+i,testPanel.nBarcodes);
	}	
	fprintf (stderr,"panel1 compare took %ld ms\n",(clock()-start)*1000/CLOCKS_PER_SEC);

	start=clock();
 for(int i=0;i<nSeqs;i++){
  dists[i]=phammingEval.bestMatch(testPanel.sequences,badBarcodes+i*nBases,testPanel.nBarcodes);
	}	
	fprintf (stderr,"direct compare took %ld ms\n",(clock()-start)*1000/CLOCKS_PER_SEC);

	start=clock();
 for(int i=0;i<nSeqs;i++){
  dists[i]=testPanel.bestMatch(badBarcodes+i*nBases);;
	}	
	fprintf (stderr,"panel1 compare took %ld ms\n",(clock()-start)*1000/CLOCKS_PER_SEC);
	
	//check
	int nCorrected=0;
	for(int i=0;i<nSeqs;i++){
		if(dists[i] > 0){
			nCorrected++;
		}	
	}
	fprintf (stderr,"%ld out %ld corrected\n",nCorrected,nSeqs);	 
}
