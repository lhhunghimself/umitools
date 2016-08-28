#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

static unsigned char countLUT[65536];
//Encoding is 
//0000 N  0
//0001 A  1
//0010 C  2
//0100 G  4
//1000 T  8

//property being that hamming bitwise distance between any two bases ACGT is 2 an is 1 to N

void initLUT(){
	//check if it is already inited
	for(unsigned int i=0;i<65536;i++){
		countLUT[i]=__builtin_popcount(i);
	}	
}

void encodeSequence(char *sequence, unsigned int length,unsigned int *encodedSequence,unsigned int chunkSize){
	//do it by char
	int i;
	char *byChar=(char*)encodedSequence;
	memset(encodedSequence,0,chunkSize*sizeof(unsigned int));
	for (i=0;i<length;i++){
		unsigned char k=(i%2)? 1 : 1<<4;
		if(sequence[i]=='A'){
			byChar[i/2]+=k;
			continue;
		}
		k=k<<1;	
		if(sequence[i]=='C'){
		 byChar[i/2]+=k;
			continue;
		}
		k=k<<1;	
		if(sequence[i]=='G'){
			byChar[i/2]+=k;
			continue;
		}
		k=k<<1;
		if(sequence[i]=='T'){
			byChar[i/2]+=k;
		}
	}
}
void encodeSequence_ll(char *sequence, long long length,long long *encodedSequence,unsigned int chunkSize){
	//do it by char
	int i;
	char *byChar=(char*)encodedSequence;
	memset(encodedSequence,0,chunkSize*sizeof(long long));
	for (i=0;i<length;i++){
		unsigned char k=(i%2)? 1 : 1<<4;
		if(sequence[i]=='A'){
			byChar[i/2]+=k;
			continue;
		}
		k=k<<1;	
		if(sequence[i]=='C'){
		 byChar[i/2]+=k;
			continue;
		}
		k=k<<1;	
		if(sequence[i]=='G'){
			byChar[i/2]+=k;
			continue;
		}
		k=k<<1;
		if(sequence[i]=='T'){
			byChar[i/2]+=k;
		}
	}
}	
unsigned int hBitWiseDistance(unsigned int *a, unsigned int *b,unsigned int seqLength,unsigned int maxdist){
	unsigned int i,dist=0;
	int n=(seqLength%(sizeof(unsigned int)*2) ==0 )? seqLength/(sizeof(unsigned int)*2): seqLength/(sizeof(unsigned int)*2)+1;
	for(int i=0;i<n && dist< maxdist;i++)
  dist+=(__builtin_popcount(a[i]^b[i]));
 return dist;
}
unsigned int hBitWiseDistance_ll(long long *a, long long *b,unsigned int seqLength,unsigned int maxdist){
	long long i,dist=0;
	int n=(seqLength%(sizeof(long long)*2) ==0 )? seqLength/(sizeof(long long)*2): seqLength/(sizeof(long long)*2)+1;
	for(int i=0;i<n && dist<maxdist;i++)
  dist+=(__builtin_popcountll(a[i]^b[i]));
 return dist;
}
unsigned int LUTDistance(unsigned int *a, unsigned int *b,unsigned int seqLength,int maxdist){
	unsigned int i,dist=0;
	int n=(seqLength%(sizeof(unsigned int)*2) ==0 )? seqLength/(sizeof(unsigned int)*2): seqLength/(sizeof(unsigned int)*2)+1;
	for(int i=0;i<n && dist <maxdist;i++){
	 unsigned int diff=a[i]^b[i];
  dist+=countLUT[diff & 0xFFFF] + countLUT[diff>>16];
	}
 return dist;
}
unsigned int slowHammingDistance(char *a, char *b, unsigned int len,int maxdist){
	unsigned int i,dist=0;
 for (i = 0; i < len && dist<maxdist; ++i){
		if(a[i] !=b[i]){
			dist+=2;
			if(a[i] == 'N' || b[i] == 'N')dist--;
		}	
	}
 return(dist);
}	
