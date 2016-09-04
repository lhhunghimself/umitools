#include <zlib.h>  
#include <stdio.h>
#include <iostream>
#include <string.h>  
#include <omp.h>
#include "kseq.h"
#include "umitools.hpp"
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#ifndef NWELLS
#define NWELLS 96
#endif
//this will divide the folders with the different wells
//this can be done in parallel

extern "C" {
 #include "optparse.h"  
}
bool Ncheck(const char *seq, const int size);

#define R1_LENGTH 16 //default barcode length

KSEQ_INIT(gzFile, gzread);  
  
int main(int argc, char *argv[]){

 char verbose=0,extract=0,append=0,barcode=0;  
 int l1,l2,emptySeqs=0,nameMismatch=0,minQual=10,UMILength=R1_LENGTH,nArgs=0;
 unsigned int barcodeSize=0,minDist=4;
 int adjustedQ=minQual+33;
 char *arg=0,*outputFileName=0;
 string barcodeFileName;
 struct optparse options;
 int opt;
 int mismatchTol=1,NTol=2;
 optparse_init(&options, argv);
 int nThreads=1;
 int filter=0;
 vector<string> inputFiles;
 

 
 //parse flags
 while ((opt = optparse(&options, "aievft:m:N:d:b:l:q:")) != -1) {
  switch (opt){
			case 'a':
			 append=1;
			break;
			case 'v':
			 verbose=1;
			break;
			case 'e':
    extract=1;
   break;
			case 'f':
    filter=1;
   break;     
   case 'd':
    minDist=atoi(options.optarg);
   break;   
   case 'm':
    mismatchTol=atoi(options.optarg);
   break;
   case 't':
    nThreads=atoi(options.optarg);
   break;         
   case 'N':
    NTol=atoi(options.optarg);
   break;
   case 'b':
    barcodeFileName=options.optarg;
   break;
   case 'l':
    UMILength=atoi(options.optarg);//not the well barcode
   break;
   case 'q':
    minQual=atoi(options.optarg);
    adjustedQ=33+minQual;
   break;      
   case '?':
    fprintf(stderr, "%s: %s\n", argv[0], options.errmsg);
    exit(0);
  }
 }
	
 //parse file arguments
 //these will be R1 followed by R2 arguments
 while ((arg = optparse_arg(&options))){
	 inputFiles.push_back(string(arg));
	}
	if((append && extract && ! outputFileName)){
		fprintf(stderr, "must give outputfileName if appending and extracting\n");  
  exit(0);
	}
 //extract fileName
	if(verbose){
		//print out the parameters
		fprintf(stderr,"Verbose mode on\n");
		fprintf(stderr,"UMI Length %d\n",UMILength);
		fprintf(stderr,"Barcode fileName %s\n",barcodeFileName.c_str());
		fprintf(stderr,"Minimum quality %d\n",minQual);
		int i=0;
		while (i<inputFiles.size()){
		 fprintf(stderr,"R1 file %s\n",inputFiles[i++].c_str());
		 if(i == inputFiles.size()){
				fprintf(stderr,"missing corresponding R2 file\n");
				exit(EXIT_FAILURE);
			}
		fprintf(stderr,"R2 file %s\n",inputFiles[i++].c_str());		
		}		
	}
 fs::create_directory(fs::system_complete("./output"));
  //change this if using 384 wells
 umipanel<uint32_t,unsigned char> **barcodePanel=new umipanel<uint32_t,unsigned char>*[nThreads];
 for(int i=0;i<nThreads;i++)
	 barcodePanel[i]=new umipanel<uint32_t,unsigned char> (barcodeFileName,mismatchTol,NTol);
	for(int i=0;i<NWELLS;i++){
		auto p=fs::path("./output/"+barcodePanel[0]->wells[i]);
		fs::create_directory(fs::system_complete(p));
	}
	//create bad directory	
	if(!filter){
	 auto p=fs::path("./output/X");
	 fs::create_directory(fs::system_complete(p));
	}
	string outputDir=fs::system_complete("output/").string();
#pragma omp parallel for num_threads (nThreads) schedule (dynamic)
	for (int i=0;i<inputFiles.size();i+=2){
		const int tid=omp_get_thread_num();
		const int wellSequenceSize=barcodePanel[tid]->barcodeSize;
		gzFile fp1=0, fp2=0;  
  kseq_t *seq1,*seq2;
  int l1,l2;
		if(!(fp1 = gzopen(inputFiles[i].c_str(), "r")))exit(EXIT_FAILURE);
		if(!(fp2 = gzopen(inputFiles[i+1].c_str(), "r")))exit(EXIT_FAILURE);
  seq1 = kseq_init(fp1);
  seq2 = kseq_init(fp2);
  
  fs::path R1(inputFiles[i]);
  fs::path R2(inputFiles[i+1]);
  string R1stem;
  if (R1.extension().string() == "gz" && (R1.stem().extension().string() == "fastq" || (R1.stem().extension().string() == ".fq")))   
   R1stem=R1.stem().stem().string();
  else 
   R1stem=R1.stem().string();
  FILE *ofp,*ofps[NWELLS+1];	
  for(int j=0;j<NWELLS+1;j++){
			ofps[j]=0;
			string file;
			if(!filter && !j){
				file=outputDir+"X"+"/"+R1stem+"R2_"+"X"+".fq";
	   ofps[0]=fopen(file.c_str(),"w");			
	   if(!ofps[j]){
				 fprintf(stderr,"%s unable to open file %s\n",inputFiles[i].c_str(),file.c_str());
			 }	
			}	
			else if(j){
    file =outputDir+barcodePanel[tid]->wells[j-1]+"/"+R1stem+"R2_"+barcodePanel[tid]->wells[j-1]+".fq";
    ofps[j]=fopen(file.c_str(),"w");
			 if(!ofps[j]){
				 fprintf(stderr,"%s unable to open file %s\n",inputFiles[i].c_str(),file.c_str());
			 }				
			}
		}

	 
  while ((l1 = kseq_read(seq1)) >= 0 && (l2 = kseq_read(seq2)) >=0) {
	  //check for errors 
		 if(l1==0 || l2 == 0 || strcmp(seq1->name.s,seq2->name.s)){
		 	if(l1==0 || l2 == 0) emptySeqs++;
		 	if(strcmp(seq1->name.s,seq2->name.s)){
		 		nameMismatch++;
		 		if(verbose)fprintf(stderr,"Warning - mismatch of names in R1 and R2\n");
		 	}	
		 	continue;
		 }
		 //adjust quality and find the well

			char *cptr=seq1->seq.s;
		 char *qptr=seq1->qual.s;
		 int k=0; 
			while(*cptr && k < UMILength){
		  if(*qptr<adjustedQ)*cptr='N';
		  cptr++;qptr++;k++;
		 }
		 //check if there is a N
		 if(filter && Ncheck((seq1->seq.s)+wellSequenceSize,UMILength-wellSequenceSize)) continue;
		 const unsigned int barcodeIndex=barcodePanel[tid]->bestMatch(seq1->seq.s);
		 if(filter && !barcodeIndex)continue;
		 ofp=ofps[barcodeIndex];
			fputc('@',ofp);
		 fputs(seq2->name.s,ofp);
		 fputc(':',ofp);				
		 fwrite(seq1->seq.s,UMILength,1,ofp);
			fputc('\n',ofp);	
			fputs(seq2->seq.s,ofp);
		 fputs("\n+\n",ofp);		 
   fputs(seq2->qual.s,ofp);
   fputc('\n',ofp);
						 
		}
		if(!filter)
		 fclose(ofps[0]);		   
		for(int j=1;j<NWELLS+1;j++)
		 fclose(ofps[j]);
		  
  kseq_destroy(seq1);
  kseq_destroy(seq2);
  gzclose(fp1);
  gzclose(fp2);		 	 	
	}
	for(int i=0;i<nThreads;i++)
  if(barcodePanel[i])delete barcodePanel[i];
 delete[] barcodePanel;
 return 0;  
}

bool Ncheck(const char *seq, const int size){
	for(int i=0;i<size;i++)
	 if(seq[i] == 'N')return 1;
	return 0;
}
