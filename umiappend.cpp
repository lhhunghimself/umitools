#include <zlib.h>  
#include <stdio.h>
#include <string.h>  
#include "kseq.h"
#include "umitools.hpp"
extern "C" {
 #include "optparse.h"  
}

#define R1_LENGTH 16 //default barcode length

KSEQ_INIT(gzFile, gzread);  
  
int main(int argc, char *argv[]){
 gzFile fp1=0, fp2=0;  
 kseq_t *seq1,*seq2;
 char *extractFileName=0;
 FILE *outputfp=0,*extractfp=0;
 char verbose=0,extract=0,append=0,barcode=0;  
 int l1,l2,emptySeqs=0,nameMismatch=0,minQual=10,UMILength=R1_LENGTH,nArgs=0;
 unsigned int barcodeSize=0,minDist=4;
 int adjustedQ=minQual+33;
 char *arg=0,*outputFileName=0,*inputFileNames[2]={0,0};
 string barcodeFileName;
 struct optparse options;
 int opt;
 int hashTolerance=1;
 optparse_init(&options, argv);
 
 //change this if using 384 wells
 umipanel<uint32_t,unsigned char> *barcodePanel=0;
 
 //parse flags
 while ((opt = optparse(&options, "aievd:b:o:l:q:")) != -1) {
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
   case 'o':
    outputFileName=options.optarg;
   break;
   case 'd':
    minDist=atoi(options.optarg);
   break;
   case 'b':
    barcode=1;
    barcodeFileName=options.optarg;
    barcodePanel=new umipanel<uint32_t,unsigned char> (barcodeFileName,hashTolerance);
   break;
   case 'l':
    UMILength=atoi(options.optarg);
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
 while ((arg = optparse_arg(&options))){
		inputFileNames[nArgs++]=arg;
  if(nArgs >2){
		 fprintf(stderr, "Usage: %s -lf <in.r1> <in.r2>\n", argv[0]);  
   exit(0); 
		}	
	}
	if(nArgs !=2){
		fprintf(stderr, "Usage: %s -lf <in.r1> <in.r2>\n", argv[0]);  
  exit(0);
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
		if(outputFileName)fprintf(stderr,"Outputs to %s\n",outputFileName);
				if(outputFileName && extract)fprintf(stderr,"Outputs umis to %s\n",extractFileName);
		else fprintf(stderr,"outputs to stdout\n");
		fprintf(stderr,"R1 file %s\n",inputFileNames[0]);		
		fprintf(stderr,"R2 file %s\n",inputFileNames[1]);
	}
	//check if the 		       
 if((fp1 = gzopen(inputFileNames[0], "r"))){
		if(verbose)fprintf(stderr,"opened R1 file %s\n",inputFileNames[0]);
	}
	else{
	 fprintf(stderr,"unable to open R1 file %s\n",inputFileNames[0]),	 
		exit(0);
	}	
 if((fp2 = gzopen(inputFileNames[1], "r"))){
		if(verbose)fprintf(stderr,"opened R2 file %s\n",inputFileNames[1]);
	}
	else{
	 fprintf(stderr,"unable to open R2 file %s\n",inputFileNames[1]),	 
		exit(0);
	}
 seq1 = kseq_init(fp1);
 seq2 = kseq_init(fp2);
 if(outputFileName){
		if((outputfp=fopen(outputFileName,"w"))){
			if(verbose) fprintf(stderr,"opening output file %s\n",outputFileName);
		}	
		else{	
		 fprintf(stderr,"unable to open output file %s\n",outputFileName),	 
		 exit(0);
		}
		if(extract){
			extractFileName=new char[strlen(outputFileName)+6];
			strcpy(extractFileName,outputFileName);
			strcat(extractFileName,".umis");
	  if((extractfp=fopen(extractFileName,"w"))){
			 if(verbose) fprintf(stderr,"opening umi file %s\n",extractFileName);
		 }	
		 else{	
		  fprintf(stderr,"unable to open output file %s\n",outputFileName),	 
		  exit(0);
		 }
		}
	}
	else{
		if(extract)extractfp=stdout;
		else outputfp=stdout;
	}	
 while ((l1 = kseq_read(seq1)) >= 0 && (l2 = kseq_read(seq2)) >=0) {
		int i=0;
		unsigned int nBest=0;
		string bestWell,bestSequence;
	 char *cptr,*qptr; 
		if(l1==0 || l2 == 0 || strcmp(seq1->name.s,seq2->name.s)){
			if(l1==0 || l2 == 0) emptySeqs++;
			if(strcmp(seq1->name.s,seq2->name.s)){
				nameMismatch++;
				if(verbose)fprintf(stderr,"Warning - mismatch of names in R1 and R2\n");
			}	
			continue;
		}
		//print orig name
		if(outputfp){
		 fputc('@',outputfp);
		 fputs(seq2->name.s,outputfp);
		 fputc(':',outputfp);
		}
		
  //append sequence
		cptr=seq1->seq.s;
		qptr=seq1->qual.s;
		//adjust for quality
		while(*cptr && i < UMILength){
		 if(*qptr<adjustedQ)*cptr='N';
		 cptr++;qptr++;i++;
		}
		//barcode correct
		int wellIndex=barcodePanel->bestMatch(seq1->seq.s);		
  if(wellIndex--){ //wellIndex returns 0 when no match index+1 otherwise
			cptr=seq1->seq.s+barcodeSize;
		 if(extractfp){
				fwrite(barcodePanel->wells[wellIndex].c_str(),barcodePanel->wells[wellIndex].length(),1,extractfp);
				fputc(':',extractfp);
				fwrite(cptr,UMILength-barcodeSize,1,extractfp);;
				fputc('\n',extractfp);
			}
		 if(outputfp){
    fwrite(barcodePanel->wells[wellIndex].c_str(),barcodePanel->wells[wellIndex].length(),1,outputfp);
				fputc(':',outputfp);
				fwrite(cptr,UMILength-barcodeSize,1,outputfp);;
				fputc('\n',outputfp);
			}
  }
		else{
			cptr=seq1->seq.s;
		 if(extractfp){
				fputs("X",extractfp);
				fputc(':',extractfp);
				fwrite(cptr,UMILength-barcodeSize,1,extractfp);
				fputc('\n',extractfp);
			}				
		 if(outputfp){
				fputs("XX",outputfp);
				fputc(':',outputfp);
				fwrite(cptr,UMILength-barcodeSize,1,outputfp);;
				fputc('\n',outputfp);
			}					
		}
		if(outputfp){
		 fputs(seq2->seq.s,outputfp);
		 fputs("\n+\n",outputfp);		 
   fputs(seq2->qual.s,outputfp);
   fputc('\n',outputfp);
		}
	}
	if(outputfp && outputfp != stdout)fclose(outputfp);		 
 kseq_destroy(seq1);
 kseq_destroy(seq2);
 gzclose(fp1);
 gzclose(fp2);
 if(extract){
		if(extractfp !=stdout)fclose(extractfp);
		if(extractFileName)delete [] extractFileName;
	}	
 if(emptySeqs || nameMismatch){
		fprintf(stderr,"WARNING %d empty sequences and %d name mismatches encountered\n",emptySeqs,nameMismatch);
	}
	//if(outputFileName) free(outputFileName);	
	if(barcodePanel)delete barcodePanel;
 return 0;  
}

