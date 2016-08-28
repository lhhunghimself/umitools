#include <zlib.h>  
#include <stdio.h>
#include <string.h>  
#include "kseq.h"
#include "optparse.h"  
#include "fasthamming.h"
#define R1_LENGTH 16 //default barcode length

KSEQ_INIT(gzFile, gzread);  
  
int main(int argc, char *argv[]){
 gzFile fp1=0, fp2=0;  
 kseq_t *seq1,*seq2;
 char *extractFilename=0;
 FILE *outputfp=0,*extractfp=0;
 char verbose=0,extract=0,append=0;  
 int l1,l2,emptySeqs=0,nameMismatch=0,minQual=10,barcodeLength=R1_LENGTH,nArgs=0;
 int adjustedQ=minQual+33;
 char *arg=0,*outputFilename=0,*inputFilenames[2]={0,0};
 struct optparse options;
 int opt;
 optparse_init(&options, argv);

 //parse flags
 while ((opt = optparse(&options, "aiev:o:l:q:")) != -1) {
  switch (opt){
			case 'a'
			 append=1;
			break;
			case 'v':
			 verbose=1;
			break;
			case 'e':
    extract=1;
   break; 
   case 'o':
    outputFilename=options.optarg;
   break;
   case 'l':
    barcodeLength=atoi(options.optarg);
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
		inputFilenames[nArgs++]=arg;
  if(nArgs >2){
		 fprintf(stderr, "Usage: %s -lf <in.r1> <in.r2>\n", argv[0]);  
   exit(0); 
		}	
	}
	if(nArgs !=2){
		fprintf(stderr, "Usage: %s -lf <in.r1> <in.r2>\n", argv[0]);  
  exit(0);
	}
	if((append && extract && ! outputfileName)){
		fprintf(stderr, "must give outputfileName if appending and extracting\n");  
  exit(0);
	}
 //extract fileName
 if(extract){
  if(!(extractFilename=malloc(strlen(outputFilename)+5))){
	 	exit(0);
	 }
	}
	if(verbose){
		//print out the parameters
		fprintf(stderr,"Verbose mode on\n");
		fprintf(stderr,"Barcode Length %d\n",barcodeLength);
		fprintf(stderr,"Minimum quality %d\n",minQual);
		if(outputFilename)fprintf(stderr,"Outputs to %s\n",outputFilename);
				if(outputFilename && extract)fprintf(stderr,"Outputs umis to %s\n",extractFilename);
		else fprintf(stderr,"outputs to stdout\n");
		fprintf(stderr,"R1 file %s\n",inputFilenames[0]);		
		fprintf(stderr,"R2 file %s\n",inputFilenames[1]);
	}
	//check if the 		       
 if(fp1 = gzopen(inputFilenames[0], "r")){
		if(verbose)fprintf(stderr,"opened R1 file %s\n",inputFilenames[0]);
	}
	else{
	 fprintf(stderr,"unable to open R1 file %s\n",inputFilenames[0]),	 
		exit(0);
	}	
 if(fp2 = gzopen(inputFilenames[1], "r")){
		if(verbose)fprintf(stderr,"opened R2 file %s\n",inputFilenames[1]);
	}
	else{
	 fprintf(stderr,"unable to open R2 file %s\n",inputFilenames[1]),	 
		exit(0);
	}
 seq1 = kseq_init(fp1);
 seq2 = kseq_init(fp2);
 if(outputFilename){
		if(outputfp=fopen(outputFilename,"w")){
			if(verbose) fprintf(stderr,"opening output file %s\n",outputFilename);
		}	
		else{	
		 fprintf(stderr,"unable to open output file %s\n",outputFilename),	 
		 exit(0);
		}
		if(extract){
			strcpy(extractFilename,outputFilename);
			strcat(extractFilename,".umis");
	  if(extractfp=fopen(extractFilename,"w")){
			 if(verbose) fprintf(stderr,"opening umi file %s\n",extractFilename);
		 }	
		 else{	
		  fprintf(stderr,"unable to open output file %s\n",outputFilename),	 
		  exit(0);
		 }
		}
	}
 while ((l1 = kseq_read(seq1)) >= 0 && (l2 = kseq_read(seq2)) >=0) {
		int i=0;
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
		fputc('@',outputfp);
		fputs(seq2->name.s,outputfp);
		fputc(':',outputfp);
		//append sequence
		cptr=seq1->seq.s;
		qptr=seq1->qual.s;
		if(extract){
			while(*cptr && i < barcodeLength){
			 if(*qptr<adjustedQ){
					if(outputfp)fputc('N',outputfp);
					fputc('N',extractfp);
				}
			 else{
					if(outputfp)fputc(*cptr,outputfp);
					fputc(*cptr,extractfp);
				}
			 cptr++;qptr++;i++;
		 }
		 fputc('\n',extractfp);	
		}
		else{	
		 while(*cptr && i < barcodeLength){
		 	if(*qptr<adjustedQ)fputc('N',outputfp);
		 	else fputc(*cptr,outputfp);
		 	cptr++;qptr++;i++;
		 }
		}
		fputc('\n',outputfp);
		fputs(seq2->seq.s,outputfp);
		fputs("\n+\n",outputfp);		 
  fputs(seq2->qual.s,outputfp);
  fputc('\n',outputfp);
	}
	if(outputfp != stdout)fclose(outputfp);		 
 kseq_destroy(seq1);
 kseq_destroy(seq2);
 gzclose(fp1);
 gzclose(fp2);
 if(extract){
		fclose(extractfp);
		if(extractFileName)free(extractFileName)
	}	
 if(emptySeqs || nameMismatch){
		fprintf(stderr,"WARNING %d empty sequences and %d name mismatches encountered\n",emptySeqs,nameMismatch);
	}
	if(outputFileName) free(outputFileName);	
 return 0;  
}

