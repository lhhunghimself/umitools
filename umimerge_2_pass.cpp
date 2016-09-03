#include <zlib.h>  
#include <stdio.h>
#include <string.h>
#include <string> 
#include <iostream> 
#include <unordered_map> 
#include <unordered_set> 
#include <set> 
#include <vector> 
#include "kseq.h"
#include <glob.h>
#include "umitools.hpp"

#define MAX_EDIT_DISTANCE 1
#define MAX_BEST 20
#ifndef NWELLS
#define NWELLS 96
#endif

using namespace std;   

bool polyACheck(string &sequence);
void readRefseq(char *filename, unordered_map<string,string> &refseq_to_gene,vector<string> &geneList);
string readWells(char *filename, unordered_map<string,unsigned int> &well_to_index,vector<string> &wellList);
void readERCC(char *filename, vector<string> &erccList);
void splitStr(char *cstr,const char *delim, vector<string> &items);
void splitStr(string str,const char *delim, vector<string> &items);
string splitStrIndex(string str,const char *delim, int index);
bool multiGeneHit(vector<string> &best_list, string gene, unordered_map<string,string> &refseq_to_gene);
unsigned int sumCounts(unsigned int **counts,unsigned int size1, unsigned int size2,int offset);
unsigned int sumCountsi(unsigned int **counts,unsigned int i, unsigned int size2,int offset);
int main(int argc, char *argv[]){
	if(argc != 8){
		fprintf(stderr,"%d parms requires params: sample_id, sym2ref, ercc_fasta, barcodes, aligned_dir, dge_dir, loose_barcodes\n",argc);
		for(int i=1;i<argc;i++){
			fprintf(stderr,"%d\t%s\n",i,argv[i]);
		}	
		exit(EXIT_FAILURE);
	}
	int k=1;	
 char *sample_id=argv[k++];
 char *sym2ref=argv[k++];
 char *ercc_fasta=argv[k++];
 char *barcodes=argv[k++];
 char *aligned_dir=argv[k++];
 char *dge_dir=argv[k++];
 char *loose_barcodes=argv[k++];
 
 unordered_map<string,string> refseq_to_gene;
 vector<string>erccList, geneList,wellList;
 unordered_map<string,unsigned int>well_to_index;
 unordered_map<string,unsigned int>ercc_to_index;
 unordered_map<string,unsigned int>gene_to_index;
 vector<string> unknown_list;
 
 string plateID=readWells(barcodes,well_to_index,wellList);
 readERCC(ercc_fasta,erccList);
 readRefseq(sym2ref,refseq_to_gene ,geneList);
 
 for(int i=0;i<geneList.size();i++)
  gene_to_index[geneList[i]]=i;
 for(int i=0;i<erccList.size();i++)
  ercc_to_index[erccList[i]]=i;
 
 //outputfile names
 string logFile[2]={string(dge_dir)+"/"+string(sample_id)+".unq.log.dat",string(dge_dir)+"/"+string(sample_id)+".all.log.dat"};
 string unknownFile[2]={string(dge_dir)+"/"+string(sample_id)+".unq.unknown_list",string(dge_dir)+"/"+string(sample_id)+".all.unknown_list"};
 string totalAlignFile[2]={string(dge_dir)+"/"+string(sample_id)+".unq.refseq.total.dat",string(dge_dir)+"/"+string(sample_id)+".all.refseq.total.dat"};
 string umiFile[2]={string(dge_dir)+"/"+string(sample_id)+".unq.refseq.umi.dat",string(dge_dir)+"/"+string(sample_id)+".all.refseq.umi.dat"};
 string erccFile[2]={string(dge_dir)+"/"+string(sample_id)+".unq.spike.total.dat",string(dge_dir)+"/"+string(sample_id)+".all.spike.total.dat"};
 string erccUMIFile[2]={string(dge_dir)+"/"+string(sample_id)+".unq.spike.umi.dat",string(dge_dir)+"/"+string(sample_id)+".all.spike.umi.dat"};
 string wellTotalFile[2]={string(dge_dir)+"/"+string(sample_id)+".unq.well_summary.dat",string(dge_dir)+"/"+string(sample_id)+".all.well_summary.dat"};
 
 
 //get all the sam files
 glob_t glob_result;
 char glob_str[1024];
 sprintf(glob_str,"%s/*.sam",aligned_dir);
 glob(glob_str,GLOB_TILDE,NULL,&glob_result);
 
 int mismatchTol=0;
 int NTol=0;
 
 //umipanel<uint32_t,unsigned char> *barcodePanel=0; //change to <uint32_t,uint16_t> for 384 wells
 umipanel<uint32_t,unsigned char> barcodePanel(string(barcodes),mismatchTol,NTol);  		
	//initialize counters
	unsigned int 
	 //counts
	 count_total_reads[2] = {0,0},
	 count_assigned_reads[2] = {0,0},
	 count_assigned_aligned_reads[2] = {0,0},	
	 count_assigned_mito_reads[2] = {0,0},
	 count_assigned_mito_umi[2] = {0,0},
 	count_assigned_unknown_reads[2] = {0,0},
	 count_assigned_unknown_umi[2] = {0,0},
  
  //3D arrays for spike ins
  //arrays for pointers
  
  
		*count_spike_total[NWELLS],
		*count_spike_umi[NWELLS],   
		*count_total[NWELLS],
		*count_umi[NWELLS];   

  //allocate memory for total
  //first half of memory is to store unq - second for multigene counts
		count_spike_total[0]=new unsigned int[erccList.size()*NWELLS*2];
		count_spike_umi[0]=new unsigned int[erccList.size()*NWELLS*2];      		
		count_total[0]=new unsigned int[geneList.size()*NWELLS*2];
		count_umi[0]=new unsigned int[geneList.size()*NWELLS*2];

		
	//set to zero
		
		memset(count_spike_total[0],0,erccList.size()*NWELLS*2*sizeof(unsigned int));
		memset(count_spike_umi[0],0,erccList.size()*NWELLS*2*sizeof(unsigned int));    		
		memset(count_total[0],0,geneList.size()*NWELLS*2*sizeof(unsigned int));
		memset(count_umi[0],0,geneList.size()*NWELLS*2*sizeof(unsigned int));

  
		//initialize pointers
		for(int i=1;i<NWELLS;i++){
			count_spike_total[i]=count_spike_total[i-1]+erccList.size()*2;
			count_spike_umi[i]=count_spike_umi[i-1]+erccList.size()*2;      
			count_total[i]=count_total[i-1]+geneList.size()*2;
			count_umi[i]=count_umi[i-1]+geneList.size()*2;
		}
  
		
	//text hash - may change trie later

	unordered_set<string>unknown_set;
	unordered_set<string>unknown_set_mm;


	//for hashvalue - gives rise to 96*1024*1024*28000 ~ 42 bits - faster than 20 byte text to hash

	const uint64_t barcodeBase=geneList.size()+erccList.size()+1; //add1 for chrM
	const uint64_t wellBase=(1<<20)*barcodeBase;
 for(unsigned int i=0; i<glob_result.gl_pathc; ++i) {
  FILE *fp =fopen(glob_result.gl_pathv[i],"r");
  if(!fp){
			fprintf(stderr,"unable to open file %s\n", glob_result.gl_pathv[i]);
			continue;
		}
		char line[1024];
  for(int mm=0;mm<2;mm++){
			const unsigned int erccOffset=erccList.size()*mm;
			const unsigned int geneOffset=geneList.size()*mm;
   unordered_set<uint64_t>umi_seen;
  	unordered_set<string>unknown_umi_seen;
			while(fgets(line,sizeof(line),fp)){
		  if(line[0] != '@'){
					vector <string> items;
					vector <string> tempItems;
					count_total_reads[mm]++;
					splitStr(line," \t",items);
		   //get rid of all the bad counts
	
		   //parse for barcode
	
	 	  splitStr(items[0],":",tempItems);
	 	  string fullBarcode=tempItems[tempItems.size()-1];
	 	  
	 	  //string well=tempItems[tempItems.size()-2];
	 	  int wellIndex=barcodePanel.bestMatch(fullBarcode.c_str());
	 	  if(!wellIndex)continue;
	 	  wellIndex--;
		   string barcode=fullBarcode.substr(6,10); //change if barcode size changes
		   
		   //skip if ambiguous barcode
		   unsigned int barcodeIndex=0;
		   if(ambigCheck(barcode,barcodeIndex))continue;
	
	    
		   count_assigned_reads[mm]++;
		   string aligned_id=items[2];
		   if(aligned_id == "*")continue;
		   string read=items[9];
		   int edit_dist=stoi(splitStrIndex(items[12],":",-1));
		   int best_hits=stoi(splitStrIndex(items[13],":",-1));
		    
	    if(edit_dist > MAX_EDIT_DISTANCE || best_hits > MAX_BEST || polyACheck(read)) continue;
	    
	    vector<string> best_hits_list;
	    if(items.size() > 19){				
						string best_hits_loc = splitStrIndex(items[19],":",-1);
		    vector<string> split_loc;
		    splitStr(best_hits_loc,";",split_loc);
		    for (int k=0;k<split_loc.size()-1;k++){//extra ; at end so need the -1
							best_hits_list.push_back(splitStrIndex(split_loc[k],",",0));			
						} 
					}

					count_assigned_aligned_reads[mm]++;
					if (aligned_id.substr(0,4) == "ERCC"){
						if(!mm && best_hits_list.size()) continue;
						const int erccIndex=ercc_to_index[aligned_id];
						uint64_t umi=wellIndex*wellBase+barcodeIndex*barcodeBase+erccIndex+geneList.size();
						count_spike_total[wellIndex][erccIndex+erccOffset]++;
						if(!umi_seen.count(umi)){
							umi_seen.insert(umi);
							count_spike_umi[wellIndex][erccIndex+erccOffset]++;
						}
					}
					else if 	(aligned_id.substr(0,4) == "chrM"){
						if(!mm && best_hits_list.size()) continue;
						uint64_t umi=wellIndex*wellBase+barcodeIndex*barcodeBase+barcodeBase-2;
						count_assigned_mito_reads[mm]++;
						if(!umi_seen.count(umi)){
							count_assigned_mito_umi[mm]++;
							umi_seen.insert(umi);
						}
					}
					else if (refseq_to_gene.count(aligned_id)){
						string gene = refseq_to_gene[aligned_id];
						const unsigned int geneIndex =gene_to_index[gene];
						uint64_t umi=wellIndex*wellBase+barcodeIndex*barcodeBase+geneIndex;
						bool multiGene=multiGeneHit(best_hits_list,gene,refseq_to_gene);
						if(!mm && multiGene)continue;
						count_total[wellIndex][geneIndex+geneOffset]++;
						if(!umi_seen.count(umi)){
							count_umi[wellIndex][geneIndex+geneOffset]++;
							umi_seen.insert(umi);
						}					
					}
					else{
						if(!mm && best_hits_list.size()) continue;
						string umi=barcodePanel.wells[wellIndex] + barcode + aligned_id;
						if(!unknown_set.count(aligned_id)){
						 unknown_set.insert(aligned_id);
						}
						count_assigned_unknown_reads[mm]++;
						if(!unknown_umi_seen.count(umi)){
						 count_assigned_unknown_umi[mm]++;
							unknown_umi_seen.insert(umi);
						}				
					}		 
			 }
			}
			rewind(fp);
		}
	 fclose(fp);
	 
	}

 for(int mm=0;mm<2;mm++){
					const unsigned int erccOffset=erccList.size()*mm;
			const unsigned int geneOffset=geneList.size()*mm;
	 FILE *fp=fopen(logFile[mm].c_str(),"w");
	 if(!fp)exit(EXIT_FAILURE);
	 fprintf(fp,"Sample_ID\tTotal\tAssigned\tAligned\tSpike_Total\tSpike_UMI\tMito_Total\tMito_UMI\tRefseq_Total\tRefseq_UMI\tUnknown_Total\tUnknown_UMI\n");
	 fprintf(fp,"%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
	           sample_id,
	           count_total_reads[mm],
												count_assigned_reads[mm],
												count_assigned_aligned_reads[mm],
												sumCounts(count_spike_total,NWELLS,erccList.size(),erccOffset),
												sumCounts(count_spike_umi,NWELLS,erccList.size(),erccOffset),
												count_assigned_mito_reads[mm],
												count_assigned_mito_umi[mm],
												sumCounts(count_total,NWELLS,geneList.size(),geneOffset),
												sumCounts(count_umi,NWELLS,geneList.size(),geneOffset),
												count_assigned_unknown_reads[mm],
												count_assigned_unknown_umi[mm]
										);
	 fclose(fp);
	 if(!(fp=fopen(unknownFile[mm].c_str(),"w")))exit(EXIT_FAILURE);
	 for (auto itr = unknown_set.begin(); itr != unknown_set.end(); ++itr) {
	  fprintf(fp,"%s\n",itr->c_str());
	 }
	 fclose(fp);
	
	 if(!(fp=fopen(totalAlignFile[mm].c_str(),"w")))exit(EXIT_FAILURE);
	 for(int i=0;i<wellList.size();i++){
		 fprintf(fp,"\t%s",wellList[i].c_str());
		}
		fprintf(fp,"\n");
	 for(int i=0;i<geneList.size();i++){
	  fprintf(fp,"%s\t",geneList[i].c_str());
	  for(int j=0;j<NWELLS-1;j++){
				fprintf(fp,"%d\t",count_total[j][erccOffset]);
			}
			fprintf(fp,"%d\n",count_total[NWELLS-1][i+erccOffset]);
		} 
	 fclose(fp);
	 
	 if(!(fp=fopen(umiFile[mm].c_str(),"w")))exit(EXIT_FAILURE);
	 for(int i=0;i<wellList.size();i++){
		 fprintf(fp,"\t%s",wellList[i].c_str());
		}
		fprintf(fp,"\n");  
	 for(int i=0;i<geneList.size();i++){
	  fprintf(fp,"%s\t",geneList[i].c_str());
	  for(int j=0;j<NWELLS-1;j++){
				fprintf(fp,"%d\t",count_umi[j][i+geneOffset]);
			}
			fprintf(fp,"%d\n",count_umi[NWELLS-1][i+geneOffset]);
		}  
	 fclose(fp);
	 if(!(fp=fopen(erccFile[mm].c_str(),"w")))exit(EXIT_FAILURE);
	 for(int i=0;i<wellList.size();i++){
		 fprintf(fp,"\t%s",wellList[i].c_str());
		}
		fprintf(fp,"\n"); 
	 for(int i=0;i<erccList.size();i++){
	  fprintf(fp,"%s",erccList[i].c_str());
	  for(int j=0;j<NWELLS;j++){
				fprintf(fp,"\t%d",count_spike_umi[j][i+erccOffset]);
			}
			fprintf(fp,"\n");
		}
	 fclose(fp);
	 if(!(fp=fopen(erccUMIFile[mm].c_str(),"w")))exit(EXIT_FAILURE);
	 for(int i=0;i<wellList.size();i++){
		 fprintf(fp,"\t%s",wellList[i].c_str());
		}
		fprintf(fp,"\n");
	 for(int i=0;i<erccList.size();i++){
	  fprintf(fp,"%s",erccList[i].c_str());
	  for(int j=0;j<NWELLS;j++){
				fprintf(fp,"\t%d",count_spike_total[j][i+erccOffset]);
			}
			fprintf(fp,"\n");
		}
	 fclose(fp);
	
	 if(!(fp=fopen(wellTotalFile[mm].c_str(),"w")))exit(EXIT_FAILURE);
	 for(int i=0;i<wellList.size();i++){
		 fprintf(fp,"\t%s",wellList[i].c_str());
		}
		fprintf(fp,"\n");
	 
	 fprintf(fp,"Refseq_Total"); 
	 for(int i=0;i<NWELLS;i++){
			fprintf(fp,"\t%d",sumCountsi(count_total,i,geneList.size(),geneOffset));
		}
		fprintf(fp,"\n");
		
	 fprintf(fp,"Refseq_UMI"); 
	 for(int i=0;i<NWELLS;i++){
			fprintf(fp,"\t%d",sumCountsi(count_umi,i,geneList.size(),geneOffset));
		}
		fprintf(fp,"\n");
		
	 fprintf(fp,"Spike_Total"); 
	 for(int i=0;i<NWELLS;i++){
			fprintf(fp,"\t%d",sumCountsi(count_spike_total,i,erccList.size(),erccOffset));
		}
		fprintf(fp,"\n");
		
	 fprintf(fp,"Spike_UMI"); 
	 for(int i=0;i<NWELLS;i++){
			fprintf(fp,"\t%d",sumCountsi(count_spike_umi,i,erccList.size(),erccOffset));
		}
		fprintf(fp,"\n");
	 fclose(fp);
	}
 //clean up
	delete[]count_spike_total[0];
	delete[]count_spike_umi[0];   
	delete[]count_total[0];
	delete[]count_umi[0];   

 return 1;		
	 
}

bool multiGeneHit(vector<string> &best_list, string gene, unordered_map<string,string> &refseq_to_gene){
	for(int i=0;i<best_list.size();i++)
	 if(!refseq_to_gene.count(best_list[i]) || gene != refseq_to_gene[best_list[i]])return 1;
	return 0; 
}	
void splitStr(char *cstr,const char *delim, vector<string> &items){
	char *p=strtok(cstr,delim);
	items.resize(0);		
	while(p){
		items.push_back(string(p));
		p=strtok(0,delim);
	}	
}
void splitStr(string str,const char *delim, vector<string> &items){
	//have to make a copy of str in this case
	if(str.size()<1024){
		char cstr[1024];
		strcpy(cstr,str.c_str());
	 char *p=strtok(cstr,delim);
	 items.resize(0);
	 while(p){
	 	items.push_back(string(p));
	 	p=strtok(0,delim);
	 }
	}
	else{
		char *cstr=(char*) malloc(str.size()+1);
		strcpy(cstr,str.c_str());
	 char *p=strtok(cstr,delim);
	 items.resize(0);
	 while(p){
	 	items.push_back(string(p));
	 	p=strtok(0,delim);
	 }
	 free(cstr);	
	}	 	
}
string splitStrIndex(string str,const char *delim, int index){
	//have to make a copy of str in this case
	vector <string> temp;
	int n=0;
	if(str.size()<1024){
		char cstr[1024];
		strcpy(cstr,str.c_str());
	 char *p=strtok(cstr,delim);
	 while(p){
			if(index < 0) temp.push_back(string(p));
			else if(n==index){
				return string(p);
			}	 
			n++;
	 	p=strtok(0,delim);
	 }
	 if(index <0 && index+temp.size() >=0){
			return temp[index+temp.size()];
		}	
	}
	else{
		char *cstr=(char*) malloc(str.size()+1);
		strcpy(cstr,str.c_str());
	 char *p=strtok(cstr,delim);
	 while(p){
			if(index < 0) temp.push_back(string(p));
			else if(n==index){
				return string(p);
			}	 
			n++;
	 	p=strtok(0,delim);
	 }
	 free(cstr);	
	 if(index <0 && index+temp.size() >=0){
			return temp[index+temp.size()];
		}	
	}	 	
}

string readWells(char *filename, unordered_map<string,unsigned int> &well_to_index,vector<string> &wellList){
	FILE *fp=fopen(filename,"r");
	if(!fp)exit(EXIT_FAILURE);
	char line[64],id[64],well[64],seq[64];
	unsigned int k=0;
	while(fgets(line,sizeof(line),fp)){
	 sscanf(line,"%s %s %s",id,well,seq);
	 wellList.push_back(string(id)+"_"+string(well));
	 well_to_index[string(well)]=k;
	 well_to_index[string(id)+"_"+string(well)]=k++;
	}
	fclose(fp);
	return string(id);	
}
	
void readERCC(char *filename, vector<string> &erccList){
	FILE *fp=fopen(filename,"r");
	if(!fp)exit(EXIT_FAILURE);
	char line[64]; //max line width is 50
	while(fgets(line,sizeof(line),fp)){
	 if(line[0] == '>'){ //fasta file - just want comment line without carat
			char temp[64];
			sscanf(line+1,"%s",temp);
	  erccList.push_back(string(temp));
		}
	}
	fclose(fp);	
}

void readRefseq(char *filename, unordered_map<string,string> &refseq_to_gene, vector<string> &geneList){
	FILE *fp=fopen(filename,"r");
	if(!fp)exit(EXIT_FAILURE);
	char line[1024]; //max line width is 843
	char gene[256],refseq[1024];
	while(fgets(line,sizeof(line),fp)){
	 sscanf(line,"%s %s",gene,refseq);
	 geneList.push_back(string(gene));
	 //replace first comma with 0
	 char *p=strtok(refseq,",");
	 while(p){
	  refseq_to_gene[string(p)]=string(gene);
	  p=strtok(0,",");
		}
	}
	fclose(fp);	
}	

bool polyACheck(string &sequence){
	if(sequence.size() < 20) return 0;
	const char *c=sequence.c_str()+sequence.size()-20;
	while(*c){
		if(*c != 'A') return 0;
		c++;
	}
	return 1;	
}
	
		
unsigned int sumCounts(unsigned int **counts,unsigned int size1, unsigned int size2, int offset){
	unsigned int sum=0;
	for(int i=0;i<size1;i++)
	 for(int j=offset;j<size2+offset;j++)
	  sum+=counts[i][j];
	return sum;
}
		
unsigned int sumCountsi(unsigned int **counts,unsigned int i, unsigned int size2, int offset){
	unsigned int sum=0;
	for(int j=offset;j<size2+offset;j++)
	 sum+=counts[i][j];
	return sum;
}
