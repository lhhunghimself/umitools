
#include <zlib.h>  
#include <omp.h>
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
unsigned int sumCounts(unsigned int *counts,unsigned int size1);
void readRefseq(char *filename, unordered_map<string,string> &refseq_to_gene,vector<string> &geneList);
string readWells(char *filename, unordered_map<string,unsigned int> &well_to_index,vector<string> &wellList);
void readERCC(char *filename, vector<string> &erccList);
void splitStr(char *cstr,const char *delim, vector<string> &items);
void splitStr(string str,const char *delim, vector<string> &items);
string splitStrIndex(string str,const char *delim, int index);
bool multiGeneHit(vector<string> &best_list, string gene, unordered_map<string,string> &refseq_to_gene);
unsigned int sumCounts(unsigned int **counts,unsigned int size1, unsigned int size2);
unsigned int sumCountsi(unsigned int **counts,unsigned int i, unsigned int size2);
int main(int argc, char *argv[]){
	if(argc < 8 || argc > 9){
		fprintf(stderr,"%d parms requires params: sample_id, sym2ref, ercc_fasta, barcodes, aligned_dir, dge_dir, nThreads\n",argc);
		for(int i=1;i<argc;i++){
			fprintf(stderr,"%d\t%s\n",i,argv[i]);
		}	
		exit(EXIT_FAILURE);
	}
	int unipos=0;
	int k=1;	
 char *sample_id=argv[k++];
 char *sym2ref=argv[k++];
 char *ercc_fasta=argv[k++];
 char *barcodes=argv[k++];
 char *aligned_dir=argv[k++];
 char *dge_dir=argv[k++];
 int nThreads=atoi(argv[k++]);
	if(argc ==9 )unipos=atoi(argv[k++]);
 
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
 string logFile=string(dge_dir)+"/"+string(sample_id)+".unq.log.dat";
 string unknownFile=string(dge_dir)+"/"+string(sample_id)+".unq.unknown_list";
 string totalAlignFile=string(dge_dir)+"/"+string(sample_id)+".unq.refseq.total.dat";
 string umiFile=string(dge_dir)+"/"+string(sample_id)+".unq.refseq.umi.dat";
 string erccFile=string(dge_dir)+"/"+string(sample_id)+".unq.spike.total.dat";
 string erccUMIFile=string(dge_dir)+"/"+string(sample_id)+".unq.spike.umi.dat";
 string wellTotalFile=string(dge_dir)+"/"+string(sample_id)+".unq.well_summary.dat";
 
 string logFile_mm=string(dge_dir)+"/"+string(sample_id)+".all.log.dat";
 string unknownFile_mm=string(dge_dir)+"/"+string(sample_id)+".all.unknown_list";
 string totalAlignFile_mm=string(dge_dir)+"/"+string(sample_id)+".all.refseq.total.dat";
 string umiFile_mm=string(dge_dir)+"/"+string(sample_id)+".all.refseq.umi.dat";
 string erccFile_mm=string(dge_dir)+"/"+string(sample_id)+".all.spike.total.dat";
 string erccUMIFile_mm=string(dge_dir)+"/"+string(sample_id)+".all.spike.umi.dat";
 string wellTotalFile_mm=string(dge_dir)+"/"+string(sample_id)+".all.well_summary.dat";
 
 int mismatchTol=0;
 int NTol=0;
 
 //umipanel<uint32_t,unsigned char> *barcodePanel=0; //change to <uint32_t,uint16_t> for 384 wells
 umipanel<uint32_t,unsigned char> barcodePanel(string(barcodes),mismatchTol,NTol);  		
	//initialize counters
	unsigned int 
	 //unq counts
	 *count_total_reads = new unsigned int [nThreads],
	 *count_assigned_reads = new unsigned int [nThreads],
	 *count_assigned_aligned_reads = new unsigned int [nThreads],	
	 *count_assigned_mito_reads = new unsigned int [nThreads],
	 *count_assigned_mito_umi = new unsigned int [nThreads],
 	*count_assigned_unknown_reads = new unsigned int [nThreads],
	 *count_assigned_unknown_umi = new unsigned int [nThreads],
  
  //mm *counts
  *count_total_reads_mm = new unsigned int [nThreads],
	 *count_assigned_reads_mm = new unsigned int [nThreads],
	 *count_assigned_aligned_reads_mm = new unsigned int [nThreads],	
	 *count_assigned_mito_reads_mm = new unsigned int [nThreads],
	 *count_assigned_mito_umi_mm = new unsigned int [nThreads],
 	*count_assigned_unknown_reads_mm = new unsigned int [nThreads],
	 *count_assigned_unknown_umi_mm = new unsigned int [nThreads],
  
  //2D arrays for spike ins
  //arrays for pointers
		*count_spike_total[NWELLS],
		*count_spike_umi[NWELLS],   
		*count_total[NWELLS],
		*count_umi[NWELLS],   
		*count_spike_total_mm[NWELLS],
		*count_spike_umi_mm[NWELLS],   
		*count_total_mm[NWELLS],
		*count_umi_mm[NWELLS];
  
  //allocate memory for total
		count_spike_total[0]=new unsigned int[erccList.size()*NWELLS];
		count_spike_umi[0]=new unsigned int[erccList.size()*NWELLS];      
		count_spike_total_mm[0]=new unsigned int[erccList.size()*NWELLS];
		count_spike_umi_mm[0]=new unsigned int[erccList.size()*NWELLS];   
		
		count_total[0]=new unsigned int[geneList.size()*NWELLS];
		count_umi[0]=new unsigned int[geneList.size()*NWELLS];
		count_total_mm[0]=new unsigned int[geneList.size()*NWELLS];
		count_umi_mm[0]=new unsigned int[geneList.size()*NWELLS];
		
	//set to zero
		
		memset(count_spike_total[0],0,erccList.size()*NWELLS*sizeof(unsigned int));
		memset(count_spike_umi[0],0,erccList.size()*NWELLS*sizeof(unsigned int));      
		memset(count_spike_total_mm[0],0,erccList.size()*NWELLS*sizeof(unsigned int));
		memset(count_spike_umi_mm[0],0,erccList.size()*NWELLS*sizeof(unsigned int));   
		
		memset(count_total[0],0,geneList.size()*NWELLS*sizeof(unsigned int));
		memset(count_umi[0],0,geneList.size()*NWELLS*sizeof(unsigned int));
		memset(count_total_mm[0],0,geneList.size()*NWELLS*sizeof(unsigned int));
		memset(count_umi_mm[0],0,geneList.size()*NWELLS*sizeof(unsigned int));
		
		memset(count_total_reads,0,nThreads*sizeof(unsigned int));
	 memset(count_assigned_reads,0,nThreads*sizeof(unsigned int));
	 memset(count_assigned_aligned_reads,0,nThreads*sizeof(unsigned int));	
	 memset(count_assigned_mito_reads,0,nThreads*sizeof(unsigned int));
	 memset(count_assigned_mito_umi,0,nThreads*sizeof(unsigned int));
 	memset(count_assigned_unknown_reads,0,nThreads*sizeof(unsigned int));
	 memset(count_assigned_unknown_umi,0,nThreads*sizeof(unsigned int));
  memset(count_total_reads_mm,0,nThreads*sizeof(unsigned int));
	 memset(count_assigned_reads_mm,0,nThreads*sizeof(unsigned int));
	 memset(count_assigned_aligned_reads_mm,0,nThreads*sizeof(unsigned int));	
	 memset(count_assigned_mito_reads_mm,0,nThreads*sizeof(unsigned int));
	 memset(count_assigned_mito_umi_mm,0,nThreads*sizeof(unsigned int));
 	memset(count_assigned_unknown_reads_mm,0,nThreads*sizeof(unsigned int));
	 memset(count_assigned_unknown_umi_mm,0,nThreads*sizeof(unsigned int));
	
		
  
		//initialize pointers
		for(int i=1;i<NWELLS;i++){
		 count_spike_total[i]=count_spike_total[i-1]+erccList.size();
			count_spike_umi[i]=count_spike_umi[i-1]+erccList.size();      
			count_spike_total_mm[i]=count_spike_total_mm[i-1]+erccList.size();
			count_spike_umi_mm[i]=count_spike_umi_mm[i-1]+erccList.size();   
			count_total[i]=count_total[i-1]+geneList.size();
			count_umi[i]=count_umi[i-1]+geneList.size();
			count_total_mm[i]=count_total_mm[i-1]+geneList.size();
			count_umi_mm[i]=count_umi_mm[i-1]+geneList.size();
		}
		
	//text hash - may change trie later

	//for hashvalue - gives rise to 96*1024*1024*28000 ~ 42 bits - faster than 20 byte text to hash

	const uint64_t barcodeBase=geneList.size()+erccList.size()+1; //add1 for chrM
	const uint64_t wellBase=(1<<20)*barcodeBase;


 //this hash is necessary to join the unknown hashes	 
	unordered_set<string>unknown_set;
 unordered_set<string>unknown_set_w[NWELLS];
 unordered_set<string>unknown_set_mm;
 unordered_set<string>unknown_set_mm_w[NWELLS];
 
#pragma omp parallel for num_threads (nThreads) schedule (dynamic)
 for	(int wellIndex=0; wellIndex<=NWELLS; wellIndex++){
		//hashes
		unordered_set<uint64_t>umi_seen;
  unordered_set<uint64_t>umi_seen_mm;
	 unordered_set<string>unknown_umi_seen;
	 unordered_set<string>unknown_umi_seen_mm;
		
		//find matching files in the separated wells directory
		vector<string> inputFiles;
		char glob_str[1024];
	 glob_t glob_result;
	 string well=(wellIndex >=NWELLS) ? "X" :barcodePanel.wells[wellIndex];
		sprintf(glob_str,"%s/%s/*.sam",aligned_dir,well.c_str());
		glob(glob_str,GLOB_TILDE,NULL,&glob_result);
	 for(int j=0; j<glob_result.gl_pathc; j++){
			inputFiles.push_back(string(glob_result.gl_pathv[j]));
	 }
		for(int i=0;i<inputFiles.size();i++){
			const unsigned int tid=omp_get_thread_num();
			fprintf(stderr,"working on %s\n",inputFiles[i].c_str());
	  FILE *fp =fopen(inputFiles[i].c_str(),"r");
	  if(!fp){
				fprintf(stderr,"unable to open file %s\n", glob_result.gl_pathv[i]);
				continue;
			}
			char line[1024];
	
			while(fgets(line,sizeof(line),fp)){
		  if(line[0] != '@'){
	
					count_total_reads[tid]++;
					count_total_reads_mm[tid]++;;
	 	  if(wellIndex >= NWELLS)continue; //bad counts for statistics are placed at the end

		    
		   //get barcode
					vector <string> items;
					vector <string> tempItems;
					splitStr(line," \t",items);
	 	  splitStr(items[0],":",tempItems);
	 	  string fullBarcode=tempItems[tempItems.size()-1];
	 	  
	 	  //string well=tempItems[tempItems.size()-2];
	
		   string barcode=fullBarcode.substr(6,10); //change if barcode size changes
		   
		   //skip if ambiguous barcode and get unique barcode index from sequence
		   unsigned int barcodeIndex=0;
		   if(ambigCheck(barcode,barcodeIndex))continue;
	     
		   count_assigned_reads[tid]++;
		   count_assigned_reads_mm[tid]++;
		   string aligned_id=items[2];
     unsigned pos = stoi(items[3]);
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
					count_assigned_aligned_reads[tid]++;
					count_assigned_aligned_reads_mm[tid]++;
					
					if (aligned_id.substr(0,4) == "ERCC"){
						const int erccIndex=ercc_to_index[aligned_id];
						uint64_t umi=barcodeIndex*barcodeBase+erccIndex+geneList.size();
						if(unipos) umi = (umi << 16) | (pos >> 6 & 0x03FF);
						if(!best_hits_list.size()){ 
						 count_spike_total[wellIndex][erccIndex]++;						
						 if(!umi_seen.count(umi)){
							 umi_seen.insert(umi);
							 count_spike_umi[wellIndex][erccIndex]++;
						 }
						} 
						count_spike_total_mm[wellIndex][erccIndex]++;
						if(!umi_seen_mm.count(umi)){
							umi_seen_mm.insert(umi);
							count_spike_umi_mm[wellIndex][erccIndex]++;
						}
					}
					else if 	(aligned_id.substr(0,4) == "chrM"){
						uint64_t umi=barcodeIndex*barcodeBase+barcodeBase-1;
						if(unipos) umi = (umi << 16) | (pos >> 6 & 0x03FF);
						if(!best_hits_list.size()){
							count_assigned_mito_reads[tid]++;
						 if(!umi_seen.count(umi)){
						 	count_assigned_mito_umi[tid]++;
						 	umi_seen.insert(umi);
						 }						
						} 
						count_assigned_mito_reads_mm[tid]++;
						if(!umi_seen_mm.count(umi)){
							count_assigned_mito_umi_mm[tid]++;
							umi_seen_mm.insert(umi);
						}
					}
					else if (refseq_to_gene.count(aligned_id)){
						string gene = refseq_to_gene[aligned_id];
						const unsigned int geneIndex =gene_to_index[gene];
						uint64_t umi=barcodeIndex*barcodeBase+geneIndex;
						if(unipos) umi = (umi << 16) | (pos >> 6 & 0x03FF);
						bool multiGene=multiGeneHit(best_hits_list,gene,refseq_to_gene);
						if(!multiGene){
							count_total[wellIndex][geneIndex]++;
						 if(!umi_seen.count(umi)){
							 count_umi[wellIndex][geneIndex]++;
							 umi_seen.insert(umi);
						 }	
						}
						count_total_mm[wellIndex][geneIndex]++;
						if(!umi_seen_mm.count(umi)){
							count_umi_mm[wellIndex][geneIndex]++;
							umi_seen_mm.insert(umi);
						}					
					}
					else{
						string umi=barcodePanel.wells[wellIndex] + barcode + aligned_id;
						if(!best_hits_list.size()){
							count_assigned_unknown_reads[tid]++;
						 if(!unknown_set_w[wellIndex].count(aligned_id)){
						 	unknown_set.insert(aligned_id);
						 }
						 if(!unknown_umi_seen.count(umi)){
						  count_assigned_unknown_umi[tid]++;	
							 unknown_umi_seen.insert(umi);
						 }
						}
						if(!unknown_set_mm.count(aligned_id)){
						 unknown_set_mm_w[wellIndex].insert(aligned_id);
						}
						count_assigned_unknown_reads_mm[tid]++;
						if(!unknown_umi_seen_mm.count(umi)){
						 count_assigned_unknown_umi_mm[tid]++;
							unknown_umi_seen_mm.insert(umi);
						}				
					}		 
			 }
			}
		 fclose(fp);
		 
		}
		//update counts
	
	}//end wells loop

 //combine unknown sets - assume that bad wells are ignored
 for(int w=0; w< NWELLS; w++){
		for (auto itr = unknown_set_w[w].begin(); itr != unknown_set_w[w].end(); ++itr){
		 unknown_set.insert(*itr);
		}
	}
	for(int w=0; w< NWELLS; w++){
		for (auto itr = unknown_set_mm_w[w].begin(); itr != unknown_set_mm_w[w].end(); ++itr){
		 unknown_set_mm.insert(*itr);
		}
	}
 FILE *fp=fopen(logFile.c_str(),"w");
 if(!fp)exit(EXIT_FAILURE);
 fprintf(fp,"Sample_ID\tTotal\tAssigned\tAligned\tSpike_Total\tSpike_UMI\tMito_Total\tMito_UMI\tRefseq_Total\tRefseq_UMI\tUnknown_Total\tUnknown_UMI\n");
 fprintf(fp,"%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
           sample_id,
           sumCounts(count_total_reads,nThreads),
											sumCounts(count_assigned_reads,nThreads),
											sumCounts(count_assigned_aligned_reads,nThreads),
											sumCounts(count_spike_total,NWELLS,erccList.size()),
											sumCounts(count_spike_umi,NWELLS,erccList.size()),
											sumCounts(count_assigned_mito_reads,nThreads),
											sumCounts(count_assigned_mito_umi,nThreads),
											sumCounts(count_total,NWELLS,geneList.size()),
											sumCounts(count_umi,NWELLS,geneList.size()),
											sumCounts(count_assigned_unknown_reads,nThreads),
											sumCounts(count_assigned_unknown_umi,nThreads)
									);
 fclose(fp);
 if(!(fp=fopen(logFile_mm.c_str(),"w")))exit(EXIT_FAILURE);
 fprintf(fp,"Sample_ID\tTotal\tAssigned\tAligned\tSpike_Total\tSpike_UMI\tMito_Total\tMito_UMI\tRefseq_Total\tRefseq_UMI\tUnknown_Total\tUnknown_UMI\n");
 fprintf(fp,"%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
           sample_id,
           sumCounts(count_total_reads_mm,nThreads),
											sumCounts(count_assigned_reads_mm,nThreads),
											sumCounts(count_assigned_aligned_reads_mm,nThreads),		
           sumCounts(count_spike_total_mm,NWELLS,erccList.size()),
											sumCounts(count_spike_umi_mm,NWELLS,erccList.size()),
											sumCounts(count_assigned_mito_reads_mm,nThreads),
											sumCounts(count_assigned_mito_umi_mm,nThreads),
											sumCounts(count_total_mm,NWELLS,geneList.size()),
											sumCounts(count_umi_mm,NWELLS,geneList.size()),
											sumCounts(count_assigned_unknown_reads_mm,nThreads),
											sumCounts(count_assigned_unknown_umi_mm,nThreads)
									);
 fclose(fp);
 
 if(!(fp=fopen(unknownFile.c_str(),"w")))exit(EXIT_FAILURE);
 for (auto itr = unknown_set.begin(); itr != unknown_set.end(); ++itr) {
  fprintf(fp,"%s\n",itr->c_str());
 }
 fclose(fp);
 if(!(fp=fopen(unknownFile_mm.c_str(),"w")))exit(EXIT_FAILURE);
 for (auto itr = unknown_set_mm.begin(); itr != unknown_set_mm.end(); ++itr) {
  fprintf(fp,"%s\n",itr->c_str());
 }
 fclose(fp);
 
 if(!(fp=fopen(totalAlignFile.c_str(),"w")))exit(EXIT_FAILURE);
 for(int i=0;i<wellList.size();i++){
	 fprintf(fp,"\t%s",wellList[i].c_str());
	}
	fprintf(fp,"\n");
 for(int i=0;i<geneList.size();i++){
  fprintf(fp,"%s\t",geneList[i].c_str());
  for(int j=0;j<NWELLS-1;j++){
			fprintf(fp,"%d\t",count_total[j][i]);
		}
		fprintf(fp,"%d\n",count_total[NWELLS-1][i]);
	} 
 fclose(fp);


 if(!(fp=fopen(totalAlignFile_mm.c_str(),"w")))exit(EXIT_FAILURE);
 for(int i=0;i<wellList.size();i++){
	 fprintf(fp,"\t%s",wellList[i].c_str());
	}
	fprintf(fp,"\n");
 for(int i=0;i<geneList.size();i++){
  fprintf(fp,"%s\t",geneList[i].c_str());
  for(int j=0;j<NWELLS-1;j++){
			fprintf(fp,"%d\t",count_total_mm[j][i]);
		}
		fprintf(fp,"%d\n",count_total_mm[NWELLS-1][i]);
	} 
 fclose(fp);
 
 if(!(fp=fopen(umiFile.c_str(),"w")))exit(EXIT_FAILURE);
 for(int i=0;i<wellList.size();i++){
	 fprintf(fp,"\t%s",wellList[i].c_str());
	}
	fprintf(fp,"\n");  
 for(int i=0;i<geneList.size();i++){
  fprintf(fp,"%s\t",geneList[i].c_str());
  for(int j=0;j<NWELLS-1;j++){
			fprintf(fp,"%d\t",count_umi[j][i]);
		}
		fprintf(fp,"%d\n",count_umi[NWELLS-1][i]);
	}  
 fclose(fp);
 if(!(fp=fopen(umiFile_mm.c_str(),"w")))exit(EXIT_FAILURE);
 for(int i=0;i<wellList.size();i++){
	 fprintf(fp,"\t%s",wellList[i].c_str());
	}
	fprintf(fp,"\n"); 
 for(int i=0;i<geneList.size();i++){
  fprintf(fp,"%s\t",geneList[i].c_str());
  for(int j=0;j<NWELLS-1;j++){
			fprintf(fp,"%d\t",count_umi_mm[j][i]);
		}
		fprintf(fp,"%d\n",count_umi_mm[NWELLS-1][i]);
	}
 fclose(fp);
 if(!(fp=fopen(erccFile.c_str(),"w")))exit(EXIT_FAILURE);
 for(int i=0;i<wellList.size();i++){
	 fprintf(fp,"\t%s",wellList[i].c_str());
	}
	fprintf(fp,"\n"); 
 for(int i=0;i<erccList.size();i++){
  fprintf(fp,"%s",erccList[i].c_str());
  for(int j=0;j<NWELLS;j++){
			fprintf(fp,"\t%d",count_spike_umi[j][i]);
		}
		fprintf(fp,"\n");
	}
 fclose(fp);
 if(!(fp=fopen(erccFile_mm.c_str(),"w")))exit(EXIT_FAILURE);
 for(int i=0;i<wellList.size();i++){
	 fprintf(fp,"\t%s",wellList[i].c_str());
	}
	fprintf(fp,"\n"); 
 for(int i=0;i<erccList.size();i++){
  fprintf(fp,"%s",erccList[i].c_str());
  for(int j=0;j<NWELLS;j++){
			fprintf(fp,"\t%d",count_spike_umi_mm[j][i]);
		}
		fprintf(fp,"\n");
	}
 fclose(fp);   
 if(!(fp=fopen(erccUMIFile.c_str(),"w")))exit(EXIT_FAILURE);
 for(int i=0;i<wellList.size();i++){
	 fprintf(fp,"\t%s",wellList[i].c_str());
	}
	fprintf(fp,"\n");
 for(int i=0;i<erccList.size();i++){
  fprintf(fp,"%s",erccList[i].c_str());
  for(int j=0;j<NWELLS;j++){
			fprintf(fp,"\t%d",count_spike_total[j][i]);
		}
		fprintf(fp,"\n");
	}
 fclose(fp);
 if(!(fp=fopen(erccUMIFile_mm.c_str(),"w")))exit(EXIT_FAILURE);
 for(int i=0;i<wellList.size();i++){
	 fprintf(fp,"\t%s",wellList[i].c_str());
	}
	fprintf(fp,"\n"); 
 for(int i=0;i<erccList.size();i++){
  fprintf(fp,"%s",erccList[i].c_str());
  for(int j=0;j<NWELLS;j++){
			fprintf(fp,"\t%d",count_spike_total_mm[j][i]);
		}
		fprintf(fp,"\n");
	}
 fclose(fp);
 if(!(fp=fopen(wellTotalFile.c_str(),"w")))exit(EXIT_FAILURE);
 for(int i=0;i<wellList.size();i++){
	 fprintf(fp,"\t%s",wellList[i].c_str());
	}
	fprintf(fp,"\n");
 
 fprintf(fp,"Refseq_Total"); 
 for(int i=0;i<NWELLS;i++){
		fprintf(fp,"\t%d",sumCountsi(count_total,i,geneList.size()));
	}
	fprintf(fp,"\n");
	
 fprintf(fp,"Refseq_UMI"); 
 for(int i=0;i<NWELLS;i++){
		fprintf(fp,"\t%d",sumCountsi(count_umi,i,geneList.size()));
	}
	fprintf(fp,"\n");
	
 fprintf(fp,"Spike_Total"); 
 for(int i=0;i<NWELLS;i++){
		fprintf(fp,"\t%d",sumCountsi(count_spike_total,i,erccList.size()));
	}
	fprintf(fp,"\n");
	
 fprintf(fp,"Spike_UMI"); 
 for(int i=0;i<NWELLS;i++){
		fprintf(fp,"\t%d",sumCountsi(count_spike_umi,i,erccList.size()));
	}
	fprintf(fp,"\n");
 fclose(fp);
 
 if(!(fp=fopen(wellTotalFile_mm.c_str(),"w")))exit(EXIT_FAILURE); 
 for(int i=0;i<wellList.size();i++){
	 fprintf(fp,"\t%s",wellList[i].c_str());
	}
	fprintf(fp,"\n");
 
 fprintf(fp,"Refseq_Total"); 
 for(int i=0;i<NWELLS;i++){
		fprintf(fp,"\t%d",sumCountsi(count_total_mm,i,geneList.size()));
	}
	fprintf(fp,"\n",sumCountsi(count_total_mm,NWELLS-1,geneList.size()));
	
 fprintf(fp,"Refseq_UMI"); 
 for(int i=0;i<NWELLS;i++){
		fprintf(fp,"\t%d",sumCountsi(count_umi_mm,i,geneList.size()));
	}
	fprintf(fp,"\n");
	
 fprintf(fp,"Spike_Total"); 
 for(int i=0;i<NWELLS;i++){
		fprintf(fp,"\t%d",sumCountsi(count_spike_total_mm,i,erccList.size()));
	}
	fprintf(fp,"\n");
	
 fprintf(fp,"Spike_UMI"); 
 for(int i=0;i<NWELLS;i++){
		fprintf(fp,"\t%d",sumCountsi(count_spike_umi_mm,i,erccList.size()));
	}
	fprintf(fp,"\n");
 fclose(fp);
 	 //clean up
	delete[]count_spike_total[0];
	delete[]count_spike_umi[0];   
	delete[]count_total[0];
	delete[]count_umi[0];   
	delete[]count_spike_total_mm[0];
	delete[]count_spike_umi_mm[0];   
	delete[]count_total_mm[0];
	delete[]count_umi_mm[0];
		delete[] count_total_reads;
	 delete[] count_assigned_reads;
	 delete[] count_assigned_aligned_reads;	
	 delete[] count_assigned_mito_reads;
	 delete[] count_assigned_mito_umi;
 	delete[] count_assigned_unknown_reads;
	 delete[] count_assigned_unknown_umi;
  delete[] count_total_reads_mm;
	 delete[] count_assigned_reads_mm;
	 delete[] count_assigned_aligned_reads_mm;	
	 delete[] count_assigned_mito_reads_mm;
	 delete[] count_assigned_mito_umi_mm;
 	delete[] count_assigned_unknown_reads_mm;
	 delete[] count_assigned_unknown_umi_mm;
	
 return 1;		
	 
}

bool multiGeneHit(vector<string> &best_list, string gene, unordered_map<string,string> &refseq_to_gene){
	for(int i=0;i<best_list.size();i++)
	 if(!refseq_to_gene.count(best_list[i]) || gene != refseq_to_gene[best_list[i]])return 1;
	return 0; 
}	
void splitStr(char *cstr,const char *delim, vector<string> &items){
	char *save;
	char *p=strtok_r(cstr,delim,&save);
	items.resize(0);		
	while(p){
		items.push_back(string(p));
		p=strtok_r(0,delim,&save);
	}	
}
void splitStr(string str,const char *delim, vector<string> &items){
	//have to make a copy of str in this case
	if(str.size()<1024){
		char cstr[1024];
		strcpy(cstr,str.c_str());
	char *save;
	char *p=strtok_r(cstr,delim,&save);
	 items.resize(0);
	 while(p){
	 	items.push_back(string(p));
		 p=strtok_r(0,delim,&save);
	 }
	}
	else{
		char *cstr=(char*) malloc(str.size()+1);
		strcpy(cstr,str.c_str());
		char *save;
	 char *p=strtok_r(cstr,delim,&save);
	 items.resize(0);
	 while(p){
	 	items.push_back(string(p));
	 	p=strtok_r(0,delim,&save);
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
		char *save;
	 char *p=strtok_r(cstr,delim,&save);
	 while(p){
			if(index < 0) temp.push_back(string(p));
			else if(n==index){
				return string(p);
			}	 
			n++;
	 	p=strtok_r(0,delim,&save);
	 }
	 if(index <0 && index+temp.size() >=0){
			return temp[index+temp.size()];
		}	
	}
	else{
		char *cstr=(char*) malloc(str.size()+1);
		strcpy(cstr,str.c_str());
			char *save;
	 char *p=strtok_r(cstr,delim,&save);
	 while(p){
			if(index < 0) temp.push_back(string(p));
			else if(n==index){
				return string(p);
			}	 
			n++;
	 	p=strtok_r(0,delim,&save);
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
	 char *save;
	 char *p=strtok_r(refseq,",",&save);
	 while(p){
	  refseq_to_gene[string(p)]=string(gene);
	  p=strtok_r(0,",",&save);
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
unsigned int sumCounts(unsigned int *counts,unsigned int size1){
	unsigned int sum=0;
	for(int i=0;i<size1;i++)
	  sum+=counts[i];
	return sum;
}	
		
unsigned int sumCounts(unsigned int **counts,unsigned int size1, unsigned int size2){
	unsigned int sum=0;
	for(int i=0;i<size1;i++)
	 for(int j=0;j<size2;j++)
	  sum+=counts[i][j];
	return sum;
}
		
unsigned int sumCountsi(unsigned int **counts,unsigned int i, unsigned int size2){
	unsigned int sum=0;
	for(int j=0;j<size2;j++)
	 sum+=counts[i][j];
	return sum;
}
