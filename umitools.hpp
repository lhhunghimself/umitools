#include <string>      
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include "fasthamming.hpp"
using namespace std;
template <class T1, class T2>class umipanel{
	//T1 is used to choose between 32 and 64 bit uint
	//T2 is used to choose between unsigned char for 96 weill and uint16_t for 384 wells
	public:
	 T2 *hash;  //for barcode
	 vector<string> sequences; //contains panel of barcodes 
	 vector<string> wells;     //contains well of barcodes 
  vector <acgt<T1>> Aseqs;
 	string panelID;
	 unsigned int nBarcodes=0;
	 unsigned int barcodeSize=0;
	 unsigned int hashSize=0;
	 unsigned int mismatchTol; //how many basepairs difference for matching panel
	 unsigned int NTol; //how many basepairs in query sequence
	 hamming<T1> hammingEval;
	 umipanel(){
			nBarcodes=0;
			barcodeSize=0;
		}	
 	umipanel (string fileName,int _mismatchTol,int _NTol){
			//read in panel
			mismatchTol=_mismatchTol;
			NTol=_NTol;
 		string line,name,well,sequence;
   ifstream inFile(fileName,ifstream::in);
   if(!inFile){
				fprintf(stderr,"unable to read in barcode file %s\n",fileName.c_str());
				exit(EXIT_FAILURE);
			}
   while(getline(inFile,line)){
    istringstream iss(line);
    iss >> panelID; iss >> well;iss >> sequence;
    wells.push_back(well);
    sequences.push_back(sequence);
    if(!barcodeSize) barcodeSize=sequences[0].size();
    Aseqs.push_back(acgt<T1>(sequence.c_str(),barcodeSize)); 
    nBarcodes++;
 	 }
   inFile.close();
   hammingEval=hamming<T1>(barcodeSize);
   hashSize=5;
   for(int i=1;i<barcodeSize;i++){
				hashSize*=5;
			}
   hash=new T2[hashSize];
			memset(hash,0,hashSize*sizeof(T2));
   fillHash();
	 }
	 ~umipanel(){
			if(hash)delete[] hash;
		}	
	 void fillHash(){
			char bp[5]={'N','A','C','G','T'};
			char *seq=new char[barcodeSize];
			char *indices=new char[barcodeSize];
			memset(indices,0,barcodeSize);
			hash[0]=0;
			for(int i=0;i<barcodeSize;i++)
			 seq[i]='N';
   for(int i=1;i<hashSize;i++){
				int numberofNs=0;
				int divisor=5;
				seq[0]=bp[i%divisor];
				int j=1;
				while(i%divisor == 0 && j<barcodeSize){
					indices[j]=(indices[j]+1)%5;
					seq[j]=bp[indices[j]];
					if(!seq[j])numberofNs++;
				 divisor*=5;
				 j++;
				}
				hash[i]=hammingEval.bestMatch(sequences,seq,nBarcodes,mismatchTol,NTol)+1;		
			}
			delete[] seq;	
			delete[] indices;	
		}
		unsigned int bestMatch(const char *query){
			return (unsigned int) hash[hashCode(query)];
		}		 
		unsigned int hashCode(const char *sequence){
			unsigned int code =0;
			int k=1;
			for (int i=0;i<barcodeSize;i++){
				switch (sequence[i]){
					case 'A':
					 code+=k;
					break;
					case 'C':
					 code+=2*k;
					 break;
					case 'G':
					 code+=3*k;
					 break;
				 case 'T':
				  code+=4*k;
				  break;  
				}
				k*=5;				
			}
			return code;	
		}	
};	

