#include <string>      
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include "fasthamming.hpp"
using namespace std;
template <class T>class umipanel{
	//T is used to choose between 32 and 64 bit uint
	public:
	 vector<string> sequences; //contains panel of barcodes 
	 vector<string> wells;     //contains well of barcodes 
  vector <acgt<T>> Aseqs;
 	string panelID;
	 unsigned int nBarCodes=0;
	 unsigned int barcodeSize=0;
	 hamming<T> hammingEval;
	 umipanel(){
			nBarCodes=0;
			barcodeSize=0;
		}	
 	umipanel(string fileName){
			//read in panel
 		string line,name,well,sequence;
   ifstream inFile(fileName,ifstream::in);
   while(getline(inFile,line)){
    istringstream iss(line);
    iss >> panelID; iss >> well;iss >> sequence;
    wells.push_back(well);
    sequences.push_back(sequence);
    if(!barcodeSize) barcodeSize=sequences[0].size();
    Aseqs.push_back(acgt<T>(sequence.c_str(),barcodeSize)); 
    nBarCodes++;
 	 }
   inFile.close();

   hammingEval=hamming<uint64_t>(barcodeSize);
	 }
	 unsigned int findClosest(char *sequence,unsigned int maxDist, string &bestWell,string &bestSequence, unsigned int &nBest){
		 unsigned int bestDist=maxDist;
		 unsigned int bestIndex=0;
		 
		 //hammingEval.encodeSequence(sequence,querySequence);
		 nBest=0;
		 for(unsigned int i=0;i<nBarCodes;i++){
				unsigned int dist=0;
		 	//unsigned int dist=hammingEval.hpopcount_l(encodedSequences+i*encodeOffset, querySequence,bestDist);
    if(dist <bestDist){
					nBest=1;
		 		bestDist=dist;
		 		bestIndex=i;
		 	}
		 	else if(nBest && dist == bestDist){
					nBest++;
				}		
		 }
   bestWell=wells[bestIndex];
   bestSequence=sequences[bestIndex];
   return bestDist;
		}
};	

