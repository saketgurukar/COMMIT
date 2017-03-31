#ifndef _SPAN_H_
#define _SPAN_H_

#define DEBUGG 0
#define DEBUGTMP 0

#define traverse(c,it) for(typeof(c.begin()) it= c.begin(); it!=c.end(); it++)

#define traverseReverse(c,it) for(typeof(c.rbegin()) it= c.rbegin(); it!=c.rend(); it++)

#define beginToEnd(c) (c).begin(),(c).end()

#include<iostream>
using namespace std;

#include <vector>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <string>
#include <sstream>
#include <map>
#include <strstream>
#include <fstream>
#include <utility>



#include <time.h>
#include <sys/time.h>
#include <list>
#include <unistd.h>



#include "Motifs.h"


//typedef pair<int,int> edge;
typedef pair<unsigned int, vector<int> > Transaction;
typedef map<int, vector< vector<int> > > ClosedTable;



struct Pairdata {
	vector<Transaction> database;
	vector<unsigned int> indeces;
	vector<unsigned int> lasttimes;
	map<int,int >  tidToDBIndex;

	void clear() {
		database.clear();
		indeces.clear();
	}

	void operator=(const Pairdata& obj){
		database=obj.database;
		indeces=obj.indeces;
		lasttimes=obj.lasttimes;
		tidToDBIndex=obj.tidToDBIndex;
	}
};

struct closeness{
	vector<bool> in;
	vector<int > ids;
};

class storeLandmarks{
public:

	multimap<int,vector<int> > I;
	vector<int > pattern;

};



class Embeddings{
public:
	int counts, sizePDB;
	multimap<int,vector<int> > I;
	vector<int > pattern,Ivec;  // PATTERN + SUBSEQUENCE LIKE 1,3 1,3 1,3
	friend bool operator<(const Embeddings& l,const Embeddings& r){ return l.counts > r.counts; }
};


typedef map<int, vector< multimap<int, storeLandmarks >::iterator > > ClosedTableStore;

enum relation {
	in12,  // 1st pattern is sub pattern, 2nd pattern is super pattern.
	in21, // 2nd pattern is sub pattern, 1st pattern is super pattern.
	notrel
};

class SpanFvm {
	vector<int> pattern;
	int  min_sup;
	double fvm;
	int ppid;
	int Sptime;
	int k,length, PatID;
	fstream fpsp,fout;
	ClosedTable closedtab;
	map<vector<int>, int > patternToId;
	map<int,bool> idSkip ;
	map< int, multimap<int,vector<int> > > idToEmbedding;

	map<int, vector<int> > pidToPat;  // Dictionary
	map<vector<int>, int > PatToPid;
	map< int, vector<int>  > idToEmbed;
	multimap <int, storeLandmarks > topkResult;
	vector<Embeddings> topkResultV;
	map<int,int> mainTidToIndex;
	map< int, map<int,int>  > SeqIdToInd;
	int nextSup;
public:
	Motifs *in;

	SpanFvm(){

	}
	SpanFvm(int m,int p,int min,int len,double f,Motifs *ind){
		min_sup=min;		k=m;
		ppid=1;
		Sptime=p;
		length=len;
		PatID=0;
		fvm=f;
		in= ind;
		nextSup=0;
	}
	vector<Embeddings>& gettopkResultV(){
		return topkResultV;
	}
	multimap <int, storeLandmarks >& getTopkResult(){
		return topkResult;
	}
	ClosedTable& getClosedTable(){
		return closedtab;
	}
	map<vector<int>, int >& getPatternToId(){
		return patternToId;
	}

	map<int,bool>& getIdexists(){
		return idSkip;
	}

	map<int, vector<int> >& getpidToPat(){
		return pidToPat;
	}
	map<vector<int>, int >& getPatToPid(){
		return PatToPid;
	}
	map< int, vector<int>  >& getidToEmbed(){
		return idToEmbed;
	}

	map< int, map<int,int>  >& getSeqIdToInd(){
		return SeqIdToInd;
	}
	map<int,int>& getMainTidToIndex(){
		return mainTidToIndex;
	}


	map< int, multimap<int,vector<int> > >& getIdToEmbedding(){
		return idToEmbedding;
	}
	/*	int featureVMmodi(vector<int> minVector, vector<int>& db,vector<int>& items, int index) ;*/

	int minFreVec(Pairdata& pairdata,multimap<int,vector<int> >& I,vector<int>& Ivec,vector<int> pattern,vector<Embeddings>& globalK);
	int patternExtensionVec(vector<Embeddings>& topk,Pairdata& pairdata);
	int getClosedVec(vector<Embeddings> & global,vector<Embeddings> & topk);
	int getTopkVec(vector<Embeddings> & global,vector<Embeddings> & topk);

	int getItems(multimap<int,vector<int> >& I,vector<int>& Ivec,Pairdata& pairdata, vector<int >& newItems);

	void getRelation(vector<int> pattern1,vector<int> pattern2,relation& rel );
	int printClosedPatterns(Pairdata& pairdata);
	void run(fstream& fp,Pairdata& pairdata);

	int getTopk(multimap<int, storeLandmarks >  & global,	multimap<int, storeLandmarks > & topk);
	int getClosed(multimap<int, storeLandmarks > & global);
	int patternExtension(multimap<int, storeLandmarks >& topk,Pairdata& pairdata);

	int INSGrow(Pairdata& pairdata,multimap<int,vector<int> >& I,vector<int>& Ivec,vector<int> Pattern,
			int a,int b, multimap<int,vector<int> >& Istar,vector<int>& IvecStar);
	int countfrequencyN(vector<int> minVector, vector<int> db);
	int minFre(Pairdata& pairdata,multimap<int,vector<int> >& I,vector<int> Pattern,multimap<int, storeLandmarks >& globalK);
	int printDB(fstream& f,vector<Transaction> &new_database);
	int getFloorN(vector<int>& minVector, vector<int> db);
	int featureVM(vector<int> minVector, vector<int>& db,vector<int>& items, int index);
	int checkAdditionOfItem(unsigned int pos,const vector<int>& sequence,vector<unsigned int>& nodeSet,int t,int CheckS,int CheckD );
	int generatePatterns(Pairdata& pairdata);
	void print_pattern(Pairdata &projected);
	void analyze(fstream& fp);



};
#endif
