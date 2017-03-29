#ifndef _INDUCE_H_
#define _INDUCE_H_


#include <iostream>
using namespace std;
#include "Snap.h"

#undef min
#undef max

#include "GraphToSeqDB.h"
#include <cmath>
#include <sstream>
#include "fnv.hh"
#include "nauty.h"
#include <algorithm>
#include <iterator>

#define DEBUGTM 0
#define DEBUG 0
#define PRINT 0
typedef multimap<int,pair<int,int > > timeEdgeList;  // <time,<Srcnid,Destnid >  >
typedef multimap<int,int > timeNewNodeList;
typedef std::map<unsigned long,int > StoreCounts;
typedef std::map< vector<int> , unsigned long > StorepatternCounts;


#define tr(c,it) for(typeof(c.begin()) it= c.begin(); it!=c.end(); it++)

#define trr(c,it) for(typeof(c.rbegin()) it= c.rbegin(); it!=c.rend(); it++)

#define all(c) (c).begin(),(c).end()



class InducedGraph{

	StoreCounts storeCounts;

	StorepatternCounts storepatternCounts;
	fstream finGraph;
	int Stime;
	fstream frd;
	int k;

public:
	GraphToSeqDB *tm;
	InducedGraph(){

	}
	InducedGraph(int t,int vk){
		Stime=t;
		k=vk;
	}
	InducedGraph(int t,int vk,GraphToSeqDB *ttm){
		Stime=t;
		tm=ttm;
		k=vk;
	}


	int setTime(int t){
		Stime=t;

	}
	bool isGrow(){
		return storeCounts.size()>=2*k ?  true : false;
	}
	int getTopKClosedGraphsVonline(vector<int > Ivec,int counts, int k,
			map< int, map<int,int>  > &SeqIdToInd, map<int, vector<int> > & pidToPat);

	int readGraph(int index,PTimeNENet& subgraph);
	bool checkTemporal(PTimeNENet& subgraph,int eid);
	void getSubgraph(const PTimeNENet& Graph, const TIntV& NIdV, PTimeNENet& subgraph);
	unsigned long countGraph(PTimeNENet& subgraph);
	unsigned long blissCall(PTimeNENet& validSubgraph);
	unsigned long nautyCal(PTimeNENet& validSubgraph);
	int writeResult();
	void getTimeConstrainedSubgraph(PTimeNENet& subgraph,PTimeNENet& validSubgraph);
	void addColorNodetheEdges(PTimeNENet& validSubgraph,
			timeNewNodeList timeListNode);
};
#endif
