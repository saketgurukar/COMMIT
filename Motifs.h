#include <iostream>
using namespace std;
#include <cmath>
#include <sstream>
#include <algorithm>
#include <iterator>

#include "SmartPtr.cpp"
#include "GraphIsomorphism.h"
#include "Snap.h"
#undef min
#undef max

typedef multimap<int, pair<int, int> > timeEdgeList;
typedef multimap<int, int> timeNewNodeList;
typedef std::map<unsigned long, int> StoreCounts;
typedef std::map<vector<int>, unsigned long> StorepatternCounts;

#define traverse(c,it) for(typeof(c.begin()) it= c.begin(); it!=c.end(); it++)
#define traverseReverse(c,it) for(typeof(c.rbegin()) it= c.rbegin(); it!=c.rend(); it++)
#define beginToEnd(c) (c).begin(),(c).end()

class Motifs {
	map<int, PTimeNENet> seqId_to_Graph_Map;
	int deltaTime, top_k;

	fstream frd;
	fstream finGraph;
	std::map<unsigned long, int> storeCounts;
	std::map<vector<int>, unsigned long> storepatternCounts;
public:

	Motifs(int t, int k, map<int, PTimeNENet> seqGraphMappings) {
		deltaTime = t;
		seqId_to_Graph_Map = seqGraphMappings;
		top_k = k;
	}

	int getTopKClosedGraphsVonline(vector<int> Ivec, int counts, int k, map<int, map<int, int> > &SeqIdToInd,
			map<int, vector<int> > & pidToPat);

	bool checkTemporal(PTimeNENet& subgraph, int eid);
	unsigned long countGraph(PTimeNENet& subgraph);

	int writeResult();
	void getTimeConstrainedSubgraph(PTimeNENet& subgraph, PTimeNENet& validSubgraph);
	void addColorNodetheEdges(PTimeNENet& validSubgraph, timeNewNodeList timeListNode);
};

