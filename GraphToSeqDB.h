#include <iostream>
using namespace std;
#include "CommitUtil.h"
#include "Snap.h" // For Snap library
#undef min
#undef max
#include <vector>
#include <map>
#include <algorithm>
#include <iterator>
#include <list>
#include <fstream>
#include <time.h>
#include <sys/time.h>
#include <cstdlib>
#include <sstream>
#include <utility>


typedef multimap<unsigned int, unsigned int> TimeStampEdgeIDMap;
typedef map<unsigned int, bool> IsEdgeVisitedMap;
typedef list<unsigned long> EdgeTimeStamps;
typedef map<pair<int, int>, EdgeTimeStamps> EdgeTimesMap;

class GraphToSeqDB {
	map<int, PTimeNENet> seqId_to_Graph_Map;
	TimeStampEdgeIDMap edgeTimeStampEdgeIDMap; // MAP of Edge Time n Edge ID
	IsEdgeVisitedMap isEdgeVisitedMap;
	PTimeNENet graph;
	fstream fseq;
	int id, deltaTime;
	string file;

public:
	GraphToSeqDB(int deltaTm) {
		id = 0;
		graph = TTimeNENet::New();
		fseq.open(FILENAME, fstream::out);
		deltaTime = deltaTm;
		file = NULL;
	}
	GraphToSeqDB(string fle, int deltaTm) {
		id = 0;
		graph = TTimeNENet::New();
		fseq.open(FILENAME, fstream::out);
		deltaTime = deltaTm;
		file = fle;
	}
	void convertGraphToSequences();
	void ReadTimeGraph(string file);
	void DBScan();
	PTimeNENet& getMaximalGraph(int edgeId, bool deleteEdges);
	void printDegreeSequence(PTimeNENet& graph);
	void getNeighborEdges(list<int>& edgeids, int nodeId, long edgeTime, map<int, int> scannedEdge);
	map<int, PTimeNENet>& getSeqIdToGraphMap() {
		return seqId_to_Graph_Map;
	}

//	friend class SpanFvm;
};

