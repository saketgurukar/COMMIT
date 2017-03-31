#include "Motifs.h"

int Motifs::getTopKClosedGraphsVonline(vector<int> Ivec, int counts, int k, map<int, map<int, int> > &SeqIdToInd,
		map<int, vector<int> > & pidToPat) {
	int lk = 0;
	vector<int> checkPattern;
	PTimeNENet graph;
	lk++;
	traverse(Ivec, iit)
	{
		checkPattern.clear();
		checkPattern = pidToPat[(*iit)];
		int traverse = checkPattern[0];
		typeof(storepatternCounts.begin()) spcit = storepatternCounts.find(checkPattern);
		if (spcit == storepatternCounts.end()) {
			graph = seqId_to_Graph_Map[traverse];
			map<int, int> EdgeID = SeqIdToInd[traverse];
			TIntV subgraphnodes;
			for (typeof(checkPattern.begin()) eit = checkPattern.begin() + 1; eit != checkPattern.end(); eit++) {
				subgraphnodes.Add(EdgeID[(*eit)]);
			}
			PTimeNENet subgraph = graph->GetSubGraph(subgraphnodes);
			if (subgraphnodes.Len() < 1) {
				exit(0);
			}
			int eid = subgraphnodes[0];
			unsigned long name;
			if (checkTemporal(subgraph, eid)) {
				name = countGraph(subgraph);
				storepatternCounts[checkPattern] = name;
			}
		}
	}
}

void Motifs::getTimeConstrainedSubgraph(PTimeNENet& subgraph, PTimeNENet& validSubgraph) {

	int maxNodeid = 0, nextNodeId, maxEdgeId = 0, nextEdgeId;
	timeEdgeList timeListEdge;
	timeNewNodeList timeListNode;
	for (TNodeEdgeNet<TSecTm, TSecTm>::TNodeI ni = subgraph->BegNI(); ni != subgraph->EndNI(); ni++) {
		validSubgraph->AddNode(ni.GetId());
		if (maxNodeid < ni.GetId())
			maxNodeid = ni.GetId();
	}

	nextNodeId = maxNodeid + 1;
	for (TNodeEdgeNet<TSecTm, TSecTm>::TEdgeI ei = subgraph->BegEI(); ei != subgraph->EndEI(); ei++) {
		validSubgraph->AddEdge(ei);
		if (maxEdgeId < ei.GetId())
			maxEdgeId = ei.GetId();
	}
	nextEdgeId = maxEdgeId + 1;

	for (TNodeEdgeNet<TSecTm, TSecTm>::TEdgeI ei = validSubgraph->BegEI(); ei != validSubgraph->EndEI(); ei++) {
		validSubgraph->AddNode(nextNodeId, TSecTm(1)); //Indicates coloured node.
		timeListEdge.insert(timeEdgeList::value_type(nextNodeId, pair<int, int>(ei.GetSrcNId(), ei.GetDstNId())));
		timeListNode.insert(pair<int, int>(ei.GetDat(), nextNodeId));
		validSubgraph->DelEdge(ei.GetId());
		nextNodeId++;
	}
	for (timeEdgeList::iterator it = timeListEdge.begin(); it != timeListEdge.end(); it++) {
		validSubgraph->AddEdge(((*it).second).first, (*it).first, -1); // No need to time info.. still
		validSubgraph->AddEdge((*it).first, ((*it).second).second);
	}
	addColorNodetheEdges(validSubgraph, timeListNode);

}

void Motifs::addColorNodetheEdges(PTimeNENet& validSubgraph, timeNewNodeList timeListNode) {
	timeNewNodeList::iterator current;
	int currentTime;
	current = timeListNode.begin();
	int i = 0;
	while (current != timeListNode.end()) {

		timeNewNodeList::iterator first, second, Pointer;
		i++;
		first = current;
		second = current;
		Pointer = current;
		currentTime = (*current).first;
		do {
			second++;
		} while (second != timeListNode.end() && (*second).first == currentTime);
		first = second;
		if (first == timeListNode.end()) {
			break;
		}
		do {
			second++;
			if (second == timeListNode.end()) {
				break;
			}
		} while ((*second).first == (*first).first);

		while (current != first) {
			Pointer = first;
			while (Pointer != second) {
				validSubgraph->AddEdge((*current).second, (*Pointer).second);
				++Pointer;
			}
			current++;
		}
	}
}

unsigned long Motifs::countGraph(PTimeNENet& subgraph) {

	unsigned long name;
	SmartPtr<GraphIsomorphism> graphIso(new GraphIsomorphism());

	PTimeNENet validSubgraph = TTimeNENet::New();
	getTimeConstrainedSubgraph(subgraph, validSubgraph);   // COLORED GRAPH / TEMPORAL GRAPH
	name = graphIso->nautyCal(validSubgraph);
	StoreCounts::iterator sit;

	sit = storeCounts.find(name);
	if (sit != storeCounts.end()) {
		(*sit).second += 1;
	} else {
		FILE *F = CommitUtil::writeMotifToFile(name);
		for (TNodeEdgeNet<TSecTm, TSecTm>::TEdgeI it = subgraph->BegEI(); it != subgraph->EndEI(); it++) {
			fprintf(F, " %d\t%d\t%d \n", it.GetSrcNId(), it.GetDstNId(), it.GetDat().GetAbsSecs());
		}
		fclose(F);
		storeCounts.insert(StoreCounts::value_type(name, 1));
	}
	return name;

}


bool Motifs::checkTemporal(PTimeNENet& graph, int eid) {

	SmartPtr<GraphToSeqDB> graphToSeqDb(new GraphToSeqDB(deltaTime));
	PTimeNENet smallGraph = graphToSeqDb->getMaximalGraph(eid, true);
	if ((smallGraph->GetNodes() == graph->GetNodes()) && (smallGraph->GetEdges() == graph->GetEdges()))
		return true;
	else
		return false;

}

int Motifs::writeResult() {
	StoreCounts::iterator sit;
	int c = 0;
	multimap<int, unsigned long> countsCanon;
	for (StoreCounts::iterator sit = storeCounts.begin(); sit != storeCounts.end(); sit++) {
		countsCanon.insert(make_pair((*sit).second, (*sit).first));
	}


	bool del = false;
	for (multimap<int, unsigned long>::reverse_iterator mit = countsCanon.rbegin(); mit != countsCanon.rend(); mit++) {
		if (c >= top_k) {
			continue;
		}
		ostringstream s;
		s << "Mined_Motifs/" << (*mit).second << ".txt";
		c++;
		sit = storeCounts.find((*mit).second);
		FILE *F = fopen(s.str().c_str(), "a");
		if (F == NULL) {
			cout << " Didnt open " << s.str();
			exit(0);
		}
		fseek(F, 0, 0);
		fprintf(F, "Count : %d \n", (*mit).first);
		fclose(F);
	}
}
