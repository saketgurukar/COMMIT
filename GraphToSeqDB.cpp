#include "GraphToSeqDB.h"

void GraphToSeqDB::convertGraphToSequences() {
	ReadTimeGraph(file);
	cout << "\n Read graph Graph Nodes " << graph->GetNodes() << " edges " << graph->GetEdges()
			<< "\n Running Temporal Motif Mining ..\n";
	DBScan();
}

long getId(map<long, long> renumGraph, long Id, long &nodeID) {
	if (!renumGraph.count(Id)) {
		renumGraph.insert(make_pair(Id, nodeID++));
	}
	return renumGraph.find(Id)->second;
}

void GraphToSeqDB::ReadTimeGraph(string file) {

	fstream fp = CommitUtil::openFile(file);
	long SrcNId, DstNId, nodeID = 1, edgeID = 1, duration, contactTime;
	map<long, long> renumGraph; // NORMALIZE Node numbers

	while (true) {
		fp >> contactTime;
		if (fp.eof())
			break;
		fp >> duration;
		fp >> SrcNId;
		fp >> DstNId;
		if (SrcNId == DstNId)
			continue;
		long normSrcId = getId(renumGraph, SrcNId, nodeID);
		long normDstId = getId(renumGraph, DstNId, nodeID);

		if (!graph->IsNode(normSrcId))
			graph->AddNode(normSrcId);
		if (!graph->IsNode(normDstId))
			graph->AddNode(normDstId);

		edgeTimeStampEdgeIDMap.insert(TimeStampEdgeIDMap::value_type(contactTime, edgeID));
		isEdgeVisitedMap.insert(IsEdgeVisitedMap::value_type(edgeID, 0));
		graph->AddEdge(normSrcId, normDstId, edgeID, TSecTm(contactTime));
		edgeID++;

	}
	fp.close();
}

void GraphToSeqDB::DBScan() {
	for (TimeStampEdgeIDMap::iterator it = edgeTimeStampEdgeIDMap.begin(); it != edgeTimeStampEdgeIDMap.end(); it++) {
		unsigned int edgeId = (*it).second;
		if (!(*isEdgeVisitedMap.find(edgeId)).second) {
			PTimeNENet smallGraph = getMaximalGraph(edgeId, true);
			if (smallGraph->GetEdges() >= 3) {
				printDegreeSequence (smallGraph);
				seqId_to_Graph_Map[id] = smallGraph;
				id++;
			}
		}
	}
}

void GraphToSeqDB::getNeighborEdges(list<int>& edgeids, int nodeId, long edgeTime, map<int, int> scannedEdge) {

	TNodeEdgeNet<TSecTm, TSecTm>::TNodeI NI = graph->GetNI(nodeId);
	for (int e = 0; e < NI.GetDeg(); e++) {
		if (scannedEdge.count(NI.GetNbrEId(e)))
			continue;
		TNodeEdgeNet<TSecTm, TSecTm>::TEdgeI EI = graph->GetEI(NI.GetNbrEId(e));
		int dat = EI.GetDat();
		int diff = abs(dat - edgeTime);
		if (diff <= deltaTime) {
			edgeids.push_back(EI.GetId());

		}
	}
}

PTimeNENet& GraphToSeqDB::getMaximalGraph(int edgeId, bool deleteEdges) {

	TimeStampEdgeIDMap edgetimeEdgeID;
	list<int> edgeids;
	map<int, int> scannedEdge;
	EdgeTimesMap timelistMap;
	PTimeNENet smallGraph = TTimeNENet::New();

	edgeids.push_front(edgeId);

	while (!edgeids.empty()) {

		int edgeID = edgeids.front();
		edgeids.pop_front();
		if ((*isEdgeVisitedMap.find(edgeID)).second) {
			continue;
		}

		TNodeEdgeNet<TSecTm, TSecTm>::TEdgeI ei = graph->GetEI(edgeID);
		int edgeTime = ei.GetDat();

		if (!smallGraph->IsNode(ei.GetSrcNId())) {
			smallGraph->AddNode(ei.GetSrcNId());
		}
		if (!smallGraph->IsNode(ei.GetDstNId())) {
			smallGraph->AddNode(ei.GetDstNId());
		}

		if (!smallGraph->IsEdge(ei.GetId())) {
			edgetimeEdgeID.insert(TimeStampEdgeIDMap::value_type(edgeTime, edgeID));
			smallGraph->AddEdge(ei.GetSrcNId(), ei.GetDstNId(), ei.GetId(), ei.GetDat());
			EdgeTimesMap::iterator cit = timelistMap.find(pair<int, int>(ei.GetSrcNId(), ei.GetDstNId()));

			if (cit != timelistMap.end()) {
				(*cit).second.push_back(ei.GetDat());
			} else {
				EdgeTimeStamps timeMap;
				timeMap.push_back(ei.GetDat());
				timelistMap.insert(EdgeTimesMap::value_type(pair<int, int>(ei.GetSrcNId(), ei.GetDstNId()), timeMap));
			}
		} else {
			continue;
		}

		scannedEdge.insert(map<int, int>::value_type(pair<int, int>(ei.GetId(), ei.GetDat())));
		getNeighborEdges(edgeids, ei.GetSrcNId(), edgeTime, scannedEdge);
		getNeighborEdges(edgeids, ei.GetDstNId(), edgeTime, scannedEdge);
	}

	if (deleteEdges) {
		for (map<int, int>::iterator mit = scannedEdge.begin(); mit != scannedEdge.end(); mit++) {
			long unsigned edgeid = (*mit).first;
			(*isEdgeVisitedMap.find(edgeid)).second = 1;
			graph->DelEdge(edgeid);
		}
	}

}

void GraphToSeqDB::printDegreeSequence(PTimeNENet& graph) {
	fseq << "";
	for (TimeStampEdgeIDMap::iterator it = edgeTimeStampEdgeIDMap.begin(); it != edgeTimeStampEdgeIDMap.end(); it++) {
		int edgeId = (*it).second;
		int srcDegree = graph->GetNI(graph->GetEI(edgeId).GetSrcNId()).GetDeg();
		int destDegree = graph->GetNI(graph->GetEI(edgeId).GetDstNId()).GetDeg();
		if (srcDegree > destDegree)
			fseq << "" << destDegree << "," << srcDegree << ",";
		else
			fseq << "" << srcDegree << "," << destDegree << ",";

		fseq << edgeId << "," << graph->GetEI(edgeId).GetDat() << " ";
	}
	fseq << "\n";
}
