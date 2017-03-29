#include "InducedGraph.h"

template <typename T>
int printVec1(vector<T> v){
	cout<<"";
	tr(v, vit){
		cout<<(*vit)<<",";
		vit++;
		cout<<(*vit)<<" ";
	}

}

int InducedGraph::getTopKClosedGraphsVonline(vector<int > Ivec, int counts, int k,
		map< int, map<int,int>  > &SeqIdToInd, map<int, vector<int> > & pidToPat) {

/*
	if(storeCounts.size()>=2*k)
		return 100;
*/

	int lk=0;

	vector<int> checkPattern;
	PTimeNENet graph;
	map<int,PTimeNENet>& idToGraph=(*tm).getSeqIdToGraphMap();
	lk++;

	tr(Ivec, iit){

		checkPattern.clear();

		checkPattern=pidToPat[(*iit)];

		int tr=checkPattern[0];


		typeof(storepatternCounts.begin()) spcit=storepatternCounts.find(checkPattern);


		if(spcit==storepatternCounts.end()){

			//	int in=mainTidToIndex[tr];

			graph=idToGraph[tr];
			map<int,int> EdgeID=SeqIdToInd[tr];


			TIntV subgraphnodes;
			if(DEBUG)
				cout<<"\n"<<tr<<"  ";
			for(typeof(checkPattern.begin()) eit=checkPattern.begin()+1;eit!=checkPattern.end();eit++){

				subgraphnodes.Add(EdgeID[(*eit)]);
				if(DEBUG)
					cout<<EdgeID[(*eit)]<<" ";
			}

			if(DEBUGTM)
				cout<<"\n";
			PTimeNENet subgraph=TTimeNENet::New();

			if(DEBUGTM)
				cout<<" getting subgraph \n";
			getSubgraph(graph,subgraphnodes,subgraph);
			if(DEBUGTM)
				cout<<" got subgraph \n";
			//
			if(subgraphnodes.Len()<1){
				cout<<" CHECK ";
				exit(0);
			}

			int eid=subgraphnodes[0];

			if(DEBUGTM){
				cout<<" \n -------------------------- \n";
				for (TNodeEdgeNet<TSecTm, TSecTm>::TEdgeI it = subgraph->BegEI();
						it != subgraph->EndEI(); it++) {
					cout<<" src "<<it.GetSrcNId()<<" dst "<< it.GetDstNId()<<"  id "<<it.GetId()<<" time "<<it.GetDat()<<"\n";
				}
				cout<<" \n -------------------------- \n";
			}

			//		cout<<" TIME IS "<<time;

			unsigned long name;
			if(checkTemporal(subgraph,eid) ) {

				name=countGraph(subgraph);
				storepatternCounts[checkPattern]=name;
				//					cout<<" counting";
			}


		}

	}


}

unsigned long InducedGraph::nautyCal(PTimeNENet& validSubgraph) {
	mhash::fnv<64> fnv_64;
	std::ostringstream s;

	DYNALLSTAT(int,lab1,lab1_sz);
	DYNALLSTAT(int,ptn,ptn_sz);
	DYNALLSTAT(int,orbits,orbits_sz);
	DYNALLSTAT(int,mapp,map_sz);
	DYNALLSTAT(graph,g1,g1_sz);
	DYNALLSTAT(graph,cg1,cg1_sz);

	static DEFAULTOPTIONS_GRAPH(options);
	statsblk stats;

	int n,m,i;
	size_t k;



	options.getcanon = TRUE;
	options.digraph = FALSE;

	n=validSubgraph->GetNodes()*2;


	m = SETWORDSNEEDED(n);
	nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);

	DYNALLOC1(int,lab1,lab1_sz,n,"malloc");
	DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
	DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
	DYNALLOC1(int,mapp,map_sz,n,"malloc");
	DYNALLOC2(graph,g1,g1_sz,n,m,"malloc");
	DYNALLOC1(graph,cg1,cg1_sz,m*(size_t)n,"malloc");

	EMPTYGRAPH(g1,m,n);



	for(TNodeEdgeNet<TSecTm,TSecTm>::TEdgeI vit=validSubgraph->BegEI();vit!=validSubgraph->EndEI();vit++){
		ADDONEEDGE1(g1, (int) vit.GetSrcNId(), (int)vit.GetDstNId(), (int) m);
	}


	densenauty(g1,lab1,ptn,orbits,&options,&stats,m,n,cg1);




	for (long unsigned k = 0; k < m*(size_t)n; ++k){
		s<<cg1[k];
	}

	return (fnv_64(s.str()));


}
int InducedGraph::readGraph(int index,PTimeNENet& graph){
	stringstream s;
	s<<"graphs/";
	s<<index;
	s<<".txt";

	//	cout<<"File is "<<s.str()<<"\n";

	finGraph.open(s.str().c_str(),ios::in);
	if(!finGraph){
		cout<<" Unable to open file";
		exit(0);
	}

	int src,dst,eid;
	long unsigned time;

	while (true) {

		finGraph >> src;
		if (finGraph.eof())
			break;
		finGraph >> dst;
		finGraph >> eid;
		finGraph >> time;

		if(!graph->IsNode(src))
			graph->AddNode(src);

		if(!graph->IsNode(dst) )
			graph->AddNode(dst);


		graph->AddEdge(src,dst, eid , TSecTm(time));

	}
	finGraph.close();

}

void InducedGraph::getTimeConstrainedSubgraph(PTimeNENet& subgraph,	PTimeNENet& validSubgraph) {

	int maxNodeid = 0, nextNodeId, maxEdgeId = 0, nextEdgeId;

	timeEdgeList timeListEdge;
	timeNewNodeList timeListNode;

	for (TNodeEdgeNet<TSecTm, TSecTm>::TNodeI ni = subgraph->BegNI();
			ni != subgraph->EndNI(); ni++) {
		validSubgraph->AddNode(ni.GetId());
		if (maxNodeid < ni.GetId())
			maxNodeid = ni.GetId();
	}

	//	std::cout<<"\n Maximum node "<<maxNodeid<<endl;
	nextNodeId = maxNodeid + 1;

	//	cout<<"\n Nodes and Edges "<<subgraph->GetNodes()<<"  "<<subgraph->GetEdges();

	//valid subgraph contains graph with edges in time period start to end
	for (TNodeEdgeNet<TSecTm, TSecTm>::TEdgeI ei = subgraph->BegEI();
			ei != subgraph->EndEI(); ei++) {


		validSubgraph->AddEdge(ei);
		//cout<<"\n Edge Data "<<ei.GetDat();

		if (maxEdgeId < ei.GetId())
			maxEdgeId = ei.GetId();

	}

	nextEdgeId = maxEdgeId + 1;

	//cout<<"\n Nodes and Edges "<<validSubgraph->GetNodes()<<"   "<<validSubgraph->GetEdges();

	for (TNodeEdgeNet<TSecTm, TSecTm>::TEdgeI ei = validSubgraph->BegEI();
			ei != validSubgraph->EndEI(); ei++) {

		validSubgraph->AddNode(nextNodeId, TSecTm(1)); //Indicates coloured node.
		timeListEdge.insert(
				timeEdgeList::value_type(nextNodeId,
						pair<int, int>(ei.GetSrcNId(), ei.GetDstNId())));
		timeListNode.insert(pair<int, int>(ei.GetDat(), nextNodeId));
		validSubgraph->DelEdge(ei.GetId());
		nextNodeId++;
	}

	for (timeEdgeList::iterator it = timeListEdge.begin();
			it != timeListEdge.end(); it++) {
		validSubgraph->AddEdge(((*it).second).first, (*it).first, -1); // No need to time info.. still
		validSubgraph->AddEdge((*it).first, ((*it).second).second);
	}

	//cout<<"\n Nodes and Edges "<<validSubgraph->GetNodes()<<"   "<<validSubgraph->GetEdges();
	/*
		 cout<<"\n\n\n time \n\n";
		 for(timeNewNodeList::iterator it=timeListNode.begin();it!=timeListNode.end();it++){
		 cout<<"\n Time  "<<(*it).first<<" New Node ID "<<(*it).second;
		 }
	 */
	addColorNodetheEdges(validSubgraph, timeListNode);

}

void InducedGraph::addColorNodetheEdges(PTimeNENet& validSubgraph,
		timeNewNodeList timeListNode) {
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
				//			cout<<"\nTime "<<(*current).first<<" Id "<<(*current).second<<"  Time "<<(*Pointer).first<<" id "<<(*Pointer).second;
				validSubgraph->AddEdge((*current).second, (*Pointer).second);

				++Pointer;
			}
			current++;
		}

	}

}

unsigned long InducedGraph::countGraph(PTimeNENet& subgraph){
	std::ostringstream s;
	unsigned long name;

	/*for (TNodeEdgeNet<TSecTm, TSecTm>::TNodeI it = validSubgraph->BegNI();
			it != validSubgraph->EndNI(); it++) {
		validSubgraph->DelNode(it.GetId());
	}

	map<int,int> renum;
	int lid=0;
	for (TNodeEdgeNet<TSecTm, TSecTm>::TNodeI ni = subgraph->BegNI();
			ni != subgraph->EndNI(); ni++) {
		if(renum.find(ni.GetId())==renum.end()){
			renum[ni.GetId()]=lid;
			validSubgraph->AddNode(lid);
			lid++;
		}
	}
	//	cout<<"---------------------------------------\n ";
	for (TNodeEdgeNet<TSecTm, TSecTm>::TEdgeI ei = subgraph->BegEI();
			ei != subgraph->EndEI(); ei++) {
		validSubgraph->AddEdge(renum[ei.GetSrcNId()],renum[ei.GetDstNId()],ei.GetId(),ei.GetDat());
		//		cout<<renum[ei.GetSrcNId()]<<" "<<renum[ei.GetDstNId()]<<" "<<ei.GetId()<<" "<<ei.GetDat()<<"\n";
	}*/
	//	cout<<"---------------------------------------\n ";



	PTimeNENet validSubgraph=TTimeNENet::New();
	getTimeConstrainedSubgraph(subgraph,validSubgraph);   // COLORED GRAPH / TEMPORAL GRAPH

	name =nautyCal(validSubgraph);

	//cout<<"\t "<<name<<" \t";
	StoreCounts::iterator sit;
	sit = storeCounts.find(name);

	if (sit != storeCounts.end()) {
		(*sit).second += 1;
	} else {

		ostringstream s;
		s << "Mined_Motifs/";
		s << name;
		s << ".txt";

		FILE *F = fopen(s.str().c_str(), "wt");
		if (F == NULL) {
			cout << " Didnt open " << s.str();
			exit(0);
		}

		for (TNodeEdgeNet<TSecTm, TSecTm>::TEdgeI it = subgraph->BegEI();
				it != subgraph->EndEI(); it++) {
			fprintf(F, " %d\t%d\t%d \n", it.GetSrcNId(), it.GetDstNId(),
					it.GetDat().GetAbsSecs());
		}
		fclose(F);

		storeCounts.insert(StoreCounts::value_type(name, 1));
	}

	return name;

}

void InducedGraph::getSubgraph(const PTimeNENet& Graph, const TIntV& NIdV,
		PTimeNENet& subgraph) {


	map<int,int> renum;
	int lid=0;

	for (int n = 0; n < NIdV.Len(); n++) {

		TNodeEdgeNet<TSecTm, TSecTm>::TEdgeI EI=Graph->GetEI(NIdV[n]);

		if(renum.find(EI.GetSrcNId())==renum.end()){
			renum[EI.GetSrcNId()]=lid;
			subgraph->AddNode(lid);
			lid++;
		}
		if(renum.find(EI.GetDstNId())==renum.end()){
			renum[EI.GetDstNId()]=lid;
			subgraph->AddNode(lid);
			lid++;
		}

		subgraph->AddEdge(renum[EI.GetSrcNId()],renum[EI.GetDstNId()], EI.GetId() ,EI.GetDat());

	}
}

bool InducedGraph::checkTemporal(PTimeNENet& graph,int eid){

	list<int> edgeids;

	map<int, int> scannedEdge;
	unsigned int start = 0, end = 0;
	bool first = true;

	THash<TInt, TInt> nodeRenum;
	int nId = 1;

	edgeids.push_front(eid);
	int edgeTime;

	TNodeEdgeNet<TSecTm, TSecTm>::TNodeI NI;
	TNodeEdgeNet<TSecTm, TSecTm>::TEdgeI EI;
	map<int, int>::iterator mit;

	EdgeTimesMap timelistMap;

	PTimeNENet smallGraph = TTimeNENet::New();


	while (!edgeids.empty()) {

		int edgeID = edgeids.front();
		edgeids.pop_front();


		TNodeEdgeNet<TSecTm, TSecTm>::TEdgeI ei = graph->GetEI(edgeID);

		edgeTime = ei.GetDat();
		if(DEBUGTM){
			cout<<"ID "<<edgeID;
			cout<<" t "<<edgeTime<<"\n";
		}

		if (!nodeRenum.IsKey(ei.GetSrcNId())) {
			nodeRenum.AddDat(ei.GetSrcNId(), nId);
			smallGraph->AddNode(nId);
			nId++;

		}
		if (!nodeRenum.IsKey(ei.GetDstNId())) {
			nodeRenum.AddDat(ei.GetDstNId(), nId);
			smallGraph->AddNode(nId);
			nId++;
		}

		if (!smallGraph->IsEdge(ei.GetId())) {

			smallGraph->AddEdge(nodeRenum.GetDat(ei.GetSrcNId()),
					nodeRenum.GetDat(ei.GetDstNId()), ei.GetId(), ei.GetDat());
		} else {
			continue;
		}

		scannedEdge.insert(
				map<int, int>::value_type(
						pair<int, int>(ei.GetId(), ei.GetDat())));

		NI = graph->GetNI(ei.GetSrcNId());
		for (int e = 0; e < NI.GetDeg(); e++) {

			if (scannedEdge.count(NI.GetNbrEId(e)))
				continue;

			EI = graph->GetEI(NI.GetNbrEId(e));
			int edat=EI.GetDat();

			if(DEBUGTM){
				cout<<"\n IN  "<<NI.GetNbrEId(e);
				cout<<" Time "<<EI.GetDat();
				cout<<" EdgeTime "<<edgeTime;
				cout<<" Diff "<<abs(edgeTime - edat);
			}

			if (abs(edgeTime - edat) <= Stime) {
				if(DEBUGTM)
					cout<<" ADDING "<<NI.GetNbrEId(e);
				edgeids.push_back(EI.GetId());

			}else{
				if(DEBUGTM)
					cout<<" NOT ADDING "<<NI.GetNbrEId(e);
			}
		}

		NI = graph->GetNI(ei.GetDstNId());
		for (int e = 0; e < NI.GetDeg(); e++) {
			if (scannedEdge.count(NI.GetNbrEId(e)))
				continue;

			EI = graph->GetEI(NI.GetNbrEId(e));
			int edat=EI.GetDat();

			if(DEBUGTM){
				cout<<"\n IN  "<<NI.GetNbrEId(e);
				cout<<" Time "<<EI.GetDat();
				cout<<" EdgeTime "<<edgeTime;
				cout<<" Diff "<<abs(edgeTime - edat);
			}

			int diff=abs(edgeTime-edat);

			if (diff <= Stime) {
				if(DEBUGTM)
					cout<<" ADDING "<<NI.GetNbrEId(e);
				edgeids.push_back(EI.GetId());
			}else{
				if(DEBUGTM)
					cout<<" NOT ADDING "<<NI.GetNbrEId(e);
			}
		}
	}

	//	cout<<" new graph nodes "<<smallGraph->GetNodes()<<" graph edges "<<graph->GetNodes()<<"\n" ;

	if(DEBUGTM){
		cout<<" \n -------------------------- \n";
		for (TNodeEdgeNet<TSecTm, TSecTm>::TEdgeI it = smallGraph->BegEI();
				it != smallGraph->EndEI(); it++) {
			cout<<" src "<<it.GetSrcNId()<<" dst "<< it.GetDstNId()<<"  id "<<it.GetId()<<" time "<<it.GetDat()<<"\n";
		}
		cout<<" \n -------------------------- \n";
	}

	if((smallGraph->GetNodes()==graph->GetNodes()) && (smallGraph->GetEdges()==graph->GetEdges())){

		//		cout<<" TRUE \n";
		return true;
	}else{
		//		cout<<" False \n";
		return false;
	}



}

int InducedGraph::writeResult(){
	StoreCounts::iterator sit;
	int c=0;
	multimap<int, unsigned long > countsCanon;
//	cout<<" size "<<storeCounts.size();
	for (StoreCounts::iterator sit = storeCounts.begin();
			sit != storeCounts.end(); sit++) {

		countsCanon.insert(make_pair((*sit).second,(*sit).first));
//		cout<<(*sit).second<<"\n";
	}

	//cout<<" k is "<<k;
	bool del=false;

	for(multimap<int, unsigned long >::reverse_iterator mit=countsCanon.rbegin();mit!=countsCanon.rend();mit++){

		if(c>=k){
			del=true;
		}


		ostringstream s;
		s << "Mined_Motifs/";
		s << (*mit).second;
		s << ".txt";

		if(del){

			remove(s.str().c_str());
			continue;
		}


		//cout<<" c is "<<c;
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
