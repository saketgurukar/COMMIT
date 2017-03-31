#include <unistd.h>

#include "SmartPtr.cpp"
#include "CommitUtil.h"
#include "GraphToSeqDB.h"
#include "Motifs.h"
#include "spanNfvm.h"

#define NDEBUG 0

int main(int argc, char *argv[]) {

	if (argc != 5) {
		cout << " Usuage: ./CommitTopk <filepath> <motifSizeUpto> <top_k> <deltaTime> <fvm> \n";
		exit(0);
	}
	string file = argv[1]; int	motifSize = atoi(argv[2]);	int top_k = atoi(argv[3]);
	int deltaTime = atoi(argv[4]);	double fvm = atof(argv[5]);

	unsigned long startTime = CommitUtil::getTime();

	SmartPtr<GraphToSeqDB> graphToSeqDb(new GraphToSeqDB(file, deltaTime));
	graphToSeqDb->convertGraphToSequences();

	Motifs *inGraph = new Motifs(deltaTime, top_k, graphToSeqDb->getSeqIdToGraphMap());
	SpanFvm span(top_k, time, 1, motifSize, fvm, inGraph);

	Pairdata pairdata;
	fstream fseq = CommitUtil::openFile(FILENAME);
	span.run(fseq, pairdata);
	unsigned long endTime = CommitUtil::getTime();

	cout << "\n Run time \n " << (endTime - startTime) / 1000000.0 << " seconds \n ";

}
