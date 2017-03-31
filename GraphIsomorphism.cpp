#include "GraphIsomorphism.h"


unsigned long GraphIsomorphism::nautyCal(PTimeNENet& validSubgraph) {

	mhash::fnv<64> fnv_64;
	std::ostringstream s;

	DYNALLSTAT(int,lab1,lab1_sz);
	DYNALLSTAT(int,ptn,ptn_sz);
	DYNALLSTAT(int,orbits,orbits_sz);
	DYNALLSTAT(int,mapp,map_sz);
	DYNALLSTAT(graph, g1, g1_sz);
	DYNALLSTAT(graph, cg1, cg1_sz);

	static DEFAULTOPTIONS_GRAPH (options);
	statsblk stats;

	int n, m, i;
	size_t k;

	options.getcanon = TRUE;
	options.digraph = FALSE;

	n = validSubgraph->GetNodes() * 2;

	m = SETWORDSNEEDED(n);
	nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);

	DYNALLOC1(int,lab1,lab1_sz,n,"malloc");
	DYNALLOC1(int,ptn,ptn_sz,n,"malloc");
	DYNALLOC1(int,orbits,orbits_sz,n,"malloc");
	DYNALLOC1(int,mapp,map_sz,n,"malloc");
	DYNALLOC2(graph, g1, g1_sz, n, m, "malloc");
	DYNALLOC1(graph, cg1, cg1_sz, m * (size_t) n, "malloc");

	EMPTYGRAPH(g1, m, n);

	for (TNodeEdgeNet<TSecTm, TSecTm>::TEdgeI vit = validSubgraph->BegEI(); vit != validSubgraph->EndEI(); vit++) {
		ADDONEEDGE1(g1, (int) vit.GetSrcNId(), (int) vit.GetDstNId(), (int) m);
	}

	densenauty(g1, lab1, ptn, orbits, &options, &stats, m, n, cg1);

	for (long unsigned k = 0; k < m * (size_t) n; ++k) {
		s << cg1[k];
	}

	return (fnv_64(s.str()));

}
