#include "spanNfvm.h"

unsigned long getticks2() {
	struct timeval t;
	gettimeofday(&t, 0);
	return t.tv_sec * 1000000ULL + t.tv_usec;
}

unsigned long FVMstart,FVMend,getTopkStart,getTopkEnd,minfreStart,minfreEnd,ftot,gtot,mtot;

int PatternID=0;
int patlen;

void stringToInt(string s, int&a, int &b,int &c,int &d){
	unsigned p,p1,p2,p3; int t;
	p = s.find(','); 	stringstream val1(s.substr(0, p));
	p1=s.find(',',p+1); stringstream val2(s.substr(p+1,p1));
	p2=s.find(',',p1+1);	stringstream val3(s.substr(p1+1,p2));

	//	p3=s.find(',',p2+1);	stringstream val4(s.substr(p2+1,p3));
	stringstream val4(s.substr(p2+1));
	val1 >> a; 	val2 >> b; 	val3 >> c; 	val4 >> d; //	val5 >> e;
	if(a>b){		t=a; a=b; b=t;	}
	if(DEBUGTMP)
		cout<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<" \n";
}

int printVec(vector<int> v){
	cout<<"\t";
	traverse(v, vit){
		cout<<(*vit)<<",";
		vit++;
		cout<<(*vit)<<" ";
	}
	cout<<"\n";
}

int SpanFvm::countfrequencyN(vector<int> minVector, vector<int> db) {

	int v1, v2;
	int count=0;
	bool flag=1;
	for (vector<int>::iterator dit = db.begin(); dit != db.end();
			dit++) {

		int i = 0;
		int num = 0;

		flag=1;
		num = minVector.at(i);

		v1 = (*dit);		dit++;		v2 = (*dit);

		if (num > v1)
			flag=0;

		i = i + 1;

		num = minVector.at(i);

		if (num > v2)
			flag=0;

		if(flag)
			count++;
	}
	return count;
}

int SpanFvm::getFloorN(vector<int>& minVector, vector<int> db) {

	minVector.push_back(999999); // Hopefully no single person will interact with 0.9 million people in a network.
	minVector.push_back(999999);

	int v1, v2;
	for (vector<int>::iterator dit = db.begin(); dit != db.end();
			dit++) {

		int i = 0;
		int num = 0;

		num = minVector.at(i);
		v1 = (*dit);
		dit++;
		v2 = (*dit);

		if (num > v1)
			minVector.at(i) = v1;

		i = i + 1;

		num = minVector.at(i);
		if (num > v2)
			minVector.at(i) = v2;

	}

	return 0;
}

/*
int SpanFvm::featureVMmodi(vector<int> minVector, vector<int>& db,vector<int>& items, int index) {
	if(db.size()==0)
		return 0;

	items.push_back(minVector[0]);
	items.push_back(minVector[1]);


	int a,b,c,d;
	while (index < minVector.size()) {
		vector<int> newdb;
		minVector[index]+=1;
		a=minVector[0];
		b=minVector[1];

		for(vector<int>::iterator it=db.begin();it!=db.end();it++){
			c=(*it);
			it++;
			d=(*it);
			if(a<=c && b<=d){
				newdb.push_back(c);
				newdb.push_back(d);
			}
		}

		if(newdb.size()==0){
			minVector[index]-=1;
			index++;
			continue;
		}

		vector<int> localMinVec;
		getFloorN(localMinVec,newdb);

		if(index==1){
			if(localMinVec[0]>minVector[0]){
				minVector[index]-=1;
				index++;
				continue;
			}
		}

		featureVMmodi(localMinVec,newdb,items,index);
		minVector[index]-=1;
		index++;

	}
	if(minVector[0]==999999)
		items.clear();

}*/

int SpanFvm::featureVM(vector<int> minVector, vector<int>& db,vector<int>& items, int index) {

	//	cout<<" minimum support is "<<min_sup;

	int count=countfrequencyN(minVector,db);
	if(count < min_sup)
		return 0;
	items.push_back(minVector[0]);
	items.push_back(minVector[1]);


	int a,b,c,d;
	while (index < minVector.size()) {
		vector<int> newdb;
		minVector[index]+=1;
		a=minVector[0];
		b=minVector[1];

		for(vector<int>::iterator it=db.begin();it!=db.end();it++){
			c=(*it);
			it++;
			d=(*it);
			if(a<=c && b<=d){
				newdb.push_back(c);
				newdb.push_back(d);
			}
		}

		if(newdb.size()<min_sup){
			minVector[index]-=1;
			index++;
			continue;
		}

		vector<int> localMinVec;
		getFloorN(localMinVec,newdb);

		if(index==1){
			if(localMinVec[0]>minVector[0]){
				minVector[index]-=1;
				index++;
				continue;
			}
		}

		featureVM(localMinVec,newdb,items,index);
		minVector[index]-=1;
		index++;
	}
	if(minVector[0]==999999)
		items.clear();

}

int mapToVec(map<string,bool>& OnceAdd, vector<int>& items){
	int a1,b1;
	traverse(OnceAdd,it){
		unsigned p,p1;
		p = (it)->first.find(',');
		stringstream val1((it)->first.substr(0, p));
		stringstream val2((it)->first.substr(p+1));
		val1>>a1;
		val2>>b1;
		items.push_back(a1);
		items.push_back(b1);
	}

}

int SpanFvm::INSGrow(Pairdata& pairdata,multimap<int,vector<int> >& I,vector<int>& Ivec,vector<int> pattern,
		int a,int b, multimap<int,vector<int> >& Istar,vector<int>& IvecStar){

	int index;
	int lastVistedTran=-1;
	int last_position=0;
	int d1,d2;
	int sizePDb=0; 	vector<int> patEmbed;

	if(DEBUGTMP){
		cout<<"\n Landmarks \n";
		traverse(I, iit){

			cout<<(*iit).first<<"  ";
			copy(beginToEnd((*iit).second),ostream_iterator<int>(cout," "));
			cout<<"\n";
		}
		cout<<"\n";
	}


	traverse(Ivec , mit){

		patEmbed.clear();
		patEmbed=pidToPat[(*mit)];


		if(lastVistedTran==-1){
			//lastVistedTran=(*mit).first;
			lastVistedTran=patEmbed[0];
		}

		if(lastVistedTran != patEmbed[0] ){
			last_position=0;
			lastVistedTran=patEmbed[0];
		}



		//	index=pairdata.tidToDBIndex[((*mit).first)];


		const Transaction &transaction=pairdata.database[patEmbed[0]];

		if(DEBUGTMP){
			cout<<"in "<<transaction.first;
			cout<<"\n";
		}

		int tsize=transaction.second.size();
		const vector<int> &itemset = transaction.second;
		int le=patEmbed.back();
		bool flag=false;
		int lasttime=itemset[le-1];
		int currtime=0;



		for(int iter=le;iter<tsize;){
			d1=itemset[iter++];
			d2=itemset[iter++];
			iter++; // edge id
			currtime=itemset[iter++];

			if(abs(currtime-lasttime)>Sptime){
				if(DEBUGTMP)
					cout<<" NOT ADDING main \n";

				break;
			}

			if(a<=d1 && b<=d2 ){

				/*			if(DEBUGTMP)
					cout<<"\n Index "<<index<<" le "<<le<<" iter "<<iter<<" a= "<<a<<" b="<<b<<" tid "<<(*mit).first<<"\n";
				 */
				if(last_position<iter){
					last_position=iter;
					flag=true;
					sizePDb+=(tsize- iter)+transaction.first;
					break;

				}
			}

		}
		if(flag){
			vector<int> tmp(patEmbed.begin(),patEmbed.end());
			tmp.push_back(last_position);



			typeof(PatToPid.begin()) pit=PatToPid.find(tmp);

			if(pit==PatToPid.end()){

				PatToPid[tmp]=PatID;
				pidToPat[PatID]=tmp;
				IvecStar.push_back(PatID);
				PatID++;
			}else{
				IvecStar.push_back(pit->second);

			}



			//			Istar.insert(pair<int,vector<int> >((*mit).first,tmp));
		}
	}


	if(DEBUGTMP){
		cout<<"\n **** Istar ***** \n";
		traverse(Istar, iit){

			cout<<(*iit).first<<"  ";
			copy(beginToEnd((*iit).second),ostream_iterator<int>(cout," "));
			cout<<"\n";
		}
		cout<<"\n ********* ";
		cout<<"\n";
	}
	//exit(0);
	return sizePDb;

}

int SpanFvm::printClosedPatterns(Pairdata& pairdata){
	/*

	int id;
	cout<<"\n ------- Closed Pattern ---------------- \n";
	tr(closedtab, hit){


		tr((*hit).second,pt){
			cout<<"Support => "<<(*hit).first<<"  ";

			//copy(all((*pt)),ostream_iterator<int>(cout," "));
			printVec((*pt));

			cout<<"\n";
			id=patternToId[(*pt)];
			multimap<int,vector<int> > Embeddings = idToEmbedding[id];

			if(DEBUGTMP){

				cout<<"\n **** Present in this edges ***** \n";
				tr(Embeddings, iit){

					//cout<<(*iit).first<<"  .. ";
					int in=SeqIdToIndex[(*iit).first];
					cout<<in<<" :: ";
					map<int,int> EdgeID = pairdata.indexToEdgeID[in];
					//copy(all((*iit).second),ostream_iterator<int>(cout," "));

					tr((*iit).second,vit){
						cout<<EdgeID[(*vit)]<<" ";
					}

					cout<<"\n";
				}
				cout<<"\n ********* ";
				cout<<"\n";
				cout<<"\n";
			}





		}
		cout<<"\n";
	}
	cout<<"\n -----------------------\n";
	 */

}

void SpanFvm::getRelation(vector<int> pattern1,vector<int> pattern2,relation& rel ){
	rel=notrel;
	int a,b,c,d,j,i;

	if(pattern1.size()<pattern2.size()){

		bool incI=true;
		j=0;
		for( i=0;i<pattern1.size() && j<pattern2.size();){
			a=pattern1[i];			c=pattern2[j];
			if(incI)
				i++;
			j++;
			b=pattern1[i];			d=pattern2[j];

			//			cout<<a<<","<<b<<" "<<c<<","<<d<<"\n";

			if(a<=c && b <=d){
				incI=true;
				i++; j++;
			}else{
				j++;
				incI=false;
			}
		}
		if(i==pattern1.size())
			rel=in12;
	}else if(pattern1.size()>pattern2.size()){


		bool incI=true;
		j=0;
		for( i=0;i<pattern2.size() && j<pattern1.size();){
			a=pattern2[i];			c=pattern1[j];
			if(incI)
				i++;
			j++;
			b=pattern2[i];			d=pattern1[j];

			if(a<=c && b <=d){
				incI=true;
				i++; j++;
			}else{
				j++;
				incI=false;
			}
		}
		if(i==pattern2.size())
			rel=in21;
	}
	else if(pattern1.size()==pattern2.size()){
		bool contain12=true;
		for(int i=0;i<pattern1.size();i++){
			a=pattern1[i];			c=pattern2[i];
			i++;
			b=pattern1[i];			d=pattern2[i];
			if(a>c || b>d)
				contain12=false;
		}
		if(contain12){
			rel=in12;

		}else{

			bool contain21=true;
			for(int i=0;i<pattern2.size();i++){
				a=pattern1[i];			c=pattern2[i];
				i++;
				b=pattern1[i];			d=pattern2[i];
				if(a<c || b<d)
					contain21=false;
			}
			if(contain21)
				rel=in21;
		}
	}
}

int SpanFvm::minFreVec(Pairdata& pairdata,multimap<int,vector<int> >& I,vector<int>& Ivec,vector<int> pattern,vector<Embeddings>& globalK){

	if(pattern.size()>=2*length){
		cout<<" Size reached \n";
		return 0;
	}

	multimap<int,vector<int> > Istar;
	int d1,d2,a,b;
	bool closed=true;
	Embeddings emb;
	vector<int> newitems;
	vector<int> IvecStar;

	getItems(I,Ivec,pairdata,newitems); // VALID EDGE EXTENSION


	ftot+=FVMend-FVMstart;

	int sz=newitems.size()-1;

	if(DEBUGTMP){
		cout<<" %%% ";
		printVec(newitems);
	}
	int sizePDB=0;


	for(int i=0;i<=sz;i++){

		//		b=newitems[i--];
		//		a=newitems[i];

		a=newitems[i++];
		b=newitems[i];


		sizePDB=0;
		pattern.push_back(a);
		pattern.push_back(b);

		Istar.clear();
		IvecStar.clear();

		if(DEBUGTMP){
			cout<<"\n Growing ";
			printVec(pattern);
		}

		sizePDB=INSGrow(pairdata,I,Ivec,pattern,a,b,Istar,IvecStar);  // GETSUP

		emb.I=Istar;
		emb.pattern=pattern;
		emb.counts=IvecStar.size();
		emb.sizePDB=sizePDB;
		emb.Ivec=IvecStar;

		globalK.push_back(emb);



		pattern.pop_back();
		pattern.pop_back();

	}


	return -1;

}


int SpanFvm::getItems(multimap<int,vector<int> >& I,vector<int>& Ivec,Pairdata& pairdata, vector<int >& newItems){

	int index,d1,d2;
	vector<int > Candidate_items;
	int last=-1;
	int distinct=0;
	int prevTime=0;

	vector<int> patEmbed;
	traverse(Ivec, mit){
		patEmbed.clear();
		patEmbed=pidToPat[(*mit)];


		//		index=pairdata.tidToDBIndex[((*mit).first)];

		//		const Transaction &transaction=pairdata.database[((*mit).first)];

		const Transaction &transaction=pairdata.database[patEmbed[0]];

		if(last==-1){
			last=transaction.first;

		}else if(last==transaction.first){
			continue;
		}else{
			last=transaction.first;
		}
		distinct++;

		if(DEBUGTMP)
			cout<<"\n Transaction is "<<transaction.first;


		if(DEBUGTMP){
			cout<<"in "<<transaction.first;
			cout<<"\n";
		}

		int tsize=transaction.second.size();
		const vector<int> &itemset = transaction.second;

		//	int le=(*mit).second[(*mit).second.size()-1];

		int le=patEmbed.back();

		int lasttime=itemset[le-1];
		int currtime=0;


		for(int iter=le;iter<tsize;){
			d1=itemset[iter++];
			d2=itemset[iter++];
			iter++; // edge id
			currtime=itemset[iter++];

			if(abs(currtime-lasttime)>Sptime){
				if(DEBUGTMP)
					cout<<" NOT ADDING  \n";
				break;
			}

			Candidate_items.push_back(d1);
			Candidate_items.push_back(d2);

		}

	}

	if(distinct>nextSup)
		nextSup=distinct;


	vector<int> minVec;
	getFloorN(minVec, Candidate_items);
	featureVM(minVec, Candidate_items,newItems,0); // 0 is index. 0, projected

}


int SpanFvm::patternExtensionVec(vector<Embeddings>& topk, Pairdata& pairdata){


	nextSup=0;

	vector<Embeddings> globalK;
	vector<Embeddings> projTopK;
	int psize=0;
	traverse(topk,tpit){
		minFreVec(pairdata,(*tpit).I,(*tpit).Ivec,(*tpit).pattern,globalK);  // EXTENSION OF SUBSEQUENCE
		psize=(*tpit).pattern.size()/2 + 1;
	}

	//	cout<<"\n\n PDB size -> "<<nextSup<<"\n";



	double sub=min_sup*0.4;
	min_sup=min_sup-sub;


	cout<<"\n\n Support "<<min_sup<<"\n";

	if(DEBUGTMP){
		cout<<"\n";
		cout<<" ******************************* \n";
		traverse(globalK,lit){
			cout<<"\n -------------- \n ";
			cout<<(*lit).counts<<" ";
			printVec((*lit).pattern);

			/*	tr((*lit).I, iit){
				cout<<(*iit).first<<"  ";
				copy(all((*iit).second),ostream_iterator<int>(cout," "));
				cout<<"\n";
			}*/

			traverse((*lit).Ivec,ivit){
				copy(beginToEnd(pidToPat[(*ivit)]),ostream_iterator<int>(cout," "));
				cout<<"\n";
			}

			cout<<" PDB Size "<<(*lit).sizePDB<<"\n";

			cout<<" -------------- \n ";
		}
		//	cout<<" ******************************* \n";
	}

	cout<<" Pattern size "<<psize<<"\n";

	if(DEBUGG)
		cout<<" GSize "<<globalK.size();


	getTopkVec(globalK,projTopK);


	if(DEBUGG)
		cout<<" PSize "<<projTopK.size()<<"\n";



	if(DEBUGTMP){
		cout<<"\n  After topk .. \n ";
		traverse(projTopK,lit){
			cout<<"\n -------------- \n ";
			cout<<(*lit).counts<<" ";
			printVec((*lit).pattern);

			traverse((*lit).I, iit){
				cout<<(*iit).first<<"  ";
				copy(beginToEnd((*iit).second),ostream_iterator<int>(cout," "));
				cout<<"\n";
			}
			cout<<" -------------- \n ";
		}
		cout<<" ******************************* \n";
	}

	bool ch=false;
	int ret=0;

	int index=0,del=0;

	traverse(projTopK,lit){
		if((*lit).pattern.size()>=6){
			//	in.getTopKClosedGraphsVonline((*lit).Ivec,(*lit).counts,k,SeqIdToInd,pidToPat,idToGraphS);


			ret=(*in).getTopKClosedGraphsVonline((*lit).Ivec,(*lit).counts,k,SeqIdToInd,pidToPat);
			ch=true;
			if(ret==100){ //100 is return code
				cout<<" reached k patterns ";
				return 1;
			}

			/*		if(ret==-1){
			//	cout<<" Inside ";
				impPatterns.erase(impPatterns.begin()+index-del);
				del++;
			}*/

		}

		index++;

		//	printVec((*lit).pattern);
	}
/*
	if((*in).isGrow() && ch){
		cout<<" Full patterns traveresed \n";
		return 1;
	}*/


	typeof(projTopK.begin()) it= projTopK.begin();

	if(it!=projTopK.end()){
		int patternSize=(*it).pattern.size()/2;



		if(patternSize<length)
			patternExtensionVec(projTopK,pairdata);
		else
			cout<<" pattern length reached";

	}

}



int SpanFvm::getClosedVec(vector<Embeddings> & global,vector<Embeddings> & topk){
}

int SpanFvm::getTopkVec(vector<Embeddings> & global,vector<Embeddings> & topk){

	vector<int> chPattern,insidePattern;
	map< int, closeness > closedtable;

	traverse(global,git){


		typeof(closedtable.begin()) it=closedtable.find((*git).counts);

		if(it==closedtable.end()){
			closeness c;
			vector< bool > in;
			vector< int > patterns;
			in.push_back(true);
			int location = git-global.begin();
			patterns.push_back(location);

			c.in=in;
			c.ids=patterns;
			closedtable[(*git).counts]=c;



		}else{


			chPattern=(*git).pattern;
			int i=0;
			bool flag=true;


			if(flag){
				int location = git-global.begin();
				(*it).second.in.push_back(true);
				(*it).second.ids.push_back(location);
			}
		}
	}

	int i=0;

	traverseReverse(closedtable,cit){
		int m=0;
		traverse((*cit).second.ids, idit){
			if((*cit).second.in[m++]==true){
				if(i<k)
					topk.push_back(global[*idit]);
				else
					break;
				i++;
			}
		}
	}
	patlen=0;
	global.clear();
}



void SpanFvm::run(fstream& fp,Pairdata& pairdata ) {

	Pairdata projected; // DATABASE
	Transaction transaction;

	string value;
	int id = 0,ind=0;
	vector<int> Candidate_items,pattern;
	vector<int> itemset;
	map<int,int> indexToEdge;
	int inputIndex;
	bool first=true;
	bool Add=true;


	while (getline(fp, value)) {
		itemset.clear();
		indexToEdge.clear();

		int a,b,c,d,e,lasttime,firsttime;


		first=true;
		Add=true;

		inputIndex=0;

		transaction.second.clear();

		std::stringstream strstr(value);
		transaction.first = id;
		string item;


		while(strstr>>item){

			stringToInt(item,a,b,c,d);
			itemset.push_back(a);			itemset.push_back(b);
			itemset.push_back(c);			itemset.push_back(d);

			if(first){
				firsttime=d;
				first=false;
			}

			if(abs(d-firsttime)>Sptime)
				Add=false;

			inputIndex+=4;

			if(Add){
				Candidate_items.push_back(a);
				Candidate_items.push_back(b);
			}
			indexToEdge[inputIndex]=c;

		}
		transaction.second=itemset;
		pairdata.database.push_back(transaction);
		//	pairdata.indeces.push_back(0);
		//	mainTidToIndex[id]=ind;


		SeqIdToInd[id]=indexToEdge;
		id++;	//ind++;
	}

	if(DEBUGG){
		cout<<" Read Patterns \n";
	}


	multimap<int,vector<int> > I;

	//	multimap<int, storeLandmarks > LandmarksGlobal,topkLandmarks;

	vector<Embeddings> LandmarksGlobalV;


	vector<int> minVec;
	getFloorN(minVec, Candidate_items);




	min_sup=id * fvm;


	nextSup=INT32_MAX;

//	min_sup=1;

	vector<int> items;

	FVMstart=getticks2();

	featureVM(minVec, Candidate_items,items,0); // 0 is index. 0, projected EXTENSIONMINER

	FVMend=getticks2();

	ftot+=FVMend-FVMstart;

	if(DEBUGG){
		//	cout<<"\n --------------";
		//printVec(items);
		cout<<" Size "<<items.size()/2<<"\n";
		//s	cout<<"\n --------------";
	}


	vector<int> patterns, Ivec, pos;;
	vector<int> Landmark;
	int d1,d2,a,b;
	int sizePD=0;

	//	for(int i=items.size()-1;i>=0;i--){

	// ITEM CONTAINS ALL POSSIBLE CLOSED EDGE EXTENSIONS

	for(int i=0;i<items.size();i++){



		a=items[i++];
		b=items[i];


		sizePD=0;
		//	a=items[i++];
		//	b=items[i];
		//storeLandmarks storeObj;

		Embeddings emb;
		projected.clear();
		patterns.clear();
		Ivec.clear();


		patterns.push_back(a);
		patterns.push_back(b);

		if(DEBUGTMP)
			cout<<"\n \t %%%% "<<a<<", "<<b<<" %%%% \n";

		int index=0;
		int dbsize=pairdata.database.size();

		for (unsigned int i = 0; i < dbsize; i++) {

			const Transaction &transaction = pairdata.database[i];
			const vector<int> &itemset = transaction.second;
			int tid=transaction.first;

			int itemsetSize= itemset.size();

			for (int iter = 0; iter <itemsetSize;) {

				d1=itemset[iter++];
				d2=itemset[iter++];
				iter+=2;
				if(a<=d1 && b<=d2){
					pos.clear();

					pos.push_back(tid); // POS ==> INSTANCE OF SEQDB`
					pos.push_back(iter);

					typeof(PatToPid.begin()) pit=PatToPid.find(pos);

					if(pit==PatToPid.end()){

						PatToPid[pos]=PatID;
						pidToPat[PatID]=pos;
						Ivec.push_back(PatID);
						PatID++;
					}else{
						Ivec.push_back(pit->second);
					}
					sizePD+=(itemsetSize-iter)+tid;
				}
			}
		}

		emb.pattern=patterns;
		emb.counts=Ivec.size();
		emb.sizePDB=sizePD;
		emb.Ivec=Ivec;

		LandmarksGlobalV.push_back(emb);

		//		LandmarksGlobal.insert(make_pair(I.size(),storeObj));


		//	cout<<" Inserted ";


		//	minFre(projected,I,pattern,a,b,items,1);

		if(DEBUGTMP){
			cout<<"\n -------------- \n ";
			traverse(I, iit){
				cout<<(*iit).first<<"  ";
				copy(beginToEnd((*iit).second),ostream_iterator<int>(cout," "));
				cout<<"\n";
			}
			cout<<" Compressed \n ";
			traverse(Ivec,ivit){
				copy(beginToEnd(pidToPat[(*ivit)]),ostream_iterator<int>(cout," "));
				cout<<"\n";
			}

			cout<<" -------------- \n ";
		}

		I.clear();
	}



	if(DEBUGTMP){
		cout<<"\n Global size "<<LandmarksGlobalV.size();
	}


	vector<Embeddings> topkV;
	patlen=0;

	getTopkVec(LandmarksGlobalV,topkV);



	if(DEBUGTMP){
		cout<<"\n";
		traverse(topkV,lit){
			cout<<"\n -------------- \n ";
			cout<<(*lit).counts<<" ";
			printVec((*lit).pattern);


			traverse((*lit).Ivec,ivit){
				copy(beginToEnd(pidToPat[(*ivit)]),ostream_iterator<int>(cout," "));
				cout<<"\n";
			}

			cout<<" -------------- \n ";
		}
		//	cout<<" ******************************* \n";
	}



	/*	getTopkStart=getticks2();

	getTopk(LandmarksGlobal,topkLandmarks);

	getTopkEnd=getticks2();

	gtot+=getTopkEnd-getTopkStart;


	if(DEBUGTMP){
		cout<<"\n local size "<<topkLandmarks.size()<<"\n";

		tr(topkLandmarks,lit){
			printVec((*lit).second.pattern);
		}

	}
	 */




	patternExtensionVec(topkV,pairdata); // SEGGROW

	if(DEBUGTMP){
		cout<<" TOPK RESULT \n";
		cout<<"\n";
		traverse(topkResultV,lit){

			cout<<"\n -------------- \n ";
			cout<<(*lit).counts<<" ";
			printVec((*lit).pattern);

			traverse((*lit).Ivec,ivit){
				copy(beginToEnd(pidToPat[(*ivit)]),ostream_iterator<int>(cout," "));
				cout<<"\n";
			}
			cout<<" -------------- \n ";
		}
		//	cout<<" ******************************* \n";
	}



	(*in).writeResult();

	fout.close();
}
