#include "CommitUtil.h"

unsigned long static CommitUtil::getTime(){
	struct timeval t;
	gettimeofday(&t, 0);
	return t.tv_sec * 1000000ULL + t.tv_usec;
}

static fstream& CommitUtil::openFile(string file) {
	fstream fp;
	fp.setf(std::ios_base::unitbuf);
	fp.open(file.c_str(), ios::in | ios::out);
	if (!fp.is_open()) {
		cout << " Input File Not Found \n";
		exit(0);
	}
	return fp;
}

static FILE& CommitUtil::writeMotifToFile(unsigned long name) {
	ostringstream s;
	s << "Mined_Motifs/" << name << ".txt";
	FILE *F = fopen(s.str().c_str(), "wt");
	if (F == NULL) {
		cout << " Didnt open " << s.str();
		exit(0);
	}
	return F;
}
