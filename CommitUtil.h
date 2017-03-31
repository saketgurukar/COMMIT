#include <iostream>
using namespace std;

#define FILENAME "GraphSequence.txt"

class CommitUtil {
	unsigned long static getTime();
	static fstream& openFile(string file);
	static FILE& writeMotifToFile(unsigned long file);
};
