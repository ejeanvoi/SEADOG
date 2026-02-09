//#define DEBUG

#include <cstdlib>
#include <ctime>
#include <cctype>

#include <iostream>
#include <istream>
#include <ostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <typeinfo>
#include <list>
#include <ctime>

using namespace std;

// checks if object is of a certain data types
#ifdef DEBUG
	#define EXCEPTIONDATATYPE(function,src,dest) \
		if (typeid(src)!=typeid(dest)) { \
			cerr << "ERROR: internal error in function " << function << " - unexpected data type (" \
			<< typeid(src).name() << " != " << typeid(dest).name() << ')' << endl; \
			exit(1); \
		}

	// checks if object is of one of 2 data types
	#define EXCEPTIONDATATYPE2(function,src,dest1,dest2) \
		if ( (typeid(src)!=typeid(dest1)) && (typeid(src)!=typeid(dest2)) ) { \
			cerr << "ERROR: internal error in function " << function << " - unexpected data type (" \
			<< typeid(src).name() << " != " << typeid(dest1).name() \
			<< " and " \
			<< typeid(src).name() << " != " << typeid(dest2).name() << ')' << endl; \
			exit(1); \
		}
#else
	#define EXCEPTIONDATATYPE(function,src,dest)
	#define EXCEPTIONDATATYPE2(function,src,dest1,dest2)
#endif

// report an error
#define EXCEPTION(msg) {\
		cerr << "ERROR: " << msg << endl; \
		exit(1); \
	}

// report a warning
#define WARNING(msg) {\
		cerr << "WARNING: " << msg << endl; \
	}

// output messages
bool quiet = true;
#define msgout if (!quiet) cout

// version number
string version = "1.0";

// legal characters
string legalChars4Name = "abcdefghijklmnopqrtsuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_";
#define legalChar4Name(c) (\
		((c>='a') && (c<='z')) || \
		((c>='A') && (c<='Z')) || \
		((c>='0') && (c<='9')) | \
		(c=='_') \
	)
#define legalChar4Number(c) (\
		((c>='0') && (c<='9')) || \
		(c=='.') || \
		(c=='-') || \
		(c=='+') \
	)

// encapsulate in '' if nessesary
#define name2newick(name) ((name.find_first_not_of(legalChars4Name) == string::npos) ? name : string('\''+name+'\''))

// called in case of an interrupt signel
void interruptFunction(int signal){
	msgout << endl << "User interrupt... stop immediately, bye!" << endl;
	exit(1);
}

// common initializations
void initialization(){
	// initialize random generator
// 	srand(unsigned(time(NULL)));
}

enum ReRoot {ALL=1, OPT=0};


