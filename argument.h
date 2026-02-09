#ifndef ARGUMENT_H
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include "tree.h"
#include <fstream>
#include <cmath>
#include "libs.h"
#include <map>


using namespace std;


class Argument {
public:
	static map<string, Argument> args;
	string key;
	string value;
	bool visited;

	Argument() {
		key = "";
		value = "";
		visited = false;
	}
	Argument(string key) {
		this->key = key;
		value = "";
		visited = false;
	}
	Argument(string key, string value) {
		this->key = key;
		this->value = value;
		visited = false;
	}

	inline bool operator < (const Argument& ref) const {
		return key < ref.key;
	}
	inline bool operator > (const Argument& ref) const {
		return key > ref.key;
	}
	inline bool operator == (const Argument& ref) const {
		return key == ref.key;
	}

	// adds arguments from the command line
	static void add(int argc, char *argsv[]) {
		if (argc <= 1) return;
		string key = "";
		for (int i=1; i<argc; i++) {
			string str(argsv[i]);
			if (str[0] == '-') {
				if (!key.empty())
				{
					args[key] = Argument(key, "");
					key = "";
				}
				const int pos = str.find("=");
				if (pos != string::npos) {
					key = str.substr(0,pos);
					string value = str.substr(pos+1);
					args[key] = Argument(key, value);
					key = "";
				} else {
					key = str;
				}
			} else {
				string value = str;
				args[key] = Argument(key, value);
				key = "";
			}
		}
		if (!key.empty()) 
		{
			args[key] = Argument(key, "");
			key = "";
		}
		
	}

	// searches for a particular argument
	inline static Argument *find(const string &key) {
		map<string, Argument>::iterator itr = args.find(key);
		if (itr == args.end()) return NULL;
		Argument &arg = itr->second;
		arg.visited = true;
		return &arg;
	}

	// searches for any argument
	inline static Argument *findAny(vector<string> &list) {
		for (vector<string>::iterator itr = list.begin(); itr != list.end(); itr++) {
			Argument *arg = find(*itr);
			if (arg != NULL) return arg;
		}
		return NULL;
	}

	// convert a character string into another datatype
	template<class T> void convert(T &var) const;

	// true if argument has a value
	bool hasValue() const {
		return (!value.empty());
	}

	// return unused arguments
	static vector<Argument*> unusedArgs() {
		vector<Argument*> unused;
		for (map<string, Argument>::iterator itr = args.begin(); itr != args.end(); itr++) {
			if (!itr->second.visited) unused.push_back(&itr->second);
		}
		return unused;
	}
};
map<string, Argument> Argument::args;

// convert a character string into another datatype
template<class T>
void Argument::convert(T &var) const {
	istringstream ist(value);
}
// specialization for string
template<>
void Argument::convert(string &var) const {
	var = value;
}

#endif
