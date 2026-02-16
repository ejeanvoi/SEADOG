#ifndef tree_H
#define tree_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>

using namespace std;
// Data structure for a tree node
struct node
{
	int DoubleOrder1;
	int DoubleOrder2;
	int key;
	int localkey;	//used in gene trees
	int depth;
	string GeneID;
	int visited;	//used in domain composition.
	int visited2;	//To compare domain composition

	bool valid;	//in the case of deleting nodes.
	bool isfly;	//in the case of deleting nodes.
	bool isleaf;
	bool isroot;
	bool isleft;
	bool istransfer;
	bool isreceiver;
	bool issensitive;
	bool isfakeleaf;

	node *itsroot;
	node *left;
	node *right;
	node *parent;
	node *sibling;
	//node *receiver;	// for the gene node

	char *name;
	string treename;

	int eventtype;
	char *mappedSpecies;	//prefix of the genename;
	//char *Uniprotname;
	int mappedGeneFamily;	//used in domain tree.
	int mappedNode;
	int mappedtoSpecies;	// used in domain trees for domain contents.

	int LCAscore;	// used in gene trees

	int tryMappedNode;	// Used in Local Search for iterations
	int tryEventtype;
	vector<node*> linked;

	// Constructor to initialize fields
	node() : DoubleOrder1(0), DoubleOrder2(0), key(0), localkey(0), depth(0),
	         visited(-1), visited2(-1), valid(false), isfly(false), isleaf(false),
	         isroot(false), isleft(false), istransfer(false), isreceiver(false),
	         issensitive(false), isfakeleaf(false), itsroot(NULL), left(NULL),
	         right(NULL), parent(NULL), sibling(NULL), name(NULL), eventtype(0),
	         mappedSpecies(NULL), mappedGeneFamily(-1), mappedNode(-1),
	         mappedtoSpecies(-1), LCAscore(0), tryMappedNode(-1), tryEventtype(0) {}
};

struct tree
{
	node* root;
	int index;
	vector<node*> leafnodes;
	vector<string> mappedTrees;
	int lowerbound;
	int type; //0 is domain tree, 1: gene tree, 2: species Tree
};

node* GetSibling(node* thisnode)
{
	if (thisnode==thisnode->parent->left)
		return thisnode->parent->right;
	else return thisnode->parent->left;
}

void DoubleTraverse(node *root, int *start)
{
	root->DoubleOrder1=*start;
	*start=*start+1;
	//root->visited=0;
	if(root->isleaf==0)
	{
		DoubleTraverse(root->left,start);
		DoubleTraverse(root->right,start);
	}
	root->DoubleOrder2=*start;
	*start=*start+1;
}

void DoubleTraverse(node *root, int *start, node **pointers)
{
	root->DoubleOrder1=*start;
	pointers[*start]=root;
	*start=*start+1;
	//root->visited=0;
	if(root->isleaf==0)
	{
		DoubleTraverse(root->left,start,pointers);
		DoubleTraverse(root->right,start,pointers);
		root->DoubleOrder2=*start;
		pointers[*start]=root;
		*start=*start+1;
	}
}

#endif
