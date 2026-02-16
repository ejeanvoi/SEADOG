#ifndef DL_H
#define DL_H

#include <string.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include "tree.h"

#ifndef MAX
#define MAX 5000
#endif

int GStraceback(int *events, node *genenode, int mappedindex, int *cleft, int *cright, int speciestreesize)
{
	//cout<<genenode->name<<" : "<<mappedindex<<endl;
	genenode->mappedNode=mappedindex;
	genenode->tryMappedNode=mappedindex;

	genenode->eventtype=events[genenode->localkey*speciestreesize+mappedindex];
	genenode->tryEventtype=events[genenode->localkey*speciestreesize+mappedindex];

	if (genenode->isleaf==0)
	{
		GStraceback(events,genenode->left,cleft[genenode->localkey*speciestreesize+mappedindex],cleft,cright,speciestreesize);
		GStraceback(events,genenode->right,cright[genenode->localkey*speciestreesize+mappedindex],cleft,cright,speciestreesize);
	}
	return 0;
}

int postorderSpeciesNode(int *events, int *GSc1, int *GSc2, int *GSc, int *cleft, int *cright, tree *genetree, tree *speciestree, int genetreesize, int speciestreesize, node *genenode, node *currentnode, node **pointers, int GSDuplicatecost, int GSLosscost)
{
	if (genenode->isleaf) {
		return 0;
	}

	int leftlow=MAX;
	int rightlow=MAX;
	int low=MAX;
	int i,j,k;
	int loss=0;

	int left1,right1,left2,right2,left3,right3;

	left1=88;
	left2=88;
	left3=88;
	right1=88;
	right2=88;
	right3=88;

	if(currentnode->isleaf==0)
	{
		postorderSpeciesNode(events,GSc1,GSc2,GSc,cleft,cright,genetree,speciestree,genetreesize,speciestreesize,genenode,currentnode->left,pointers,GSDuplicatecost,GSLosscost);
		postorderSpeciesNode(events,GSc1,GSc2,GSc,cleft,cright,genetree,speciestree,genetreesize,speciestreesize,genenode,currentnode->right,pointers,GSDuplicatecost,GSLosscost);


		for (i=lowestkey(currentnode->left);i<=currentnode->left->key;i++)
		{
			loss=pointers[i]->depth-currentnode->left->depth;
			if (leftlow>GSc[genenode->left->localkey*speciestreesize+i]+loss*GSLosscost)
			{
				leftlow=GSc[genenode->left->localkey*speciestreesize+i]+loss*GSLosscost;
				left1=i;
			}
		}

		for (i=lowestkey(currentnode->right);i<=currentnode->right->key;i++)
		{
			loss=pointers[i]->depth-currentnode->right->depth;
			if (rightlow>GSc[genenode->right->localkey*speciestreesize+i]+loss*GSLosscost)
			{
				rightlow=GSc[genenode->right->localkey*speciestreesize+i]+loss*GSLosscost;
				right1=i;
			}
		}
		low = leftlow+rightlow;
		leftlow=MAX;
		rightlow=MAX;
		for (i=lowestkey(currentnode->left);i<=currentnode->left->key;i++)
		{
			loss=pointers[i]->depth-currentnode->left->depth;
			if (leftlow>GSc[genenode->right->localkey*speciestreesize+i]+loss*GSLosscost)
			{
				leftlow=GSc[genenode->right->localkey*speciestreesize+i]+loss*GSLosscost;
				right2=i;
			}
		}

		for (i=lowestkey(currentnode->right);i<=currentnode->right->key;i++)
		{
			loss=pointers[i]->depth-currentnode->right->depth;
			if (rightlow>GSc[genenode->left->localkey*speciestreesize+i]+loss*GSLosscost)
			{
				rightlow=GSc[genenode->left->localkey*speciestreesize+i]+loss*GSLosscost;
				left2=i;
			}
		}
		if (low>leftlow+rightlow)
		{
			low=leftlow+rightlow;
			left1=left2;
			right1=right2;
		}

		GSc1[genenode->localkey*speciestreesize+currentnode->key]=low;
	}

	leftlow=MAX;
	rightlow=MAX;
	low=MAX;

	for (i=lowestkey(currentnode);i<=currentnode->key;i++)
	{
		loss=pointers[i]->depth-currentnode->depth;
		if (leftlow>GSc[genenode->left->localkey*speciestreesize+i]+loss*GSLosscost)
		{
			left3=i;
			leftlow=GSc[genenode->left->localkey*speciestreesize+i]+loss*GSLosscost;
		}
	}

	for (i=lowestkey(currentnode);i<=currentnode->key;i++)
	{
		loss=pointers[i]->depth-currentnode->depth;
		//cout<<loss<<endl;
		if (rightlow>GSc[genenode->right->localkey*speciestreesize+i]+loss*GSLosscost)
		{
			right3=i;
			rightlow=GSc[genenode->right->localkey*speciestreesize+i]+loss*GSLosscost;
		}
	}


	GSc2[genenode->localkey*speciestreesize+currentnode->key]=rightlow+leftlow+GSDuplicatecost;

	if (GSc2[genenode->localkey*speciestreesize+currentnode->key]>GSc1[genenode->localkey*speciestreesize+currentnode->key])
	{
		cleft[genenode->localkey*speciestreesize+currentnode->key]=left1;
		cright[genenode->localkey*speciestreesize+currentnode->key]=right1;
		GSc[genenode->localkey*speciestreesize+currentnode->key]=GSc1[genenode->localkey*speciestreesize+currentnode->key];
		events[genenode->localkey*speciestreesize+currentnode->key]=1;
	}
	else
	{
		cleft[genenode->localkey*speciestreesize+currentnode->key]=left3;
		cright[genenode->localkey*speciestreesize+currentnode->key]=right3;
		GSc[genenode->localkey*speciestreesize+currentnode->key]=GSc2[genenode->localkey*speciestreesize+currentnode->key];
		events[genenode->localkey*speciestreesize+currentnode->key]=2;
	}
	//cout<<"Cost "<<GSc[genenode->localkey*speciestreesize+currentnode->key]<<" Left "<<cleft[genenode->localkey*speciestreesize+currentnode->key]<<" Right "<<cright[genenode->localkey*speciestreesize+currentnode->key]<<endl;
	return 0;
}

int postorderGeneNode(int *events, int *GSc1, int *GSc2, int *GSc, int *cleft, int *cright, tree *genetree, tree *speciestree, int genetreesize, int speciestreesize, node *currentnode, node **pointers, int GSDuplicatecost, int GSLosscost)
{
	//cout<<"Now Gene node "<<currentnode->key<<endl<<endl<<endl;
   if(currentnode->isleaf==0)
	{
		postorderGeneNode(events,GSc1,GSc2,GSc,cleft,cright,genetree,speciestree,genetreesize,speciestreesize,currentnode->left,pointers,GSDuplicatecost,GSLosscost);
		postorderGeneNode(events,GSc1,GSc2,GSc,cleft,cright,genetree,speciestree,genetreesize,speciestreesize,currentnode->right,pointers,GSDuplicatecost,GSLosscost);
		//cout<<"Return to Gene node "<<currentnode->key<<endl<<endl;
		postorderSpeciesNode(events,GSc1,GSc2,GSc,cleft,cright,genetree,speciestree,genetreesize,speciestreesize,currentnode,speciestree->root,pointers,GSDuplicatecost,GSLosscost);
                return 0;
	}
	else return 0;
}

int DLdynamicalgorithm(tree *genetree, tree *speciestree, int genetreesize, int speciestreesize, node **pointers, int GSDuplicatecost, int GSLosscost)
{
	int i,j;
	node *domainiter;
	node *genenodeiter;
	int *GSc1 = new int[genetreesize*speciestreesize];
	int *GSc2 = new int[genetreesize*speciestreesize];
	int *GSc = new int[genetreesize*speciestreesize];
	int *cleft = new int[genetreesize*speciestreesize];
	int *cright = new int[genetreesize*speciestreesize];
	int *events = new int[genetreesize*speciestreesize];

	for (i=0;i<genetreesize*speciestreesize;i++)
	{
		GSc1[i]=MAX;
		GSc2[i]=MAX;
		GSc[i]=MAX;
		cleft[i]=-1;
		cright[i]=-1;
		events[i]=0;
	}

	for (i=0;i<genetree->leafnodes.size();i++)
	{
		if (genetree->leafnodes[i]->mappedNode < 0 || genetree->leafnodes[i]->mappedNode >= speciestreesize)
		{
			cout<<"ERROR: Gene tree leaf '"<<genetree->leafnodes[i]->name<<"' was not properly mapped to species tree!"<<endl;
			cout<<"       mappedNode="<<genetree->leafnodes[i]->mappedNode<<" (expected 0-"<<speciestreesize-1<<")"<<endl;
			delete[] GSc1;
			delete[] GSc2;
			delete[] GSc;
			delete[] cleft;
			delete[] cright;
			delete[] events;
			return -1;
		}
		GSc[genetree->leafnodes[i]->localkey*speciestreesize+genetree->leafnodes[i]->mappedNode]=0;
	}

	postorderGeneNode(events,GSc1,GSc2,GSc,cleft,cright,genetree,speciestree,genetreesize,speciestreesize,genetree->root,pointers,GSDuplicatecost,GSLosscost);

	int rootlow=MAX;
	int rootmapindex;
	for (i=0;i<speciestreesize;i++)
	{
		if (rootlow>GSc[genetree->root->localkey*speciestreesize+i])
		{
			rootlow = GSc[genetree->root->localkey*speciestreesize+i];
			rootmapindex=i;
		}
	}
	GStraceback(events,genetree->root,rootmapindex,cleft,cright,speciestreesize);
	genetree->root->LCAscore=rootlow;

	delete[] events;
	delete[] GSc1;
	delete[] GSc2;
	delete[] GSc;
	delete[] cleft;
	delete[] cright;
	return rootlow;
}

#endif
