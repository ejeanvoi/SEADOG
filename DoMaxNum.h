#ifndef Max_H
#define Max_H

#include <string.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include "tree.h"
#include <fstream>
#include <cmath>
#include "libs.h"

#define MAX 5000

// Helper function to safely add costs, preventing integer overflow
// Allows costs to grow naturally, only capping at a very high threshold
// to prevent overflow. The algorithm will determine what costs are "valid".
inline int SafeAdd(int a, int b) {
	long long result = (long long)a + (long long)b;
	// Cap at a very high value to prevent overflow, but allow costs to grow
	if (result > 1000000) return 1000000;
	return (int)result;
}

// Overload for three operands
inline int SafeAdd(int a, int b, int c) {
	long long result = (long long)a + (long long)b + (long long)c;
	if (result > 1000000) return 1000000;
	return (int)result;
}

int CheckGeneCost(int* geneflag, int* genestep, node *node1, node *node2,node **speciespointers)
{
	node *lca=LCA(speciespointers[node1->mappedNode],speciespointers[node2->mappedNode]);
	node *iter;
	iter=node1;
	int flag=0;
	while (iter!=NULL)
	{
		if(iter->mappedNode<lca->key || (iter->mappedNode==lca->key && iter!=node1))
		{
			if(iter->isroot==0)
				genestep[iter->key]=max(genestep[iter->key],(speciespointers[iter->mappedNode]->depth-lca->depth));
			else
				genestep[iter->key]=max(genestep[iter->key],2*(speciespointers[iter->mappedNode]->depth-lca->depth));
			
			if(iter->eventtype==1)
			{
				geneflag[iter->key]=1;
				flag=1;
			}
			else if (iter->mappedNode<lca->key)
				flag=1;
			//cout<<iter->name<<" : "<<genestep[iter->key]<<" : "<<geneflag[iter->key]<<endl;
		}
		else break;
		iter=iter->parent;
	}
	iter=node2;
	while (iter!=NULL && iter!=LCA(node1,node2))
	{
		if(iter->mappedNode<lca->key || (iter->mappedNode==lca->key && iter!=node2))
		{
			if(iter->isroot==0)
				genestep[iter->key]=max(genestep[iter->key],(speciespointers[iter->mappedNode]->depth-lca->depth));
			else
				genestep[iter->key]=max(genestep[iter->key],2*(speciespointers[iter->mappedNode]->depth-lca->depth));
			//cout<<iter->name<<" : "<<genestep[iter->key]<<endl;
			if(iter->eventtype==1)
			{
				geneflag[iter->key]=1;
				flag=1;
			}
			else if (iter->mappedNode<lca->key)
				flag=1;
			//cout<<iter->name<<" : "<<genestep[iter->key]<<" : "<<geneflag[iter->key]<<endl;
		}
		else break;
		iter=iter->parent;
	}
	return flag;
}

int ExtraLoss(node *node1, node *node2,node **speciespointers)
{
	if(node1->isleaf==1 || node2->isleaf==1) return MAX;
	int total=0;
	node *lca=LCA(speciespointers[node1->mappedNode],speciespointers[node2->mappedNode]);
	//cout<<"Changing "<<node1->mappedNode<<" "<<node2->mappedNode<<" to "<<lca->key<<endl;
	//cout<<speciespointers[node1->mappedNode]->isleaf<<endl;
	node *iter;
	iter=node1;
	while (iter!=NULL)
	{
		if(iter->mappedNode<lca->key || (iter->mappedNode==lca->key && iter!=node1))
		{
			if(iter->isroot==0)
				total+=speciespointers[iter->mappedNode]->depth-lca->depth;
			else
				total+=2*(speciespointers[iter->mappedNode]->depth-lca->depth);
			if(iter->eventtype==1)
				total+=4;
		}
		else break;
		iter=iter->parent;
	}
	iter=node2;
	while (iter!=NULL && iter!=LCA(node1,node2))
	{
		if(iter->mappedNode<lca->key || (iter->mappedNode==lca->key && iter!=node2))
		{
			if(iter->isroot==0)
				total+=speciespointers[iter->mappedNode]->depth-lca->depth;
			else
				total+=2*(speciespointers[iter->mappedNode]->depth-lca->depth);
			if(iter->eventtype==1)
				total+=4;
		}
		else break;
		iter=iter->parent;
	}
	return total;
}


int MaxPostorderGene(int *inindex, node **re, int genetreeindex, int *events, int *c, int *in, int* cleft, int* cright, tree* domainTree, tree** geneTrees, tree* speciesTree, int genesize, node* domainnode, node* currentnode, node** domainpointers, node** genepointers, node** speciespointers, int twoTreeTransferCost, int OneTreeTransferCost, int domainDuplicationcost, int geneDuplicationcost, int domainLosscost, int geneLosscost)
{
	static int debugCount2850 = 0;
	bool debug2850 = (domainnode->key == 2850 && debugCount2850 < 5);
	if (debug2850) {
	}
	int low=MAX;
	int leftlow=MAX;
	int rightlow=MAX;

	//int loss=0;
	int left1,left2,left3,left4,left5,right1,right2,right3,right4,right5;
	//cout<<"Entering tree "<<genetreeindex<<endl;
	int transferType;
	int j,k;
	int LocalCost=0;
	node *lca;
	if(currentnode->isleaf==0)
	{
		MaxPostorderGene(inindex,re,genetreeindex,events,c,in,cleft,cright,domainTree,geneTrees,speciesTree,genesize,domainnode,currentnode->left,domainpointers,genepointers,speciespointers,twoTreeTransferCost,OneTreeTransferCost,domainDuplicationcost,geneDuplicationcost,domainLosscost,geneLosscost);
		MaxPostorderGene(inindex,re,genetreeindex,events,c,in,cleft,cright,domainTree,geneTrees,speciesTree,genesize,domainnode,currentnode->right,domainpointers,genepointers,speciespointers,twoTreeTransferCost,OneTreeTransferCost,domainDuplicationcost,geneDuplicationcost,domainLosscost,geneLosscost);

		//	Variation Case
		leftlow=SafeAdd(in[domainnode->left->key*genesize+currentnode->left->key], in[domainnode->right->key*genesize+currentnode->right->key]);
		rightlow=SafeAdd(in[domainnode->right->key*genesize+currentnode->left->key], in[domainnode->left->key*genesize+currentnode->right->key]);
		if (debug2850 && debugCount2850 < 3) {
		}

		c[domainnode->key*genesize+currentnode->key]=GetMin(leftlow,rightlow);
		
		if(GetMin(leftlow,rightlow)==leftlow)
		{
			cleft[domainnode->key*genesize+currentnode->key]=inindex[domainnode->left->key*genesize+currentnode->left->key];
			cright[domainnode->key*genesize+currentnode->key]=inindex[domainnode->right->key*genesize+currentnode->right->key];
		}
		else
		{
			cleft[domainnode->key*genesize+currentnode->key]=inindex[domainnode->left->key*genesize+currentnode->right->key];
			cright[domainnode->key*genesize+currentnode->key]=inindex[domainnode->right->key*genesize+currentnode->left->key];
		}

		events[domainnode->key*genesize+currentnode->key]=1;

		leftlow=MAX;
		rightlow=MAX;

		// Duplication
		low=SafeAdd(c[domainnode->left->key*genesize+currentnode->key], in[domainnode->right->key*genesize+currentnode->left->key], 1);
		if(c[domainnode->key*genesize+currentnode->key]>SafeAdd(low, domainDuplicationcost))
		{
			c[domainnode->key*genesize+currentnode->key]=SafeAdd(low, domainDuplicationcost);
			events[domainnode->key*genesize+currentnode->key]=2;
			cleft[domainnode->key*genesize+currentnode->key]=currentnode->key;
			cright[domainnode->key*genesize+currentnode->key]=inindex[domainnode->right->key*genesize+currentnode->left->key];
		}
		low=GetMin(low,SafeAdd(c[domainnode->left->key*genesize+currentnode->key], in[domainnode->right->key*genesize+currentnode->right->key], 1));
		if(c[domainnode->key*genesize+currentnode->key]>SafeAdd(low, domainDuplicationcost))		
		{
			c[domainnode->key*genesize+currentnode->key]=SafeAdd(low, domainDuplicationcost);
			events[domainnode->key*genesize+currentnode->key]=2;
			cleft[domainnode->key*genesize+currentnode->key]=currentnode->key;
			cright[domainnode->key*genesize+currentnode->key]=inindex[domainnode->right->key*genesize+currentnode->right->key];
		}
		low=GetMin(low,SafeAdd(c[domainnode->right->key*genesize+currentnode->key], in[domainnode->left->key*genesize+currentnode->left->key], 1));
		if(c[domainnode->key*genesize+currentnode->key]>SafeAdd(low, domainDuplicationcost))		
		{
			c[domainnode->key*genesize+currentnode->key]=SafeAdd(low, domainDuplicationcost);
			events[domainnode->key*genesize+currentnode->key]=2;
			cright[domainnode->key*genesize+currentnode->key]=currentnode->key;
			cleft[domainnode->key*genesize+currentnode->key]=inindex[domainnode->left->key*genesize+currentnode->left->key];
		}
		low=GetMin(low,SafeAdd(c[domainnode->right->key*genesize+currentnode->key], in[domainnode->left->key*genesize+currentnode->right->key], 1));
		if(c[domainnode->key*genesize+currentnode->key]>SafeAdd(low, domainDuplicationcost))		
		{
			c[domainnode->key*genesize+currentnode->key]=SafeAdd(low, domainDuplicationcost);
			events[domainnode->key*genesize+currentnode->key]=2;
			cright[domainnode->key*genesize+currentnode->key]=currentnode->key;
			cleft[domainnode->key*genesize+currentnode->key]=inindex[domainnode->left->key*genesize+currentnode->right->key];
		}
		low=GetMin(low,SafeAdd(in[domainnode->left->key*genesize+currentnode->left->key], in[domainnode->right->key*genesize+currentnode->right->key], 2));
		if(c[domainnode->key*genesize+currentnode->key]>SafeAdd(low, domainDuplicationcost))		
		{
			c[domainnode->key*genesize+currentnode->key]=SafeAdd(low, domainDuplicationcost);
			events[domainnode->key*genesize+currentnode->key]=2;
			cleft[domainnode->key*genesize+currentnode->key]=inindex[domainnode->left->key*genesize+currentnode->left->key];
			cright[domainnode->key*genesize+currentnode->key]=inindex[domainnode->right->key*genesize+currentnode->right->key];
		}
		low=GetMin(low,SafeAdd(in[domainnode->left->key*genesize+currentnode->right->key], in[domainnode->right->key*genesize+currentnode->left->key], 2));
		if(c[domainnode->key*genesize+currentnode->key]>SafeAdd(low, domainDuplicationcost))		
		{
			c[domainnode->key*genesize+currentnode->key]=SafeAdd(low, domainDuplicationcost);
			events[domainnode->key*genesize+currentnode->key]=2;
			cleft[domainnode->key*genesize+currentnode->key]=inindex[domainnode->left->key*genesize+currentnode->right->key];
			cright[domainnode->key*genesize+currentnode->key]=inindex[domainnode->right->key*genesize+currentnode->left->key];
		}
		low=GetMin(low,SafeAdd(in[domainnode->left->key*genesize+currentnode->left->key], in[domainnode->right->key*genesize+currentnode->left->key], 2));
		if(c[domainnode->key*genesize+currentnode->key]>SafeAdd(low, domainDuplicationcost))		
		{
			c[domainnode->key*genesize+currentnode->key]=SafeAdd(low, domainDuplicationcost);
			events[domainnode->key*genesize+currentnode->key]=2;
			cleft[domainnode->key*genesize+currentnode->key]=inindex[domainnode->left->key*genesize+currentnode->left->key];
			cright[domainnode->key*genesize+currentnode->key]=inindex[domainnode->right->key*genesize+currentnode->left->key];
		}
		low=GetMin(low,SafeAdd(in[domainnode->left->key*genesize+currentnode->right->key], in[domainnode->right->key*genesize+currentnode->right->key], 2));
		if(c[domainnode->key*genesize+currentnode->key]>SafeAdd(low, domainDuplicationcost))		
		{
			c[domainnode->key*genesize+currentnode->key]=SafeAdd(low, domainDuplicationcost);
			events[domainnode->key*genesize+currentnode->key]=2;
			cleft[domainnode->key*genesize+currentnode->key]=inindex[domainnode->left->key*genesize+currentnode->right->key];
			cright[domainnode->key*genesize+currentnode->key]=inindex[domainnode->right->key*genesize+currentnode->right->key];
		}
		low=GetMin(low,SafeAdd(c[domainnode->left->key*genesize+currentnode->key], c[domainnode->right->key*genesize+currentnode->key]));
		if(c[domainnode->key*genesize+currentnode->key]>SafeAdd(low, domainDuplicationcost))		
		{
			c[domainnode->key*genesize+currentnode->key]=SafeAdd(low, domainDuplicationcost);
			events[domainnode->key*genesize+currentnode->key]=2;
			cleft[domainnode->key*genesize+currentnode->key]=currentnode->key;
			cright[domainnode->key*genesize+currentnode->key]=currentnode->key;
		}
		//cout<<c[domainnode->key*genesize+currentnode->key]<<endl;
	}
	else
	{
		low=SafeAdd(c[domainnode->left->key*genesize+currentnode->key], c[domainnode->right->key*genesize+currentnode->key]);
		if(c[domainnode->key*genesize+currentnode->key]>SafeAdd(low, domainDuplicationcost))
		{
			c[domainnode->key*genesize+currentnode->key]=SafeAdd(low, domainDuplicationcost);
			events[domainnode->key*genesize+currentnode->key]=2;
			cleft[domainnode->key*genesize+currentnode->key]=currentnode->key;
			cright[domainnode->key*genesize+currentnode->key]=currentnode->key;
		}
		//out[domainnode->key*genesize+currentnode->key]=c[domainnode->key*genesize+currentnode->key];
	}
	int stand=0;
	// Transfer Cases
	for (int i=0; i<genesize;i++) // See where is the real receiver.
	{
		stand=MAX;
		// Check the tree range
		node* high=currentnode;
		while(high->parent!=NULL)
		{
			high=high->parent;
		}
		j=lowestkey(high);
		k=high->key;

		// Check if this is a two tree or one tree transfer.
		if (i<j || i>k) transferType=1;
		else if ((i<=k && i>=j) && !(comparable(currentnode,i,genepointers) || comparable(genepointers[i],currentnode->key,genepointers)))
			transferType=0;
		else continue;

		leftlow =SafeAdd(in[domainnode->left->key*genesize+currentnode->key], in[domainnode->right->key*genesize+i]);
		rightlow=SafeAdd(in[domainnode->right->key*genesize+currentnode->key], in[domainnode->left->key*genesize+i]);
		low=SafeAdd(GetMin(leftlow,rightlow), transferType*(twoTreeTransferCost-OneTreeTransferCost)+OneTreeTransferCost);



		low=SafeAdd(low, ExtraLoss(currentnode,genepointers[i],speciespointers));
		                                                                                 
		

		lca=LCA(speciespointers[genepointers[i]->mappedNode],speciespointers[currentnode->mappedNode]);
		

		if(low<c[domainnode->key*genesize+currentnode->key])
		{
			re[domainnode->key*genesize+currentnode->key]=genepointers[i];
			c[domainnode->key*genesize+currentnode->key]=low;
			if(currentnode->itsroot==genepointers[i]->itsroot)
				events[domainnode->key*genesize+currentnode->key]=3;
			else events[domainnode->key*genesize+currentnode->key]=4;
			if(GetMin(leftlow,rightlow)==leftlow)
			{
				cleft[domainnode->key*genesize+currentnode->key]=inindex[domainnode->left->key*genesize+currentnode->key];
				cright[domainnode->key*genesize+currentnode->key]=inindex[domainnode->right->key*genesize+i];
				//cout<<cright[domainnode->key*genesize+currentnode->key]<<endl;
			}
			else
			{
				cleft[domainnode->key*genesize+currentnode->key]=inindex[domainnode->left->key*genesize+i];
				cright[domainnode->key*genesize+currentnode->key]=inindex[domainnode->right->key*genesize+currentnode->key];
				//cout<<cright[domainnode->key*genesize+currentnode->key]<<endl;
			}
		}
	}
	

		if(currentnode->isleaf==1)
		{
			in[domainnode->key*genesize+currentnode->key]=c[domainnode->key*genesize+currentnode->key];
			inindex[domainnode->key*genesize+currentnode->key]=currentnode->key;
		}
		else
		{
			in[domainnode->key*genesize+currentnode->key]=GetMin(c[domainnode->key*genesize+currentnode->key],SafeAdd(in[domainnode->key*genesize+currentnode->left->key], 1));
			if(GetMin(c[domainnode->key*genesize+currentnode->key],SafeAdd(in[domainnode->key*genesize+currentnode->left->key], 1))==c[domainnode->key*genesize+currentnode->key])
				inindex[domainnode->key*genesize+currentnode->key]=currentnode->key;
			else
				inindex[domainnode->key*genesize+currentnode->key]=inindex[domainnode->key*genesize+currentnode->left->key];

			if(GetMin(in[domainnode->key*genesize+currentnode->key],SafeAdd(in[domainnode->key*genesize+currentnode->right->key], 1))==SafeAdd(in[domainnode->key*genesize+currentnode->right->key], 1))
				inindex[domainnode->key*genesize+currentnode->key]=inindex[domainnode->key*genesize+currentnode->right->key];
			in[domainnode->key*genesize+currentnode->key]=GetMin(in[domainnode->key*genesize+currentnode->key],SafeAdd(in[domainnode->key*genesize+currentnode->right->key], 1));
		}
	//if(c[domainnode->key*genesize+currentnode->key]<5000)
	//cout<<domainnode->key<<"\t"<<currentnode->key<<"\t"<<c[domainnode->key*genesize+currentnode->key]<<"\t"<<in[domainnode->key*genesize+currentnode->key]<<"\t"<<inindex[domainnode->key*genesize+currentnode->key]<<"\t"<<cleft[domainnode->key*genesize+currentnode->key]<<"\t"<<cright[domainnode->key*genesize+currentnode->key]<<endl;

	if (debug2850) {
		debugCount2850++;
	}
	return 0;
}

int MaxPostorderDomain(int *inindex, node **r, int *events, int* c, int* in, int* cleft, int* cright, tree *domainTree, tree **geneTrees, tree* speciesTree, int genesize, node* currentnode, node** domainpointers, node **genepointers, node** speciespointers, int twoTreeTransferCost, int OneTreeTransferCost, int domainDuplicationcost, int geneDuplicationcost, int domainLosscost, int geneLosscost)
{
	int genetreeindex=0;
	node *genetreeroot;

	// Debug output for root node
	static bool isFirstCall = true;
	if (isFirstCall && currentnode->key == domainTree->root->key) {
		if (currentnode->left != NULL) {
		} else {
		}
		if (currentnode->right != NULL) {
		} else {
		}
		isFirstCall = false;
	}

	// Special debugging for node 2852 (the problematic left child)
	if (currentnode->key == 2852) {
		if (currentnode->isleaf == 0) {
			if (currentnode->left != NULL) {
			} else {
			}
			if (currentnode->right != NULL) {
			} else {
			}
		}
	}

	// Special debugging for node 2851 (problematic right child of 2852)
	if (currentnode->key == 2851) {
		if (currentnode->isleaf == 0) {
			if (currentnode->left != NULL) {
			} else {
			}
			if (currentnode->right != NULL) {
			} else {
			}
		}
	}

	// Special debugging for node 2850 (problematic right child of 2851)
	if (currentnode->key == 2850) {
		if (currentnode->isleaf == 0) {
			if (currentnode->left != NULL) {
			} else {
			}
			if (currentnode->right != NULL) {
			} else {
			}
		}
	}

	if(currentnode->isleaf==0)
	{
		if (currentnode->key == domainTree->root->key) {
		}
		MaxPostorderDomain(inindex,r,events,c,in,cleft,cright,domainTree,geneTrees,speciesTree,genesize,currentnode->left,domainpointers,genepointers,speciespointers,twoTreeTransferCost,OneTreeTransferCost,domainDuplicationcost,geneDuplicationcost,domainLosscost,geneLosscost);
		MaxPostorderDomain(inindex,r,events,c,in,cleft,cright,domainTree,geneTrees,speciesTree,genesize,currentnode->right,domainpointers,genepointers,speciespointers,twoTreeTransferCost,OneTreeTransferCost,domainDuplicationcost,geneDuplicationcost,domainLosscost,geneLosscost);

		// Check costs after processing children (for node 2852 and 2851)
		if (currentnode->key == 2852) {
			int leftValidCosts = 0, rightValidCosts = 0;
			for (int i = 0; i < genesize; i++) {
				if (c[currentnode->left->key*genesize+i] < 5000) leftValidCosts++;
				if (c[currentnode->right->key*genesize+i] < 5000) rightValidCosts++;
			}
		}
		if (currentnode->key == 2851) {
			int leftValidCosts = 0, rightValidCosts = 0;
			for (int i = 0; i < genesize; i++) {
				if (c[currentnode->left->key*genesize+i] < 5000) leftValidCosts++;
				if (c[currentnode->right->key*genesize+i] < 5000) rightValidCosts++;
			}
		}
		if (currentnode->key == 2850) {
			int leftValidCosts = 0, rightValidCosts = 0;
			int bothValidCosts = 0;
			for (int i = 0; i < genesize; i++) {
				bool leftValid = (c[currentnode->left->key*genesize+i] < 5000);
				bool rightValid = (c[currentnode->right->key*genesize+i] < 5000);
				if (leftValid) leftValidCosts++;
				if (rightValid) rightValidCosts++;
				if (leftValid && rightValid) bothValidCosts++;
			}
		}

		// need to check if there exitsts a solution
		//cout<<"Doing domain node "<<currentnode->key<<endl;
		if (currentnode->key == domainTree->root->key) {
		}
		if (currentnode->key == 2852) {
		}
		if (currentnode->key == 2851) {
		}
		if (currentnode->key == 2850) {
		}
		for (genetreeindex=0;genetreeindex<domainTree->mappedTrees.size();genetreeindex++)
		{
			genetreeroot=geneTrees[genetreeindex]->root;
			MaxPostorderGene(inindex,r,genetreeindex,events,c,in,cleft,cright,domainTree,geneTrees,speciesTree,genesize,currentnode,genetreeroot,domainpointers,genepointers,speciespointers,twoTreeTransferCost,OneTreeTransferCost,domainDuplicationcost,geneDuplicationcost,domainLosscost,geneLosscost);
		}
		if (currentnode->key == domainTree->root->key) {
		}
		if (currentnode->key == 2852) {
			// Check how many valid costs were computed for node 2852
			int validCosts = 0;
			int minCost = 5000;
			for (int i = 0; i < genesize; i++) {
				if (c[currentnode->key*genesize+i] < 5000) {
					validCosts++;
					if (c[currentnode->key*genesize+i] < minCost) {
						minCost = c[currentnode->key*genesize+i];
					}
				}
			}
		}
		if (currentnode->key == 2851) {
			// Check how many valid costs were computed for node 2851
			int validCosts = 0;
			int minCost = 5000;
			for (int i = 0; i < genesize; i++) {
				if (c[currentnode->key*genesize+i] < 5000) {
					validCosts++;
					if (c[currentnode->key*genesize+i] < minCost) {
						minCost = c[currentnode->key*genesize+i];
					}
				}
			}
		}
		if (currentnode->key == 2850) {
			// Check how many valid costs were computed for node 2850
			int validCosts = 0;
			int minCost = 5000;
			for (int i = 0; i < genesize; i++) {
				if (c[currentnode->key*genesize+i] < 5000) {
					validCosts++;
					if (c[currentnode->key*genesize+i] < minCost) {
						minCost = c[currentnode->key*genesize+i];
					}
				}
			}
		}

	}
	else if (currentnode->key == domainTree->root->key) {
	}

	//else cout<<"skipping leaf "<<currentnode->key<<endl;
	return 0;
}

#endif