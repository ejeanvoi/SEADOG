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
		leftlow=(in[domainnode->left->key*genesize+currentnode->left->key]+ in[domainnode->right->key*genesize+currentnode->right->key]);
		rightlow=(in[domainnode->right->key*genesize+currentnode->left->key]+ in[domainnode->left->key*genesize+currentnode->right->key]);

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
		low=c[domainnode->left->key*genesize+currentnode->key]+in[domainnode->right->key*genesize+currentnode->left->key]+1;
		if(c[domainnode->key*genesize+currentnode->key]>low+domainDuplicationcost)
		{
			c[domainnode->key*genesize+currentnode->key]=low+domainDuplicationcost;
			events[domainnode->key*genesize+currentnode->key]=2;
			cleft[domainnode->key*genesize+currentnode->key]=currentnode->key;
			cright[domainnode->key*genesize+currentnode->key]=inindex[domainnode->right->key*genesize+currentnode->left->key];
		}
		low=GetMin(low,c[domainnode->left->key*genesize+currentnode->key]+in[domainnode->right->key*genesize+currentnode->right->key]+1);
		if(c[domainnode->key*genesize+currentnode->key]>low+domainDuplicationcost)		
		{
			c[domainnode->key*genesize+currentnode->key]=low+domainDuplicationcost;
			events[domainnode->key*genesize+currentnode->key]=2;
			cleft[domainnode->key*genesize+currentnode->key]=currentnode->key;
			cright[domainnode->key*genesize+currentnode->key]=inindex[domainnode->right->key*genesize+currentnode->right->key];
		}
		low=GetMin(low,c[domainnode->right->key*genesize+currentnode->key]+in[domainnode->left->key*genesize+currentnode->left->key]+1);
		if(c[domainnode->key*genesize+currentnode->key]>low+domainDuplicationcost)		
		{
			c[domainnode->key*genesize+currentnode->key]=low+domainDuplicationcost;
			events[domainnode->key*genesize+currentnode->key]=2;
			cright[domainnode->key*genesize+currentnode->key]=currentnode->key;
			cleft[domainnode->key*genesize+currentnode->key]=inindex[domainnode->left->key*genesize+currentnode->left->key];
		}
		low=GetMin(low,c[domainnode->right->key*genesize+currentnode->key]+in[domainnode->left->key*genesize+currentnode->right->key]+1);
		if(c[domainnode->key*genesize+currentnode->key]>low+domainDuplicationcost)		
		{
			c[domainnode->key*genesize+currentnode->key]=low+domainDuplicationcost;
			events[domainnode->key*genesize+currentnode->key]=2;
			cright[domainnode->key*genesize+currentnode->key]=currentnode->key;
			cleft[domainnode->key*genesize+currentnode->key]=inindex[domainnode->left->key*genesize+currentnode->right->key];
		}
		low=GetMin(low,in[domainnode->left->key*genesize+currentnode->left->key]+in[domainnode->right->key*genesize+currentnode->right->key]+2);
		if(c[domainnode->key*genesize+currentnode->key]>low+domainDuplicationcost)		
		{
			c[domainnode->key*genesize+currentnode->key]=low+domainDuplicationcost;
			events[domainnode->key*genesize+currentnode->key]=2;
			cleft[domainnode->key*genesize+currentnode->key]=inindex[domainnode->left->key*genesize+currentnode->left->key];
			cright[domainnode->key*genesize+currentnode->key]=inindex[domainnode->right->key*genesize+currentnode->right->key];
		}
		low=GetMin(low,in[domainnode->left->key*genesize+currentnode->right->key]+in[domainnode->right->key*genesize+currentnode->left->key]+2);
		if(c[domainnode->key*genesize+currentnode->key]>low+domainDuplicationcost)		
		{
			c[domainnode->key*genesize+currentnode->key]=low+domainDuplicationcost;
			events[domainnode->key*genesize+currentnode->key]=2;
			cleft[domainnode->key*genesize+currentnode->key]=inindex[domainnode->left->key*genesize+currentnode->right->key];
			cright[domainnode->key*genesize+currentnode->key]=inindex[domainnode->right->key*genesize+currentnode->left->key];
		}
		low=GetMin(low,in[domainnode->left->key*genesize+currentnode->left->key]+in[domainnode->right->key*genesize+currentnode->left->key]+2);
		if(c[domainnode->key*genesize+currentnode->key]>low+domainDuplicationcost)		
		{
			c[domainnode->key*genesize+currentnode->key]=low+domainDuplicationcost;
			events[domainnode->key*genesize+currentnode->key]=2;
			cleft[domainnode->key*genesize+currentnode->key]=inindex[domainnode->left->key*genesize+currentnode->left->key];
			cright[domainnode->key*genesize+currentnode->key]=inindex[domainnode->right->key*genesize+currentnode->left->key];
		}
		low=GetMin(low,in[domainnode->left->key*genesize+currentnode->right->key]+in[domainnode->right->key*genesize+currentnode->right->key]+2);
		if(c[domainnode->key*genesize+currentnode->key]>low+domainDuplicationcost)		
		{
			c[domainnode->key*genesize+currentnode->key]=low+domainDuplicationcost;
			events[domainnode->key*genesize+currentnode->key]=2;
			cleft[domainnode->key*genesize+currentnode->key]=inindex[domainnode->left->key*genesize+currentnode->right->key];
			cright[domainnode->key*genesize+currentnode->key]=inindex[domainnode->right->key*genesize+currentnode->right->key];
		}
		low=GetMin(low,c[domainnode->left->key*genesize+currentnode->key]+c[domainnode->right->key*genesize+currentnode->key]);
		if(c[domainnode->key*genesize+currentnode->key]>low+domainDuplicationcost)		
		{
			c[domainnode->key*genesize+currentnode->key]=low+domainDuplicationcost;
			events[domainnode->key*genesize+currentnode->key]=2;
			cleft[domainnode->key*genesize+currentnode->key]=currentnode->key;
			cright[domainnode->key*genesize+currentnode->key]=currentnode->key;
		}
		//cout<<c[domainnode->key*genesize+currentnode->key]<<endl;
	}
	else
	{
		low=c[domainnode->left->key*genesize+currentnode->key]+c[domainnode->right->key*genesize+currentnode->key];
		if(c[domainnode->key*genesize+currentnode->key]>low+domainDuplicationcost)
		{
			c[domainnode->key*genesize+currentnode->key]=low+domainDuplicationcost;
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

		leftlow =in[domainnode->left->key*genesize+currentnode->key]+in[domainnode->right->key*genesize+i];
		rightlow=in[domainnode->right->key*genesize+currentnode->key]+in[domainnode->left->key*genesize+i];
		low=GetMin(leftlow,rightlow)+transferType*(twoTreeTransferCost-OneTreeTransferCost)+OneTreeTransferCost;

		
		
		low+=ExtraLoss(currentnode,genepointers[i],speciespointers);
		                                                                                 
		

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
			in[domainnode->key*genesize+currentnode->key]=GetMin(c[domainnode->key*genesize+currentnode->key],in[domainnode->key*genesize+currentnode->left->key]+1);
			if(GetMin(c[domainnode->key*genesize+currentnode->key],in[domainnode->key*genesize+currentnode->left->key]+1)==c[domainnode->key*genesize+currentnode->key])
				inindex[domainnode->key*genesize+currentnode->key]=currentnode->key;
			else
				inindex[domainnode->key*genesize+currentnode->key]=inindex[domainnode->key*genesize+currentnode->left->key];

			if(GetMin(in[domainnode->key*genesize+currentnode->key],in[domainnode->key*genesize+currentnode->right->key]+1)==(in[domainnode->key*genesize+currentnode->right->key]+1))
				inindex[domainnode->key*genesize+currentnode->key]=inindex[domainnode->key*genesize+currentnode->right->key];
			in[domainnode->key*genesize+currentnode->key]=GetMin(in[domainnode->key*genesize+currentnode->key],in[domainnode->key*genesize+currentnode->right->key]+1);
		}
	//if(c[domainnode->key*genesize+currentnode->key]<5000)
	//cout<<domainnode->key<<"\t"<<currentnode->key<<"\t"<<c[domainnode->key*genesize+currentnode->key]<<"\t"<<in[domainnode->key*genesize+currentnode->key]<<"\t"<<inindex[domainnode->key*genesize+currentnode->key]<<"\t"<<cleft[domainnode->key*genesize+currentnode->key]<<"\t"<<cright[domainnode->key*genesize+currentnode->key]<<endl;
	return 0;
}

int MaxPostorderDomain(int *inindex, node **r, int *events, int* c, int* in, int* cleft, int* cright, tree *domainTree, tree **geneTrees, tree* speciesTree, int genesize, node* currentnode, node** domainpointers, node **genepointers, node** speciespointers, int twoTreeTransferCost, int OneTreeTransferCost, int domainDuplicationcost, int geneDuplicationcost, int domainLosscost, int geneLosscost)
{
	int genetreeindex=0;
	node *genetreeroot;

	// Debug output for root node
	static bool isFirstCall = true;
	if (isFirstCall && currentnode->key == domainTree->root->key) {
		cout<<"[MPD-1] Processing domain tree root: key="<<currentnode->key<<", name="<<currentnode->name<<", isleaf="<<currentnode->isleaf<<endl;
		if (currentnode->left != NULL) {
			cout<<"[MPD-2] Root left child: key="<<currentnode->left->key<<", name="<<currentnode->left->name<<", isleaf="<<currentnode->left->isleaf<<endl;
		} else {
			cout<<"[MPD-2] Root left child is NULL!"<<endl;
		}
		if (currentnode->right != NULL) {
			cout<<"[MPD-3] Root right child: key="<<currentnode->right->key<<", name="<<currentnode->right->name<<", isleaf="<<currentnode->right->isleaf<<endl;
		} else {
			cout<<"[MPD-3] Root right child is NULL!"<<endl;
		}
		isFirstCall = false;
	}

	if(currentnode->isleaf==0)
	{
		if (currentnode->key == domainTree->root->key) {
			cout<<"[MPD-4] Processing root's children..."<<endl;
		}
		MaxPostorderDomain(inindex,r,events,c,in,cleft,cright,domainTree,geneTrees,speciesTree,genesize,currentnode->left,domainpointers,genepointers,speciespointers,twoTreeTransferCost,OneTreeTransferCost,domainDuplicationcost,geneDuplicationcost,domainLosscost,geneLosscost);
		MaxPostorderDomain(inindex,r,events,c,in,cleft,cright,domainTree,geneTrees,speciesTree,genesize,currentnode->right,domainpointers,genepointers,speciespointers,twoTreeTransferCost,OneTreeTransferCost,domainDuplicationcost,geneDuplicationcost,domainLosscost,geneLosscost);
		// need to check if there exitsts a solution
		//cout<<"Doing domain node "<<currentnode->key<<endl;
		if (currentnode->key == domainTree->root->key) {
			cout<<"[MPD-5] Processing root with "<<domainTree->mappedTrees.size()<<" gene trees..."<<endl;
		}
		for (genetreeindex=0;genetreeindex<domainTree->mappedTrees.size();genetreeindex++)
		{
			genetreeroot=geneTrees[genetreeindex]->root;
			MaxPostorderGene(inindex,r,genetreeindex,events,c,in,cleft,cright,domainTree,geneTrees,speciesTree,genesize,currentnode,genetreeroot,domainpointers,genepointers,speciespointers,twoTreeTransferCost,OneTreeTransferCost,domainDuplicationcost,geneDuplicationcost,domainLosscost,geneLosscost);
		}
		if (currentnode->key == domainTree->root->key) {
			cout<<"[MPD-6] Finished processing root gene trees"<<endl;
		}

	}
	else if (currentnode->key == domainTree->root->key) {
		cout<<"[MPD-ERROR] Domain tree root is marked as a leaf (isleaf="<<currentnode->isleaf<<")! This should never happen."<<endl;
	}

	//else cout<<"skipping leaf "<<currentnode->key<<endl;
	return 0;
}

#endif