#ifndef GeneUpper_H
#define GeneUpper_H

#include <string.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include "tree.h"
#include <fstream>
#include <cmath>
#include "libs.h"
#include "DoMaxNum.h"

#ifndef MAX
#define MAX 5000
#endif

int UpLocaltraceback(int *dup, int *trans, int *trans2, int flag, int *low, int* geneflag, int* genestep, node **r, node *domainnode, int mappedindex, int *events, int *cleft, int *cright, int genesize, node **speciespointers, node **genepointers, string domainFileName, int twoTreeTransferCost, int OneTreeTransferCost, int domainDuplicationcost, int geneDuplicationcost, int domainLosscost, int geneLosscost)
{
	// Validate parameters
	if (domainnode == NULL) {
		cout<<"ERROR: UpLocaltraceback called with NULL domainnode!"<<endl;
		return -1;
	}

	// Check array index that will be used
	int arrayIndex = domainnode->key * genesize + mappedindex;
	cout<<"[ULT] Processing domain node '"<<domainnode->name<<"' (key="<<domainnode->key<<", isleaf="<<domainnode->isleaf<<") mapped to gene index "<<mappedindex<<", arrayIndex="<<arrayIndex<<endl;

	if (mappedindex < 0 || mappedindex >= genesize) {
		cout<<"ERROR: UpLocaltraceback called with invalid mappedindex="<<mappedindex<<" (genesize="<<genesize<<")"<<endl;
		cout<<"       Domain node: "<<domainnode->name<<" (key="<<domainnode->key<<")"<<endl;
		return -1;
	}
	if (genepointers[mappedindex] == NULL) {
		cout<<"ERROR: genepointers["<<mappedindex<<"] is NULL!"<<endl;
		cout<<"       Domain node: "<<domainnode->name<<endl;
		return -1;
	}

	if(flag)
	{
		fout.open(domainFileName.c_str(),ios::app);
		fout<<domainnode->name<<": ";

		if (events[domainnode->key*genesize+mappedindex]==0)
			fout<<"Leaf";
		if (events[domainnode->key*genesize+mappedindex]==1)
			fout<<"Co-divergence";
		if (events[domainnode->key*genesize+mappedindex]==2)
		{
			fout<<"Domain duplication";
			*dup=*dup+1;
		}
		if (events[domainnode->key*genesize+mappedindex]==3)
		{
			fout<<"Intra-gene-tree domain transfer";
			*trans=*trans+1;
		}

		if (events[domainnode->key*genesize+mappedindex]==4)
		{
			fout<<"Inter-gene-tree domain transfer";
			*trans2=*trans2+1;
		}

		fout<<", Mapping -> "<<genepointers[mappedindex]->name;

		if(events[domainnode->key*genesize+mappedindex]==3 || events[domainnode->key*genesize+mappedindex]==4)
		{
			// Validate r pointer before accessing
			if (r[domainnode->key*genesize+mappedindex] == NULL) {
				cout<<"ERROR: r["<<domainnode->key*genesize+mappedindex<<"] is NULL for transfer event!"<<endl;
				cout<<"       Domain node: "<<domainnode->name<<", mappedindex="<<mappedindex<<endl;
				fout.close();
				return -1;
			}
			fout<<", Recipient -> "<<r[domainnode->key*genesize+mappedindex]->name;
		}



		fout<<endl;
		fout.close();
	}

	//cout<<*low<<endl;
	//cout<<"Domain node "<<domainnode->name<<" mapped to "<<mappedindex<<" with event type "<<events[domainnode->key*genesize+mappedindex]<<endl;
	domainnode->mappedNode=mappedindex;
	domainnode->eventtype=events[domainnode->key*genesize+mappedindex];
	domainnode->mappedtoSpecies=genepointers[mappedindex]->mappedNode;
	//cout<<cleft[domainnode->key*genesize+mappedindex]<<" "<<cright[domainnode->key*genesize+mappedindex]<<endl;
	if (domainnode->isleaf==0)
	{
		// Validate child mapping indices before recursive calls
		int leftIndex = cleft[domainnode->key*genesize+mappedindex];
		int rightIndex = cright[domainnode->key*genesize+mappedindex];

		if (leftIndex < 0 || leftIndex >= genesize) {
			cout<<"ERROR: Invalid left child mapping index="<<leftIndex<<" (genesize="<<genesize<<")"<<endl;
			cout<<"       Domain node: "<<domainnode->name<<" (key="<<domainnode->key<<"), mappedindex="<<mappedindex<<endl;
			return -1;
		}
		if (rightIndex < 0 || rightIndex >= genesize) {
			cout<<"ERROR: Invalid right child mapping index="<<rightIndex<<" (genesize="<<genesize<<")"<<endl;
			cout<<"       Domain node: "<<domainnode->name<<" (key="<<domainnode->key<<"), mappedindex="<<mappedindex<<endl;
			return -1;
		}

		UpLocaltraceback(dup,trans,trans2,flag,low,geneflag,genestep,r,domainnode->left,leftIndex,events,cleft,cright,genesize,speciespointers,genepointers,domainFileName,twoTreeTransferCost,OneTreeTransferCost,domainDuplicationcost,geneDuplicationcost,domainLosscost,geneLosscost);
		UpLocaltraceback(dup,trans,trans2,flag,low,geneflag,genestep,r,domainnode->right,rightIndex,events,cleft,cright,genesize,speciespointers,genepointers,domainFileName,twoTreeTransferCost,OneTreeTransferCost,domainDuplicationcost,geneDuplicationcost,domainLosscost,geneLosscost);

	}
	if (domainnode->isleaf==1)
		return 0;
	if(events[domainnode->key*genesize+mappedindex]==1)
	{
		*low=*low+genepointers[cleft[domainnode->key*genesize+mappedindex]]->depth-genepointers[mappedindex]->depth-1;
		*low=*low+genepointers[cright[domainnode->key*genesize+mappedindex]]->depth-genepointers[mappedindex]->depth-1;

	}

	if(events[domainnode->key*genesize+mappedindex]==2)
	{
		*low=*low+genepointers[cleft[domainnode->key*genesize+mappedindex]]->depth-genepointers[mappedindex]->depth;
		*low=*low+genepointers[cright[domainnode->key*genesize+mappedindex]]->depth-genepointers[mappedindex]->depth;
		*low=*low+domainDuplicationcost;

	}

	if(events[domainnode->key*genesize+mappedindex]==3|| events[domainnode->key*genesize+mappedindex]==4)
	{
		//cout<<cleft[domainnode->key*genesize+mappedindex]<<"\t"<<cright[domainnode->key*genesize+mappedindex]<<endl;
		*low=*low+genepointers[cleft[domainnode->key*genesize+mappedindex]]->depth-genepointers[mappedindex]->depth;
		*low=*low+genepointers[cright[domainnode->key*genesize+mappedindex]]->depth-r[domainnode->key*genesize+mappedindex]->depth;
		if(r[domainnode->key*genesize+mappedindex]->itsroot==genepointers[mappedindex]->itsroot)
		{*low=*low+OneTreeTransferCost;}
		else {*low=*low+twoTreeTransferCost;}
	}

	return 0;
}

int UptracebackReal(int flag, tree *speciesTree, int *added, node **r, tree *domaintree, tree** geneTrees, int *c, int *events, int *cleft, int *cright,  int genesize, node**domainpointers, node **speciespointers, node **genepointers, string domainFileName, int twoTreeTransferCost, int OneTreeTransferCost, int domainDuplicationcost, int geneDuplicationcost, int domainLosscost, int geneLosscost)
{
	cout<<"[UTR-1] Entered UptracebackReal, flag="<<flag<<endl;
	if(flag){
		cout<<"[UTR-2] Writing trees to output file..."<<endl;
		fout.open(domainFileName.c_str());
		fout<<"Domain Tree: "<<endl;
		fout.close();

		WritIntoFile(domaintree->root,domainFileName);
		fout.open(domainFileName.c_str(),ios::app);
		fout.close();
		cout<<"[UTR-3] Writing "<<domaintree->mappedTrees.size()<<" gene trees..."<<endl;
		for (int genetreeindex=0;genetreeindex<domaintree->mappedTrees.size();genetreeindex++)
		{
			fout.open(domainFileName.c_str(),ios::app);
			fout<<"Gene Tree "<<genetreeindex+1<<": "<<endl;
			fout<<geneTrees[genetreeindex]->root->treename<<endl;

			fout.close();
			WritIntoFile(geneTrees[genetreeindex]->root,domainFileName);
		}

		fout<<endl;
		fout.open(domainFileName.c_str(),ios::app);
		fout<<"Species Tree: "<<endl;
		fout.close();
		WritIntoFile(speciesTree->root,domainFileName);

		fout.open(domainFileName.c_str(),ios::app);
		fout<<endl<<"Reconciliation between domain tree and gene trees: "<<endl;
		fout.close();
		cout<<"[UTR-4] Tree writing complete"<<endl;
	}
	cout<<"[UTR-5] Finding root mapping..."<<endl;
	int low=MAX;
	int rootmapindex = -1;
	int validCount = 0;
	const int COST_THRESHOLD = MAX * 10;
	int minCostWithValidChildren = COST_THRESHOLD;  // Initialize to threshold, not MAX
	int bestValidIndex = -1;

	// First pass: find minimum cost and count valid mappings
	// Accept costs up to 10*MAX as potentially valid (allowing for complex reconciliations)
	for (int i=0;i<genesize;i++){
		if (c[domaintree->root->key*(genesize)+i] < COST_THRESHOLD) {
			int arrayIdx = domaintree->root->key*genesize+i;
			bool hasValidChildren = (cleft[arrayIdx] != -1 && cright[arrayIdx] != -1);
			if (hasValidChildren) {
				validCount++;
				if (c[domaintree->root->key*(genesize)+i] < minCostWithValidChildren) {
					minCostWithValidChildren = c[domaintree->root->key*(genesize)+i];
					bestValidIndex = i;
				}
			}
		}
		if (low>c[domaintree->root->key*(genesize)+i])
		{
			low = c[domaintree->root->key*(genesize)+i];
			rootmapindex=i;
		}
	}

	cout<<"[UTR-6] Found "<<validCount<<" gene nodes with valid child mappings for domain root"<<endl;
	cout<<"[UTR-6a] Original minimum cost: "<<low<<" at gene index "<<rootmapindex<<endl;
	if (bestValidIndex != -1) {
		cout<<"[UTR-6b] Best VALID mapping: cost="<<minCostWithValidChildren<<" at gene index "<<bestValidIndex<<" (name="<<genepointers[bestValidIndex]->name<<")"<<endl;
		// Use the best valid mapping instead
		rootmapindex = bestValidIndex;
		low = minCostWithValidChildren;
		cout<<"[UTR-6c] Switching to use best valid mapping"<<endl;
	}
	cout<<"[UTR-6d] Final choice: Root maps to index "<<rootmapindex<<" with cost "<<low<<endl;

	// Validate that we found at least one valid mapping
	if (bestValidIndex == -1 || rootmapindex == -1) {
		cout<<"ERROR: No valid reconciliation found for domain tree root!"<<endl;
		cout<<"       All gene nodes either have cost=MAX or invalid child indices"<<endl;
		cout<<"       This indicates a fundamental problem with the input data or algorithm"<<endl;
		return -1;
	}

	// Final validation
	int rootArrayIndex = domaintree->root->key*genesize+rootmapindex;
	cout<<"[UTR-6e] Final validation: cleft["<<rootArrayIndex<<"]="<<cleft[rootArrayIndex]<<", cright["<<rootArrayIndex<<"]="<<cright[rootArrayIndex]<<endl;

	int *geneflag =	new int[genesize];
	int *genestep =	new int[genesize];
	for (int i=0;i<genesize;i++)
	{
		geneflag[i]=0;
		genestep[i]=0;
	}
	cout<<"[UTR-7] Calling UpLocaltraceback..."<<endl;

	low=0;
	int dup=0, trans=0, trans2=0;

	UpLocaltraceback(&dup, &trans, &trans2, flag,&low,geneflag,genestep,r,domaintree->root,rootmapindex,events,cleft,cright,genesize,speciespointers,genepointers,domainFileName,twoTreeTransferCost,OneTreeTransferCost,domainDuplicationcost,geneDuplicationcost,domainLosscost,geneLosscost);
	cout<<"[UTR-8] UpLocaltraceback complete"<<endl;


	// Update Gene to Species
	cout<<"[UTR-9] Processing domain nodes for gene-to-species mapping..."<<endl;
	node *domainnode;
	int allsetflag=1;
	node *Iter;

	//while (allsetflag)
	{
		allsetflag=0;
		cout<<"[UTR-10] Checking "<<domaintree->leafnodes.size()*2-1<<" domain nodes..."<<endl;
		for(int i=0;i<domaintree->leafnodes.size()*2-1;i++)
		{
			domainnode=domainpointers[i];
			if(events[domainnode->key*genesize+domainnode->mappedNode]==3 || events[domainnode->key*genesize+domainnode->mappedNode]==4)
			{
				if(CheckGeneCost(geneflag,genestep, r[domainnode->key*genesize+domainnode->mappedNode], genepointers[domainnode->mappedNode],speciespointers)==1)
					allsetflag=1;
			}
		}
		cout<<"[UTR-11] Processing gene nodes for event types..."<<endl;

		for (int i=0;i<genesize;i++)
		{
			if(geneflag[i])
				genepointers[i]->eventtype=2;
			Iter=speciespointers[genepointers[i]->mappedNode];
			for (int j=0;j<genestep[i];j++)
			{
				if(genepointers[i]->isroot && i%2==0) continue;
				if (Iter == NULL || Iter->parent == NULL) break;
				Iter=Iter->parent;
			}
			genepointers[i]->mappedNode=Iter->key;
		}
	}



	//cout<<"Cost on the domain -> gene tree: "<<low<<endl;

	int domainloss=(low-dup*domainDuplicationcost-trans*OneTreeTransferCost-trans2*twoTreeTransferCost);


	for (int i=0;i<genesize;i++)
	{
		//cout<<*added<<endl;
		*added=*added+geneflag[i]*(geneDuplicationcost+2);
		*added=*added+genestep[i];
	}

	int Lcascore=0;
	for (int genetreeindex=0;genetreeindex<domaintree->mappedTrees.size();genetreeindex++)
	{
		Lcascore+=geneTrees[genetreeindex]->root->LCAscore;
	}

	//for (int i=0;i<genesize;i++)
	//{
		//cout<<i<<endl;
		//for (int j=0;j<genestep[i];j++)
		//{
			//genepointers[i]->mappedNode=speciespointers[genepointers[i]->mappedNode]->parent->key;
			//genepointers[i]->eventtype=2;
		//}
	//}

	int genedup=0;
	int geneloss=0;

	for (int i=0;i<genesize;i++)
	{
		if(genepointers[i]->eventtype==2 || geneflag[i])
			genedup++;
	}

	geneloss=Lcascore+*added-genedup*geneDuplicationcost;
	int last=0;

	if(flag){
		fout.open(domainFileName.c_str(),ios::app);
		fout<<endl<<"Reconciliation between Gene trees and Species tree: "<<endl;
		for (int genetreeindex=0;genetreeindex<domaintree->mappedTrees.size();genetreeindex++)
		{
			fout<<"Gene Tree "<<genetreeindex+1<<endl;

			for (int i=last;i<last+geneTrees[genetreeindex]->leafnodes.size()*2-1;i++)
			{
				fout<<genepointers[i]->name<<": ";
				if (genepointers[i]->eventtype==0)
					fout<<" Leaf";
				if (genepointers[i]->eventtype==1)
					fout<<" Speciation";
				if (genepointers[i]->eventtype==2)
					fout<<" Gene duplication";
				fout<<", Mapping -> "<<speciespointers[genepointers[i]->mappedNode]->name;
				fout<<endl;
			}
			fout<<endl;
			last=last+geneTrees[genetreeindex]->leafnodes.size()*2-1;
		}

		fout<<endl;
		fout<<"Minimal DGS reconciliation cost: "<<low+*added+Lcascore<<endl;

		fout<<"Domain-Gene reconciliation cost: "<<low<<" [domain duplications: "<<dup<<", Intra-gene-tree domain transfers: "<<trans<<", Inter-gene-tree domain transfers: "<<trans2<<", domain losses: "<<domainloss<<"]"<<endl;

		fout<<"Gene-Species reconciliation cost: "<<Lcascore+*added<<" [gene duplications: "<<genedup<<", gene loss: "<<geneloss<<"]"<<endl;

		fout.close();
	}
	low=low+*added;


	delete geneflag;
	delete genestep;
	return low;
}



int DynamicHeuristic(int flag, tree *domainTree, tree **geneTrees, tree* speciesTree, int genesize, node **domainpointers, node **genepointers, node** speciespointers, int twoTreeTransferCost, int OneTreeTransferCost, int domainDuplicationcost, int geneDuplicationcost, int domainLosscost, int geneLosscost, string domainFileName)
{
	cout<<"[DH-1] Entered DynamicHeuristic"<<endl;
	string filename = domainFileName;
	int domainsize=domainTree->leafnodes.size()*2-1;
	int speciessize=speciesTree->leafnodes.size()*2-1;
	cout<<"[DH-2] domainsize="<<domainsize<<", genesize="<<genesize<<", speciessize="<<speciessize<<endl;

	int matrixsize=domainsize*genesize;
	cout<<"[DH-3] matrixsize="<<matrixsize<<", allocating arrays..."<<endl;

	int *c		= new int[matrixsize];
	int *cleft		= new int[matrixsize];
	int *cright		= new int[matrixsize];
	int *events		= new int[matrixsize];
	int *in		= new int[matrixsize];
	node **r		= new node*[matrixsize];
	int *inindex	= new int[matrixsize];
	cout<<"[DH-4] Arrays allocated, initializing..."<<endl;
	//memset
	for(int i=0;i<matrixsize;i++)
	{

		cleft[i]=-1;
		cright[i]=-1;
		c[i]=MAX;
		events[i]=-1;
		in[i]=MAX;
		inindex[i]=-1;
	}
	cout<<"[DH-5] Arrays initialized, processing "<<domainTree->leafnodes.size()<<" domain tree leaves..."<<endl;

	node* pa;

	for (int i=0;i<domainTree->leafnodes.size();i++){
		// Validate that domain leaf is properly mapped to a gene tree node
		if (domainTree->leafnodes[i]->mappedNode < 0 || domainTree->leafnodes[i]->mappedNode >= genesize)
		{
			cout<<"ERROR: Domain tree leaf '"<<domainTree->leafnodes[i]->name<<"' was not properly mapped to gene tree!"<<endl;
			cout<<"       mappedNode="<<domainTree->leafnodes[i]->mappedNode<<" (expected 0-"<<genesize-1<<")"<<endl;
			cout<<"       Check that domain names match gene tree node names in format: domainname_genename_familyid"<<endl;
			delete[] c;
			delete[] cleft;
			delete[] cright;
			delete[] events;
			delete[] in;
			delete[] r;
			delete[] inindex;
			return -1;
		}
		c[domainTree->leafnodes[i]->key*(genesize)+ domainTree->leafnodes[i]->mappedNode]=0;
		events[domainTree->leafnodes[i]->key*(genesize)+ domainTree->leafnodes[i]->mappedNode]=0;
		in[domainTree->leafnodes[i]->key*(genesize)+ domainTree->leafnodes[i]->mappedNode]=0;
		inindex[domainTree->leafnodes[i]->key*(genesize)+ domainTree->leafnodes[i]->mappedNode]=domainTree->leafnodes[i]->mappedNode;
		pa=genepointers[domainTree->leafnodes[i]->mappedNode];
		while(pa!=NULL)
		{
			in[domainTree->leafnodes[i]->key*(genesize)+ pa->key]=genepointers[domainTree->leafnodes[i]->mappedNode]->depth-pa->depth;
			inindex[domainTree->leafnodes[i]->key*(genesize)+ pa->key]=domainTree->leafnodes[i]->mappedNode;
			pa=pa->parent;
		}
	}
	cout<<"[DH-6] Domain leaf processing complete, calling MaxPostorderDomain..."<<endl;
	MaxPostorderDomain(inindex,r,events,c,in,cleft,cright,domainTree,geneTrees,speciesTree,genesize,domainTree->root,domainpointers,genepointers,speciespointers,twoTreeTransferCost,OneTreeTransferCost,domainDuplicationcost,geneDuplicationcost,domainLosscost,geneLosscost);
	cout<<"[DH-7] MaxPostorderDomain complete"<<endl;

	// Check if any valid costs were computed for the root and its children
	int rootKey = domainTree->root->key;
	int leftKey = domainTree->root->left->key;
	int rightKey = domainTree->root->right->key;

	int rootValidCosts = 0, leftValidCosts = 0, rightValidCosts = 0;
	int rootMinCost = MAX, leftMinCost = MAX, rightMinCost = MAX;

	for (int i = 0; i < genesize; i++) {
		int rootCost = c[rootKey*genesize+i];
		int leftCost = c[leftKey*genesize+i];
		int rightCost = c[rightKey*genesize+i];

		if (rootCost < MAX) {
			rootValidCosts++;
			if (rootCost < rootMinCost) rootMinCost = rootCost;
		}
		if (leftCost < MAX) {
			leftValidCosts++;
			if (leftCost < leftMinCost) leftMinCost = leftCost;
		}
		if (rightCost < MAX) {
			rightValidCosts++;
			if (rightCost < rightMinCost) rightMinCost = rightCost;
		}
	}

	cout<<"[DH-7b] Domain root (key="<<rootKey<<") has "<<rootValidCosts<<" gene mappings with cost < MAX (min cost="<<rootMinCost<<")"<<endl;
	cout<<"[DH-7c] Root left child (key="<<leftKey<<") has "<<leftValidCosts<<" gene mappings with cost < MAX (min cost="<<leftMinCost<<")"<<endl;
	cout<<"[DH-7d] Root right child (key="<<rightKey<<") has "<<rightValidCosts<<" gene mappings with cost < MAX (min cost="<<rightMinCost<<")"<<endl;

	if (leftValidCosts == 0 || rightValidCosts == 0) {
		cout<<"[DH-7e] ERROR: Root's children have no valid mappings! This propagates up to root."<<endl;
		cout<<"[DH-7f] Checking a few leaf nodes to see if THEY have valid mappings..."<<endl;
		for (int leafIdx = 0; leafIdx < min(5, (int)domainTree->leafnodes.size()); leafIdx++) {
			int leafKey = domainTree->leafnodes[leafIdx]->key;
			int leafValidCosts = 0;
			for (int i = 0; i < genesize; i++) {
				if (c[leafKey*genesize+i] < MAX) {
					leafValidCosts++;
				}
			}
			cout<<"[DH-7g] Leaf "<<leafIdx<<" (key="<<leafKey<<", name="<<domainTree->leafnodes[leafIdx]->name<<") has "<<leafValidCosts<<" valid mappings"<<endl;
		}
	}

	//UpperPostorderDomain(2,MaxNum,r,events,upper,in,out,cleft,cright,domainTree,geneTrees,speciesTree,genesize,domainTree->root,domainpointers,genepointers,speciespointers,twoTreeTransferCost,OneTreeTransferCost,domainDuplicationcost,geneDuplicationcost,domainLosscost,geneLosscost);

	int movenum=0;
	int movestep=0;
	int transnum=0;
	int added=0;
	int Lcascore=0;
	int diffnum=0;
	int moreduplicationnumber=0;
	//cout<<"Traversal finished "<<endl;
	//cout<<"Fine"<<endl;
	cout<<"[DH-8] Calling UptracebackReal..."<<endl;
	int Upresult=UptracebackReal(flag,speciesTree, &added, r, domainTree, geneTrees, c, events, cleft, cright, genesize, domainpointers, speciespointers, genepointers,domainFileName,twoTreeTransferCost,OneTreeTransferCost,domainDuplicationcost,geneDuplicationcost,domainLosscost,geneLosscost);
	cout<<"[DH-9] UptracebackReal returned"<<endl;


	delete c;
	delete events;
	delete cleft;
	delete cright;
	delete in;
	return Upresult;
}

#endif
