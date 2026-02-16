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
			fout<<", Recipient -> "<<r[domainnode->key*genesize+mappedindex]->name;



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
		UpLocaltraceback(dup,trans,trans2,flag,low,geneflag,genestep,r,domainnode->left,cleft[domainnode->key*genesize+mappedindex],events,cleft,cright,genesize,speciespointers,genepointers,domainFileName,twoTreeTransferCost,OneTreeTransferCost,domainDuplicationcost,geneDuplicationcost,domainLosscost,geneLosscost);
		UpLocaltraceback(dup,trans,trans2,flag,low,geneflag,genestep,r,domainnode->right,cright[domainnode->key*genesize+mappedindex],events,cleft,cright,genesize,speciespointers,genepointers,domainFileName,twoTreeTransferCost,OneTreeTransferCost,domainDuplicationcost,geneDuplicationcost,domainLosscost,geneLosscost);

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
	if(flag){
		fout.open(domainFileName.c_str());
		fout<<"Domain Tree: "<<endl;
		fout.close();

		WritIntoFile(domaintree->root,domainFileName);
		fout.open(domainFileName.c_str(),ios::app);
		fout.close();
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
	}
	int low=MAX;
	int rootmapindex;

	for (int i=0;i<genesize;i++){
		if (low>c[domaintree->root->key*(genesize)+i])
		{
			low = c[domaintree->root->key*(genesize)+i];
			rootmapindex=i;
		}
	}
	//cout<<"Root maps to "<<rootmapindex<<" cost "<<low<<endl;

	int *geneflag =	new int[genesize];
	int *genestep =	new int[genesize];
	for (int i=0;i<genesize;i++)
	{
		geneflag[i]=0;
		genestep[i]=0;
	}

	low=0;
	int dup=0, trans=0, trans2=0;

	UpLocaltraceback(&dup, &trans, &trans2, flag,&low,geneflag,genestep,r,domaintree->root,rootmapindex,events,cleft,cright,genesize,speciespointers,genepointers,domainFileName,twoTreeTransferCost,OneTreeTransferCost,domainDuplicationcost,geneDuplicationcost,domainLosscost,geneLosscost);


	// Update Gene to Species
	node *domainnode;
	int allsetflag=1;
	node *Iter;

	//while (allsetflag)
	{
		allsetflag=0;
		for(int i=0;i<domaintree->leafnodes.size()*2-1;i++)
		{
			domainnode=domainpointers[i];
			if(events[domainnode->key*genesize+domainnode->mappedNode]==3 || events[domainnode->key*genesize+domainnode->mappedNode]==4)
			{
				if(CheckGeneCost(geneflag,genestep, r[domainnode->key*genesize+domainnode->mappedNode], genepointers[domainnode->mappedNode],speciespointers)==1)
					allsetflag=1;
			}
		}

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
	string filename = domainFileName;
	int domainsize=domainTree->leafnodes.size()*2-1;
	int speciessize=speciesTree->leafnodes.size()*2-1;

	int matrixsize=domainsize*genesize;

	int *c		= new int[matrixsize];
	int *cleft		= new int[matrixsize];
	int *cright		= new int[matrixsize];
	int *events		= new int[matrixsize];
	int *in		= new int[matrixsize];
	node **r		= new node*[matrixsize];
	int *inindex	= new int[matrixsize];
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
	MaxPostorderDomain(inindex,r,events,c,in,cleft,cright,domainTree,geneTrees,speciesTree,genesize,domainTree->root,domainpointers,genepointers,speciespointers,twoTreeTransferCost,OneTreeTransferCost,domainDuplicationcost,geneDuplicationcost,domainLosscost,geneLosscost);
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
	int Upresult=UptracebackReal(flag,speciesTree, &added, r, domainTree, geneTrees, c, events, cleft, cright, genesize, domainpointers, speciespointers, genepointers,domainFileName,twoTreeTransferCost,OneTreeTransferCost,domainDuplicationcost,geneDuplicationcost,domainLosscost,geneLosscost);


	delete c;
	delete events;
	delete cleft;
	delete cright;
	delete in;
	return Upresult;
}

#endif
