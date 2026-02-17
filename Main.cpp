/*
 * *   Copyright (C) 2017 Lei Li (lei.li@uconn.edu) and Mukul S. Bansal (mukul.bansal@uconn.edu).
 * *
 * *   This program is free software: you can redistribute it and/or modify
 * *   it under the terms of the GNU General Public License as published by
 * *   the Free Software Foundation, either version 3 of the License, or
 * *   (at your option) any later version.
 * *
 * *   This program is distributed in the hope that it will be useful,
 * *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 * *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * *   GNU General Public License for more details.
 * *
 * *   You should have received a copy of the GNU General Public License
 * *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * */


#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "tree.h"
#include "libs.h"
#include "Parser.h"
#include "DL.h"
//#include "SelectGene.h"
#include "DynamicHeuristic.h"
#include "argument.h"

int main(int argc,char *argv[])
{
	/*************************************************
	The event costs
	*************************************************/

	int TCostOneTree;		// cost for domains transfer within one gene family.
	int TCostTwoTree;		// cost for domains transfer between gene families.
	int domainDCost;
	int domainLCost;
	int geneLCost;
	int geneDCost;


	Argument::add(argc, argv);

	bool isUnrooted = (Argument::find("-u") != NULL) || (Argument::find("-U") != NULL);

	const Argument *arg1 = Argument::find("-d");
	string domainFileName;
	if (arg1 == NULL) {cout<<"No Input Domain Tree"; return 0;}
	else {
		arg1->convert(domainFileName);
	}

	string geneFileName;
	const Argument *arg2 = Argument::find("-g");
	if (arg2 == NULL) {cout<<"No Input Gene Trees Folder"; return 0;}
	else {
		arg2->convert(geneFileName);
	}

	string speciesFileName;
	const Argument *arg3 = Argument::find("-s");
	if (arg3 == NULL) {cout<<"No Input Species Tree"; return 0;}
	else {
		arg3->convert(speciesFileName);
	}


	string filename;

	const Argument *arg4 = Argument::find("-o");
	if (arg4 == NULL) {filename=domainFileName+".output";}
	else {
		arg4->convert(filename);
	}
	cout<<"Output file is: "<<filename<<endl;
	const Argument *argDD = Argument::find("-DD");
	if (argDD == NULL) domainDCost=2;
	else
	{
		argDD->convert(domainDCost);
	}

	const Argument *argDL = Argument::find("-DL");
	if (argDL == NULL) domainLCost=1;
	else
	{
		argDL->convert(domainLCost);
	}
	const Argument *argDTA = Argument::find("-DTA");
	if (argDTA == NULL) TCostOneTree=4;
	else
	{
		argDTA->convert(TCostOneTree);
	}
	const Argument *argGD = Argument::find("-GD");
	if (argGD == NULL) geneDCost=2;
	else
	{
		argGD->convert(geneDCost);
	}
	const Argument *argGL = Argument::find("-GL");
	if (argGL == NULL) geneLCost=1;
	else
	{
		argGL->convert(geneLCost);
	}
	const Argument *argDTB = Argument::find("-DTB");
	if (argDTB == NULL) TCostTwoTree=6;
	else

	{
		argDTB->convert(TCostTwoTree);
	}

	/************************************************/
	//int iflarge;

	//domain tree parser
	//string domainFileName="./1452.treefix.tree.FBGeneID";
	//string domainFileName=argv[1];
	cout<<"Input Domain file is: "<<domainFileName<<endl;
	tree *origindomainTree = ParserToTree(domainFileName, 0, 1);	// build the tree
	if (origindomainTree == NULL) { cout<<"Failed to parse domain tree"<<endl; return 1; }
	cout<<"Domain Tree built, size: "<<origindomainTree->leafnodes.size()*2-1<<endl;
	node **origindomainpointers = new node* [origindomainTree->leafnodes.size()*2-1];	// A matrix recording all nodes
	int start=0;
	setkey(origindomainTree,origindomainTree->root,&start,0,origindomainpointers);	// Domain tree nodes only have "key"
	start=0;
	DoubleTraverse(origindomainTree->root,&start);

	cout<<"Mapped Trees are : "<<endl;
	for (int ite=0;ite<origindomainTree->mappedTrees.size();ite++)
		cout<<origindomainTree->mappedTrees[ite]<<endl;

	//if(domainTree->leafnodes.size()>100) return 0;
	//gene tree part
	int totalGeneNodeNum=0;	// Number of all nodes in all gene families.
	tree **geneTrees = new tree*[origindomainTree->mappedTrees.size()];	//define the valid gene trees.
	string ThisgeneFileName;
	for (int geneTreeIter=0;geneTreeIter<origindomainTree->mappedTrees.size();geneTreeIter++)	//for each valid gene tree
	{
		ThisgeneFileName=geneFileName;
		ThisgeneFileName.append(origindomainTree->mappedTrees[geneTreeIter]).append(".tree");
		geneTrees[geneTreeIter]=ParserToTree(ThisgeneFileName, 1, 0);
		if(geneTrees[geneTreeIter]==NULL) { cout<<"Failed to parse gene tree: "<<ThisgeneFileName<<endl; return 1; }
		geneTrees[geneTreeIter]->root->treename=ThisgeneFileName;
		//SelectFly(geneTrees[geneTreeIter]);	//delete non-fly nodes in each gene tree
		totalGeneNodeNum+=2*geneTrees[geneTreeIter]->leafnodes.size()-1;
		//cout<<"Gene tree total size add up to: "<<totalGeneNodeNum<<endl;
	}
	node **genepointers = new node* [totalGeneNodeNum];
	setgenekey(geneTrees,origindomainTree->mappedTrees.size(),genepointers);
	cout<<"Gene Trees built, total size: "<<totalGeneNodeNum<<endl;
	/*
	cout<<"Gene Tree Topology: index left right parent for internal nodes: "<<endl;
	for(int index=0;index<geneTrees[0]->leafnodes.size()*2-1;index++)
	{
		if(genepointers[index]->isleaf==0 && genepointers[index]->isroot==0)
		cout<<genepointers[index]->key<<" "<<genepointers[index]->left->key<<" "<<genepointers[index]->right->key<<" "<<genepointers[index]->parent->key<<endl;
		else if (genepointers[index]->isroot==1)
			cout<<genepointers[index]->key<<endl;
		else
			cout<<genepointers[index]->key<<" "<<genepointers[index]->parent->key<<endl;
	}
	*/

	//if(totalGeneNodeNum>400) return 0;
	// Species tree part, similar to domain tree
	tree *speciesTree = ParserToTree(speciesFileName, 0, 0);
	if (speciesTree == NULL) { cout<<"Failed to parse species tree"<<endl; return 1; }
	cout<<"Species Tree built, size: "<<speciesTree->leafnodes.size()*2-1<<endl;
	node **speciespointers = new node* [2*speciesTree->leafnodes.size()-1];	//postorder
	start=0;
	setkey(speciesTree,speciesTree->root,&start,0,speciespointers);
	start=0;
	DoubleTraverse(speciesTree->root,&start);
	start=0;

	MapLeafNodes(geneTrees,speciesTree,origindomainTree->mappedTrees.size());
	MapLeafNodes(origindomainTree,geneTrees,origindomainTree->mappedTrees.size());

	// Find LCA mapping between gene trees and the species tree
	for (int i=0;i<origindomainTree->mappedTrees.size();i++){
		DLdynamicalgorithm(geneTrees[i],speciesTree,2*geneTrees[i]->leafnodes.size()-1,2*speciesTree->leafnodes.size()-1,speciespointers,geneDCost,geneLCost);
	}
	//cout<<geneTrees[0]->root->key<<endl;
	// Start the algorithm

	tree **DomainTrees = new tree*[origindomainTree->leafnodes.size()*2-2];
	node **domainpointers = new node* [origindomainTree->leafnodes.size()*2-1];
	node **BestIndex = new node* [origindomainTree->leafnodes.size()*2-1];
	node* iter;
	node* savedChild;
	node* savedParent;
	node* Parent;

	int loops=isUnrooted?origindomainTree->leafnodes.size()*2-2:0;

	int bestIndex=0;
	int bestScore=MAX;
	int score;
	for(int rootindex=0;rootindex<loops;rootindex++)	// For each rooting position
	{

		DomainTrees[rootindex] = ParserToTree(domainFileName, 0, 1);	// build the tree
		// Now the magic: change the root
		if(origindomainpointers[rootindex]->isleft==0 && origindomainpointers[rootindex]->parent->isroot==1)
		{
			cout<<"Skipping root index: "<<rootindex<<endl;
			continue;
		}
		start=0;
		setkey(DomainTrees[rootindex],DomainTrees[rootindex]->root,&start,0,domainpointers);

		//cout<<"Doing the root index: "<<rootindex<<" Parent "<<domainpointers[rootindex]->parent->key<<endl;

		if(rootindex<=DomainTrees[rootindex]->root->left->key)
			savedChild=DomainTrees[rootindex]->root->right;
		else
			savedChild=DomainTrees[rootindex]->root->left;

		iter=domainpointers[rootindex]->parent;
		DomainTrees[rootindex]->root->left=domainpointers[rootindex];
		domainpointers[rootindex]->parent=DomainTrees[rootindex]->root;

		if(iter->isroot==0)
			DomainTrees[rootindex]->root->right=iter;

		//iter=domainpointers[rootindex]->parent;
		savedParent=DomainTrees[rootindex]->root;



		while(iter->isroot==0)
		{
			Parent=iter->parent;

			if(iter->left->parent==iter)
			{
				if(Parent->isroot==1)
					iter->right=savedChild;
				else iter->right=Parent;
			}
			else if(iter->right->parent==iter)
			{
				if(Parent->isroot==1)
					iter->left=savedChild;
				else iter->left=Parent;
			}
			iter->parent=savedParent;
			savedParent=iter;
			iter=Parent;
		}
		savedChild->parent=savedParent;

		start=0;
		setkey(DomainTrees[rootindex],DomainTrees[rootindex]->root,&start,0,domainpointers);
		/*
		cout<<"Tree Topology: index, left, right, parent, (children for internal nodes only): "<<endl;
		for(int index=0;index<origindomainTree->leafnodes.size()*2-1;index++)
		{
			if(domainpointers[index]->isleaf==0 && domainpointers[index]->isroot==0)
			cout<<domainpointers[index]->key<<" "<<domainpointers[index]->left->key<<" "<<domainpointers[index]->right->key<<" "<<domainpointers[index]->parent->key<<endl;
			else if (domainpointers[index]->isroot==1)
				cout<<domainpointers[index]->key<<endl;
			else
				cout<<domainpointers[index]->key<<" "<<domainpointers[index]->parent->key<<endl;
		}
		*/
		MapLeafNodes(DomainTrees[rootindex],geneTrees,origindomainTree->mappedTrees.size());
		score=DynamicHeuristic(0,DomainTrees[rootindex],geneTrees,speciesTree,totalGeneNodeNum,domainpointers,genepointers,speciespointers,TCostTwoTree,TCostOneTree,domainDCost,geneDCost,domainLCost,geneLCost,filename);
		if(score<bestScore)
		{
			bestScore=score;
			bestIndex=rootindex;
			for(int i=0;i<origindomainTree->leafnodes.size()*2-1;i++)
			{
				BestIndex[i]=domainpointers[i];
			}
		}
		//JessicaSimulator(DomainTrees[rootindex],geneTrees,speciesTree,totalGeneNodeNum,domainpointers,genepointers,speciespointers,TCostTwoTree,TCostOneTree,domainDCost,geneDCost,domainLCost,geneLCost,filename);

		//GeneUpperBound(DomainTrees[rootindex],geneTrees,speciesTree,totalGeneNodeNum,domainpointers,genepointers,speciespointers,TCostTwoTree,TCostOneTree,domainDCost,geneDCost,domainLCost,geneLCost,filename);
		//cout<<"Computing Domain Content"<<endl;
		//OurDomainContentSpecies(DomainTrees[rootindex], speciesTree, filename, domainpointers, speciespointers);
	}

	//GeneUpperBound(MaxNum,domainTree,geneTrees,speciesTree,totalGeneNodeNum,domainpointers,genepointers,speciespointers,TCostTwoTree,TCostOneTree,domainDCost,geneDCost,domainLCost,geneLCost,filename);
	if(loops==0)
		DynamicHeuristic(1,origindomainTree,geneTrees,speciesTree,totalGeneNodeNum,origindomainpointers,genepointers,speciespointers,TCostTwoTree,TCostOneTree,domainDCost,geneDCost,domainLCost,geneLCost,filename);
	else
		DynamicHeuristic(1,DomainTrees[bestIndex],geneTrees,speciesTree,totalGeneNodeNum,BestIndex,genepointers,speciespointers,TCostTwoTree,TCostOneTree,domainDCost,geneDCost,domainLCost,geneLCost,filename);
	//OldOne(domainTree,geneTrees,speciesTree,totalGeneNodeNum,domainpointers,genepointers,speciespointers,TCostTwoTree,TCostOneTree,domainDCost,geneDCost,domainLCost,geneLCost,filename);

	//JessicaSimulator(MaxNum,domainTree,geneTrees,speciesTree,totalGeneNodeNum,domainpointers,genepointers,speciespointers,TCostTwoTree,TCostOneTree,domainDCost,geneDCost,domainLCost,geneLCost,domainFileName);
	cout<<"Done"<<endl;
	delete[] speciespointers;
	delete[] genepointers;
	delete[] domainpointers;
	delete[] origindomainpointers;
	delete[] BestIndex;
	delete[] geneTrees;
	delete[] DomainTrees;

	return 0;
}
