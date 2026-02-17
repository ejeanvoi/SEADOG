#ifndef filefunctions_H
#define filefunctions_H

#include <string.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "tree.h"
#include <sstream>
#include <string>
ofstream testfout;
ofstream fout;

void writetree(node *node,int dir);
void WritIntoFile(node *node, string filename);
void setgenekey(tree **genetrees, int treenumber, node** pointers);
int setLocalkey(tree *currenttree, node* root, int *start);
int setkey(tree *currenttree, node* root, int *start, int depth, node** pointers);



node* LCA(node* node1, node* node2)
{
	if (node1->key==node2->key)
		return node1;
	node* lca = node1;
	if (lca->isleaf)
		lca=lca->parent;
	while (!lca->isroot)
	{
		if (node2->DoubleOrder1>=lca->DoubleOrder1 && node2->DoubleOrder2<=lca->DoubleOrder2)
			break;
		lca=lca->parent;
	}
	return lca;
}

float GetMin(float a, float b)
{
	return (a<=b?a:b);
}

int GetMin(int a, int b)
{
	return (a<=b?a:b);
}

int GetMax(int a, int b)
{
	return (a>=b?a:b);
}

int comparable(node *des, int key, node** pointers)
{
	if (des->key==key) return 1;
	if (des->DoubleOrder2<pointers[key]->DoubleOrder2 && des->DoubleOrder1>pointers[key]->DoubleOrder1)
		return 1;
	else return 0;

}

int lowestkey(node* node)
{
	if (node->isleaf)
		return node->key;
	else return lowestkey(node->left);
}

int highestkey(node* node)
{
	if (node->isleaf)
		return node->key;
	else return highestkey(node->right);
}


void setgenekey(tree **genetrees, int treenumber, node** pointers)	// a combination of key, depth and pointers.
{
	int totalindex=0;
	int localindex=0;
	int j=0;
	int depth=0;
	for (j=0;j<treenumber;j++)
	{
		localindex=0;
		setkey(genetrees[j],genetrees[j]->root,&totalindex,depth,pointers);
		setLocalkey(genetrees[j],genetrees[j]->root,&localindex);
		localindex=0;
		DoubleTraverse(genetrees[j]->root,&localindex);
	}
}

int setLocalkey(tree *currenttree, node* root, int *start)
{
	//*start=0;
	//cout<<root->isleaf<<endl;
	if(root->isleaf==0)
	{
		setLocalkey(currenttree,root->left,start);
		//cout<<"Left ends "<<endl;
		setLocalkey(currenttree,root->right,start);
		//cout<<"right ends "<<endl;
	}
	root->localkey=*start;
	*start=*start+1;
	//cout<<"root : This node "<<root->key<<endl;
	return 0;
}

int setkey(tree *currenttree, node* root, int *start, int depth, node** pointers)
{
	//int localstart=0;
	//cout<<root->isleaf<<endl;
	if(root->isleaf==0)
	{
		setkey(currenttree,root->left,start,depth+1,pointers);
		//cout<<"Left ends "<<endl;
		setkey(currenttree,root->right,start,depth+1,pointers);
		//cout<<"right ends "<<endl;
	}
	root->key=*start;
	root->depth=depth;
	if(root->isroot==0)
		root->sibling=GetSibling(root);
	root->isfakeleaf=root->isleaf;
	root->visited=-1;
	root->visited2=-1;
	//root->group=0;
	root->itsroot=currenttree->root;
	root->isreceiver=0;
	root->istransfer=0;
	pointers[*start]=root;
	root->issensitive=0;
	*start=*start+1;

	if(root->isleaf==0)
	{
		string Intername;
		if(currenttree->type==0)
			Intername="D";
		if(currenttree->type==1)
			Intername="G";
		if(currenttree->type==2)
			Intername="S";
		stringstream ss;
		ss << root->key;
		Intername+=ss.str()+"_";
		ss.str("");
		ss.clear();
		ss << root->left->key;


		Intername+=ss.str()+"_";
		ss.str("");
		ss.clear();
		ss << root->right->key;
		Intername+=ss.str();
		char *Inter = new char[Intername.size()+1];
		for(int i=0;i<Intername.size();i++)
		{
			Inter[i]=Intername[i];
		}
		Inter[Intername.size()]='\0';
		root->name=Inter;
		//cout<<root->name<<endl;
		Inter=NULL;
	}
	//cout<<"root : This node "<<root->key<<endl;
	return 0;
}


void WritIntoFile(node *node, string filename)
{
	testfout.open(filename.c_str(), ios::app);
	writetree(node,0);
	testfout.close();
}

void writetree(node *node, int dir)
{
	if (node->isroot==1)
	{
		testfout<<"(";
		if(node->isleaf==0)
		{
			writetree(node->left,0);
			writetree(node->right,1);
		}
		testfout<<")";
		 testfout<<node->name;
		testfout<<";"<<endl;
	}
	else if (node->isleaf==1)
	{
		if (dir==0)
		{
			if (node->isleaf)
			{

				for (int i=0;i<strlen(node->name);i++)
					testfout<<node->name[i];
				if(node->mappedSpecies)
					testfout<<"_"<<node->mappedSpecies;
			}
			testfout<<"_"<<node->key;
			testfout<<",";
		}
		else if (dir==1)
		{
			if (node->isleaf)
			{
				for (int i=0;i<strlen(node->name);i++)
					testfout<<node->name[i];
				if(node->mappedSpecies)
					testfout<<"_"<<node->mappedSpecies;
			}
			testfout<<"_"<<node->key;
		}

	}
	else if (node->isleaf==0)
	{
		testfout<<"(";
		writetree(node->left,0);
		writetree(node->right,1);
		testfout<<")";
        testfout<<node->name;
		if (dir==0)
		{
			testfout<<",";
		}
	}
}

int getlength(string filename)
{
	int linenumber=0;
	char c;
	ifstream thisfin;
	thisfin.open(filename.c_str());
	if (!thisfin)
	{
		//cout<<"cannot open "<<filename<<endl;
		return 0;
	}
	for(c=thisfin.get();c!=EOF;c=thisfin.get())
	{
		if (c=='\n')
			linenumber++;
	}
	thisfin.close();
	return linenumber;
}

void MapLeafNodes(tree **genetrees, tree* speciestree, int genetreeNum)
{
	ofstream fout;
	fout.open("Maplog.txt", ios::app);
	int matched=0;
	int flag=0;
	for (int geneIter=0;geneIter<genetreeNum;geneIter++)
	{

		for (int nodeIter=0;nodeIter<genetrees[geneIter]->leafnodes.size();nodeIter++)
		{
			//cout<<geneIter<<" : ";
			matched=0;
			for (int i=0;i<speciestree->leafnodes.size();i++)
			{
				flag=1;
				// FIX: Check for NULL and length mismatch FIRST
				if (genetrees[geneIter]->leafnodes[nodeIter]->mappedSpecies == NULL)
				{
					flag=0;
				}
				else if (strlen(speciestree->leafnodes[i]->name) != strlen(genetrees[geneIter]->leafnodes[nodeIter]->mappedSpecies))
				{
					flag=0;
				}
				else
				{
					for(int k=0;k<strlen(speciestree->leafnodes[i]->name);k++)
					{
						if (speciestree->leafnodes[i]->name[k]!=genetrees[geneIter]->leafnodes[nodeIter]->mappedSpecies[k])
						{
							flag=0;
							break;
						}
					}
				}
				if(flag)
				{
					genetrees[geneIter]->leafnodes[nodeIter]->mappedNode=speciestree->leafnodes[i]->key;
					/*
					for(int k=0;k<strlen(genetrees[geneIter]->leafnodes[nodeIter]->name);k++)
					{
						cout<<genetrees[geneIter]->leafnodes[nodeIter]->name[k];
					}
					cout<<"  mapped to ";
					for(int k=0;k<strlen(speciestree->leafnodes[i]->name);k++)
					{
						cout<<speciestree->leafnodes[i]->name[k];
					}
					cout<<endl;
					*/
					matched=1;
					break;
				}
			}

			if (!matched)
			{
				// FIX: Check for NULL before accessing
				if (genetrees[geneIter]->leafnodes[nodeIter]->mappedSpecies != NULL)
				{
					for(int k=0;k<strlen(genetrees[geneIter]->leafnodes[nodeIter]->mappedSpecies);k++)
					{
						fout<<genetrees[geneIter]->leafnodes[nodeIter]->mappedSpecies[k];
					}
				}
				fout<<" failed to match a node !"<<endl;
			}
		}
	}
	fout.close();
}

void MapLeafNodes(tree *domaintree, tree** genetrees, int genetreeNum)
{
	ofstream fout;
	fout.open("Maplog.txt");
	int matched=0;
	int flag=0;

	int pos;
	int domainflag=0;
	for (int i=0;i<domaintree->leafnodes.size();i++)
	{
		//cout<<i<<" : ";

		//cout<<domaintree->leafnodes[i]->name;

		//cout<<endl;

		matched=0;
		for (int geneIter=0;geneIter<genetreeNum;geneIter++)
		{
			for (int nodeIter=0;nodeIter<genetrees[geneIter]->leafnodes.size();nodeIter++)
			{
				flag=1;
				domainflag=0;
				for(int k=0;k<strlen(domaintree->leafnodes[i]->name);k++)
				{
					if (domaintree->leafnodes[i]->name[k]=='_') {domainflag=1;pos=k+1; continue;}

					if(domainflag==1)
					{
						// FIX: Check bounds before accessing gene name
						if (k-pos >= strlen(genetrees[geneIter]->leafnodes[nodeIter]->name) ||
						    domaintree->leafnodes[i]->name[k]!=genetrees[geneIter]->leafnodes[nodeIter]->name[k-pos])
						{
							flag=0;
							break;
						}
					}

				}
				if(flag)
				{
					domaintree->leafnodes[i]->mappedGeneFamily=geneIter;
					domaintree->leafnodes[i]->mappedNode=genetrees[geneIter]->leafnodes[nodeIter]->key;
					genetrees[geneIter]->leafnodes[nodeIter]->valid=1;
					/*
					cout<<i<<" ";
					for(int k=0;k<strlen(domaintree->leafnodes[i]->name);k++)
					{
						cout<<domaintree->leafnodes[i]->name[k];
					}
					cout<<"  mapped to ";
					for(int k=0;k<strlen(genetrees[geneIter]->leafnodes[nodeIter]->name);k++)
					{
						cout<<genetrees[geneIter]->leafnodes[nodeIter]->name[k];
					}
					cout<<endl;
					*/

					matched=1;
					break;
				}
			}
			if (matched)
				break;
		}
		if (!matched)
		{
			for(int k=0;k<strlen(domaintree->leafnodes[i]->name);k++)
			{
				fout<<domaintree->leafnodes[i]->name[k];
			}
			fout<<" failed to match a node !"<<endl;
		}
	}
	fout.close();
}

#endif
