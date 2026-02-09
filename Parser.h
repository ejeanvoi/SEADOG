#ifndef Parser_H
#define Parser_H

#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include "tree.h"
#include "libs.h"


tree* ParserToTree(string filename, int hasPrefix, int hasSuffix)
{
	// hasPrfix==1 means gene node in format: genename_speciesname
	// hasSuffix==1 means domain tree in format: domain name_genename_geneTreeIndex

	ifstream Pfin;
	int i,j,k;
	Pfin.open(filename.c_str());
	if (!Pfin)
	{
		cout<<filename<<" Not Found"<<endl;
		return NULL;
	}

	char c;
	string wholeFile;

	Pfin.seekg(0, Pfin.end);
	int file_size = Pfin.tellg();
	char* file_buf = new char [file_size+1];
	//memset(file_buf, 0, file_size+1);
	//cout<<file_size<<endl;
	Pfin.seekg(0, ios::beg);
	j=0;
	int ignore=0;
	for (i=0;i<file_size;i++)
	{
		//cout<<c;
		c=Pfin.get();
		if(c==' ' || c=='\n' || c=='\t' ||c==EOF)
			continue;
		if(c==':') {
			ignore=1;
		}
		if(ignore==1 && (c==','||c==')'))
			ignore=0;
		if(!ignore)
		{
			//cout<<c;
			file_buf[j++]=c;
		}
	}
	//cout<<endl;
	file_buf[j]='\0';
	//cout<<"j = "<<j<<endl;
	wholeFile.append(file_buf);
	//cout<<wholeFile<<endl;
	delete []file_buf;

	vector<int> lastLeft;
	vector<int> lastComma;
	int indirection=0;
	int outdirection=0;

	tree *newTree = new tree();

	vector<node*> savedLeftNode;
	vector<node*> savedRightNode;

	int nameLengthLeft;
	int nameLengthRight;

	int start,end;
	//cout<<wholeFile.size()<<endl;
	for (i=0;i<wholeFile.size();i++)
	{
		if (wholeFile[i]=='('){
			lastLeft.push_back(i);
		}
		else if(wholeFile[i]==','){
			lastComma.push_back(i);
		}
		else if(wholeFile[i]==')'){    // When there is a node

			for (j=i+1;j<wholeFile.size();j++){    // Check if this node is left or right
				if (wholeFile[j]==','){
					outdirection=0;
					break;
				}
				else if(wholeFile[j]==')'){
					outdirection=1;
					break;
				}
			}
			if (wholeFile[lastLeft.back()+1]=='(' && wholeFile[lastComma.back()+1]=='(')
				indirection=3;
			else if (wholeFile[lastLeft.back()+1]=='(' && wholeFile[lastComma.back()+1]!='(')
				indirection=2;
			else if (wholeFile[lastLeft.back()+1]!='(' && wholeFile[lastComma.back()+1]=='(')
				indirection=1;
			else if (wholeFile[lastLeft.back()+1]!='(' && wholeFile[lastComma.back()+1]!='(')
				indirection=0;

			node *newNode= new node();
			//char interName[] = "Inter";
			//newNode->name=interName;
			newNode->isroot=0;
			newTree->root=newNode;

			if (indirection==0) // two leaf nodes
			{
				node *newLeft = new node();
				node *newRight = new node();
				newLeft->isleft = 1;
				newRight->isleft=0;
				newLeft->isleaf=1;
				newRight->isleaf=1;
				newLeft->isroot=0;
				newRight->isroot=0;
				newLeft->parent=newNode;
				newRight->parent=newNode;
				newNode->isleaf=0;
				newNode->left=newLeft;
				newNode->right=newRight;
				newTree->leafnodes.push_back(newLeft);
				newTree->leafnodes.push_back(newRight);

				if (outdirection==0)
				{
					newNode->isleft=1;
					savedLeftNode.push_back(newNode);
				}
				else if (outdirection==1)
				{
					newNode->isleft=0;
					savedRightNode.push_back(newNode);
				}

				start=lastLeft.back()+1;
				end=lastComma.back();

				for (j=start;j<lastComma.back();j++){
						if (wholeFile[j]==':'){
							end=j;
							break;}
				}


				if (hasPrefix){
					for (j=lastLeft.back()+1;j<lastComma.back();j++){
							if (wholeFile[j]=='_'||wholeFile[j]=='-'){
								end=j;
								//start=j+1;
								break;}
					}
					char *SpeciesID = new char[lastComma.back()-end];
					for (j=end+1;j<lastComma.back();j++){
							SpeciesID[j-end-1]=wholeFile[j];
							//cout<<SpeciesID[j-end-1];
					}
					//cout<<endl;
					SpeciesID[lastComma.back()-end-1]='\0';
					newLeft->mappedSpecies=SpeciesID;
				}

				if (hasSuffix){
					for (j=lastComma.back()-1;j>=start;j--){
							if (wholeFile[j]=='_'){
								end=j;
								break;}
					}
					char *GeneID = new char[lastComma.back()-end];
					for (j=end+1;j<lastComma.back();j++){
							GeneID[j-end-1]=wholeFile[j];
					}
					GeneID[lastComma.back()-end-1]='\0';
					//int NewID;
					//sscanf(GeneID,"%d",&NewID);
					newLeft->GeneID=GeneID;
					newTree->mappedTrees.push_back(GeneID);
				}

				nameLengthLeft=end-start;
				char *newNameLeft = new char[nameLengthLeft+1];
				
				for (j=start;j<end;j++)
				{
					newNameLeft[j-start]=wholeFile[j];
					//cout<<wholeFile[j];
				}
				//cout<<" built "<<endl;
				newNameLeft[end-start]='\0';
				newLeft->name=newNameLeft;


				//	do right node
				start=lastComma.back()+1;
				end=i;

				for (j=start;j<i;j++){
						if (wholeFile[j]==':'){
							end=j;
							break;}
				}

				if (hasPrefix){
					for (j=lastComma.back()+1;j<i;j++){
							if (wholeFile[j]=='_'||wholeFile[j]=='-'){
								end=j;
								//start=j+1;
								break;}
					}

					char *SpeciesID = new char[i-end];
					for (j=end+1;j<i;j++){
							SpeciesID[j-end-1]=wholeFile[j];
						//cout<<SpeciesID[j-end-1];
						}
					//cout<<endl;
					SpeciesID[i-end-1]='\0';
					newRight->mappedSpecies=SpeciesID;
				}

				if (hasSuffix){
					for (j=i-1;j>=start;j--){
							if (wholeFile[j]=='_'){
								end=j;
								break;}
					}
					char *GeneID = new char[i-end];
					for (j=end+1;j<i;j++){
							GeneID[j-end-1]=wholeFile[j];
					}
					GeneID[i-end-1]='\0';
					int NewID;
					sscanf(GeneID,"%d",&NewID);
					newRight->GeneID=GeneID;
					newTree->mappedTrees.push_back(GeneID);
				}


				nameLengthLeft=end-start;
				char *newNameRight = new char[nameLengthRight+1];
				for (j=start;j<end;j++)
				{
					newNameRight[j-start]=wholeFile[j];
					//cout<<wholeFile[j];
				}
				//cout<<" built ";
				newNameRight[end-start]='\0';
				newRight->name=newNameRight;
				//cout<<newRight->name<<endl;
				lastComma.pop_back();
				lastLeft.pop_back();

			}
			else if (indirection==1)
			{
				node *newLeft = new node();
				newLeft->isleft = 1;
				newLeft->isleaf=1;
				newLeft->isroot=0;
				newLeft->parent=newNode;
				newNode->isleaf=0;
				newNode->left=newLeft;
				newNode->right=savedRightNode.back();
				savedRightNode.back()->parent=newNode;
				savedRightNode.pop_back();
				newTree->leafnodes.push_back(newLeft);

				if (outdirection==0)
				{
					newNode->isleft=1;
					savedLeftNode.push_back(newNode);
				}
				else if (outdirection==1)
				{
					newNode->isleft=0;
					savedRightNode.push_back(newNode);
				}
				start=lastLeft.back()+1;
				end=lastComma.back();

				for (j=start;j<lastComma.back();j++){
						if (wholeFile[j]==':'){
							end=j;
							break;}
				}

				if (hasPrefix){
					for (j=lastLeft.back()+1;j<lastComma.back();j++){
							if (wholeFile[j]=='_'||wholeFile[j]=='-'){
								end=j;
								//start=j+1;
								break;}
					}
					char *SpeciesID = new char[lastComma.back()-end];
					for (j=end+1;j<lastComma.back();j++){
							SpeciesID[j-end-1]=wholeFile[j];
							//cout<<SpeciesID[j-end-1];
					}
					//cout<<endl;
					SpeciesID[lastComma.back()-end-1]='\0';
					newLeft->mappedSpecies=SpeciesID;
				}

				
				if (hasSuffix){
					for (j=lastComma.back()-1;j>=start;j--){
							if (wholeFile[j]=='_'){
								end=j;
								break;}
					}
					char *GeneID = new char[lastComma.back()-end];
					for (j=end+1;j<lastComma.back();j++){
							GeneID[j-end-1]=wholeFile[j];
					}
					GeneID[lastComma.back()-end-1]='\0';
					int NewID;
					sscanf(GeneID,"%d",&NewID);
					newLeft->GeneID=GeneID;
					newTree->mappedTrees.push_back(GeneID);
				}

				nameLengthLeft=end-start;
				char *newNameLeft = new char[nameLengthLeft+1];
				
				for (j=start;j<end;j++)
				{
					newNameLeft[j-start]=wholeFile[j];
					//cout<<wholeFile[j];
				}
				//cout<<" built "<<endl;

				newNameLeft[end-start]='\0';
				
				newLeft->name=newNameLeft;

				lastComma.pop_back();
				lastLeft.pop_back();
			}
			else if (indirection==2)
			{
				node *newRight = new node();
				newRight->isleft=0;
				newRight->isleaf=1;
				newRight->isroot=0;
				newRight->parent=newNode;
				newNode->isleaf=0;
				newNode->right=newRight;
				newNode->left=savedLeftNode.back();
				savedLeftNode.back()->parent=newNode;
				savedLeftNode.pop_back();
				newTree->leafnodes.push_back(newRight);

				if (outdirection==0)
				{
					newNode->isleft=1;
					savedLeftNode.push_back(newNode);
				}
				else if (outdirection==1)
				{
					newNode->isleft=0;
					savedRightNode.push_back(newNode);
				}

				start=lastComma.back()+1;
				end=i;

				for (j=start;j<i;j++){
						if (wholeFile[j]==':'){
							end=j;
							break;}
				}

				if (hasPrefix){
					for (j=lastComma.back()+1;j<i;j++){
							if (wholeFile[j]=='_'||wholeFile[j]=='-'){
								end=j;
								//start=j+1;
								break;}
					}

					char *SpeciesID = new char[i-end];
					for (j=end+1;j<i;j++){
							SpeciesID[j-end-1]=wholeFile[j];
						//cout<<SpeciesID[j-end-1];
						}
					//cout<<endl;
					SpeciesID[i-end-1]='\0';
					newRight->mappedSpecies=SpeciesID;
				}

				if (hasSuffix){
					for (j=i-1;j>=start;j--){
							if (wholeFile[j]=='_'){
								end=j;
								break;}
					}
					char *GeneID = new char[i-end];
					for (j=end+1;j<i;j++){
							GeneID[j-end-1]=wholeFile[j];
					}
					GeneID[i-end-1]='\0';
					int NewID;
					sscanf(GeneID,"%d",&NewID);
					newRight->GeneID=GeneID;
					newTree->mappedTrees.push_back(GeneID);
				}

				nameLengthLeft=end-start;
				char *newNameRight = new char[nameLengthRight+1];
				for (j=start;j<end;j++)
				{
					newNameRight[j-start]=wholeFile[j];
					//cout<<wholeFile[j];
				}
				//cout<<" built "<<endl;
				newNameRight[end-start]='\0';
				newRight->name=newNameRight;

				lastComma.pop_back();
				lastLeft.pop_back();
			}
			else if (indirection==3)
			{
				newNode->left=savedLeftNode.back();
				newNode->right=savedRightNode.back();
				savedLeftNode.back()->parent=newNode;
				savedRightNode.back()->parent=newNode;
				savedLeftNode.pop_back();
				savedRightNode.pop_back();
				newNode->isleaf=0;
				if (outdirection==0)
				{
					newNode->isleft=1;
					savedLeftNode.push_back(newNode);
				}
				else if (outdirection==1)
				{
					newNode->isleft=0;
					savedRightNode.push_back(newNode);
				}
				lastComma.pop_back();
				lastLeft.pop_back();
			}
		}
	}
	newTree->root->isroot=1;
	if(hasSuffix)
	{
		sort(newTree->mappedTrees.begin(),newTree->mappedTrees.end());
		vector<string> ::iterator ix=unique(newTree->mappedTrees.begin(),newTree->mappedTrees.end());
		newTree->mappedTrees.erase(ix,newTree->mappedTrees.end());
	}
	/*
	int Keystart=0;
	setkey(newTree,newTree->root,&Keystart);
	Keystart=0;
	setLocalkey(newTree,newTree->root,&Keystart);
	*/
	if(hasSuffix==1)
		newTree->type=0;
	else if(hasPrefix==1)
		newTree->type=1;
	else 
		newTree->type=2;
	return newTree;
}

#endif
