#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
#include <iterator>
#include "gfa.h"
#include "common.h"

std::map<std::string, gfaNode*> read_gfa(parameters* params)
{	
	// Read the S lines in the GFA file
	std::cout << "Reading the GFA file"<< std::endl;
	
	std::string line;	
	std::vector <std::string> tokens;	
	std::ifstream fp(params->ref_graph);
	map<string, gfaNode*> ref;
	
	while(fp)
	{
		getline(fp, line);
			
		std::string tmp_str;
		std::stringstream s(line);
		tokens.clear();
		
		while(getline(s, tmp_str, '\t'))
        	tokens.push_back(tmp_str);
		
		if (tokens[0] != "S")
			continue;
		
		gfaNode *g = new gfaNode();
		g->name = tokens[1];
		g->sequence = tokens[2];
		g->len = stoi(tokens[3].substr(5));
		g->contig = tokens[4].substr(5);
		g->offset = stoi(tokens[5].substr(5));

		ref.insert(std::pair<std::string, gfaNode*>(g->name, g));
	}
	//cout <<ref["s10"]->contig<<endl;
	/*map<std::string, gfaNode*>::iterator itr;
    cout << "\nThe map is :"<<cnt<< "\n";
    cout << "\tKEY\tELEMENT\n";
    for (itr = ref.begin(); itr != ref.end(); ++itr) {
        cout << '\t' << itr->first << '\t' << itr->second->offset<< '\n';
    }*/
	
	return ref;
}

/*
void init_reference(reference** ref)
{
	*ref = (reference*) getMem(sizeof(reference));
	(*ref)->contig = NULL;
	(*ref)->sequence = NULL;
	(*ref)->len = -1;
	(*ref)->offset = -1;
}


void gfa_traverse(gfa** hashtab)
{
	gfa *ptr;
	for(int i = 0; i < HASHSIZE; i++)
	{
		for(ptr = hashtab[i]; ptr != NULL; ptr = ptr->next)
			printf("%s - %s\n", ptr->node->name, ptr->node->contig);
	}
}

gfa *gfa_lookup(gfa **hashtab, char *s)
{
    gfa *np;
	//printf("lookup %s - %u", s, hash(s));
    for (np = hashtab[hash(s)]; np != NULL; np = np->next)
	{
        if (strcmp(s, np->node->name) == 0)
        	return np;
	}
	return NULL;
}

gfa *gfa_insert(gfa **hashtab, char *node_name, reference *node)
{
    gfa *np;
    unsigned hashval;
	
	if ((np = gfa_lookup(hashtab, node_name)) == NULL)
	{
        np = (gfa*) getMem(sizeof(np));
        if (np == NULL)
        	return NULL;
		//np->node_name = strdup(node_name);
		np->node = node;
        hashval = hash(node_name); 
		np->next = hashtab[hashval];
        hashtab[hashval] = np;
    } 
	//else
    //    free((gfa *) np->node_name);
    if ((np->node = node) == NULL)
    	return NULL;
   
	//printf("\n");
   	return np;
}*/
