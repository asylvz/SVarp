#include <iostream>
#include <map>
#include "cmdline.h"
#include "gfa.h"
#include "reference.h" 
#include "alignment.h"
#include "sv.h"

//ctags *.c

int main(int argc, char** argv)
{
	std::vector <variant*> insertions;	
	
	parameters* params;
	init_params(&params);
	
	int return_value = parse_command_line(argc, argv, params);			
	printf("\nThe input files are: \n\t%s\n\t%s\n", params->gaf, params->ref_graph);
	

	std::map<std::string, gfaNode*> ref = read_gfa(params);
	
	int a = read_alignments(params, ref, insertions);	
	
	return 0;
}


	/*gfa *g = new gfa("s1","chr19","ACGT", 100, 39);
	g->display();
	
	gfa *g2 = new gfa("s2","chr1","ACGT", 100, 39);
	g2->display();
	
	List <gfa*> list3;	
	list3.add(g2);
	list3.add(g);
	gfa *t = list3.remove();
  	cout << t->contig;	
	*/
	
	/*gfa *gfa_table[HASHSIZE];
	for(int i = 0; i < HASHSIZE; i++)
		gfa_table[i] = NULL;
	
	gaf *gaf_table[HASHSIZE];
	for(int i = 0; i < HASHSIZE; i++)
		gaf_table[i] = NULL;
	*/
