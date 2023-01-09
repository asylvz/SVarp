#include <iostream>
#include <string>
#include <cstring>
#include "sv.h"

variant* generate_sv_node(std::map<std::string, gfaNode*> gfa, string path, int path_start, int path_end, int ref_pos, int cigar_len, char sv_type)
{
	variant *v = new variant();
	char *path_copy = (char *) path.c_str();
	
	char *mytoken = strtok(path_copy,"><");
	
	int offset = 0, node_count = 0, total_so_far = 0;
	int total_path_length = path_end - path_start;
	
	while(mytoken) 
	{
		v->node = mytoken;
		char strand = path[offset];
		node_count += 1;
		offset += strlen(mytoken) + 1;
		mytoken = strtok(NULL, "><");
		
		if ((node_count == 1) && (mytoken == NULL)) //means there is only a single node
		{
			if(strand == '>')
				v->ref_start = gfa[v->node]->offset + (path_start + ref_pos);
			else
				v->ref_start = gfa[v->node]->offset + (gfa[v->node]->len - (path_end + ref_pos));
			
			v->contig = gfa[v->node]->contig;
			v->node_strand = strand;
			if (sv_type == INSERTION)
			{
				v->ref_end = v->ref_start + 0;
				v->sv_type = INSERTION;			 
			}
			
			return v;	
		}
		else if((node_count == 1) && (mytoken != NULL)) //First node
		{
			int node_map_size = gfa[v->node]->len - path_start;
			if (node_map_size >= ref_pos)
			{
				if(strand == '>')
					v->ref_start = gfa[v->node]->offset + (path_start + ref_pos);
				else
					v->ref_start = gfa[v->node]->offset + (gfa[v->node]->len - (path_start + ref_pos));
			
				v->contig = gfa[v->node]->contig;
				v->node_strand = strand;
				if (sv_type == INSERTION)
				{
					v->ref_end = v->ref_start + 0;
					v->sv_type = INSERTION;			 
				}
				return v;	
			}
			else
				total_so_far += node_map_size;	
		}
		else if(mytoken == NULL) //Last node
		{
			if (total_path_length >= ref_pos)
			{
				if(strand == '>')
					v->ref_start = gfa[v->node]->offset + (ref_pos - total_so_far);
				else
					v->ref_start = gfa[v->node]->offset + (gfa[v->node]->len - (ref_pos - total_so_far));
			
				v->contig = gfa[v->node]->contig;
				v->node_strand = strand;
				if (sv_type == INSERTION)
				{
					v->ref_end = v->ref_start + 0;
					v->sv_type = INSERTION;			 
				}
				return v;	
			}
			else
				cout<<"Size problem in adding SV"<<endl;
		}
		else //middle node
		{
			int node_map_size = total_so_far + gfa[v->node]->len;
			
			if (node_map_size >= ref_pos)
			{
				if(strand == '>')
					v->ref_start = gfa[v->node]->offset + (ref_pos - total_so_far);
				else
					v->ref_start = gfa[v->node]->offset + (gfa[v->node]->len - (ref_pos - total_so_far));
			
				v->contig = gfa[v->node]->contig;
				v->node_strand = strand;
				if (sv_type == INSERTION)
				{
					v->ref_end = v->ref_start + 0;
					v->sv_type = INSERTION;			 
				}
				return v;	
			}
			else
				total_so_far += node_map_size;	
		}
		
	}
	cout<<"ERRROORRRRRR"<< total_path_length<<" "<< total_so_far<< " "<<ref_pos<<endl;
	return NULL;
}
