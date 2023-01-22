#include <iostream>
#include <string>
#include <cstring>
#include "sv.h"


variant* generate_sv_node(std::map<std::string, gfaNode*> gfa, std::string path, int path_start, int path_end, int ref_pos, int cigar_len, char sv_type)
{
	variant *v = new variant();
	char *path_copy = (char *) path.c_str();
	
	char *mytoken = strtok(path_copy,"><");
	
	int offset = 0, node_count = 0, total_so_far = 0;
	//int total_path_length = path_end - path_start;
	
	while(mytoken) 
	{
		//if (strcmp (mytoken,"s18891") == 0)
		//	cout<<"node count "<< node_count<<" - mytoken: "<<mytoken <<endl;
		
		v->node = mytoken;
		v->phased = false;
		char strand = path[offset];
		node_count += 1;
		offset += strlen(mytoken) + 1;
		mytoken = strtok(NULL, "><");
		

		if ((node_count == 1) && (!mytoken)) //means there is only a single node
		{

			//if(path_start == 20405)
			//	cout<<"Single node* mytoken "<< mytoken<<endl;

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
		else if((node_count == 1) && (mytoken)) //First node
		{
			int	node_map_size = gfa[v->node]->len - path_start;
			
			//if (v->node == "s18892")
			//	cout<<"First node* node_map_size "<< node_map_size <<" ref_pos="<<ref_pos<<endl;
			
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
			{
				total_so_far = node_map_size;	
				//if (path_start == 20405)
				//	cout<<"First node* total_so_far: "<< total_so_far<<endl;
				
			}
		}
		else if(!mytoken) //Last node
		{
			//if (v->node == "s18891")
			//	cout<<"Last node* total_so_far "<< total_so_far <<" ref_pos="<<ref_pos<<endl;
			
			if (ref_pos >= total_so_far)
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
				std::cout<<"Size problem in adding SV"<<std::endl;
		}
		else //middle node
		{
			int node_map_size = total_so_far + gfa[v->node]->len;
			
			//if (v->node == "s18891")
			//	cout<<"Middle node* node_map_size "<< node_map_size <<" ref_pos="<<ref_pos<<"mytoken= "<<mytoken <<endl;
			
			if (node_map_size >= ref_pos)
			{
				//	cout<<"contig offset: "<< gfa[v->node]->offset <<" ref_pos: "<<ref_pos<<" total_so_far: "<<total_so_far <<"diff: "<<ref_pos-total_so_far<<endl;
				
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
			{
				
				total_so_far = node_map_size;
			}
		}
	}
	std::cout<<"ERRROORRRRRR -"<< total_so_far<< " "<<ref_pos<<std::endl;
	return NULL;
}
