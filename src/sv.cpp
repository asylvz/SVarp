#include <iostream>
#include <string>
#include <vector>
#include <cstring>
#include <algorithm>
#include "sv.h"


// Compares two intervals according to ending times in descending order.
bool cmp(variant* i1, variant* i2)
{
	if (i1->ref_start != i2->ref_start)
    	return (i1->ref_start < i2->ref_start);
	else
		return (i1->ref_end > i2->ref_end);
}


std::map<std::string, std::vector<variant*>> arrange_variants(std::map<std::string, variant*>& initial_insertions)
{

	std::map<std::string, variant*>::iterator itr;
	
	std::map<std::string, std::vector<variant*>>::iterator it;
	std::map<std::string, std::vector<variant*>> ins;
	
	for (itr=initial_insertions.begin(); itr != initial_insertions.end(); ++itr)
	{
		int pos = (itr->first).find(':');
		std::string contig_name = itr->first.substr(0, pos);
	
		it = ins.find(contig_name);
		if (it != ins.end())
			it->second.push_back(itr->second);
		else
		{
			std::vector<variant*> v;
			v.push_back(itr->second);
			ins.insert(std::pair<std::string, std::vector<variant*>>(contig_name, v));
		}
	}

	//Now sort the vector of variants for each contig
	for (it=ins.begin(); it != ins.end(); ++it)
		std::sort((it->second).begin(), (it->second).end(), cmp);

	return ins;
}


//Put the variants into a map of vectors with contig name as the key
//We don't do this while iterating the gaf file in alignment.cpp because we want to 
//find the sv later in O(1) using "contig_name:start_end" in order to add reads to the
//read set of the variant
std::map<std::string, std::vector<svtig*>> refine_svs(std::map<std::string, variant*> initial_insertions)
{
	
	std::map<std::string, std::vector<variant*>>::iterator it;
	std::map<std::string, std::vector<variant*>> insertions;
	
	std::map<std::string, std::vector<svtig*>> var;
	std::vector<svtig*> var_vector;
	
	insertions = arrange_variants(initial_insertions);
	
	int merged_sv_count = 0;
	for (it=insertions.begin(); it != insertions.end(); ++it)
	{
		//std::cout<<it->first<< " " <<(it->second).size() <<std::endl;
		var_vector.clear();
		
		svtig* svtig_tmp = nullptr;
		int start_pos = -1000;
		bool first = true;
		
		//If the distance between any two sv is smaller than MIN_SV_DISTANCE, 
		//we merge the reads of these svs
		for (auto &sv : it->second) 
		{
			//std::cout<<"\t"<<sv->ref_start<<std::endl;
			if((!first) && (start_pos + MIN_SV_DISTANCE >= sv->ref_start))
			{
				//std::cout<<"\t-------------IN1 start_pos:"<<start_pos<<" - sv->ref_start"<<sv->ref_start<<std::endl;
				(svtig_tmp->reads_h1).insert((sv->reads_h1).begin(), (sv->reads_h1).end());
				start_pos = sv->ref_start;
				merged_sv_count++;
			}
			else
			{
				//if((!first) && (svtig_tmp->reads_h1.size() > MIN_READ_SUPPORT))
	
				if(!first)
					var_vector.push_back(svtig_tmp);
				//std::cout<<"\tIN2"<<std::endl;
				svtig_tmp = new svtig();
				svtig_tmp->sv_type = INSERTION;
				svtig_tmp->contig = it->first;
				svtig_tmp->phased = false;
				svtig_tmp->sv_size = sv->sv_size;
				start_pos = sv->ref_start;
				
				(svtig_tmp->reads_h1).insert((sv->reads_h1).begin(), (sv->reads_h1).end());
				first = false;
			}
		}
		//Add the last one as well
		if (!first)
			var_vector.push_back(svtig_tmp);
		
		//if (var_vector.size() >	MIN_READ_SUPPORT)
		var.insert(std::pair<std::string, std::vector<svtig*>>(it->first, var_vector));
	}
	
	std::cout<<"Merged "<<merged_sv_count<<" SVs\n\n";
	std::vector<int> read_size;
/*
	//std::cout<<"-------SVTIGS------\n";
	std::map<std::string, std::vector<svtig*>>::iterator it2;
	int cnt = 1;
	for (it2=var.begin(); it2 != var.end(); ++it2)
	{
		//std::cout<<it2->first;
		for (auto &sv : it2->second) 
		{
			std::cout<<"\n"<<it2->first<<" svtig "<<cnt++<<"\n";
			for (auto &read : sv->reads_h1)
        		std::cout << ' ' << read;
		}
        std::cout <<"\n\n";
	}
	std::cout<<"Merged "<<merged_sv_count<<" SVs\n\n";
*/	
	return var;
}


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



