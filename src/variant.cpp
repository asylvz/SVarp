#include <iostream>
#include <string>
#include <cstring>
#include <algorithm>
#include <filesystem>
#include "variant.h"



// Compares two intervals according to ending times in descending order.
inline bool cmp(Variant* i1, Variant* i2)
{
	if (i1->pos_in_node != i2->pos_in_node)
    	return (i1->pos_in_node < i2->pos_in_node);
	else
		return (i1->pos_in_node_end > i2->pos_in_node_end);
}


int arrange_variants(std::map<std::string, Variant*>& vars, std::map<std::string, std::vector<Variant*>>& vars_by_node)
{
	std::map<std::string, Variant*>::iterator itr;
	std::map<std::string, std::vector<Variant*>>::iterator it;
	
	for (itr=vars.begin(); itr != vars.end(); ++itr)
	{
		//std::cout<<itr->first<<"\n";
		int pos = (itr->first).find(':');
		std::string node_name = itr->first.substr(0, pos);

		it = vars_by_node.find(node_name);
		if (it != vars_by_node.end())
			it->second.push_back(itr->second);
		else
		{
			std::vector<Variant*> v;
			v.push_back(itr->second);
			vars_by_node.insert(std::pair<std::string, std::vector<Variant*>>(node_name, v));
		}
	}

	//std::cout<<"sorting.......................\n";
	//Now sort the vector of variants for each contig

	for (it=vars_by_node.begin(); it != vars_by_node.end(); ++it)
		std::sort((it->second).begin(), (it->second).end(), cmp);
	
	/*int cnt = 0;
	for (it=vars_by_node.begin(); it != vars_by_node.end(); ++it)
	{
		cnt++;
		std::cout<<"\n\n"<<it->first<<"\n";
		for (auto &sv : it->second) 
			std::cout<<cnt <<"\t"<< sv->pos_in_node<<"  ";
	}
	exit(1);*/
	
	return RETURN_SUCCESS;
}


int merge_svs_within_node(parameters& params, std::map<std::string, gfaNode*>& gfa, std::map<std::string, std::vector<Variant*>>::iterator& vars_by_node ,std::vector<Svtig*>& var_vector)
{
	Svtig* svtig_tmp;
	int start_pos = -1000;
	bool first = true;
	
	std::set <std::string> sets_intersect;

	for (auto &sv : vars_by_node->second) 
	{
		if(first)
		{
			svtig_tmp = new Svtig();
			svtig_tmp->node = vars_by_node->first;
			svtig_tmp->phased = false;
			svtig_tmp->contig = gfa[vars_by_node->first]->contig;
			svtig_tmp->ref_pos = sv->pos_in_ref;			
			start_pos = sv->pos_in_node;	
			svtig_tmp->start_pos = sv->pos_in_node;
			(svtig_tmp->reads_untagged).insert((sv->reads_untagged).begin(), (sv->reads_untagged).end());

			first = false;
		}
		else
		{
			sets_intersect.clear();
			std::set_intersection((svtig_tmp->reads_untagged).begin(), (svtig_tmp->reads_untagged).end(), (sv->reads_untagged).begin(), (sv->reads_untagged).end(), std::inserter(sets_intersect, sets_intersect.end()));
			
			if(start_pos + params.dist_threshold >= sv->pos_in_node)
			{
				(svtig_tmp->reads_untagged).insert((sv->reads_untagged).begin(), (sv->reads_untagged).end());
				start_pos = sv->pos_in_node;
			}
			else if(sets_intersect.size() == sv->reads_untagged.size())
				start_pos = sv->pos_in_node;
			else
			{
				if(svtig_tmp->reads_untagged.size() > 0)
					var_vector.push_back(svtig_tmp);
				else
					delete svtig_tmp;
				
				svtig_tmp = new Svtig();
				svtig_tmp->node = vars_by_node->first;
				svtig_tmp->phased = false;
				svtig_tmp->contig = gfa[vars_by_node->first]->contig;
				svtig_tmp->ref_pos = sv->pos_in_ref;
			
				start_pos = sv->pos_in_node;	
				svtig_tmp->start_pos = sv->pos_in_node;
				(svtig_tmp->reads_untagged).insert((sv->reads_untagged).begin(), (sv->reads_untagged).end());
				first = false;
			}
		}	
	}

	//Add the last one as well
	if(!first)
		var_vector.push_back(svtig_tmp);
	else
		delete svtig_tmp;

	return RETURN_SUCCESS;
}


int merge_neighbor_nodes(parameters& params, std::map<std::string, gfaNode*>& gfa, std::map<std::string, std::vector<Svtig*>> init_svtigs, std::map <std::string, std::vector<std::string>>& incoming, std::map <std::string, std::vector<std::string>>& outgoing)
{

	std::map<std::string, std::vector<std::string>>::iterator it_nodes;
	std::map<std::string, std::vector<Svtig*>>::iterator it_neighbors;
	for (auto &nd: init_svtigs)
	{
		if (nd.second.empty())
		{
			std::cout<<"Empty...\n";
			continue;
		}
		//For each node, get the sv at the front and at the back of the vector.
		//This can be done since they are sorted based on the position within the node
		//These may need to be merged with the SVs in the neighboring nodes
		Svtig* svtig_front = nd.second.front();
		//Svtig* svtig_back = nd.second.back();
		
		//First check the svtigs at the front of the vector of svs for this node
		if (svtig_front->start_pos < params.dist_threshold)
		{
			it_nodes = incoming.find(nd.first);
			if (it_nodes != incoming.end())
			{
				for(auto &incoming_node: it_nodes->second)
				{
					it_neighbors = init_svtigs.find(incoming_node);
					if(it_neighbors != init_svtigs.end())
					{
						Svtig* svtig_incoming = it_neighbors->second.back();
						
						int node_len = gfa[incoming_node]->len;
						if ((node_len - svtig_incoming->start_pos + svtig_front->start_pos) < params.dist_threshold)
						{
							(svtig_front->reads_untagged).insert((svtig_incoming->reads_untagged).begin(), (svtig_incoming->reads_untagged).end());	
							//std::string svtig_name = nd.first + ":" + std::to_string(svtig_front->start_pos);
							//svtig_front->name = svtig_name;
							svtig_incoming->filter = true;
							//std::cout<<"(FRONT) Adding "<<svtig_incoming->start_pos <<" (node-len = "<<node_len<< ") to "<< svtig_name<<"\n";
							//std::cout<<"Adding "<< svtig_incoming->reads_untagged.size() <<" reads to "<< svtig_front->reads_untagged.size() <<"\n\n";
						}
					}
				}
			}
		}
		//Omitted this because if there is an incoming svtig X to Y, we don't need to check Y to Z because when we process Z, we still merge Y with Z
		//Only missing check is the first and last nodes. To be added...

		/*
		int node_len = gfa[nd.first]->len;
		if (svtig_back->start_pos > node_len - MIN_SV_DISTANCE)
		{
			it_nodes = outgoing.find(nd.first);
			if (it_nodes != outgoing.end())
			{
				for(auto &outgoing_node: it_nodes->second)
				{
					it_neighbors = init_svtigs.find(outgoing_node);
					if(it_neighbors != init_svtigs.end())
					{
						Svtig* svtig_outgoing = it_neighbors->second.front();
						//If we added the reads of this svtig to the next one in the above step, skip.
						if (svtig_outgoing->filter == true)
							continue;
					
						if ((node_len + svtig_outgoing->start_pos - svtig_back->start_pos) < MIN_SV_DISTANCE)
						{
							(svtig_back->reads_untagged).insert((svtig_outgoing->reads_untagged).begin(), (svtig_outgoing->reads_untagged).end());	
							std::string svtig_added = outgoing_node + ":" + std::to_string(svtig_outgoing->start_pos);
							std::string svtig_name = nd.first + ":" + std::to_string(svtig_back->start_pos);
							svtig_front->name = svtig_name;
							svtig_outgoing->filter = true;
							//std::cout<<"(BACK) Adding "<<svtig_added<< " to "<< svtig_name<<" (node-len = "<<node_len<<")\n";
							//std::cout<<"Adding "<< svtig_outgoing->reads_h1.size() <<" reads to "<< svtig_back->reads_h1.size() <<"\n\n";
							std::cout<<"(BACK) Adding "<<svtig_outgoing->start_pos <<" (node-len = "<<node_len<< ") to "<< svtig_name<<"\n";
							std::cout<<"Adding "<< svtig_outgoing->reads_untagged.size() <<" reads to "<< svtig_back->reads_untagged.size() <<"\n\n";
						}
					}
				}
			}
		}*/
	}

	return RETURN_SUCCESS;
}


int find_final_svtigs(parameters& params, std::map<std::string, std::vector<Svtig*>>& init_svtigs, std::map<std::string, std::vector<Svtig*>>& final_svtigs, int& svtig_cnt_rp_filtered)
{	
	std::vector<Svtig*> var_vector;
	for (auto &node: init_svtigs)
	{
		if (node.second.empty())
			continue;
		
		var_vector.clear();
		for(auto &svtig: node.second)
		{
			if (svtig->filter)
				continue;
			
			if (svtig->reads_untagged.size() > static_cast<unsigned int>(params.support))
				var_vector.push_back(svtig);	
		}
		if (!var_vector.empty())
		{
			final_svtigs.insert(std::pair<std::string, std::vector<Svtig*>>(node.first, var_vector));
			svtig_cnt_rp_filtered += var_vector.size();
		}
	}

	return RETURN_SUCCESS;
}	



//Put the variants into a map of vectors with contig name as the key
//We don't do this while iterating the gaf file in alignment.cpp because we want to 
//find the sv later in O(1) using "contig_name:start_end" in order to add reads to the
//read set of the variant
int refine_svs(parameters& params, std::map<std::string, gfaNode*>& gfa, std::map<std::string, Variant*>& vars, std::map<std::string, std::vector<Svtig*>>& final_svtigs, std::map <std::string, std::vector<std::string>>& incoming, std::map <std::string, std::vector<std::string>>& outgoing)
{	
	std::map<std::string, std::vector<Variant*>>::iterator it;
	std::map<std::string, std::vector<Variant*>> vars_by_node;
	std::map<std::string, std::vector<Svtig*>> init_svtigs;

	std::vector<Svtig*> var_vector;

	arrange_variants(vars, vars_by_node);

	std::cout<<"\nMerging nearby SVs"<<std::endl;	
		
	int svtig_cnt = 0, svtig_cnt_rp_filtered = 0;
	for (it=vars_by_node.begin(); it != vars_by_node.end(); ++it)
	{
		//std::cout<<it->first<< " " <<(it->second).size() <<std::endl;
		var_vector.clear();
		
		//If the distance between any two SV is smaller than MIN_SV_DISTANCE, 
		//we merge the reads of these SVs
		merge_svs_within_node(params, gfa, it, var_vector);
		if (!var_vector.empty())
		{
			init_svtigs.insert(std::pair<std::string, std::vector<Svtig*>>(it->first, var_vector));
			svtig_cnt += var_vector.size();
		}
	}

	//merge inter nodes	
	merge_neighbor_nodes(params, gfa, init_svtigs, incoming, outgoing);
	find_final_svtigs(params, init_svtigs, final_svtigs, svtig_cnt_rp_filtered);
		
	/*for(auto &a: final_svtigs)
	{
		//std::cout<<a.first<<" "<< a.second.size() <<"\n";
		for (auto &tmp: a.second)
		{
			//std::cout<<tmp->name<<"\n";	
			if (tmp->node == "s439952" && tmp->start_pos == 359)
				std::cout<<"REFINE ---   "<<a.first<<" "<<tmp->contig<< " " <<tmp->ref_pos<<"\n";
			//std::cout<<tmp->node<< " " <<tmp->start_pos << "\t"<< tmp->contig<<"\n";
		}
		//std::cout<<"\n\n";
	}*/
	
	std::cout<<"--->there are "<<svtig_cnt<<" SVtigs after merging and "<<svtig_cnt_rp_filtered<<" after filtering based on minimum read support\n\n";

	return RETURN_SUCCESS;
}


int mapping_start_end(std::map<std::string, gfaNode*>& gfa, Gaf& line, std::map<std::string, Variant*>& variations_inter)
{
	bool skip_start = false, skip_end = false;
	char *path_copy = strdup(line.path.c_str());
	char *mytoken = strtok(path_copy,"><");
	int inserted_var_cnt = 0;

	//breakpoint1 is the position of the reads starting at that loci
	int node_count = 0, offset = 0, br1_start, br2_end, node_map_size = 0;
	
	if (line.query_start < MIN_READ_START_END_WINDOW)
		skip_start = true;
	if ((line.query_length - line.query_end) < MIN_READ_START_END_WINDOW)
		skip_end = true;
	if (skip_start && skip_end)
		return 0;
			
	std::string current_node, start_node, end_node;
	
	while(mytoken) 
	{
		current_node = mytoken;
		node_count++;
		char strand = line.path[offset];
		offset += strlen(mytoken) + 1;
		mytoken = strtok(NULL, "><");
		
		//std::cout<<"mytoken= "<<mytoken<<" current_token= "<< current_node<<"\n";

		if ((node_count == 1) && (!mytoken)) //means there is only a single node
		{
			if (!skip_start)
			{
				if(strand == '>')
					br1_start = line.path_start;
				else
					br1_start = gfa[current_node]->len - line.path_end;
				start_node = current_node;
			}
			
			if (!skip_end)
			{
				if(strand == '>')
					br2_end = line.path_end;
				else
					br2_end = gfa[current_node]->len - line.path_start;
				end_node = current_node;
			}
		}
		else if((node_count == 1) && (mytoken)) //First node
		{
			if (!skip_start)
			{
				if(strand == '>')
					br1_start = line.path_start;
				else
					br1_start = gfa[current_node]->len - line.path_start;
				start_node = current_node; 
			}

			node_map_size = gfa[current_node]->len;
		}
		else if(!mytoken) //Last node
		{
			if (!skip_end)
			{
				if(strand == '>')
					br2_end = line.path_end - node_map_size;
				else
					br2_end = gfa[current_node]->len - (line.path_end - node_map_size);
				end_node = current_node;
			}
		}
		else
			node_map_size += gfa[current_node]->len;
	}
	
	if (!skip_start)
	{
		std::string var_name = start_node + ":" + std::to_string(br1_start);
		std::map<std::string, Variant*>::iterator it = variations_inter.find(var_name);
						
		if (it != variations_inter.end())
			it->second->reads_untagged.insert(line.query_name);	
		else
		{
			Variant *v = new Variant();
			v->reads_untagged.insert(line.query_name);
			v->path = line.path;
			v->pos_in_node = br1_start;
			v->contig = gfa[start_node]->contig;
			v->node = start_node;
			v->pos_in_ref = gfa[v->node]->offset + v->pos_in_node;
			v->type = INTER;
			variations_inter.insert(std::pair<std::string, Variant*>(var_name, v));
			inserted_var_cnt++;
			//if (tokens[0] == "637c403e-43ba-4c67-9918-238263625f39")
			//	std::cout<<"HEREEE  INTRAA "<<v->pos_in_node<<"\n";
			if (v->node == "s439952" && v->pos_in_node == 359)
				std::cout<<"HEREEE   "<<gfa[v->node]->offset<<"\n";
		}

	}
	if (!skip_end)
	{
		std::string var_name = end_node + ":" + std::to_string(br2_end);
		
		std::map<std::string, Variant*>::iterator it = variations_inter.find(var_name);
						
		if (it != variations_inter.end())
			it->second->reads_untagged.insert(line.query_name);
		else
		{
			Variant *v = new Variant();
			v->reads_untagged.insert(line.query_name);
			v->path = line.path;
			v->pos_in_node = br2_end;
			v->contig = gfa[end_node]->contig;
			v->node = end_node;
			v->pos_in_ref = gfa[v->node]->offset + v->pos_in_node;
			v->type = INTER;	
			variations_inter.insert(std::pair<std::string, Variant*>(var_name, v));
			inserted_var_cnt++;
			if (v->node == "s439952" && v->pos_in_node == 359)
				std::cout<<"HEREEE   "<<gfa[v->node]->offset<<"\n";
		}
	}

	free(path_copy);
	return inserted_var_cnt;
}


Variant* generate_sv_node(std::map<std::string, gfaNode*>& gfa, Gaf& line, const int base_pos, int var_len, char sv_type)
{
	Variant *v = new Variant();

	char *path_copy = strdup(line.path.c_str());
	char *mytoken = strtok(path_copy,"><");
	
	int path_start = line.path_start;
	int offset = 0, node_count = 0, total_so_far = 0, node_map_size;
	v->phased = false;
	v->type = INTRA;	
	int pos_in_cigar = base_pos;
	bool del_incomplete = false;	

	while(mytoken) 
	{
		v->node = mytoken;
		char strand = line.path[offset];
		//v->path += strand + mytoken;
		node_count++;
		offset += strlen(mytoken) + 1;
		mytoken = strtok(NULL, "><");

		if ((node_count == 1) && (!mytoken)) //means there is only a single node
		{
			int pos_in_node = path_start + pos_in_cigar;
			if(strand == '>')
				v->pos_in_node = pos_in_node; //add +1 after offset because the sequence starts at offset +1
			else
			{
				if (sv_type == DELETION)
					v->pos_in_node = (gfa[v->node]->len - pos_in_node) - var_len;
				else
					v->pos_in_node = gfa[v->node]->len - pos_in_node;
			}
		
			v->pos_in_ref = gfa[v->node]->offset + v->pos_in_node;
			v->contig = gfa[v->node]->contig;
			v->node_strand = strand;
			
			if (sv_type == INSERTION)
			{
				v->pos_in_node_end = v->pos_in_node + 0;
				v->sv_type = INSERTION;			 
			}
			else if (sv_type == DELETION)
			{
				v->pos_in_node_end = v->pos_in_node + var_len;
				v->sv_type = DELETION;			 
			}
			free(path_copy);
			return v;	
		}
		else if((node_count == 1) && (mytoken)) //First node
		{
			node_map_size = gfa[v->node]->len - path_start;
			
			if ((sv_type != DELETION && pos_in_cigar < node_map_size) || (strand == '>' && pos_in_cigar < node_map_size) || (strand == '<' && sv_type == DELETION && pos_in_cigar + var_len < node_map_size))
			{
				int pos_in_node = path_start + pos_in_cigar;
				
				if(strand == '>')
					v->pos_in_node = pos_in_node;
				else
				{
					if (sv_type == DELETION)
						v->pos_in_node = (gfa[v->node]->len - pos_in_node) - var_len;
					else
						v->pos_in_node = gfa[v->node]->len - pos_in_node;
				}
				
				if (v->pos_in_node < 0)
				{
					std::cout<<"ERROR (first node) in addition of SV signal loci "<<v->node <<"\n";
					std::cout<<"v->pos_in_node = "<< v->pos_in_node<< " node len = "<<gfa[v->node]->len << " path_start = "<<path_start<< " pos_in_cigar = "<<pos_in_cigar<<" var_len = "<<var_len<<" sv_type = "<<sv_type<< " pos_in_node = "<<pos_in_node <<"\n";
					std::cout <<line.path<<"\t" <<line.path_length<<"\t"<<line.path_start<<"\t"<<line.path_end<<"\n";
				}

				v->pos_in_ref = gfa[v->node]->offset + v->pos_in_node;
				v->contig = gfa[v->node]->contig;
				v->node_strand = strand;
				if (sv_type == INSERTION)
				{
					v->pos_in_node_end = v->pos_in_node + 0;
					v->sv_type = INSERTION;			 
				}
				else if (sv_type == DELETION)
				{
					v->pos_in_node_end = v->pos_in_node + var_len;
					v->sv_type = DELETION;			 
				}

				free(path_copy);
				return v;	
			}
			else
			{	
				total_so_far = node_map_size;
				if (strand == '<' && sv_type == DELETION && pos_in_cigar < node_map_size)
				{
					pos_in_cigar += var_len;
					del_incomplete = true;
				}
			}
		}
		else if(!mytoken) //Last node
		{
			if ((sv_type != DELETION && pos_in_cigar >= node_map_size) || (del_incomplete && pos_in_cigar >= node_map_size) || (strand == '>' && pos_in_cigar >= node_map_size) || (strand == '<' && sv_type == DELETION && pos_in_cigar + var_len >= node_map_size && del_incomplete == false))
			{
				int pos_in_node = pos_in_cigar - total_so_far;
				
				if(strand == '>')
					v->pos_in_node = pos_in_node;
				else
				{
					if (sv_type == DELETION)
						if (del_incomplete)
							v->pos_in_node = gfa[v->node]->len - pos_in_node;
						else
							v->pos_in_node = (gfa[v->node]->len - pos_in_node) - var_len;
					else
						v->pos_in_node = gfa[v->node]->len - pos_in_node;
				}
				if (v->pos_in_node < 0)
				{
					std::cout<<"ERROR (last node) in addition of SV signal loci "<<v->node <<"\n";
					std::cout<<"v->pos_in_node = "<< v->pos_in_node<< "node len = "<<gfa[v->node]->len << " path_start = "<<path_start<< " pos_in_cigar = "<<pos_in_cigar<<" var_len = "<<var_len<<" sv_type = "<<sv_type<< " pos_in_node = "<<pos_in_node <<"\n";
					std::cout <<line.path<<"\t" <<line.path_length<<"\t"<<line.path_start<<"\t"<<line.path_end<<"\n";
				}

				v->pos_in_ref = gfa[v->node]->offset + v->pos_in_node;
				v->contig = gfa[v->node]->contig;
				v->node_strand = strand;
				if (sv_type == INSERTION)
				{
					v->pos_in_node_end = v->pos_in_node + 0;
					v->sv_type = INSERTION;			 
				}
				else if (sv_type == DELETION)
				{
					v->pos_in_node_end = v->pos_in_node + var_len;
					v->sv_type = DELETION;			 
				}

				free(path_copy);
				return v;	
			}
			else
				std::cout<<"Size problem in adding SV"<<std::endl;
		}
		else //middle node
		{
			node_map_size = total_so_far + gfa[v->node]->len;
			
			if ((sv_type != DELETION && pos_in_cigar < node_map_size) || (del_incomplete && pos_in_cigar < node_map_size) || (strand == '>' && pos_in_cigar < node_map_size) || (strand == '<' && sv_type == DELETION && pos_in_cigar + var_len < node_map_size && del_incomplete == false))
			{
				int pos_in_node = pos_in_cigar - total_so_far;
				
				if(strand == '>')
					v->pos_in_node = pos_in_node;
				else
				{
					if (sv_type == DELETION)
					{
						if (del_incomplete)
							v->pos_in_node = gfa[v->node]->len - pos_in_node;
						else
							v->pos_in_node = (gfa[v->node]->len - pos_in_node) - var_len;
					}
					else
						v->pos_in_node = gfa[v->node]->len - pos_in_node;
				}
	
				if (v->pos_in_node < 0)
				{
					std::cout<<"ERRRORRR middle "<<v->node <<"\n";
					std::cout<< v->pos_in_node<< "node len = "<<gfa[v->node]->len << " base pos = "<<base_pos<<" pos_in_cigar "<< pos_in_cigar<<" total_so_far = "<<total_so_far<<" var_len = "<<var_len<< "\n";
					std::cout <<line.path<<"\t" <<line.path_length<<"\t"<<line.path_start<<"\t"<<line.path_end<<"\n";
				}

				v->pos_in_ref = gfa[v->node]->offset + v->pos_in_node;
				v->contig = gfa[v->node]->contig;
				v->node_strand = strand;
				if (sv_type == INSERTION)
				{
					v->pos_in_node_end = v->pos_in_node + 0;
					v->sv_type = INSERTION;			 
				}
				else if (sv_type == DELETION)
				{
					v->pos_in_node_end = v->pos_in_node + var_len;
					v->sv_type = DELETION;			 
				}

				free(path_copy);
				return v;	
			}
		 	else
			{
				total_so_far = node_map_size;
				if (strand == '<' && sv_type == DELETION && pos_in_cigar < node_map_size)
				{
					pos_in_cigar += var_len;
					del_incomplete = true;
				}
			}
		}
	}

	free(path_copy);
	return NULL;
}



