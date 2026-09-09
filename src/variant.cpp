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
		auto pos = (itr->first).find(':');
		if (pos == std::string::npos)
			continue;
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

	//Now sort the vector of variants for each contig
	for (it=vars_by_node.begin(); it != vars_by_node.end(); ++it)
		std::sort((it->second).begin(), (it->second).end(), cmp);
	
	return RETURN_SUCCESS;
}


int merge_svs_within_node(parameters& params, std::map<std::string, gfaNode*>& gfa, std::map<std::string, std::vector<Variant*>>::iterator& vars_by_node ,std::vector<SVCluster*>& var_vector)
{
	if (gfa.count(vars_by_node->first) == 0)
		return RETURN_ERROR;

	SVCluster* svtig_tmp = nullptr;
	int start_pos = -1000;
	bool first = true;

	std::set <std::string> sets_intersect;

	for (auto &sv : vars_by_node->second) 
	{
		if(first)
		{
			svtig_tmp = new SVCluster();
			svtig_tmp->node = vars_by_node->first;
			svtig_tmp->phased = false;
			svtig_tmp->contig = gfa[vars_by_node->first]->contig;
			svtig_tmp->rank = gfa[vars_by_node->first]->rank;
			svtig_tmp->ref_pos = sv->pos_in_ref;			
			start_pos = sv->pos_in_node;	
			svtig_tmp->start_pos = sv->pos_in_node;
			svtig_tmp->end_pos = std::max(sv->pos_in_node, sv->pos_in_node_end);
			(svtig_tmp->reads_untagged).insert((sv->reads_untagged).begin(), (sv->reads_untagged).end());

			first = false;
		}
		else
		{
			if(start_pos + params.dist_threshold >= sv->pos_in_node)
			{
				(svtig_tmp->reads_untagged).insert((sv->reads_untagged).begin(), (sv->reads_untagged).end());
				start_pos = sv->pos_in_node;
				svtig_tmp->end_pos = std::max(svtig_tmp->end_pos, std::max(sv->pos_in_node, sv->pos_in_node_end));
			}
			else
			{
				// Cluster by shared read identity, not only by position: an SV beyond
				// dist_threshold whose reads are all already in this cluster stays
				// merged, and start_pos advances so a chain of same-read signals keeps
				// clustering. The intersection is only needed here, so it is computed
				// on this branch rather than for every SV.
				sets_intersect.clear();
				std::set_intersection((svtig_tmp->reads_untagged).begin(), (svtig_tmp->reads_untagged).end(), (sv->reads_untagged).begin(), (sv->reads_untagged).end(), std::inserter(sets_intersect, sets_intersect.end()));

				if(sets_intersect.size() == sv->reads_untagged.size())
				{
					start_pos = sv->pos_in_node;
					svtig_tmp->end_pos = std::max(svtig_tmp->end_pos, std::max(sv->pos_in_node, sv->pos_in_node_end));
				}
				else
				{
					if(svtig_tmp->reads_untagged.size() > 0)
						var_vector.push_back(svtig_tmp);
					else
						delete svtig_tmp;

					svtig_tmp = new SVCluster();
					svtig_tmp->node = vars_by_node->first;
					svtig_tmp->phased = false;
					svtig_tmp->contig = gfa[vars_by_node->first]->contig;
					svtig_tmp->rank = gfa[vars_by_node->first]->rank;
					svtig_tmp->ref_pos = sv->pos_in_ref;

					start_pos = sv->pos_in_node;
					svtig_tmp->start_pos = sv->pos_in_node;
					svtig_tmp->end_pos = std::max(sv->pos_in_node, sv->pos_in_node_end);
					(svtig_tmp->reads_untagged).insert((sv->reads_untagged).begin(), (sv->reads_untagged).end());
					first = false;
				}
			}
		}
	}

	// Add the last one as well (if allocated)
	if (svtig_tmp != nullptr && !first)
	{
		if (svtig_tmp->reads_untagged.size() > 0)
			var_vector.push_back(svtig_tmp);
		else
			delete svtig_tmp;
	}

	return RETURN_SUCCESS;
}


//Merge clusters that lie within dist_threshold of each other across a link. The link
//orientations say which node ends meet, so the gap is measured from the right end of each
//cluster. The merged cluster is dropped from its own node.
int merge_neighbor_nodes(parameters& params, std::map<std::string, gfaNode*>& gfa, std::map<std::string, std::vector<SVCluster*>>& init_svtigs, EdgeMap& incoming, EdgeMap& /*outgoing*/)
{
	for (auto &nd : init_svtigs)
	{
		EdgeMap::iterator it_edges = incoming.find(nd.first);
		if (nd.second.empty() || it_edges == incoming.end() || gfa.count(nd.first) == 0)
			continue;
		int len_to = gfa[nd.first]->len;

		for (const Edge& e : it_edges->second)
		{
			if (e.node == nd.first)
				continue;
			auto it_from = init_svtigs.find(e.node);
			if (it_from == init_svtigs.end() || it_from->second.empty() || gfa.count(e.node) == 0)
				continue;

			//the cluster of each node closest to the shared end
			SVCluster* to = (e.to == '+') ? nd.second.front() : nd.second.back();
			SVCluster* from = (e.from == '+') ? it_from->second.back() : it_from->second.front();
			if (to->filter) //already merged away, anything added now would go with it
				continue;

			int gap_to = (e.to == '+') ? to->start_pos : len_to - to->end_pos;
			int gap_from = (e.from == '+') ? gfa[e.node]->len - from->end_pos : from->start_pos;
			if (std::max(gap_to, 0) + std::max(gap_from, 0) < params.dist_threshold)
			{
				to->reads_untagged.insert(from->reads_untagged.begin(), from->reads_untagged.end());
				from->filter = true;
			}
		}
	}
	return RETURN_SUCCESS;
}


int find_final_svtigs(parameters& params, std::map<std::string, std::vector<SVCluster*>>& init_svtigs, std::map<std::string, std::vector<SVCluster*>>& final_svtigs, int& svtig_cnt_rp_filtered)
{	
	std::vector<SVCluster*> var_vector;
	for (auto &node: init_svtigs)
	{
		if (node.second.empty())
			continue;
		
		var_vector.clear();
		for(auto &svtig: node.second)
		{
			// Clusters not promoted to final_svtigs (merged away or below support)
			// are owned nowhere else, so free them here.
			if (svtig->filter)
			{
				delete svtig;
				continue;
			}

			if (svtig->reads_untagged.size() > static_cast<unsigned int>(params.support))
				var_vector.push_back(svtig);
			else
				delete svtig;
		}
		if (!var_vector.empty())
		{
			final_svtigs.insert(std::pair<std::string, std::vector<SVCluster*>>(node.first, var_vector));
			svtig_cnt_rp_filtered += var_vector.size();
		}
	}

	return RETURN_SUCCESS;
}	



//Put the variants into a map of vectors with contig name as the key
//We don't do this while iterating the GAF file in alignment.cpp because we want to 
//find the SV O(1) using "contig_name:start_end" in order to add reads to the
//read set of the variant during the GAF processing
int merge_svs(parameters& params, std::map<std::string, gfaNode*>& gfa, std::map<std::string, Variant*>& vars, std::map<std::string, std::vector<SVCluster*>>& final_svtigs, EdgeMap& incoming, EdgeMap& outgoing)
{	
	std::map<std::string, std::vector<Variant*>>::iterator it;
	std::map<std::string, std::vector<Variant*>> vars_by_node;
	std::map<std::string, std::vector<SVCluster*>> init_svtigs;

	std::vector<SVCluster*> var_vector;
	
	//store variants in node-based map s.t. each node record has sorted variant positons in a vector
	arrange_variants(vars, vars_by_node);
	std::cout<<"\nMerging nearby SV signals"<<std::endl;	
		
	int svtig_cnt = 0, svtig_cnt_rp_filtered = 0;
	for (it=vars_by_node.begin(); it != vars_by_node.end(); ++it)
	{
		//std::cout<<it->first<< " " <<(it->second).size() <<std::endl;
		var_vector.clear();
		
		//If the distance between any two SVs <MIN_SV_DISTANCE, merge the reads of these SVs
		merge_svs_within_node(params, gfa, it, var_vector);
		if (!var_vector.empty())
		{
			init_svtigs.insert(std::pair<std::string, std::vector<SVCluster*>>(it->first, var_vector));
			svtig_cnt += var_vector.size();
		}
	}

	//merge inter nodes	
	merge_neighbor_nodes(params, gfa, init_svtigs, incoming, outgoing);
	find_final_svtigs(params, init_svtigs, final_svtigs, svtig_cnt_rp_filtered);
		
	std::cout<<"--> "<<svtig_cnt<<" read clusters (putative svtigs) after merging\n";
	std::cout<<"--> "<<svtig_cnt_rp_filtered<<" read clusters after filtering based on minimum read support\n\n";

	if (params.fp_logs.is_open()) {
		params.fp_logs << "--> " << svtig_cnt << " read clusters (putative svtigs) after merging\n";
		params.fp_logs << "--> " << svtig_cnt_rp_filtered << " read clusters after filtering based on minimum read support\n";
	}

	return RETURN_SUCCESS;
}


int mapping_start_end(std::map<std::string, gfaNode*>& gfa, Gaf& line, std::map<std::string, Variant*>& variations_inter, int min_end)
{
	bool skip_start = false, skip_end = false;
	const std::string &path = line.path;
	int inserted_var_cnt = 0;

	// breakpoint1 is the position of the reads starting at that loci
	int node_count = 0, br1_start = -1, br2_end = -1, node_map_size = 0;
	
	
	if (line.query_start < min_end)
		skip_start = true;
	if ((line.query_length - line.query_end) < min_end)
		skip_end = true;
	if (skip_start && skip_end)
		return 0;
			
	std::string current_node, start_node, end_node;
    
	// iterate through path string, reading strand and node name
	size_t p = 0;
	while (p < path.size())
	{
		char strand = path[p];
		++p;
		size_t q = p;
		while (q < path.size() && path[q] != '>' && path[q] != '<') ++q;
		current_node = path.substr(p, q - p);
		p = q;
		node_count++;

		if (gfa.count(current_node) == 0)
			continue;

		if ((node_count == 1) && (p == path.size())) //means there is only a single node
		{
			if (!skip_start)
			{
				if(strand == '>')
					br1_start = line.path_start;
				else
					br1_start = gfa[current_node]->len - line.path_start;
				start_node = current_node;
			}

			if (!skip_end)
			{
				if(strand == '>')
					br2_end = line.path_end;
				else
					br2_end = gfa[current_node]->len - line.path_end;
				end_node = current_node;
			}
		}
		else if((node_count == 1) && (p < path.size())) //First node
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
		else if(p == path.size()) //Last node
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
	
	// If start_node/end_node is empty the boundary node was absent from the
	// graph; skip it (gfa[""] would null-dereference).
	if (!skip_start && !start_node.empty())
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
		}

	}
	if (!skip_end && !end_node.empty())
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
		}
	}

	return inserted_var_cnt;
}


//Node and forward-strand coordinate of an intra-alignment SV; 0-based boundary index as in mapping_start_end
Variant* generate_sv_node(std::map<std::string, gfaNode*>& gfa, Gaf& line, const int base_pos, int var_len, char sv_type)
{
	const std::string &path = line.path;
	int p0 = line.path_start + base_pos - 1; //first deleted base, or the base after the insertion
	int p1 = (sv_type == DELETION) ? p0 + var_len - 1 : p0; //last deleted base
	if (p0 < 0)
		return nullptr;

	//Nodes holding p0 and p1
	std::string node0, node1;
	char strand0 = 0, strand1 = 0;
	int idx0 = -1, idx1 = -1, off0 = 0, off1 = 0, acc = 0, idx = 0;
	size_t p = 0;
	while (p < path.size() && idx1 < 0)
	{
		char strand = path[p++];
		size_t q = p;
		while (q < path.size() && path[q] != '>' && path[q] != '<') ++q;
		std::string node = path.substr(p, q - p);
		p = q;

		std::map<std::string, gfaNode*>::iterator it = gfa.find(node);
		if (it == gfa.end())
			return nullptr;
		int len = it->second->len;
		if (idx0 < 0 && p0 < acc + len)
		{
			node0 = node; strand0 = strand; idx0 = idx; off0 = p0 - acc;
		}
		if (idx0 >= 0 && p1 < acc + len)
		{
			node1 = node; strand1 = strand; idx1 = idx; off1 = p1 - acc;
		}
		acc += len;
		idx++;
	}
	if (idx1 < 0)
		return nullptr;

	Variant *v = new Variant();
	v->phased = false;
	v->type = INTRA;
	v->sv_type = sv_type;

	if (sv_type == INSERTION)
	{
		v->node = node0;
		v->node_strand = strand0;
		v->pos_in_node = (strand0 == '>') ? off0 : gfa[node0]->len - off0;
	}
	else if (idx0 == idx1)
	{
		v->node = node0;
		v->node_strand = strand0;
		v->pos_in_node = (strand0 == '>') ? off0 : gfa[node0]->len - 1 - off1;
	}
	else
	{
		//Deletion across nodes: report the node whose forward end lies inside it, ties by name
		bool first = (strand0 == '>'), last = (strand1 == '<');
		if (first == last)
			first = (node0 <= node1);
		if (first)
		{
			v->node = node0;
			v->node_strand = strand0;
			v->pos_in_node = (strand0 == '>') ? off0 : 0;
		}
		else
		{
			v->node = node1;
			v->node_strand = strand1;
			v->pos_in_node = (strand1 == '<') ? gfa[node1]->len - 1 - off1 : 0;
		}
	}
	v->pos_in_node_end = v->pos_in_node + ((sv_type == DELETION) ? var_len : 0);
	v->pos_in_ref = gfa[v->node]->offset + v->pos_in_node;
	v->contig = gfa[v->node]->contig;
	return v;
}



