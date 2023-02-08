#include <iostream>
#include <string>
#include <vector>
#include <cstring>
#include <algorithm>
#include <filesystem>
#include "sv.h"



int find_deletions(parameters* params, std::map<std::string, std::vector<svtig*>> deletions)
{

	std::cout<<"\nFinding final deletions"<<std::endl;	
	std::map<std::string, std::vector<svtig*>>::iterator itr;
	int read_support_filter = 0, final_sv = 0;

	std::string cwd = std::filesystem::current_path().string();
	std::string log_path = cwd + "/log/";	
	
	std::ofstream outFile(log_path + "deletions.bed");
	
	for (itr=deletions.begin(); itr != deletions.end(); ++itr)
	{
		//Generate fastq files
		for (auto &sv : itr->second) 
		{
			int read_count = sv->reads_h1.size() + sv->reads_h2.size();
			if(read_count < MIN_READ_SUPPORT)	
			{
				read_support_filter++;
				continue;
			}
			else
			{
				if (sv->phased == true)
					outFile<< itr->first<<"\t"<< sv->start_pos<<" "<< sv->end_pos <<" "<<sv->reads_h1.size()<<" "<<sv->reads_h2.size()<<" phased" <<std::endl;		
				else		
					outFile<< itr->first<<"\t"<< sv->start_pos<<" "<< sv->end_pos <<" "<<read_count <<std::endl;		

				final_sv++;
			}
		}
	}
	std::cout<<"--->there are "<<final_sv<<" deletions (written to "<<log_path + "deletions.bed)\n";

	return RETURN_SUCCESS;
}


// Compares two intervals according to ending times in descending order.
bool cmp(variant* i1, variant* i2)
{
	if (i1->ref_start != i2->ref_start)
    	return (i1->ref_start < i2->ref_start);
	else
		return (i1->ref_end > i2->ref_end);
}


int arrange_variants(std::map<std::string, variant*>& initial_vars, std::map<std::string, std::vector<variant*>>& final_ins, std::map<std::string, std::vector<variant*>>& final_dels)
{
	std::map<std::string, variant*>::iterator itr;
	std::map<std::string, std::vector<variant*>>::iterator it;
		
	for (itr=initial_vars.begin(); itr != initial_vars.end(); ++itr)
	{
		int pos = (itr->first).find(':');
		std::string contig_name = itr->first.substr(0, pos);
	
		if (itr->second->sv_type == INSERTION)
		{
			it = final_ins.find(contig_name);
			if (it != final_ins.end())
				it->second.push_back(itr->second);
			else
			{
				std::vector<variant*> v;
				v.push_back(itr->second);
				final_ins.insert(std::pair<std::string, std::vector<variant*>>(contig_name, v));
			}
		}
		else if(itr->second->sv_type == DELETION)
		{
			it = final_dels.find(contig_name);
			if (it != final_dels.end())
				it->second.push_back(itr->second);
			else
			{
				std::vector<variant*> v;
				v.push_back(itr->second);
				final_dels.insert(std::pair<std::string, std::vector<variant*>>(contig_name, v));
			}
		}
	}

	//Now sort the vector of variants for each contig
	for (it=final_ins.begin(); it != final_ins.end(); ++it)
		std::sort((it->second).begin(), (it->second).end(), cmp);

	for (it=final_dels.begin(); it != final_dels.end(); ++it)
		std::sort((it->second).begin(), (it->second).end(), cmp);
	
	//std::cout<<final_ins.size()<< " "<<final_dels.size()<<std::endl;	
	return RETURN_SUCCESS;
}


//Put the variants into a map of vectors with contig name as the key
//We don't do this while iterating the gaf file in alignment.cpp because we want to 
//find the sv later in O(1) using "contig_name:start_end" in order to add reads to the
//read set of the variant
int refine_svs(std::map<std::string, variant*> initial_variations, std::map<std::string, std::vector<svtig*>>& final_ins, std::map<std::string, std::vector<svtig*>>& final_del)
{

	std::map<std::string, std::vector<variant*>>::iterator it;
	std::map<std::string, std::vector<variant*>> insertions, deletions;
	
	std::vector<svtig*> var_vector;
	
	arrange_variants(initial_variations, insertions, deletions);
	
	std::cout<<"Merging nearby SVs"<<std::endl;	
	int merged_ins_count = 0, merged_del_count = 0;
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
				merged_ins_count++;
			}
			else
			{
				if((!first) && (svtig_tmp->reads_h1.size() > MIN_READ_SUPPORT))	
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
		if((!first) && (svtig_tmp->reads_h1.size() > MIN_READ_SUPPORT))	
			var_vector.push_back(svtig_tmp);
		
		//if (var_vector.size() >	MIN_READ_SUPPORT)
		final_ins.insert(std::pair<std::string, std::vector<svtig*>>(it->first, var_vector));
	}

	for (it=deletions.begin(); it != deletions.end(); ++it)
	{
		//std::cout<<it->first<< " " <<(it->second).size() <<std::endl;
		var_vector.clear();
		
		svtig* svtig_tmp = nullptr;
		int start_pos = -1000;
		int end_pos = std::numeric_limits<int>::max();
		bool first = true;
		
		//If the distance between any two SVs is smaller than MIN_SV_DISTANCE, 
		//we merge the reads of these svs
		for (auto &sv : it->second) 
		{
			//std::cout<<"\t"<<sv->ref_start<<std::endl;
			if((!first) && (end_pos + (MIN_SV_DISTANCE / 2) >= sv->ref_start))
			{
				//std::cout<<"\t-------------IN1 start_pos:"<<start_pos<<" - sv->ref_start"<<sv->ref_start<<std::endl;
				(svtig_tmp->reads_h1).insert((sv->reads_h1).begin(), (sv->reads_h1).end());
				if (sv->ref_end > end_pos)
				{
					end_pos = sv->ref_end;
					svtig_tmp->end_pos = sv->ref_end;
				}
				merged_del_count++;
			}
			else
			{
				if((!first) && (svtig_tmp->reads_h1.size() > MIN_READ_SUPPORT))	
					var_vector.push_back(svtig_tmp);
					
				svtig_tmp = new svtig();
				svtig_tmp->sv_type = DELETION;
				svtig_tmp->contig = it->first;
				svtig_tmp->phased = false;
				svtig_tmp->sv_size = sv->sv_size;
				start_pos = sv->ref_start;
				end_pos = sv->ref_end;
				
				svtig_tmp->start_pos = sv->ref_start;
				svtig_tmp->end_pos = sv->ref_end;

				(svtig_tmp->reads_h1).insert((sv->reads_h1).begin(), (sv->reads_h1).end());
				first = false;
			}
		}
		//Add the last one as well
		if((!first) && (svtig_tmp->reads_h1.size() > MIN_READ_SUPPORT))	
			var_vector.push_back(svtig_tmp);
		
		//if (var_vector.size() >	MIN_READ_SUPPORT)
		final_del.insert(std::pair<std::string, std::vector<svtig*>>(it->first, var_vector));
	}

	
	std::cout<<"--->merged "<<merged_ins_count<<" insertions and "<<merged_del_count<<" deletions (with "<<MIN_SV_DISTANCE<<" as the SV distance threshold)\n--->there are "<<final_ins.size()<<" insertion and " <<final_del.size()<<" deletion svtigs\n\n";

	//std::vector<int> read_size;
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
	return RETURN_SUCCESS;
}


variant* generate_sv_node(std::map<std::string, gfaNode*>& gfa, std::string path, int path_start, int path_end, int ref_pos, int cigar_len, char sv_type)
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
			else if (sv_type == DELETION)
			{
				v->ref_end = v->ref_start + cigar_len;
				v->sv_type = DELETION;			 
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
				else if (sv_type == DELETION)
				{
					v->ref_end = v->ref_start + cigar_len;
					v->sv_type = DELETION;			 
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
				else if (sv_type == DELETION)
				{
					v->ref_end = v->ref_start + cigar_len;
					v->sv_type = DELETION;			 
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
				else if (sv_type == DELETION)
				{
					v->ref_end = v->ref_start + cigar_len;
					v->sv_type = DELETION;			 
				}
				return v;	
			}
		 	else		
				total_so_far = node_map_size;

		}
	}
	std::cout<<"ERRROORRRRRR -"<< total_so_far<< " "<<ref_pos<<std::endl;
	return NULL;
}



