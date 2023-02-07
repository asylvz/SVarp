#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <limits>
#include "alignment.h"


void contig_coverage(std::map <std::string, Contig*>& ref, std::map<std::string, gfaNode*>& gfa, std::vector <std::string>& tokens)
{
	char *path_copy = (char *) tokens[5].c_str();
	
	char *mytoken = strtok(path_copy,"><");
	
	int node_count = 0, total_so_far = 0;
	int path_start = stoi(tokens[7]);
	int path_end = stoi(tokens[8]);
	int total_path_length = path_end - path_start;
	while(mytoken) 
	{
		//std::cout<<mytoken<<std::endl;
		node_count += 1;
		std::string node = mytoken;

		mytoken = strtok(NULL, "><");
		//std::cout<<mytoken<<std::endl;
		
		if ((node_count == 1) && (!mytoken)) //the single node
		{
			std::string contig = gfa[node]->contig; 
			ref[contig]->mapped_bases += path_end - path_start;	
			ref[contig]->mapped_reads++;	

			ref["overall"]->mapped_reads++;	
			ref["overall"]->mapped_bases += path_end - path_start;
		}
		else if((node_count == 1) && (mytoken)) //first node
		{
			std::string contig = gfa[node]->contig; 
			ref[contig]->mapped_bases += gfa[node]->len - path_start;	
			ref[contig]->mapped_reads++;	

			ref["overall"]->mapped_reads++;	
			ref["overall"]->mapped_bases += gfa[node]->len - path_start;
			
			total_so_far += gfa[node]->len - path_start;
		}
		else if(!mytoken) //Last node
		{
			std::string contig = gfa[node]->contig; 
			ref[contig]->mapped_bases += total_path_length - total_so_far;
			ref[contig]->mapped_reads++;	

			ref["overall"]->mapped_reads++;	
			ref["overall"]->mapped_bases += total_path_length - total_so_far;

			//if(total_path_length - total_so_far < 0)
			//	std::cout<< "middles: "<<middles <<" " <<total_path_length<<" "<<total_so_far<<std::endl;
		}
		else //middle node
		{
			std::string contig = gfa[node]->contig; 
			ref[contig]->mapped_bases += gfa[node]->len;		
			ref[contig]->mapped_reads++;	

			ref["overall"]->mapped_reads++;	
			ref["overall"]->mapped_bases += gfa[node]->len;		

			total_so_far += gfa[node]->len;
		}
	}
}


int read_alignments(parameters *params, std::map <std::string, Contig*>& ref, std::map<std::string, gfaNode*> gfa, std::map<std::string, variant*>& variations)
{
	int secondary = 0, primary = 0, insertion_count = 0, deletion_count = 0, line_count = 0, perc = 0;
	
	std::cout<<"Reading the GAF file"<<std::endl;
	
	std::string line;	
	std::vector <std::string> tokens;
	std::vector<int> cigarLen;
	std::vector<char> cigarOp;
	
	std::ifstream fp(params->gaf);
	if(!fp.good())
	{
        std::cerr << "Error opening '"<<params->gaf<< std::endl;
        return RETURN_ERROR;
    }

    int total_line_count = 0;
    char endline_char = '\n';
    while (fp.ignore(std::numeric_limits<std::streamsize>::max(), fp.widen(endline_char)))
		++total_line_count;
	
	fp.clear() ; // clear the failed state of the stream
	fp.seekg(0) ; // seek to the first character in the file
	
	if (total_line_count > TEST_SAMPLE_SIZE)
		total_line_count = TEST_SAMPLE_SIZE;

	while(fp)
	{
		getline(fp, line);
		line_count++;

		std::string tmp_str;
		std::stringstream s(line);
		tokens.clear();
		
		while(getline(s, tmp_str, '\t'))
        	tokens.push_back(tmp_str);
		
		if(stoi(tokens[11]) < MINMAPQ)
			continue;
			
		bool isPrimary = true;	
		std::string cigar;
		for (auto& tok : tokens) 
		{
			if(strstr(tok.c_str(), "tp:A:"))
			{
				if (tok.substr(5, 6) != "P")
				{
					isPrimary = false;
					secondary++;
					break;
				}
				primary++;
			}
			else if(strstr(tok.c_str(), "cg:Z:"))
			{
				cigarLen.clear();
				cigarOp.clear();

				cigar = tok.substr(5);
				int cigar_cnt = decompose_cigars(cigar, cigarLen, cigarOp);
				int ref_pos = 0;
				for (int c = 0; c < cigar_cnt; c++)
				{
					if (cigarOp[c] == INSERTION && cigarLen[c] > MINSVSIZE)
					{
						variant* var = generate_sv_node(gfa, tokens[5], stoi(tokens[7]), stoi(tokens[8]), ref_pos, cigarLen[c], INSERTION);
						
						if (var)
						{
							var->sv_size = cigarLen[c];	
							insertion_count++;
								
							std::string var_name = var->contig + ":" + std::to_string(var->ref_start) + "_" + std::to_string(var->ref_end);
							std::map<std::string, variant*>::iterator it = variations.find(var_name);
							if (it != variations.end())
								it->second->reads_h1.insert(tokens[0]);	
							else
							{
								var->reads_h1.insert(tokens[0]);
								variations.insert(std::pair<std::string, variant*>(var_name, var));
							}
						}
						else
							std::cout<<"RETURNED NULL\n";
					}
					else if (cigarOp[c] == DELETION && cigarLen[c] > MINSVSIZE)
					{
						variant* var = generate_sv_node(gfa, tokens[5], stoi(tokens[7]), stoi(tokens[8]), ref_pos, cigarLen[c], DELETION);
						
						if (var)
						{
							var->sv_size = cigarLen[c];	
							deletion_count++;
								
							std::string var_name = var->contig + ":" + std::to_string(var->ref_start) + "_" + std::to_string(var->ref_end);
							std::map<std::string, variant*>::iterator it = variations.find(var_name);
							if (it != variations.end())
								it->second->reads_h1.insert(tokens[0]);	
							else
							{
								var->reads_h1.insert(tokens[0]);
								variations.insert(std::pair<std::string, variant*>(var_name, var));
							}
						}
						else
							std::cout<<"RETURNED NULL\n";

					}

					if (cigarOp[c] != 'I')
						ref_pos += cigarLen[c];
				}
			}
   		}
		//Check the read depth/coverage here
		contig_coverage(ref, gfa, tokens);	

		if(!isPrimary)
			continue;

		int perc_tmp = ((double) line_count / total_line_count) * 100;
		if (perc_tmp > perc)
		{
			perc = perc_tmp;
			fprintf(stderr, "--->%d%%\r",perc);
			fflush(stderr);
		}

		//alignment_within_gfa(gaf, ref, tokens);
		if(line_count > TEST_SAMPLE_SIZE)
			break;
	}
	std::cout<<"\n--->there are "<<primary<<" primary mappings and "<<insertion_count<<" insertion - "<<deletion_count<<" deletion loci\n\n";

	std::map<std::string, Contig*>::iterator it;
	for (it=ref.begin(); it != ref.end(); ++it)
	{
		it->second->coverage = (double) it->second->mapped_bases / it->second->contig_length;	
		logFile<< it->first<<"---> LEN= "<<it->second->contig_length<<" - Mapped bases= "<< it->second->mapped_bases<<" - Mapped reads= "<<it->second->mapped_reads << " - Cov= "<<it->second->coverage <<std::endl;		
	}			
	
	return RETURN_SUCCESS;
}


/*void alignment_within_gfa(std::multimap<string, alignment*>& gaf, std::map<string, gfaNode*> gfa, vector <std::string> tokens)
{
	char *path_copy = (char *) tokens[5].c_str();
	
	char *mytoken = strtok(path_copy,"><");
	
	int offset = 0, node_count = 0, total_so_far = 0, total_node_length = 0;
	int path_start = stoi(tokens[7]);
	int path_end = stoi(tokens[8]);
	int total_path_length = path_end - path_start;

	
	while(mytoken) 
	{
		alignment *aln = new alignment();
			
		aln->node = mytoken;
		aln->strand = tokens[5][offset];
		aln->read_name = tokens[0];
		aln->path = tokens[5];
		node_count += 1;
		offset += strlen(mytoken) + 1;
		mytoken = strtok(NULL, "><");
		
		if ((node_count == 1) && (mytoken == NULL)) //which means there is only a single node
		{
			if(aln->strand == '>')
			{
				aln->start = gfa[aln->node]->offset + path_start;
				aln->end = gfa[aln->node]->offset + path_end;
			}
			else
			{
				aln->start = gfa[aln->node]->offset + (gfa[aln->node]->len - path_end);
				aln->end = gfa[aln->node]->offset + (gfa[aln->node]->len - path_start);
			}	
			
			if(aln->end < aln->start)
				cout<<"problem1";
			
			total_so_far += (aln->end - aln->start);
			total_node_length += gfa[aln->node]->len;
		}
		else if((node_count == 1) && (mytoken != NULL))
		{
			if(aln->strand == '>')
			{
				aln->start = gfa[aln->node]->offset + path_start;
				aln->end = gfa[aln->node]->offset + gfa[aln->node]->len;
			}
			else
			{
				aln->start = gfa[aln->node]->offset;
				aln->end = gfa[aln->node]->offset + (gfa[aln->node]->len - path_start);
			}
			
			if(aln->end < aln->start)
				cout<<"problem2";

			total_so_far += (aln->end - aln->start);
			total_node_length += gfa[aln->node]->len;
		}
		else if(mytoken == NULL)
		{
			if(aln->strand == '>')
			{
				aln->start = gfa[aln->node]->offset;
				aln->end = gfa[aln->node]->offset + (total_path_length - total_so_far);
			}
			else
			{
				aln->start = gfa[aln->node]->offset + (gfa[aln->node]->len - (total_path_length - total_so_far));
				aln->end = gfa[aln->node]->offset + gfa[aln->node]->len;
			}

			if(aln->end < aln->start)
				cout<<"problem3"; 
		}
		else //middle node
		{
			aln->start = gfa[aln->node]->offset;
			aln->end = gfa[aln->node]->offset + gfa[aln->node]->len;
			
			total_so_far += gfa[aln->node]->len;
			total_node_length += gfa[aln->node]->len;

			if(aln->end < aln->start)
				cout<< "problem4";
		}
		
		//if (gfa[aln->node]->contig == "CHM13#0#chr1")
		//{
		//	cout<<"aln - \tstart: "<<aln->start<<"\tend: "<<aln->end<<endl;
		//}
		gaf.insert(std::pair<std::string, alignment*>(gfa[aln->node]->contig, aln));
		//cout<<"inserted "<<gfa[aln->node]->contig<<endl;
	}
}*/



/*void find_supporting_reads(std::map<std::string, gfaNode*> ref, std::multimap<std::string, alignment*> aln, std::set<std::string> contigs, std::multimap<std::string, variant*>& insertions)
{
	
	cout<<"Finding the supporting reads"<<endl;
	treenode* root;
	
	for(auto c: contigs)
	{
		//Insert the alignments of this contig into the interval tree		
		root = NULL;
		auto aln_range = aln.equal_range(c);
		for (auto i = aln_range.first; i != aln_range.second; ++i)
			root = insert_treenode(root, i->second);
		
		//if (find_height(root) >0 )
		
		//if (c == "CHM13#0#chr1")
		//	cout<<"Contig:"<<c<<" tree height: "<<find_height(root) <<endl;
		
		auto sv_range = insertions.equal_range(c);
		for (auto i = sv_range.first; i != sv_range.second; ++i)
		{
			std::set <alignment*> overlaps;
			find_overlaps(root, i->second, overlaps);
			
			//if (i->second->contig == "CHM13#0#chr1")
			//{
			//	cout<<"For SV = "<< i->second->ref_start<< " "<< i->second->ref_end<<" - "<< overlaps.size()<<endl;
			//}
			for (auto t:overlaps)
			{
				//cout<<t->read_name<<endl;
				i->second->reads_h1.insert(t->read_name);
				//cout<<i->second->reads.size();
			}
		}
	}
}

*/
