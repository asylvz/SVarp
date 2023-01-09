#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <map>
#include <iterator>
#include "alignment.h"

int decompose_cigars(string cigar, int (&cigarLen)[10000], char (&cigarOp)[10000])
{
	/*get the Cigar*/
	int cigar_offset = 0, str_offset = 0, cigar_cnt = 0;
	char* cigar_copy = (char*) cigar.c_str();	
	char *tmp_str = (char*) malloc(sizeof(char) * 6);
	//cout<< cigar_copy<<"\n\n";
	while(cigar_offset < cigar.length())
	{
		if (isdigit(*(cigar_copy + cigar_offset)) == 0)
		{
			cigarOp[cigar_cnt] = *(cigar_copy + cigar_offset);
			//cout<<cigarOp[cigar_cnt]<<endl;
			cigarLen[cigar_cnt] = atoi(tmp_str);
			//printf("-->(%d)%d%c\n", cigar_cnt, cigarLen[cigar_cnt], cigarOp[cigar_cnt]);
			free(tmp_str);
			tmp_str = (char*) malloc(sizeof(char) * 6);
			str_offset = 0;
			cigar_cnt++;
		}
		else 
		{
			*(tmp_str + str_offset) = *(cigar_copy + cigar_offset);
			//cout<<tmp_str<<endl;
			str_offset++;
		}
		//cout<<cigar_offset<<endl; 
		cigar_offset++;
	}
	if(tmp_str != NULL)
		free(tmp_str);
	
	return cigar_cnt;
}


int read_alignments(parameters *params, std::map<std::string, gfaNode*> ref, std::vector<variant*>& insertions)
{
	int secondary = 0, primary = 0, insertion_count = 0;
	
	std::cout << "Reading the GAF file"<< std::endl;
	
	std::string line;	
	std::vector <std::string> tokens;
	std::ifstream fp(params->gaf);
	std::map<std::string, alignment*> gaf;	
	
	while(fp)
	{
		getline(fp, line);

		std::string tmp_str;
		std::stringstream s(line);
		tokens.clear();
		
		while(getline(s, tmp_str, '\t'))
        	tokens.push_back(tmp_str);
		
		if(stoi(tokens[11]) < MINMAPQ)
			continue;
			
		bool isPrimary = true;	
		string cigar;
		for (auto& tok : tokens) {
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
				int cigarLen[10000] = {0};
				char cigarOp[10000] = {0};
				cigar = tok.substr(5);
				int cigar_cnt = decompose_cigars(cigar, cigarLen, cigarOp);
				int ref_pos = 0;
				for (int c = 0; c < cigar_cnt; c++)
				{
					if (cigarOp[c] == INSERTION && cigarLen[c] > MINSVSIZE)
					{
						variant* var = generate_sv_node(ref, tokens[5], stoi(tokens[7]), stoi(tokens[8]), ref_pos, cigarLen[c], INSERTION);
						if (var)
						{
							insertion_count++;
							insertions.push_back(var);
						}
					}
					//cout<<"Added"<<endl;
					if (cigarOp[c] != 'I')
						ref_pos += cigarLen[c];
				}	
			}
   		}
		if(!isPrimary)
			continue;

		alignment_within_gfa(gaf, ref, tokens);	
		if(primary > 1000)
			break;
}
	cout<<"There are "<<primary<<" primary mappings and "<<insertion_count<<" insertions\n"<<endl;
	
	return RETURN_SUCCESS;
}


void alignment_within_gfa(map<string, alignment*>& gaf, map<string, gfaNode*> gfa, vector <std::string> tokens)
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

		gaf.insert(std::pair<std::string, alignment*>(gfa[aln->node]->contig, aln));
		//cout<<"inserted "<<gfa[aln->node]->contig<<endl;
	}
}
