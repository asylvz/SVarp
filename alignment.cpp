#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
#include <iterator>
#include "alignment.h"


void find_cigars(int (&cigarLen)[10000], string (&cigarOp)[10000])
{
	cout<<"k";	
}


int read_alignments(parameters *params)
{
	int secondary = 0;
	int primary = 0;
	
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

		alignment *aln = new alignment();
		aln->read_name = tokens[0];
		aln->strand = tokens[4][0];
		aln->path = tokens[5];
		
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
				int cigarLen[10000];
				string cigarOp[10000];
				cigar = tok.substr(5);
				find_cigars(cigarLen, cigarOp);
				//cout<< cigar<<endl;
			}
   		}
		if(!isPrimary)
			continue;

		//aln->path_length = stoi(tokens[5].substr(5));

		//ref.insert(std::pair<std::string, gfaNode*>(g->name, g));
	}
	cout<<"There are "<<primary<<" primary mappings\n"<<endl;
	/*
				case 6:
					path_length = atoi(mytoken);	
					break;
				case 7:
					path_start = atoi(mytoken);
					break;
				case 8:
					path_end = atoi(mytoken);
					break;
				default:

					else if (strstr(mytoken, "cg:Z:") == mytoken)
					{
						//get the Cigar
						cigar = substr(mytoken, 5, strlen(mytoken) - 1);
						//printf("%s\n", cigar);
						int cigar_offset = 0, str_offset = 0; 
						cigar_cnt = 0;
						tmp_str = (char*) getMem(sizeof(char) * 6);
						memset(cigarLen, 0, 50000);
						memset(cigarOp, 0, 50000);
						
						while(cigar_offset < strlen(cigar))
						{
							if (isdigit(*(cigar + cigar_offset)) == 0)
							{
								cigarOp[cigar_cnt] = *(cigar + cigar_offset);
								cigarLen[cigar_cnt] = atoi(tmp_str);
								//printf("-->(%d)%d%c\n", cigar_cnt, cigarLen[cigar_cnt], cigarOp[cigar_cnt]);
								free(tmp_str);
								tmp_str = (char*) getMem(sizeof(char) * 6);	
								str_offset = 0;
								cigar_cnt++;
							}
							else 
							{
								*(tmp_str + str_offset) = *(cigar + cigar_offset);
								str_offset++;
							}
							cigar_offset++;
						}
						break;			
					}
			}
			field_cnt++;
			mytoken = strtok(NULL, "\t");
		}
		skip_alignment:;
		if(line_cnt == 50000)
			break;
		line_cnt++;
		//printf("Here\n");
		find_alignment_in_gfa(gfa_table, gaf_table, read_name, strand, path, path_start, path_end);
		
		//Find SVs
		
		//for (int c = 0; c < cigar_cnt; c++)
		//{
		//	if (cigarOp[c] == 'I' && cigarLen[c] > 50)
	//		{
	//			sv* var = generate_sv_node(gfa_table, ref_pos, strand, path, path_start, path_end, cigarLen[c]);
	//		}
	//		if (cigarOp[c] != 'I')
	//			ref_pos += cigarLen[c];
	//	}/
		ref_pos = 0;	
    }
    fclose(fp);
    
	if (line)
        free(line);
    
	//printf("Height is %d\n",find_height(t));
	//inorder(t);	
		*/
	return RETURN_SUCCESS;
}


/*
void find_alignment_in_gfa(gfa* gfatable[], gaf* gaftable[], char* read_name, char strand, char* path, int path_start, int path_end)
{
	alignment* aln = NULL;
	gfa* tmp = NULL;
	
	char *copy = strdup(path);
	char *delim = "><";
	char *mytoken = strtok(copy, delim);
	int offset = 0;
	int node_count = 0;
	
	int total_so_far = 0;
	int total_node_length = 0;
	int total_path_length = path_end - path_start;

	while(mytoken) 
	{
		init_alignment(&aln);
		aln->node = strdup(mytoken);
		aln->strand = path[offset];
		aln->read_name = strdup(read_name);
		aln->path = strdup(path);
		node_count += 1;
		offset += strlen(mytoken) + 1;
		mytoken = strtok(NULL, delim);
		
		if ((node_count == 1) && (mytoken == NULL)) //which means there is only a single node
		{
			if(aln->strand == '>')
			{
				tmp = gfa_lookup(gfatable, aln->node);
				aln->start = tmp->node->offset + path_start;
				aln->end = tmp->node->offset + path_end;
			}
			else
			{
				tmp = gfa_lookup(gfatable, aln->node);
				aln->start = tmp->node->offset + (tmp->node->len - path_end);
				aln->end = tmp->node->offset + (tmp->node->len - path_start);
			}	
			
			if(aln->end < aln->start)
				printf("problem1");
			total_so_far += (aln->end - aln->start);
			total_node_length += tmp->node->len;
		}
		else if((node_count == 1) && (mytoken != NULL))
		{
			if(aln->strand == '>')
			{
				tmp = gfa_lookup(gfatable, aln->node);
				aln->start = tmp->node->offset + path_start;
				aln->end = tmp->node->offset + tmp->node->len;
			}
			else
			{
				tmp = gfa_lookup(gfatable, aln->node);
				aln->start = tmp->node->offset;
				aln->end = tmp->node->offset + (tmp->node->len - path_start);
			}
			
			if(aln->end < aln->start)
				printf("\nproblem2 - node= %s, node_len = %d, path_start = %d\n", tmp->node->name, tmp->node->len, path_start);
			total_so_far += (aln->end - aln->start);
			total_node_length += tmp->node->len;
		}
		else if(mytoken == NULL)
		{
			if(aln->strand == '>')
			{
				tmp = gfa_lookup(gfatable, aln->node);
				aln->start = tmp->node->offset;
				aln->end = tmp->node->offset + (total_path_length - total_so_far);
			}
			else
			{
				tmp = gfa_lookup(gfatable, aln->node);
				aln->start = tmp->node->offset + (tmp->node->len - (total_path_length - total_so_far));
				aln->end = tmp->node->offset + tmp->node->len;
			}

			if(aln->end < aln->start)
				printf("problem3");
		}
		else //middle node
		{
			tmp = gfa_lookup(gfatable, aln->node);
			aln->start = tmp->node->offset;
			aln->end = tmp->node->offset + tmp->node->len;
			
			total_so_far += tmp->node->len;
			total_node_length += tmp->node->len;

			if(aln->end < aln->start)
				printf("problem4");
		}

		gaf_insert(gaftable, tmp->node->contig, aln);	
		//root = insert_treenode(root, aln);
		
		//printf("mytoken = %s\n", mytoken);
		//printf("%c%s ", aln->strand, aln->node);
	}
	free(copy);
	
}
*/
