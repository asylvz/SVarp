#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "alignment.h"
#include "interval_tree.h"

void init_alignment(alignment** aln)
{
	/* initialize parameters */
	*aln = (alignment*) getMem(sizeof(alignment));
	(*aln)->read_name = NULL;
	(*aln)->start = -1;
	(*aln)->end = -1;
	(*aln)->strand = -1;
	(*aln)->node = NULL;
	(*aln)->path = NULL;
}


treenode* find_alignment_in_gfa(gfa* hashtab[], treenode* root, char* read_name, char strand, char* path, int path_start, int path_end)
{
	alignment* aln = NULL;
	gfa* tmp = NULL;
	
	char *copy = strdup(path);
	char *delim = "><";
	char *mytoken = strtok(copy, delim);
	int offset = 0;
	int node_count = 1;
	
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

		offset += strlen(mytoken) + 1;
		mytoken = strtok(NULL, delim);
		
		if ((node_count == 1) && (mytoken == NULL)) //which means there is only a single node
		{
			if(aln->strand == '>')
			{
				tmp = gfa_lookup(hashtab, aln->node);
				aln->start = tmp->node->offset + path_start;
				aln->end = tmp->node->offset + path_end;
			}
			else
			{
				tmp = gfa_lookup(hashtab, aln->node);
				aln->start = tmp->node->offset + (tmp->node->len - path_end);
				aln->end = tmp->node->offset + (tmp->node->len - path_start);
			}
			total_so_far += (aln->end - aln->start);
			total_node_length += tmp->node->len;
		}
		else if((node_count == 1) && (mytoken != NULL))
		{
			if(aln->strand == '>')
			{
				tmp = gfa_lookup(hashtab, aln->node);
				aln->start = tmp->node->offset + path_start;
				aln->end = tmp->node->offset + tmp->node->len;
			}
			else
			{
				tmp = gfa_lookup(hashtab, aln->node);
				aln->start = tmp->node->offset;
				aln->end = tmp->node->offset + (tmp->node->len - path_start);
			}
			total_so_far += (aln->end - aln->start);
			total_node_length += tmp->node->len;
		}
		else if(mytoken == NULL)
		{
			if(aln->strand == '>')
			{
				tmp = gfa_lookup(hashtab, aln->node);
				aln->start = tmp->node->offset;
				aln->end = tmp->node->offset + (total_path_length - total_so_far);
			}
			else
			{
				tmp = gfa_lookup(hashtab, aln->node);
				aln->start = tmp->node->offset + (tmp->node->len - (total_path_length - total_so_far));
				aln->end = tmp->node->offset + tmp->node->len;
			}
		}
		else //middle node
		{
			tmp = gfa_lookup(hashtab, aln->node);
			aln->start = tmp->node->offset;
			aln->end = tmp->node->offset + tmp->node->len;
			total_so_far += tmp->node->len;
			total_node_length += tmp->node->len;
		}
		root = insert_treenode(root, aln);
		
		//printf("mytoken = %s\n", mytoken);
		//printf("%c%s ", aln->strand, aln->node);
	}
	free(copy);
	
	return root;
}

int read_alignments(parameters *params, gfa* gfa_table[])
{
	FILE* fp;
    char* line = NULL;
    size_t len = 0;
    ssize_t read;
	char *mytoken = NULL;
	int node_len, offset, field_cnt = 0, line_cnt = 1;
	int secondary_filter = 0;
	char *cigar = NULL, *read_name = NULL, *path = NULL;
	char strand;
	int path_length, path_start, path_end;
	int cigarLen[50000];
	char cigarOp[50000];
	char *tmp_str = NULL;
	
	treenode* t = NULL;
	
	fp = safe_fopen(params->gaf, "r");
	while ((read = getline(&line, &len, fp)) != -1) {
        //printf("%s", line);
		//init_alignment(&aln);	
		
		field_cnt = 0;
		mytoken = strtok(line, "\t" );
		while(mytoken != NULL)
		{
			//printf("%d - %s\n", field_cnt, mytoken);
      		switch(field_cnt)
			{
				case 0:
					set_str(&read_name, mytoken);
					break;
				case 4:
					strand = mytoken[0];
					break;
				case 5:
					set_str(&path, mytoken);
					break;
				case 6:
					path_length = atoi(mytoken);	
					break;
				case 7:
					path_start = atoi(mytoken);
					break;
				case 8:
					path_end = atoi(mytoken);
					break;
				case 11:
					if (atoi(mytoken) < 20)
						goto skip_alignment;
					break;
				default:

					if (strstr(mytoken, "tp:A:") == mytoken)
					{
						//primary mapping
						char* tmp = substr(mytoken, 5, 6);
						if ( strcmp(tmp, "P") == 0)
							secondary_filter++;
							break;
					}
					else if (strstr(mytoken, "cg:Z:") == mytoken)
					{
						//get the Cigar
						cigar = substr(mytoken, 5, strlen(mytoken) - 1);
						//printf("%s\n", cigar);
						int cigar_offset = 0, str_offset = 0, cigar_cnt = 0;
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
		if(line_cnt == 100)
			break;
		line_cnt++;
		t = find_alignment_in_gfa(gfa_table, t, read_name, strand, path, path_start, path_end);
    }

    fclose(fp);
    
	if (line)
        free(line);
    
	printf("Height is %d\n",find_height(t));
	inorder(t);	
		
	return RETURN_SUCCESS;
}

