#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "alignment.h"
#include "interval_tree.h"
#include "sv.h"

void gaf_traverse(gaf** hashtab)
{
	gaf *ptr = NULL;
	for(int i = 0; i < HASHSIZE - 1; i++)
	{
		printf("BUCKET %d\n", i);
		for(ptr = hashtab[i]; ptr != NULL; ptr = ptr->next)
		{
			printf("%s - %s (%d - %d)\n", ptr->contig, ptr->node->read_name, ptr->node->start, ptr->node->end);
		}
		printf("\n\n");
	}
}

gaf *gaf_lookup(gaf **hashtab, char *s)
{
    gaf *np;
	//printf("lookup %s - %u", s, hash(s));
    for (np = hashtab[hash(s)]; np != NULL; np = np->next)
	{
        if (strcmp(s, np->contig) == 0)
        	return np;
	}
	return NULL;
}

gaf *gaf_insert(gaf **hashtab, char *contig_name, alignment *node)
{
    gaf *np;
    unsigned hashval;
	
	if ((np = gaf_lookup(hashtab, contig_name)) == NULL)
	{
        np = (gaf*) getMem(sizeof(*np));
        if (np == NULL)
        	return NULL;
		np->contig = strdup(contig_name);
		np->node = node;
        hashval = hash(contig_name); 
		np->next = hashtab[hashval];
        hashtab[hashval] = np;
    } 
    if ((np->node = node) == NULL)
    	return NULL;
   
   	return np;
}


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

int read_alignments(parameters *params, gfa* gfa_table[], gaf* gaf_table[])
{
	FILE* fp;
    char* line = NULL;
    size_t len = 0;
    ssize_t read;
	char *mytoken = NULL;
	int node_len, offset, field_cnt = 0, line_cnt = 1;
	int secondary_filter = 0, cigar_cnt, ref_pos = 0;
	char *cigar = NULL, *read_name = NULL, *path = NULL;
	char strand;
	int path_length, path_start, path_end;
	int cigarLen[50000];
	char cigarOp[50000];
	char *tmp_str = NULL;
	

	fp = safe_fopen(params->gaf, "r");
	while ((read = getline(&line, &len, fp)) != -1) {
        //printf("%s", line);
		
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
		
		/*for (int c = 0; c < cigar_cnt; c++)
		{
			if (cigarOp[c] == 'I' && cigarLen[c] > 50)
			{
				sv* var = generate_sv_node(gfa_table, ref_pos, strand, path, path_start, path_end, cigarLen[c]);
			}
			if (cigarOp[c] != 'I')
				ref_pos += cigarLen[c];
		}*/
		ref_pos = 0;	
    }
    fclose(fp);
    
	if (line)
        free(line);
    
	//printf("Height is %d\n",find_height(t));
	//inorder(t);	
		
	return RETURN_SUCCESS;
}

