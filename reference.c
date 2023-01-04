#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "reference.h"


int read_reference_graph(parameters* params, gfa* hashtab[])
{
	FILE* fp;
    char* line = NULL;
    size_t len = 0;
    ssize_t read;
	char *mytoken = NULL, *str_tmp = NULL;
	int field_cnt = 0;
	bool insert_line;
	
	reference* ref;
	fp = safe_fopen(params->ref_graph, "r");
	
	while ((read = getline(&line, &len, fp)) != -1) {
        //printf("%s", line);
		
		init_reference(&ref);
		insert_line = true;
		field_cnt = 0;
		mytoken = strtok(line, "\t" );
		while(mytoken != NULL)
		{
      		switch(field_cnt)
			{
				case 0:
					if (strcmp(mytoken, "S") != 0)
					{
						insert_line = false;
						goto exit_while_loop;
					}
					break;
				case 1:
					set_str(&(ref->name), mytoken);
					break;
				case 2:
					set_str(&(ref->sequence), mytoken);
					break;
				case 3:
					str_tmp = substr(mytoken, 5, strlen(mytoken));
					ref->len = atoi(str_tmp);
					free(str_tmp);
					break;
				case 4:
					str_tmp = substr(mytoken, 5, strlen(mytoken));
					set_str(&(ref->contig), str_tmp);
					free(str_tmp);
					break;
				case 5:
					str_tmp = substr(mytoken, 5, strlen(mytoken));
					ref->offset = atoi(str_tmp);
					free(str_tmp);
					//printf("%s - %d\n",mytoken, ref->offset);
					break;
			}
			field_cnt++;
			mytoken = strtok(NULL, "\t");
		}
		exit_while_loop:;
		if (insert_line)
			gfa_insert(hashtab, ref->name, ref);	
    }

    fclose(fp);
   
	if (line)
        free(line);
    
	return RETURN_SUCCESS;
}


void init_reference(reference** ref)
{
	/* initialize parameters */
	*ref = (reference*) getMem(sizeof(reference));
	(*ref)->contig = NULL;
	(*ref)->sequence = NULL;
	(*ref)->len = -1;
	(*ref)->offset = -1;
}


void gfa_traverse(gfa** hashtab)
{
	gfa *ptr;
	for(int i = 0; i < HASHSIZE; i++)
	{
		for(ptr = hashtab[i]; ptr != NULL; ptr = ptr->next)
			printf("%s - %s\n", ptr->node->name, ptr->node->contig);
	}
}

gfa *gfa_lookup(gfa **hashtab, char *s)
{
    gfa *np;
	//printf("lookup %s - %u", s, hash(s));
    for (np = hashtab[hash(s)]; np != NULL; np = np->next)
	{
        if (strcmp(s, np->node->name) == 0)
        	return np;
	}
	return NULL;
}

gfa *gfa_insert(gfa **hashtab, char *node_name, reference *node)
{
    gfa *np;
    unsigned hashval;
	
	if ((np = gfa_lookup(hashtab, node_name)) == NULL)
	{
        np = (gfa*) getMem(sizeof(np));
        if (np == NULL)
        	return NULL;
		//np->node_name = strdup(node_name);
		np->node = node;
        hashval = hash(node_name); 
		np->next = hashtab[hashval];
        hashtab[hashval] = np;
    } 
	//else
    //    free((gfa *) np->node_name);
    if ((np->node = node) == NULL)
    	return NULL;
   
	//printf("\n");
   	return np;
}
