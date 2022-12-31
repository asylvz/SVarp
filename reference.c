#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "reference.h"


static gfa *hashtab[HASHSIZE]; /* pointer table */

/*treenode insert_treenode(treenode *root)
{
	
}
*/

int read_reference_graph(parameters* params)
{
	FILE* fp;
    char* line = NULL;
    size_t len = 0;
    ssize_t read;
	char *mytoken = NULL, *str_tmp = NULL;
	int field_cnt = 0;

	reference* ref;
	init_reference(&ref);

	fp = safe_fopen(params->ref_graph, "r");
	while ((read = getline(&line, &len, fp)) != -1) {
        //printf("%s", line);
		
		field_cnt = 0;
		mytoken = strtok(line, "\t" );
		while(mytoken != NULL)
		{
      		switch(field_cnt)
			{
				case 0:
					if (strcmp(mytoken, "S") != 0)
						goto exit_while_loop;
					break;
				case 1:
					set_str(&(ref->node), mytoken);
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

gfa *gfa_lookup(char *s)
{
    gfa *np;
    for (np = hashtab[hash(s)]; np != NULL; np = np->next)
        if (strcmp(s, np->node_name) == 0)
          return np; /* found */
    return NULL; /* not found */
}

gfa *gfa_insert(char *node_name, char *defn)
{
    gfa *np;
    unsigned hashval;
    
	if ((np = gfa_lookup(node_name)) == NULL) /* not found */
	{ 
        np = (gfa*) malloc(sizeof(*np));
        if (np == NULL || (np->node_name = strdup(node_name)) == NULL)
        	return NULL;
        hashval = hash(node_name);
        np->next = hashtab[hashval];
        hashtab[hashval] = np;
    } 
	else /* already there */
        free((void *) np->defn); /*free previous defn */
    if ((np->defn = strdup(defn)) == NULL)
    	return NULL;
    return np;
}
