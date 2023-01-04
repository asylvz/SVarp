#include <stdio.h>
#include <stdlib.h>
#include "reference.h"
#include "alignment.h"

int free_gfa(gfa** hash_table)
{
	gfa *ptr, *ptr_next;
	for (int i = 0; i < HASHSIZE; i++)
	{
		for (ptr = hash_table[i]; ptr != NULL; ptr = ptr->next)
		{
			if (ptr->node->name != NULL)	
				free(ptr->node->name);
			if (ptr->node->contig != NULL)	
				free(ptr->node->contig);
			if (ptr->node->sequence != NULL)
				free(ptr->node->sequence);
		}
		free(hash_table[i]);
	}
	return RETURN_SUCCESS;
}


int free_gaf(gaf** hash_table)
{
	gaf *ptr, *ptr_next;
	for (int i = 0; i < HASHSIZE; i++)
	{
		for (ptr = hash_table[i]; ptr != NULL; ptr = ptr->next)
		{
			if (ptr->contig != NULL)	
				free(ptr->contig);
			if (ptr->node->read_name != NULL)	
				free(ptr->node->read_name);
			if (ptr->node->node != NULL)	
				free(ptr->node->node);
			if (ptr->node->path != NULL)
				free(ptr->node->path);
		}
		free(hash_table[i]);
	}
	return RETURN_SUCCESS;
}
