#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "alignment.h"

/*int read_reference_graph(parameters *params)
{
	FILE* fp;
    char* line = NULL;
    size_t len = 0;
    ssize_t read;
	char *mytoken = NULL, *str_tmp = NULL, *node_name = NULL, *contig_name = NULL, *seq = NULL;
	int node_len, offset, field_cnt = 0;

	fp = safe_fopen(params->ref_graph, "r");
	while ((read = getline(&line, &len, fp)) != -1) {
        //printf("%s", line);
		
		field_cnt = 0;
		mytoken = strtok(line, "\t" );
		while(mytoken != NULL)
		{
			if(field_cnt == 0 && strcmp(mytoken, "S") != 0)
				break;

      		switch(field_cnt)
			{
				case 1:
					set_str(&node_name, mytoken);
					break;
				case 2:
					set_str(&seq, mytoken);
					break;
				case 3:
					str_tmp = substr(mytoken, 5, strlen(mytoken));
					node_len = atoi(str_tmp);
					free(str_tmp);
					break;
				case 4:
					str_tmp = substr(mytoken, 5, strlen(mytoken));
					set_str(&contig_name, str_tmp);
					free(str_tmp);
					break;
				case 5:
					str_tmp = substr(mytoken, 5, strlen(mytoken));
					offset = atoi(str_tmp);
					free(str_tmp);
					//printf("%s - %d\n",mytoken, offset);
					break;
			}
			field_cnt++;
			mytoken = strtok(NULL, "\t");
		}
    }

    fclose(fp);
    
	if (line)
        free(line);
    
	return RETURN_SUCCESS;
}*/
