
typedef struct _refs
{
	char* node; /* name of the node */
	char* contig; /* name of the contig that the node belongs to, e.g., SN:Z:chr1 */
	char* sequence; /* sequence of the node - second column of the GFA */
	int len; /*node length, e.g., LN:i:74 */
	int offset; /*node offset within the contig e.g., SO:i:10746 */
} reference;

typedef struct _treenode
{
	int pos;
	int max;
	int height;
	struct _treenode *left, *right;
} treenode;


typedef struct _gfa { /* table entry: */
    struct _gfa *next; /* next entry in chain */
    char *node_name; /* defined name */
    char *defn; /* replacement text */
} gfa;


gfa *gfa_insert(char *name, char *defn);
gfa *gfa_lookup(char *s);

int read_reference_graph(parameters* params);
void init_reference(reference** ref);