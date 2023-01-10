#ifndef __INTERVAL_TREE
#define __INTERVAL_TREE
#include "alignment.h"

typedef struct _treenode
{
	alignment *node;
	int max;
	int height;
	struct _treenode *left, *right;
} treenode;


//TREE
treenode* insert_treenode(treenode *root, alignment *aln);
int find_height(treenode* root);
void inorder(treenode *root);
treenode* find_overlaps(treenode* root, variant* sv, std::set <alignment*>& overlaps);

#endif
