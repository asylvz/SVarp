#include <stdio.h>
#include "interval_tree.h"
#include "sv.h"

treenode* new_node(alignment* aln)
{
	treenode *temp = (treenode*) getMem(sizeof(treenode));
    temp->node = aln;
    temp->max = aln->end;
	temp->height = 1;
    temp->left = temp->right = NULL;
	
	return temp;
}

int find_height(treenode* root)
{
	if (root == NULL)
		return 0;
	else
		return root->height;
}

treenode* insert_treenode(treenode *root, alignment* aln)
{
	// Base case: Tree is empty, new node becomes root
    if (root == NULL)
        return new_node(aln);
	else if(aln->start < root->node->start)
		root->left = insert_treenode(root->left, aln);
	else
		root->right = insert_treenode(root->right, aln);
	
   	root->height = 1 + max(find_height(root->left), find_height(root->right));	
    
	// Update the max value of this ancestor if needed
    if (root->max < aln->end)
        root->max = aln->end;
 
    return root;	
}

treenode* find_overlaps(treenode* root, variant* sv, variant* lst[1001], int n)
{
	/*This is based on exercise 13.4-3 in Introduction to Algorithms 3rd edition Cormen. Each branch returns an interval (at least), as there are k branches in the tree - O(klgn). It would ideally be a red-black tree or an avl tree
	Also benefited from:
	- https://www.geeksforgeeks.org/interval-tree/
	- http://www.davismol.net/2016/02/07/data-structures-augmented-interval-tree-to-search-for-interval-overlapping/
    - https://github.com/gzc/CLRS/blob/master/C14-Augmenting-Data-Structures/14.3.md
  */
	if ((root->node->end >= sv->ref_start >= root->node->start) || (root->node->end >= sv->ref_end >= root->node->start))
	{
		lst[n] = sv;
		n += 1;
	}
	if (root->left != NULL && root->left->max >= sv->ref_start)
		find_overlaps(root->left, sv, lst, n);
	
	if (root->right != NULL && root->right->max >= sv->ref_start)
		find_overlaps(root->right, sv, lst, n);

  	
	return root;
}

void inorder(treenode *root) {
    if (root == NULL)
        return;

    inorder(root->left);

    printf("[%d, %d] max = %d\n", root->node->start, root->node->end, root->max); 
    inorder(root->right);
}