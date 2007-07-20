#include "common.hh"

bool
NHtree_is_binary(NHNode *root)
{
    /* true if root is a leaf. */
    if (root->children == 0)
        return true;

    /* false if root does not have exactly two children. */
    if (root->children->next == 0
        || root->children->next->next != 0)
        return false;

    /* otherwise it depends on the children of root. */
    return NHtree_is_binary(root->children->node) 
        && NHtree_is_binary(root->children->next->node);
}

