#ifndef COMMON_IMPL_HH
#define COMMON_IMPL_HH

#include <vector>
#include <string>
#include <fstream>
#include <iterator>
#include <map>
#include <algorithm>
#include <iostream>
#include <stack>

template<typename T>
void 
create_binary_tree(NHNode *root, Binary_tree<T> &tree,
                   typename Binary_tree<T>::vid_t ext)
{
    if (root->children == 0)
        return;

    char *left_label = root->children->node->label;
    char *right_label = root->children->next->node->label;
    tree.add_children(ext, T(1), T(1),
                      left_label == 0 ? "" : left_label,
                      right_label == 0 ? "" : right_label);

    create_binary_tree(root->children->node, tree, tree.left(ext));
    create_binary_tree(root->children->next->node, tree, tree.right(ext));
}

template<typename T> 
std::vector<typename Binary_tree<T>::vid_t>
create_gene_species_map(const Binary_tree<T> &species_tree, 
                        const Binary_tree<T> &gene_tree, 
                        std::string map_filename)
{
    using namespace std;
    typedef typename Binary_tree<T>::vid_t vid_t;

    vector<string> map_file_content;
    
    ifstream map_file(map_filename.c_str());
    copy(istream_iterator<string>(map_file), istream_iterator<string>(),
         back_inserter(map_file_content));

    /* Make sure there are even number of strings in map file. */
    if (map_file_content.size() == 0 || map_file_content.size() % 2 != 0)
        {
            throw Some_exception("error reading map file.");
        } 

    /* Create a map from gene name to species name. */
    map<string, string> str_sigma;
    for (unsigned i = 0; i < map_file_content.size(); i += 2)
        str_sigma[map_file_content[i]] = map_file_content[i+1];
    
    /* Create the final map called sigma mapping vid_t to vid_t */
    vector<vid_t> sigma(gene_tree.size());
    for (vid_t v = 0; v <= gene_tree.last(); ++v)
        {
            if (!gene_tree.is_leaf(v))
                continue;
            
            string gene_label = gene_tree.label(v);
            string species_label = str_sigma[gene_label];
            if (species_label == "")
                {
                    string message = 
                        "gene label '" + gene_label + "' "
                        "is missing in map file.";
                    throw Some_exception(message.c_str());
                }
            if (species_tree[species_label] == Binary_tree<T>::NONE)
                {
                    string message = 
                        "species label '" + species_label + "' "
                        "which occurs in map file "
                        "does not exist in species tree.";
                    throw Some_exception(message.c_str());
                }
            sigma[v] = species_tree[species_label];
        }
    
    return sigma;
}

template<typename T>
void
output_subtree(std::ostream &out, const Binary_tree<T> &tree,
               typename Binary_tree<T>::vid_t subtree)
{
    if (tree.is_leaf(subtree))
        {
            out << subtree << "_" << tree.label(subtree);
            return;
        }
    
    out << "(";
    output_subtree(out, tree, tree.left(subtree));
    out << ",";
    output_subtree(out, tree, tree.right(subtree));
    out << ")" << subtree << "_" << tree.label(subtree);
}



template<typename T>
std::ostream &
operator<<(std::ostream &out, const Binary_tree<T> &tree)
{
    output_subtree(out, tree, 0);
    return out;
}


#endif /* COMMON_IMPL_HH */
