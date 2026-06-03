#include "function.hh"

namespace coot {

int get_level_gemmi(const tree<linked_residue_t_gemmi>& glyco_tree,
                    gemmi::Residue* residue_p) {
   int level = -1;
   tree<linked_residue_t_gemmi>::iterator it;
   for (it = glyco_tree.begin(); it != glyco_tree.end(); ++it) {
      if (it->residue == residue_p) {
         level = 0;
         tree<linked_residue_t_gemmi>::iterator this_one = it;
         bool has_parent = true;
         while (has_parent) {
            if (!this_one.node->parent) {
               has_parent = false;
            } else {
               this_one = this_one.node->parent;
               level++;
            }
         }
      }
   }
   return level;
}

} // namespace coot