#pragma once

#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/unitcell.hpp>
#include <gemmi/math.hpp>
#include <string>
#include <vector>

#include "/lmb/home/jdialpuri/Development/coot-tooling/generated-tests/coot__minimol__atom__atom/gemmi/function.hh"

namespace coot {

struct node_info {
   short int symm_trans_needed_flag;
   int index;
   node_info(int index_in) {
      symm_trans_needed_flag = 0;
      index = index_in;
   }
   node_info() : symm_trans_needed_flag(0), index(0) {}
};

class atom_graph {
public:
   std::vector<coot::minimol::atom> atoms;
   std::vector<coot::minimol::atom> atom_info;
   std::vector<std::vector<coot::node_info> > nodes;
   std::vector<double> cell;
   std::string spgr;

   static void atom_graph_gemmi(gemmi::Structure* st,
                                const std::vector<std::vector<coot::node_info>> &connection_indices,
                                const std::vector<gemmi::Vec3> &coords,
                                coot::atom_graph &ag);
};

} // namespace coot