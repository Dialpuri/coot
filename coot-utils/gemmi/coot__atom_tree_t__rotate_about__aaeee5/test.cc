#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/model.hpp>
#include <vector>
#include <string>
#include <queue>
#include "coot-utils/atom-vertex.hh"
#include "coot-utils/map-index.hh"
#include "function.hh"

// Inline copy of fill_atom_vertex_vec logic (from the existing gemmi port)
static bool fill_atom_vertex_vec_using_contacts(
    const std::vector<std::vector<int> >& contact_indices,
    int base_atom_index,
    const std::vector<gemmi::Atom>& residue_atoms,
    std::vector<coot::atom_vertex>& atom_vertex_vec) {

  bool r = false;
  coot::atom_vertex av;

  atom_vertex_vec.resize(residue_atoms.size());
  av.connection_type = coot::atom_vertex::START;
  atom_vertex_vec[base_atom_index] = av;

  if (contact_indices.size() == 0) {
    atom_vertex_vec.clear();
    return 0;
  }

  std::queue<int> q;
  q.push(base_atom_index);
  std::vector<int> done;

  while (!q.empty()) {
    int this_base_atom = q.front();
    std::vector<int> av_contacts = contact_indices[this_base_atom];

    for (size_t iav = 0; iav < av_contacts.size(); iav++) {
      int i_forward = av_contacts[iav];

      bool ifound_forward = false;
      for (size_t ifo = 0; ifo < atom_vertex_vec[this_base_atom].forward.size(); ifo++) {
        if (atom_vertex_vec[this_base_atom].forward[ifo] == static_cast<size_t>(av_contacts[iav])) {
          ifound_forward = true;
          break;
        }
      }

      if (!ifound_forward) {
        bool ifound_forward_forward = false;
        for (size_t ifo = 0; ifo < atom_vertex_vec[i_forward].forward.size(); ifo++) {
          if (atom_vertex_vec[i_forward].forward[ifo] == static_cast<size_t>(this_base_atom)) {
            ifound_forward_forward = true;
            break;
          }
        }
        if (!ifound_forward_forward) {
          atom_vertex_vec[this_base_atom].forward.push_back(av_contacts[iav]);
        }
      }

      bool in_done = false;
      for (size_t idone = 0; idone < done.size(); idone++) {
        if (done[idone] == static_cast<size_t>(av_contacts[iav])) {
          in_done = true;
          break;
        }
      }
      if (!in_done)
        q.push(av_contacts[iav]);
    }

    for (size_t iav = 0; iav < av_contacts.size(); iav++) {
      if (atom_vertex_vec[av_contacts[iav]].backward.size() == 0) {
        if (atom_vertex_vec[av_contacts[iav]].connection_type != coot::atom_vertex::START)
          atom_vertex_vec[av_contacts[iav]].backward.push_back(this_base_atom);
      }
    }
    q.pop();
    done.push_back(this_base_atom);
    r = true;
  }

  return r;
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(OracleTest, rotate_about) {
    gemmi::Structure st = gemmi::read_structure_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
    ASSERT_FALSE(st.models.empty());

    gemmi::Model& model = st.models[0];
    gemmi::Residue* res = nullptr;

    for (auto& chain : model.chains) {
        if (chain.name == "A") {
            for (auto& r : chain.residues) {
                if (r.seqid.num.value == 1299 && r.name == "LZA") {
                    res = &r;
                    break;
                }
            }
        }
    }
    ASSERT_NE(res, nullptr);

    EXPECT_EQ(std::string(res->name), "LZA");
    EXPECT_EQ(res->seqid.num.value, 1299);

    int n_atoms = static_cast<int>(res->atoms.size());
    EXPECT_EQ(n_atoms, 25);

    // Build atom_selection from residue atoms
    std::vector<const gemmi::Atom*> atom_selection;
    for (auto& atom : res->atoms) {
        atom_selection.push_back(&atom);
    }

    // Build a star contact graph: atom 0 is the base, all others connect to it
    std::vector<std::vector<int>> contact_indices(n_atoms);
    int base_atom_index = 0;
    for (int i = 1; i < n_atoms; i++) {
        contact_indices[0].push_back(i);
        contact_indices[i].push_back(0);
    }

    // Build atom_vertex_vec inline
    std::vector<coot::atom_vertex> atom_vertex_vec;
    bool filled = fill_atom_vertex_vec_using_contacts(
        contact_indices, base_atom_index, res->atoms, atom_vertex_vec);
    ASSERT_TRUE(filled);

    // Case 1: rotate using int indices (0 and 1)
    {
        int index2 = 0;
        int index3 = 1;
        double angle = 30.0;
        bool reversed = false;

        double new_torsion = coot::atom_tree_t_rotate_about_gemmi(
            atom_vertex_vec, atom_selection, res,
            index2, index3, angle, reversed);
        EXPECT_NEAR(new_torsion, 0.0, 1e-4);
    }

    // Case 2: edge case — invalid index2 == -1 triggers exception
    {
        int index2 = -1;
        int index3 = 1;
        double angle = 45.0;
        bool reversed = false;

        EXPECT_THROW({
            coot::atom_tree_t_rotate_about_gemmi(
                atom_vertex_vec, atom_selection, res,
                index2, index3, angle, reversed);
        }, std::runtime_error);
    }

    // Case 3: edge case — invalid index3 == -1 triggers exception
    {
        int index2 = 0;
        int index3 = -1;
        double angle = 45.0;
        bool reversed = false;

        EXPECT_THROW({
            coot::atom_tree_t_rotate_about_gemmi(
                atom_vertex_vec, atom_selection, res,
                index2, index3, angle, reversed);
        }, std::runtime_error);
    }
}