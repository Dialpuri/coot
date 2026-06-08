#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <vector>
#include <utility>
#include <string>
#include "function.hh"
#include "geometry/protein-geometry.hh"

// Helper: find atom by chain, residue seqnum, and atom name in a gemmi::Structure
static gemmi::Atom* find_atom(gemmi::Structure& st, const std::string& chain_name, int res_seq_num, const std::string& atom_name) {
    for (auto& model : st.models) {
        for (auto& chain : model.chains) {
            if (chain.name != chain_name) continue;
            for (auto& residue : chain.residues) {
                if (residue.seqid.num.value != res_seq_num) continue;
                for (auto& atom : residue.atoms) {
                    // Trim whitespace from atom name
                    std::string trimmed = atom.name;
                    auto start = trimmed.find_first_not_of(' ');
                    if (start != std::string::npos) {
                        trimmed = trimmed.substr(start);
                    }
                    while (!trimmed.empty() && trimmed.back() == ' ') {
                        trimmed.pop_back();
                    }
                    if (trimmed == atom_name) {
                        return &atom;
                    }
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, contact_info) {
    // Load molecule
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_GE(st.models.size(), 1u);

    gemmi::Model& model = st.models[0];
    int imol = 0;

    // Init protein geometry
    coot::protein_geometry geom;
    geom.init_standard();

    // Build atom selection with atoms from residue 10 (GLY): N, CA, C, O
    std::vector<gemmi::Atom*> atoms_vec;
    std::vector<std::string> atom_names = {"N", "CA", "C", "O"};

    for (const auto& atom_name : atom_names) {
        gemmi::Atom* atom = find_atom(st, "A", 10, atom_name);
        if (atom) {
            atoms_vec.push_back(atom);
        }
    }

    // Case 1: Valid selection with atoms from a single residue (GLY 10)
    {
        EXPECT_EQ(static_cast<int>(atoms_vec.size()), 4);
        EXPECT_EQ(imol, 0);

        std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>> link_bond_atoms;
        EXPECT_EQ(link_bond_atoms.size(), 0u);

        coot::contact_info ci(model, atoms_vec, imol, &geom, link_bond_atoms);

        EXPECT_EQ(ci.n_contacts(), 6);
    }

    // Case 2: Empty atom selection (edge case)
    {
        std::vector<gemmi::Atom*> atoms_vec_empty;

        EXPECT_EQ(atoms_vec_empty.size(), 0u);

        std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>> link_bond_atoms;

        coot::contact_info ci(model, atoms_vec_empty, imol, &geom, link_bond_atoms);

        EXPECT_EQ(ci.n_contacts(), 0);
    }

    // Case 3: Add a link bond between atoms in the selection
    {
        std::vector<std::pair<gemmi::Atom*, gemmi::Atom*>> link_bond_atoms;
        if (atoms_vec.size() >= 2) {
            // Add link bond between N and CA of residue 10
            link_bond_atoms.push_back(std::make_pair(atoms_vec[0], atoms_vec[1]));
        }

        EXPECT_EQ(link_bond_atoms.size(), 1u);

        coot::contact_info ci(model, atoms_vec, imol, &geom, link_bond_atoms);

        EXPECT_EQ(ci.n_contacts(), 7);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}