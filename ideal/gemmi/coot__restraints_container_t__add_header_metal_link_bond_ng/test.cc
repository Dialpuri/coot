#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, add_header_metal_link_bond_ng) {
    // Case 1: Valid atoms
    {
        auto st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        auto& model = st.models[0];

        // Build atom entries with parent context
        std::vector<atom_entry> atoms;
        coot::atom_spec_t spec1;
        coot::atom_spec_t spec2;

        for (auto& chain : model.chains) {
            for (auto& res : chain.residues) {
                for (auto& atom : res.atoms) {
                    atoms.push_back({&atom, &res, &chain});
                    // gemmi atom names are NOT padded: "N", "CA" (not "N  ", "CA ")
                    if (chain.name == "A" && res.seqid.num.value == 1 && atom.name == "N") {
                        spec1 = coot::atom_spec_from_gemmi(atom, res, chain);
                    }
                    if (chain.name == "A" && res.seqid.num.value == 1 && atom.name == "CA") {
                        spec2 = coot::atom_spec_from_gemmi(atom, res, chain);
                    }
                }
            }
        }

        // Initialize bonded_atom_indices and restraints_vec
        std::vector<std::set<int>> bonded_atom_indices(atoms.size());
        std::vector<coot::simple_restraint> restraints_vec;

        double dist = 2.2;

        size_t bonds_before = bonded_atom_indices.size();
        size_t restraints_before = restraints_vec.size();

        bool added = coot::add_header_metal_link_bond_ng_gemmi(
            atoms, bonded_atom_indices, restraints_vec, spec1, spec2, dist);

        size_t bonds_after = bonded_atom_indices.size();
        size_t restraints_after = restraints_vec.size();

        EXPECT_TRUE(added);
        EXPECT_EQ(bonds_after, bonds_before);
        EXPECT_TRUE(restraints_after > restraints_before);
    }

    // Case 2: Invalid atom spec — atom that does not exist (res_no 9999)
    {
        auto st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        auto& model = st.models[0];

        std::vector<atom_entry> atoms;
        coot::atom_spec_t spec1;

        for (auto& chain : model.chains) {
            for (auto& res : chain.residues) {
                for (auto& atom : res.atoms) {
                    atoms.push_back({&atom, &res, &chain});
                    if (chain.name == "A" && res.seqid.num.value == 1 && atom.name == "N") {
                        spec1 = coot::atom_spec_from_gemmi(atom, res, chain);
                    }
                }
            }
        }

        std::vector<std::set<int>> bonded_atom_indices(atoms.size());
        std::vector<coot::simple_restraint> restraints_vec;

        // First add a valid bond
        coot::atom_spec_t spec2("A", 1, "", "CA", "");
        coot::add_header_metal_link_bond_ng_gemmi(
            atoms, bonded_atom_indices, restraints_vec, spec1, spec2, 2.2);

        // Now try with a bad spec
        coot::atom_spec_t bad_spec("A", 9999, "", "CA", "");

        size_t bonds_before2 = bonded_atom_indices.size();
        size_t restraints_before2 = restraints_vec.size();

        bool added = coot::add_header_metal_link_bond_ng_gemmi(
            atoms, bonded_atom_indices, restraints_vec, spec1, bad_spec, 2.0);

        size_t bonds_after2 = bonded_atom_indices.size();
        size_t restraints_after2 = restraints_vec.size();

        EXPECT_FALSE(added);
        EXPECT_EQ(bonds_after2, bonds_before2);
        EXPECT_EQ(restraints_after2, restraints_before2);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}