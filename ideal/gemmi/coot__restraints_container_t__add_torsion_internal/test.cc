#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <iostream>
#include "function.hh"

TEST(OracleTest, add_torsion_internal) {
    // Load structure via gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    // Build global atom index map
    std::map<gemmi::Atom*, int> atom_index_map;
    int idx = 0;
    for (gemmi::Model& model : st.models) {
        for (gemmi::Chain& chain : model.chains) {
            for (gemmi::Residue& r : chain.residues) {
                for (gemmi::Atom& atom : r.atoms) {
                    atom_index_map[&atom] = idx++;
                }
            }
        }
    }

    // Get chain A, residue at index 9 (10th residue)
    gemmi::Chain* gemmi_chain = st.models[0].find_chain("A");
    ASSERT_NE(gemmi_chain, nullptr);
    EXPECT_LT(9u, gemmi_chain->residues.size());

    gemmi::Residue& gemmi_res = gemmi_chain->residues[9];

    // Debug: print atom names
    std::cerr << "Residue " << gemmi_res.name << " (seqid " << gemmi_res.seqid.num.value << "): "
              << gemmi_res.atoms.size() << " atoms\n";
    for (const auto& atom : gemmi_res.atoms) {
        std::cerr << "  atom.name = [" << atom.name << "] (len=" << atom.name.size() << ")\n";
    }

    // Verify atom count
    EXPECT_EQ(gemmi_res.atoms.size(), 19u);

    std::vector<gemmi::Atom*> gemmi_res_atoms;
    for (gemmi::Atom& atom : gemmi_res.atoms) {
        gemmi_res_atoms.push_back(&atom);
    }

    // Setup standalone restraints container
    coot::restraints_container_t restraints;

    // Case 1: Valid torsion (N-CA-C-CB)
    {
        size_t n_restraints_before = restraints.restraints_vec.size();
        EXPECT_EQ(n_restraints_before, 0u);

        coot::dict_torsion_restraint_t torsion("chi1", " N  ", " CA ", " C  ", " CB ", 0.0, 10.0, 1);

        bool result = restraints.add_torsion_internal_gemmi(torsion, gemmi_res_atoms, atom_index_map, 1.0);

        size_t n_restraints_after = restraints.restraints_vec.size();
        EXPECT_EQ(result, true);
        EXPECT_EQ(n_restraints_after, 1u);
    }

    // Case 2: Zero periodicity (should fail)
    {
        size_t n_restraints_before = restraints.restraints_vec.size();
        EXPECT_EQ(n_restraints_before, 1u);

        coot::dict_torsion_restraint_t torsion("zero_period", " N  ", " CA ", " C  ", " CB ", 0.0, 10.0, 0);

        bool result = restraints.add_torsion_internal_gemmi(torsion, gemmi_res_atoms, atom_index_map, 1.0);

        size_t n_restraints_after = restraints.restraints_vec.size();
        EXPECT_EQ(result, false);
        EXPECT_EQ(n_restraints_after, 1u);
    }

    // Case 3: Zero ESD (should fail)
    {
        size_t n_restraints_before = restraints.restraints_vec.size();
        EXPECT_EQ(n_restraints_before, 1u);

        coot::dict_torsion_restraint_t torsion("zero_esd", " N  ", " CA ", " C  ", " CB ", 0.0, 0.0, 1);

        bool result = restraints.add_torsion_internal_gemmi(torsion, gemmi_res_atoms, atom_index_map, 1.0);

        size_t n_restraints_after = restraints.restraints_vec.size();
        EXPECT_EQ(result, false);
        EXPECT_EQ(n_restraints_after, 1u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
