#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

#include "geometry/protein-geometry.hh"

TEST(OracleTest, hydrogen_placement_by_dictionary) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    gemmi::Model& model = st.models[0];
    gemmi::Chain& chain = model.chains[0];

    // Case 1: non-HOH with geom_p == nullptr (warning path, no atoms added)
    {
        gemmi::Residue& res = chain.residues[0]; // MET 1
        int atoms_before = static_cast<int>(res.atoms.size());
        std::string res_name = res.name;
        int seqnum = res.seqid.num.value;

        EXPECT_EQ(res_name, "MET");
        EXPECT_EQ(seqnum, 1);
        EXPECT_EQ(atoms_before, 19);

        // No geometry — nothing to do (the outer function would warn)
        // Atoms should remain the same

        int atoms_after = static_cast<int>(res.atoms.size());
        EXPECT_EQ(atoms_after, 19);
    }

    // Case 2: non-HOH with geom_p != nullptr but empty geometry (no restraints found, no atoms added)
    {
        coot::protein_geometry pg;

        gemmi::Residue& res = chain.residues[1]; // GLU 2
        int atoms_before = static_cast<int>(res.atoms.size());
        std::string res_name = res.name;
        int seqnum = res.seqid.num.value;

        EXPECT_EQ(res_name, "GLU");
        EXPECT_EQ(seqnum, 2);
        EXPECT_EQ(atoms_before, 15);

        // Try to get restraints — if found, place hydrogens
        std::pair<bool, coot::dictionary_residue_restraints_t> p =
            pg.get_monomer_restraints(res_name, 0);
        if (p.first) {
            gemmi::CRA cra{&chain, &res, nullptr};
            coot::reduce::hydrogen_placement_by_dictionary_gemmi(
                p.second, cra, st, 1.09, 1.08, 1.01, 0.96, 1.34);
        }

        int atoms_after = static_cast<int>(res.atoms.size());
        EXPECT_EQ(atoms_after, 15);
    }

    // Case 3: non-HOH with geometry — another residue
    {
        gemmi::Residue& res = chain.residues[2]; // ASN 3
        int atoms_before = static_cast<int>(res.atoms.size());
        std::string res_name = res.name;

        EXPECT_EQ(res_name, "ASN");
        EXPECT_FALSE(res_name == "HOH");
        EXPECT_EQ(atoms_before, 14);

        coot::protein_geometry pg;
        std::pair<bool, coot::dictionary_residue_restraints_t> p =
            pg.get_monomer_restraints(res_name, 0);
        if (p.first) {
            gemmi::CRA cra{&chain, &res, nullptr};
            coot::reduce::hydrogen_placement_by_dictionary_gemmi(
                p.second, cra, st, 1.09, 1.08, 1.01, 0.96, 1.34);
        }

        int atoms_after = static_cast<int>(res.atoms.size());
        EXPECT_EQ(atoms_after, 14);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}