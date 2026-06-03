#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, atom_index_quad_torsion) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model& model = st.models[0];
    gemmi::Chain& chain = model.chains[0];

    // Find VAL 7
    const gemmi::Residue* val7 = nullptr;
    for (gemmi::Residue& res : chain.residues) {
        if (res.seqid.num.value == 7 && res.name == "VAL") {
            val7 = &res;
            break;
        }
    }
    ASSERT_NE(val7, nullptr);

    // Build atom table
    std::vector<gemmi::Atom> atoms = val7->atoms;
    int n_atoms = static_cast<int>(atoms.size());

    int idx_N = -1, idx_CA = -1, idx_C = -1, idx_O = -1;
    for (int i = 0; i < n_atoms; i++) {
        std::string aname = atoms[i].name;
        aname.erase(0, aname.find_first_not_of(" "));
        if (!aname.empty())
            aname.erase(aname.find_last_not_of(" ") + 1);
        if (aname == "N")  idx_N = i;
        if (aname == "CA") idx_CA = i;
        if (aname == "C")  idx_C = i;
        if (aname == "O")  idx_O = i;
    }

    // Case 1: valid backbone omega-like torsion N-CA-C-O
    {
        EXPECT_EQ(idx_N, 0);
        EXPECT_EQ(idx_CA, 1);
        EXPECT_EQ(idx_C, 2);
        EXPECT_EQ(idx_O, 3);

        coot::atom_index_quad quad(idx_N, idx_CA, idx_C, idx_O);
        double torsion = quad.torsion_gemmi(atoms, n_atoms);
        EXPECT_NEAR(torsion, 134.24, 1e-4);
    }

    // Case 2: chi1 torsion CA-CB-CG1-CG2
    {
        int idx_CB = -1, idx_CG1 = -1, idx_CG2 = -1;
        for (int i = 0; i < n_atoms; i++) {
            std::string aname = atoms[i].name;
            aname.erase(0, aname.find_first_not_of(" "));
            if (!aname.empty())
                aname.erase(aname.find_last_not_of(" ") + 1);
            if (aname == "CB")  idx_CB = i;
            if (aname == "CG1") idx_CG1 = i;
            if (aname == "CG2") idx_CG2 = i;
        }
        EXPECT_EQ(idx_CA, 1);
        EXPECT_EQ(idx_CB, 4);
        EXPECT_EQ(idx_CG1, 5);
        EXPECT_EQ(idx_CG2, 6);

        coot::atom_index_quad quad(idx_CA, idx_CB, idx_CG1, idx_CG2);
        double torsion = quad.torsion_gemmi(atoms, n_atoms);
        EXPECT_NEAR(torsion, 125.753, 1e-3);
    }

    // Case 3: bad index — should throw
    {
        coot::atom_index_quad quad(0, 1, 2, 999);
        EXPECT_THROW(quad.torsion_gemmi(atoms, n_atoms), std::runtime_error);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}