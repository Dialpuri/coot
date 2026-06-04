#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "geometry/protein-geometry.hh"

// Suppress INFO output (goes to cout/cerr) during test
static class output_silencer {
public:
    output_silencer() {
        saved_cout_ = std::cout.rdbuf(null_stream_.rdbuf());
        saved_cerr_ = std::cerr.rdbuf(null_stream_.rdbuf());
    }
    ~output_silencer() {
        std::cout.rdbuf(saved_cout_);
        std::cerr.rdbuf(saved_cerr_);
    }
private:
    std::streambuf *saved_cout_;
    std::streambuf *saved_cerr_;
    std::stringstream null_stream_;
} global_silencer;

TEST(OracleTest, get_acedrg_types_for_residue) {
    coot::protein_geometry geom;
    geom.init_standard();

    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-dev/coot/canyon/test/canyon-test.pdb");
    gemmi::Model& model = st.models[0];

    int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;

    // Case 1: Valid protein residue (ASP at chain A, residue 1) — has monomer restraints
    {
        const gemmi::Chain* chain = model.find_chain("A");
        ASSERT_NE(chain, nullptr);

        const gemmi::Residue* residue_p = nullptr;
        for (auto& res : chain->residues) {
            if (res.seqid.num.value == 1 && res.name == "ASP") {
                residue_p = &res;
                break;
            }
        }
        ASSERT_NE(residue_p, nullptr);
        EXPECT_STREQ(residue_p->name.c_str(), "ASP");

        coot::acedrg_types_for_residue_t types =
            coot::get_acedrg_types_for_residue_gemmi(*residue_p, imol_enc, geom);

        EXPECT_EQ(types.bond_types.size(), 7u);

        // Bond 0: N -- CA
        {
            const auto &bt = types.bond_types[0];
            EXPECT_EQ(bt.atom_id_1, " N  ");
            EXPECT_EQ(bt.atom_id_2, " CA ");
            EXPECT_EQ(bt.atom_type_1, "N(CCCH)(H)3");
            EXPECT_EQ(bt.atom_type_2, "C(CCHH)(NH3)(COO)(H)");
            EXPECT_NEAR(bt.bond_length, 1.49381, 1e-4);
            EXPECT_EQ(bt.bond_is_between_atoms_in_the_same_ring, 0);
        }

        // Bond 1: CA -- C
        {
            const auto &bt = types.bond_types[1];
            EXPECT_EQ(bt.atom_id_1, " CA ");
            EXPECT_EQ(bt.atom_id_2, " C  ");
            EXPECT_EQ(bt.atom_type_1, "C(CCHH)(NH3)(COO)(H)");
            EXPECT_EQ(bt.atom_type_2, "C(CCHN)(O)2");
            EXPECT_NEAR(bt.bond_length, 1.4933, 1e-4);
            EXPECT_EQ(bt.bond_is_between_atoms_in_the_same_ring, 0);
        }

        // Bond 2: CA -- CB
        {
            const auto &bt = types.bond_types[2];
            EXPECT_EQ(bt.atom_id_1, " CA ");
            EXPECT_EQ(bt.atom_id_2, " CB ");
            EXPECT_EQ(bt.atom_type_1, "C(CCHH)(NH3)(COO)(H)");
            EXPECT_EQ(bt.atom_type_2, "C(CCHN)(COO)(H)2");
            EXPECT_NEAR(bt.bond_length, 1.52577, 1e-4);
            EXPECT_EQ(bt.bond_is_between_atoms_in_the_same_ring, 0);
        }

        // Bond 3: C -- O
        {
            const auto &bt = types.bond_types[3];
            EXPECT_EQ(bt.atom_id_1, " C  ");
            EXPECT_EQ(bt.atom_id_2, " O  ");
            EXPECT_EQ(bt.atom_type_1, "C(CCHN)(O)2");
            EXPECT_EQ(bt.atom_type_2, "O(CCO)");
            EXPECT_NEAR(bt.bond_length, 1.21442, 1e-4);
            EXPECT_EQ(bt.bond_is_between_atoms_in_the_same_ring, 0);
        }

        // Bond 4: CB -- CG
        {
            const auto &bt = types.bond_types[4];
            EXPECT_EQ(bt.atom_id_1, " CB ");
            EXPECT_EQ(bt.atom_id_2, " CG ");
            EXPECT_EQ(bt.atom_type_1, "C(CCHN)(COO)(H)2");
            EXPECT_EQ(bt.atom_type_2, "C(CCHH)(O)2");
            EXPECT_NEAR(bt.bond_length, 1.51298, 1e-4);
            EXPECT_EQ(bt.bond_is_between_atoms_in_the_same_ring, 0);
        }

        // Bond 5: CG -- OD1
        {
            const auto &bt = types.bond_types[5];
            EXPECT_EQ(bt.atom_id_1, " CG ");
            EXPECT_EQ(bt.atom_id_2, " OD1");
            EXPECT_EQ(bt.atom_type_1, "C(CCHH)(O)2");
            EXPECT_EQ(bt.atom_type_2, "O(CCO)");
            EXPECT_NEAR(bt.bond_length, 1.25169, 1e-4);
            EXPECT_EQ(bt.bond_is_between_atoms_in_the_same_ring, 0);
        }

        // Bond 6: CG -- OD2
        {
            const auto &bt = types.bond_types[6];
            EXPECT_EQ(bt.atom_id_1, " CG ");
            EXPECT_EQ(bt.atom_id_2, " OD2");
            EXPECT_EQ(bt.atom_type_1, "C(CCHH)(O)2");
            EXPECT_EQ(bt.atom_type_2, "O(CCO)");
            EXPECT_NEAR(bt.bond_length, 1.2606, 1e-4);
            EXPECT_EQ(bt.bond_is_between_atoms_in_the_same_ring, 0);
        }
    }

    // Case 2: Non-existent residue — verifies the nullptr guard
    {
        const gemmi::Chain* chain = model.find_chain("A");
        ASSERT_NE(chain, nullptr);

        const gemmi::Residue* residue_p = nullptr;
        for (auto& res : chain->residues) {
            if (res.seqid.num.value == 9999 && res.name == "NONE") {
                residue_p = &res;
                break;
            }
        }
        EXPECT_EQ(residue_p, nullptr);
    }

    // Case 3: Residue with no dictionary entry (unknown residue type)
    {
        const gemmi::Chain* chain = model.find_chain("A");
        ASSERT_NE(chain, nullptr);

        const gemmi::Residue* residue_p = nullptr;
        for (auto& res : chain->residues) {
            if (res.seqid.num.value == 1 && res.name == "UNX") {
                residue_p = &res;
                break;
            }
        }
        EXPECT_EQ(residue_p, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}