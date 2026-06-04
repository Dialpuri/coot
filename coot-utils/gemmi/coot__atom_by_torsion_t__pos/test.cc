#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "coot-utils/glyco-torsions.hh"
#include <utility>

TEST(OracleTest, atom_by_torsion_t_pos) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    gemmi::Model* model = &st.models[0];
    ASSERT_NE(model, nullptr);

    gemmi::Chain* chain = nullptr;
    for (auto& ch : model->chains) {
        if (ch.name == "A") {
            chain = &ch;
            break;
        }
    }
    ASSERT_NE(chain, nullptr);

    const gemmi::Residue* base_res = nullptr;
    const gemmi::Residue* ext_res = nullptr;
    for (auto& res : chain->residues) {
        if (res.seqid.num.value == 1) base_res = &res;
        if (res.seqid.num.value == 2) ext_res = &res;
    }

    ASSERT_NE(base_res, nullptr);
    ASSERT_NE(ext_res, nullptr);

    // Case 1: Valid atom_by_torsion_t — all 3 prior atoms in base residue
    {
        coot::atom_by_torsion_t atb(
            coot::atom_by_torsion_base_t(
                "OXT",          // atom_name
                "O",            // element
                std::make_pair(true, "CA"),   // prior_atom_1: from base
                std::make_pair(true, "N"),    // prior_atom_2: from base
                std::make_pair(true, "C")     // prior_atom_3: from base
            ),
            1.23,   // bond_length
            120.0,  // angle (degrees)
            180.0   // torsion (degrees)
        );

        EXPECT_EQ(atb.atom_name, "OXT");
        EXPECT_EQ(atb.element, "O");
        EXPECT_EQ(atb.prior_atom_1.first, true);
        EXPECT_EQ(atb.prior_atom_1.second, "CA");
        EXPECT_EQ(atb.prior_atom_2.first, true);
        EXPECT_EQ(atb.prior_atom_2.second, "N");
        EXPECT_EQ(atb.prior_atom_3.first, true);
        EXPECT_EQ(atb.prior_atom_3.second, "C");

        clipper::Coord_orth pos = coot::pos_gemmi(
            atb.atom_name,
            atb.prior_atom_1,
            atb.prior_atom_2,
            atb.prior_atom_3,
            *base_res, *ext_res,
            1.23,    // bond_length
            120.0,   // angle (degrees)
            180.0);  // torsion (degrees)
        EXPECT_NEAR(pos.x(), 13.6869, 1e-4);
        EXPECT_NEAR(pos.y(), 4.10212, 1e-4);
        EXPECT_NEAR(pos.z(), 73.5509, 1e-4);
    }

    // Case 2: prior_atom from ext_residue (first=false)
    {
        coot::atom_by_torsion_t atb(
            coot::atom_by_torsion_base_t(
                "CB",
                "C",
                std::make_pair(false, "CA"),  // from ext_residue
                std::make_pair(true, "N"),    // from base_residue
                std::make_pair(true, "C")     // from base_residue
            ),
            1.54,   // bond_length
            109.5,  // angle
            60.0    // torsion
        );

        EXPECT_EQ(atb.atom_name, "CB");
        EXPECT_EQ(atb.prior_atom_1.first, false);
        EXPECT_EQ(atb.prior_atom_1.second, "CA");

        clipper::Coord_orth pos = coot::pos_gemmi(
            atb.atom_name,
            atb.prior_atom_1,
            atb.prior_atom_2,
            atb.prior_atom_3,
            *base_res, *ext_res,
            1.54,   // bond_length
            109.5,  // angle (degrees)
            60.0);  // torsion (degrees)
        EXPECT_NEAR(pos.x(), 13.5573, 1e-4);
        EXPECT_NEAR(pos.y(), 6.18142, 1e-4);
        EXPECT_NEAR(pos.z(), 78.2036, 1e-4);
    }

    // Case 3: Missing atom triggers exception
    {
        coot::atom_by_torsion_t atb(
            coot::atom_by_torsion_base_t(
                "NZ",
                "N",
                std::make_pair(true, "CA"),
                std::make_pair(true, "N"),
                std::make_pair(true, "NONEXISTENT")  // this atom doesn't exist
            ),
            1.47,   // bond_length
            110.0,  // angle
            120.0   // torsion
        );

        EXPECT_EQ(atb.atom_name, "NZ");
        EXPECT_EQ(atb.prior_atom_3.first, true);
        EXPECT_EQ(atb.prior_atom_3.second, "NONEXISTENT");

        EXPECT_THROW(coot::pos_gemmi(
            atb.atom_name,
            atb.prior_atom_1,
            atb.prior_atom_2,
            atb.prior_atom_3,
            *base_res, *ext_res,
            1.47,   // bond_length
            110.0,  // angle (degrees)
            120.0), // torsion (degrees)
            std::exception);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}