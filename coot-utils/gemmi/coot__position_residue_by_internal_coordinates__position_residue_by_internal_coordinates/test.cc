#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find atom by name in a gemmi::Residue
static gemmi::Atom *find_atom(gemmi::Residue *res, const char *name) {
    if (!res) return nullptr;
    std::string target(name);
    for (gemmi::Atom& a : res->atoms) {
        std::string aname;
        for (char c : a.name) {
            if (c != ' ') aname += c;
        }
        if (aname == target) return &a;
    }
    return nullptr;
}

TEST(OracleTest, position_residue_by_internal_coordinates) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_GT(st.models.size(), 0);

    gemmi::Model &model = st.models[0];
    ASSERT_GT(model.chains.size(), 0);
    gemmi::Chain &chain = model.chains[0];
    ASSERT_GT(chain.residues.size(), 1);
    gemmi::Residue &res_ref = chain.residues[0];
    gemmi::Residue &res_mov = chain.residues[1];

    // Assert input metadata from oracle
    EXPECT_EQ(std::string(res_ref.name), "MET");
    EXPECT_EQ(res_ref.seqid.num.value, 1);
    EXPECT_EQ(std::string(res_mov.name), "GLU");
    EXPECT_EQ(res_mov.seqid.num.value, 2);

    // Record moving residue atom count
    int mov_atom_count = static_cast<int>(res_mov.atoms.size());

    // ---- Case 1: Valid quad (N,CA,C,CB) ----
    {
        coot::atom_name_quad quad("N", "CA", "C", "CB");

        // Positions BEFORE
        gemmi::Atom *cb = find_atom(&res_mov, "CB");
        double xb = cb ? cb->pos.x : 0, yb = cb ? cb->pos.y : 0, zb = cb ? cb->pos.z : 0;

        coot::position_residue_by_internal_coordinates_gemmi pric(
            &res_ref, &res_mov, quad,
            1.45, 109.5, 60.0);

        bool ok = pric.move_moving_residue();

        // Positions AFTER
        cb = find_atom(&res_mov, "CB");
        double xa = cb ? cb->pos.x : 0, ya = cb ? cb->pos.y : 0, za = cb ? cb->pos.z : 0;

        EXPECT_FALSE(ok);       // oracle: move_ok: false
        EXPECT_EQ(mov_atom_count, 9);
        EXPECT_EQ(res_mov.atoms.size(), 9);
        EXPECT_NEAR(xb, 12.357, 1e-3);
        EXPECT_NEAR(yb, 4.638, 1e-3);
        EXPECT_NEAR(zb, 79.388, 1e-3);
        EXPECT_NEAR(xa, 12.357, 1e-3);
        EXPECT_NEAR(ya, 4.638, 1e-3);
        EXPECT_NEAR(za, 79.388, 1e-3);
    }

    // ---- Case 2: Bad quad (nonexistent atom names) ----
    {
        coot::atom_name_quad bq("OOO", "XXX", "YYY", "ZZZ");
        coot::position_residue_by_internal_coordinates_gemmi pric2(
            &res_ref, &res_mov, bq, 1.45, 109.5, 60.0);
        bool ok2 = pric2.move_moving_residue();
        EXPECT_FALSE(ok2);   // oracle: bad_move_ok: false
    }

    // ---- Case 3: Null ref residue ----
    {
        coot::atom_name_quad q3("N", "CA", "C", "CB");
        coot::position_residue_by_internal_coordinates_gemmi pric3(
            nullptr, &res_mov, q3, 1.45, 109.5, 60.0);
        bool ok3 = pric3.move_moving_residue();
        EXPECT_FALSE(ok3);   // oracle: null_ref_move_ok: false
    }

    // ---- Case 4: Null moving residue ----
    {
        coot::atom_name_quad q4("N", "CA", "C", "CB");
        coot::position_residue_by_internal_coordinates_gemmi pric4(
            &res_ref, nullptr, q4, 1.45, 109.5, 60.0);
        bool ok4 = pric4.move_moving_residue();
        EXPECT_FALSE(ok4);   // oracle: null_mov_move_ok: false
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}