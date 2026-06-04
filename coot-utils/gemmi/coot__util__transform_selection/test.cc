#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, transform_selection) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    // Case 1: transform atoms in chain A with a (5,5,5) translation
    {
        // Count selected atoms in chain A
        int n_selected = 0;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& res : chain.residues) {
                    n_selected += static_cast<int>(res.atoms.size());
                }
            }
        }
        EXPECT_EQ(n_selected, 2106);

        // Find first atom in chain A
        gemmi::Atom* first_atom = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "A") {
                for (gemmi::Residue& res : chain.residues) {
                    if (!res.atoms.empty()) {
                        first_atom = &res.atoms[0];
                        break;
                    }
                }
                if (first_atom) break;
            }
        }

        double x_before = 0, y_before = 0, z_before = 0;
        if (first_atom) {
            x_before = first_atom->pos.x;
            y_before = first_atom->pos.y;
            z_before = first_atom->pos.z;
        }

        EXPECT_NEAR(x_before, 12.334, 1e-4);
        EXPECT_NEAR(y_before, 2.772, 1e-4);
        EXPECT_NEAR(z_before, 74.951, 1e-4);

        clipper::Mat33<> id_mat(1, 0, 0, 0, 1, 0, 0, 0, 1);
        clipper::Vec3<> trans(5.0, 5.0, 5.0);
        clipper::RTop_orth rtop(id_mat, trans);

        coot::util::transform_selection_gemmi(model, "A", rtop);

        double x_after = 0, y_after = 0, z_after = 0;
        if (first_atom) {
            x_after = first_atom->pos.x;
            y_after = first_atom->pos.y;
            z_after = first_atom->pos.z;
        }

        EXPECT_NEAR(x_after, 17.334, 1e-4);
        EXPECT_NEAR(y_after, 7.772, 1e-4);
        EXPECT_NEAR(z_after, 79.951, 1e-4);
        EXPECT_NEAR((x_after - x_before), 5.0, 1e-4);
        EXPECT_NEAR((y_after - y_before), 5.0, 1e-4);
        EXPECT_NEAR((z_after - z_before), 5.0, 1e-4);
    }

    // Case 2: identity transform (no-op) on chain B
    {
        int n_selected = 0;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "B") {
                for (gemmi::Residue& res : chain.residues) {
                    n_selected += static_cast<int>(res.atoms.size());
                }
            }
        }
        EXPECT_EQ(n_selected, 1);

        gemmi::Atom* first_atom = nullptr;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "B") {
                for (gemmi::Residue& res : chain.residues) {
                    if (!res.atoms.empty()) {
                        first_atom = &res.atoms[0];
                        break;
                    }
                }
                if (first_atom) break;
            }
        }

        double x_before = 0;
        if (first_atom) {
            x_before = first_atom->pos.x;
        }

        EXPECT_NEAR(x_before, 30.744, 1e-4);

        clipper::RTop_orth rtop = clipper::RTop_orth().identity();

        coot::util::transform_selection_gemmi(model, "B", rtop);

        double x_after = 0;
        if (first_atom) {
            x_after = first_atom->pos.x;
        }

        EXPECT_NEAR(x_after, 30.744, 1e-4);
        EXPECT_NEAR((x_after - x_before), 0.0, 1e-4);
    }

    // Case 3: empty selection (invalid chain)
    {
        int n_selected = 0;
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name == "Z") {
                for (gemmi::Residue& res : chain.residues) {
                    n_selected += static_cast<int>(res.atoms.size());
                }
            }
        }
        EXPECT_EQ(n_selected, 0);

        clipper::Mat33<> id_mat(1, 0, 0, 0, 1, 0, 0, 0, 1);
        clipper::Vec3<> trans(1.0, 0.0, 0.0);
        clipper::RTop_orth rtop(id_mat, trans);

        EXPECT_NO_THROW(coot::util::transform_selection_gemmi(model, "Z", rtop));

        EXPECT_EQ(n_selected, 0);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}