#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/mmread.hpp>
#include <vector>
#include <string>
#include "function.hh"

static gemmi::Residue* find_residue(gemmi::Model& model, const std::string& chain_id, int resno)
{
    for (auto& ch : model.chains) {
        if (ch.name != chain_id) continue;
        for (auto& res : ch.residues) {
            if (res.seqid.num.value == resno) return &res;
        }
    }
    return nullptr;
}

TEST(OracleTest, replace_coords) {
    gemmi::Structure st = gemmi::read_structure_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
    ASSERT_EQ(st.models.size(), 1u);

    // === Case 1: Replace coords — shift selected atoms by 10 Å ===
    {
        gemmi::Model& model = st.models[0];

        gemmi::Residue* res = find_residue(model, "A", 1299);
        ASSERT_NE(res, nullptr);

        int n_atoms_in_res = res->atoms.size();
        EXPECT_EQ(n_atoms_in_res, 25);

        // Store original coordinates
        std::vector<std::pair<size_t, double>> orig_coords;
        for (size_t ia = 0; ia < res->atoms.size(); ia++) {
            auto& at = res->atoms[ia];
            orig_coords.push_back({ia, at.pos.x});
        }

        // Verify first atom original coords
        EXPECT_EQ(res->atoms[0].name, "C26");
        EXPECT_NEAR(orig_coords[0].second, 31.8560, 1e-4);

        // Create a temporary chain/residue/atoms with shifted coords
        gemmi::Chain temp_chain;
        temp_chain.name = "A";
        gemmi::Residue temp_res;
        temp_res.seqid = res->seqid;
        for (auto& at : res->atoms) {
            gemmi::Atom a = at;
            a.pos.x += 10.0;
            a.pos.y += 10.0;
            a.pos.z += 10.0;
            temp_res.atoms.push_back(std::move(a));
        }
        temp_chain.residues.push_back(std::move(temp_res));

        // Build CRAs pointing to the temp atoms
        gemmi::Chain* t_ch = &temp_chain;
        gemmi::Residue* t_res = &temp_chain.residues[0];
        std::vector<gemmi::CRA> selection;
        for (auto& at : t_res->atoms) {
            selection.push_back(gemmi::CRA{t_ch, t_res, &at});
        }
        EXPECT_EQ(selection.size(), 25);

        // Call replace_coords_gemmi
        EXPECT_NO_THROW(coot::molecule_t::replace_coords_gemmi(model, selection, false, true));

        // Check results
        double new_x = res->atoms[0].pos.x;
        double expected_x = orig_coords[0].second + 10.0;
        EXPECT_NEAR(new_x, 41.8560, 1e-4);
        EXPECT_NEAR(expected_x, 41.8560, 1e-4);
    }

    // === Case 2: Empty ASC — no-op ===
    {
        gemmi::Model& model = st.models[0];

        gemmi::Residue* res = find_residue(model, "A", 1299);
        ASSERT_NE(res, nullptr);

        gemmi::Atom& target = res->atoms[0];
        double bx = target.pos.x;

        std::vector<gemmi::CRA> empty_selection;
        EXPECT_NO_THROW(coot::molecule_t::replace_coords_gemmi(model, empty_selection, false, true));

        double ax = target.pos.x;
        EXPECT_NEAR(bx, ax, 1e-4);
    }

    // === Case 3: Single atom replace ===
    {
        gemmi::Model& model = st.models[0];

        gemmi::Residue* res = find_residue(model, "A", 1);
        ASSERT_NE(res, nullptr);
        ASSERT_FALSE(res->atoms.empty());

        gemmi::Atom& target = res->atoms[0];
        double orig_x = target.pos.x;
        EXPECT_EQ(target.name, "N");
        EXPECT_NEAR(orig_x, 12.3420, 1e-4);

        // Create single-atom temp structure with shifted coords
        gemmi::Chain temp_chain;
        temp_chain.name = "A";
        gemmi::Residue temp_res;
        temp_res.seqid = res->seqid;
        gemmi::Atom temp_atom = target;
        temp_atom.pos.x += 5.0;
        temp_res.atoms.push_back(std::move(temp_atom));
        temp_chain.residues.push_back(std::move(temp_res));

        gemmi::Chain* t_ch = &temp_chain;
        gemmi::Residue* t_res = &temp_chain.residues[0];
        std::vector<gemmi::CRA> selection;
        selection.push_back(gemmi::CRA{t_ch, t_res, &t_res->atoms[0]});

        EXPECT_NO_THROW(coot::molecule_t::replace_coords_gemmi(model, selection, false, true));

        EXPECT_NEAR(target.pos.x, 17.3420, 1e-4);
        EXPECT_NEAR(target.pos.x, orig_x + 5.0, 1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}