#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Helper: find atom index by name in residue atoms (unpadded gemmi names)
static int find_atom_index(const gemmi::Residue& res, const std::string& name) {
    for (int i = 0; i < static_cast<int>(res.atoms.size()); i++) {
        if (res.atoms[i].name == name) return i;
    }
    return -1;
}

// Build contact_indices for standard PHE connectivity
static std::vector<std::vector<int>> build_phe_contact_indices(
    const gemmi::Residue& res)
{
    std::map<std::string, int> idx;
    // Standard PHE atoms — use unpadded gemmi names
    std::vector<std::string> names = {
        "N", "CA", "C", "O", "CB", "CG",
        "CD1", "CD2", "CE1", "CE2", "CZ"
    };
    for (auto& nm : names) {
        int i = find_atom_index(res, nm);
        idx[nm] = i;
    }

    // Bonds: pairs of atom names
    std::vector<std::pair<std::string, std::string>> bonds = {
        {"N", "CA"},
        {"CA", "C"},
        {"CA", "CB"},
        {"C", "O"},
        {"CB", "CG"},
        {"CG", "CD1"},
        {"CG", "CD2"},
        {"CD1", "CE1"},
        {"CD2", "CE2"},
        {"CE1", "CZ"},
        {"CE2", "CZ"}
    };

    int n = static_cast<int>(res.atoms.size());
    std::vector<std::vector<int>> ci(n, std::vector<int>());
    for (auto& b : bonds) {
        int i1 = idx[b.first];
        int i2 = idx[b.second];
        if (i1 >= 0 && i2 >= 0) {
            ci[i1].push_back(i2);
            ci[i2].push_back(i1);
        }
    }
    return ci;
}

TEST(OracleTest, change_by_internal) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    ASSERT_GE(st.models.size(), 1);

    gemmi::Model& model = st.models[0];

    // Find a PHE residue
    gemmi::Residue* phe_res = nullptr;
    for (gemmi::Chain& chain : model.chains) {
        for (gemmi::Residue& res : chain.residues) {
            if (res.name == "PHE") {
                phe_res = &res;
                break;
            }
        }
        if (phe_res) break;
    }
    ASSERT_NE(phe_res, nullptr);

    // Build contact indices for PHE connectivity
    std::vector<std::vector<int>> contact_indices =
        build_phe_contact_indices(*phe_res);

    // Record original CZ coordinates (before rotation)
    int cz_idx = find_atom_index(*phe_res, "CZ");
    double cz_x_orig = 0, cz_y_orig = 0, cz_z_orig = 0;
    bool found_cz = (cz_idx >= 0);
    if (found_cz) {
        cz_x_orig = phe_res->atoms[cz_idx].pos.x;
        cz_y_orig = phe_res->atoms[cz_idx].pos.y;
        cz_z_orig = phe_res->atoms[cz_idx].pos.z;
    }

    // Build atom_name_pairs for PHE chi angles — use unpadded gemmi names
    std::vector<coot::chi_angles::atom_name_pair> atom_name_pairs;
    atom_name_pairs.push_back(coot::chi_angles::atom_name_pair("CB", "CG"));
    atom_name_pairs.push_back(coot::chi_angles::atom_name_pair("CG", "CD1"));
    atom_name_pairs.push_back(coot::chi_angles::atom_name_pair("CD1", "CE1"));
    atom_name_pairs.push_back(coot::chi_angles::atom_name_pair("CE1", "CZ"));

    // Tree base atom — use special "don't use" value
    coot::atom_spec_t tree_base_atom;
    tree_base_atom.atom_name = "CA";
    tree_base_atom.int_user_data = -999; // special "don't use" value

    // === Case 1: chi4 rotation by 60 degrees (should succeed) ===
    {
        int ichi = 4;
        double diff = 60.0;

        auto result = coot::chi_angles::change_by_internal_gemmi(
            ichi, diff, atom_name_pairs, contact_indices,
            *phe_res, tree_base_atom
        );

        EXPECT_EQ(result.first, 0);
        EXPECT_NEAR(result.second, 1.0315, 1e-4);

        // CZ coordinates after rotation
        if (found_cz) {
            double cz_x_after = phe_res->atoms[cz_idx].pos.x;
            EXPECT_NEAR(cz_x_orig, 18.2510, 1e-4);
            EXPECT_NEAR(cz_x_after, 18.2510, 1e-4);
        }
    }

    // === Case 2: invalid ichi (too large, should return error code 2) ===
    {
        int ichi = 99;
        double diff = 10.0;

        auto result = coot::chi_angles::change_by_internal_gemmi(
            ichi, diff, atom_name_pairs, contact_indices,
            *phe_res, tree_base_atom
        );

        EXPECT_EQ(result.first, 2);
        EXPECT_NEAR(result.second, 0.0, 1e-4);
    }

    // === Case 3: chi1 rotation (should cause actual mutation) ===
    // Addresses coverage warning: BEFORE == AFTER for all mutation observations
    {
        int cd1_idx = find_atom_index(*phe_res, "CD1");
        double cd1_x_before = 0;
        bool found_cd1 = (cd1_idx >= 0);
        if (found_cd1) {
            cd1_x_before = phe_res->atoms[cd1_idx].pos.x;
        }

        int ichi = 1;
        double diff = 60.0;

        auto result = coot::chi_angles::change_by_internal_gemmi(
            ichi, diff, atom_name_pairs, contact_indices,
            *phe_res, tree_base_atom
        );

        EXPECT_EQ(result.first, 0);

        // CD1 should have moved after chi1 rotation
        if (found_cd1) {
            double cd1_x_after = phe_res->atoms[cd1_idx].pos.x;
            // The x coordinate should be different after rotation
            EXPECT_NE(cd1_x_after, cd1_x_before);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}