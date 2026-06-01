#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "geometry/protein-geometry.hh"
#include "coot-utils/coot-coord-extras.hh"

// Helper: find a residue by chain and sequence number in a gemmi Structure
static gemmi::Residue* find_residue(const gemmi::Structure& st, const std::string& chain_id, int seq_num) {
    for (const auto& model : st.models) {
        for (const auto& chain : model.chains) {
            if (chain.name == chain_id) {
                for (const auto& res : chain.residues) {
                    if (res.seqid.num.value == seq_num) {
                        return const_cast<gemmi::Residue*>(&res);
                    }
                }
            }
        }
    }
    return nullptr;
}

TEST(OracleTest, match_torsions_gemmi) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_FALSE(st.models.empty());

    // Case 1 (oracle): two residues — oracle did NOT call init_standard(), so restraints_found is false
    {
        coot::protein_geometry geom;
        // NOTE: intentionally NOT calling geom.init_standard() to match oracle behaviour

        gemmi::Residue *res_moving = find_residue(st, "A", 5);
        gemmi::Residue *res_ref    = find_residue(st, "A", 47);

        EXPECT_NE(res_moving, nullptr);
        EXPECT_NE(res_ref, nullptr);

        std::string res_name(res_moving->name);
        EXPECT_EQ(res_name, "GLN");

        std::pair<bool, coot::dictionary_residue_restraints_t> rest_info =
            geom.get_monomer_restraints(res_name, 0);

        EXPECT_FALSE(rest_info.first); // oracle: restraints_found: false

        // Test the ported match_torsions_gemmi function
        coot::match_torsions_data_t data =
            coot::match_torsions_gemmi(res_moving, res_ref, rest_info.second);

        EXPECT_EQ(data.res_moving, res_moving);
        EXPECT_EQ(data.res_ref, res_ref);
    }

    // Case 2 (oracle): invalid reference — different residue type
    {
        coot::protein_geometry geom;
        // NOTE: intentionally NOT calling geom.init_standard() to match oracle behaviour

        gemmi::Residue *res_moving = find_residue(st, "A", 5);
        gemmi::Residue *res_ref    = find_residue(st, "A", 10);

        EXPECT_NE(res_moving, nullptr);
        EXPECT_NE(res_ref, nullptr);

        std::string res_name(res_moving->name);
        std::string ref_name(res_ref->name);
        EXPECT_EQ(res_name, "GLN");
        EXPECT_EQ(ref_name, "ILE");

        std::pair<bool, coot::dictionary_residue_restraints_t> rest_info =
            geom.get_monomer_restraints(res_name, 0);

        EXPECT_FALSE(rest_info.first); // oracle did not reach the inner block

        coot::match_torsions_data_t data =
            coot::match_torsions_gemmi(res_moving, res_ref, rest_info.second);

        EXPECT_EQ(data.res_moving, res_moving);
        EXPECT_EQ(data.res_ref, res_ref);
    }

    // Case 3 (complementary): with init_standard() so restraints ARE found — tests non-trivial path
    {
        coot::protein_geometry geom;
        geom.init_standard();

        gemmi::Residue *res_moving = find_residue(st, "A", 5);
        gemmi::Residue *res_ref    = find_residue(st, "A", 47);

        EXPECT_NE(res_moving, nullptr);
        EXPECT_NE(res_ref, nullptr);

        std::string res_name(res_moving->name);
        EXPECT_EQ(res_name, "GLN");

        std::pair<bool, coot::dictionary_residue_restraints_t> rest_info =
            geom.get_monomer_restraints(res_name, 0);

        EXPECT_TRUE(rest_info.first); // with init_standard(), restraints should be found

        // Get torsion restraints
        std::vector<coot::dict_torsion_restraint_t> tr_ligand =
            geom.get_monomer_torsions_from_geometry(res_name, 0, 0);
        std::vector<coot::dict_torsion_restraint_t> tr_ref =
            geom.get_monomer_torsions_from_geometry(res_name, 0, 0);

        EXPECT_FALSE(tr_ligand.empty());
        EXPECT_EQ(tr_ligand.size(), tr_ref.size());

        // Construct match_torsions_data and verify storage
        coot::match_torsions_data_t data =
            coot::match_torsions_gemmi(res_moving, res_ref, rest_info.second);

        EXPECT_EQ(data.res_moving, res_moving);
        EXPECT_EQ(data.res_ref, res_ref);
        EXPECT_EQ(data.moving_residue_restraints.comp_id(), rest_info.second.comp_id());
    }

    // Case 4 (complementary): different residue types with init_standard()
    {
        coot::protein_geometry geom;
        geom.init_standard();

        gemmi::Residue *res_moving = find_residue(st, "A", 5);
        gemmi::Residue *res_ref    = find_residue(st, "A", 10);

        EXPECT_NE(res_moving, nullptr);
        EXPECT_NE(res_ref, nullptr);

        std::string res_name(res_moving->name);
        std::string ref_name(res_ref->name);
        EXPECT_EQ(res_name, "GLN");
        EXPECT_EQ(ref_name, "ILE");

        std::pair<bool, coot::dictionary_residue_restraints_t> rest_info =
            geom.get_monomer_restraints(res_name, 0);

        EXPECT_TRUE(rest_info.first);

        std::vector<coot::dict_torsion_restraint_t> tr_ligand =
            geom.get_monomer_torsions_from_geometry(res_name, 0, 0);
        std::vector<coot::dict_torsion_restraint_t> tr_ref =
            geom.get_monomer_torsions_from_geometry(ref_name, 0, 0);

        EXPECT_FALSE(tr_ligand.empty());
        // GLN and ILE have different torsions, so sizes may differ

        coot::match_torsions_data_t data =
            coot::match_torsions_gemmi(res_moving, res_ref, rest_info.second);

        EXPECT_EQ(data.res_moving, res_moving);
        EXPECT_EQ(data.res_ref, res_ref);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}