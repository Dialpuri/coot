#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static std::vector<std::string> get_atom_names(const gemmi::Residue& res) {
    std::vector<std::string> names;
    for (const auto& atom : res.atoms)
        names.push_back(atom.name);
    return names;
}

TEST(OracleTest, GraphMatchInfoMatchNames) {
    // case: original test
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::Model& model = st.models[0];
        gemmi::Chain& chain = model.chains[0];
        gemmi::Residue& res_ala = chain.residues[9]; // residue A/10 is index 9

        std::vector<std::string> original_names = get_atom_names(res_ala);
        EXPECT_EQ(original_names.size(), 8u);

        // Modify atom names to simulate a "working" residue
        std::vector<std::string> modified_names;
        for (auto& atom : res_ala.atoms) {
            std::string new_name = atom.name + "X";
            atom.name = new_name;
            modified_names.push_back(new_name);
        }

        // Create matching pairs: modified -> original
        std::vector<std::pair<std::pair<std::string, std::string>,
                              std::pair<std::string, std::string>>> matching_names;
        for (size_t i = 0; i < modified_names.size(); i++) {
            matching_names.push_back(std::make_pair(
                std::make_pair(modified_names[i], ""),
                std::make_pair(original_names[i], "")
            ));
        }

        coot::graph_match_info_t gmi;
        gmi.success = true;
        gmi.matching_atom_names = matching_names;

        gmi.match_names_gemmi(res_ala);

        // Verify - check if names match original
        std::vector<std::string> restored_names = get_atom_names(res_ala);
        bool all_correct = true;
        for (size_t i = 0; i < original_names.size(); i++) {
            if (restored_names[i] != original_names[i]) {
                all_correct = false;
                break;
            }
        }
        EXPECT_TRUE(all_correct) << "Not all atom names were restored correctly";
        EXPECT_EQ(restored_names.size(), 8u);
    }

    // case: residue with insertion code
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::Model& model = st.models[0];
        gemmi::Chain& chain = model.chains[0];
        gemmi::Residue& res_ala = chain.residues[47]; // residue A/48 is index 47

        std::vector<std::string> original_names = get_atom_names(res_ala);
        ASSERT_GT(original_names.size(), 0u);

        // Modify atom names
        std::vector<std::string> modified_names;
        for (auto& atom : res_ala.atoms) {
            std::string new_name = atom.name + "_modified";
            atom.name = new_name;
            modified_names.push_back(new_name);
        }

        // Create matching pairs
        std::vector<std::pair<std::pair<std::string, std::string>,
                              std::pair<std::string, std::string>>> matching_names;
        for (size_t i = 0; i < modified_names.size(); i++) {
            matching_names.push_back(std::make_pair(
                std::make_pair(modified_names[i], ""),
                std::make_pair(original_names[i], "")
            ));
        }

        coot::graph_match_info_t gmi;
        gmi.success = true;
        gmi.matching_atom_names = matching_names;

        gmi.match_names_gemmi(res_ala);

        // Verify - check if names match original
        std::vector<std::string> restored_names = get_atom_names(res_ala);
        bool all_correct = true;
        for (size_t i = 0; i < original_names.size(); i++) {
            if (restored_names[i] != original_names[i]) {
                all_correct = false;
                break;
            }
        }
        EXPECT_TRUE(all_correct) << "Atom names not restored correctly for residue with insertion code";
    }

    // case: graph match failed (success = false)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::Model& model = st.models[0];
        gemmi::Chain& chain = model.chains[0];
        gemmi::Residue& res_ala = chain.residues[9];

        std::vector<std::string> original_names = get_atom_names(res_ala);

        coot::graph_match_info_t gmi;
        gmi.success = false;

        std::vector<std::pair<std::pair<std::string, std::string>,
                              std::pair<std::string, std::string>>> matching_names;
        for (size_t i = 0; i < original_names.size(); i++) {
            matching_names.push_back(std::make_pair(
                std::make_pair("NONEXISTENT", ""),
                std::make_pair(original_names[i], "")
            ));
        }
        gmi.matching_atom_names = matching_names;

        gmi.match_names_gemmi(res_ala);

        std::vector<std::string> current_names = get_atom_names(res_ala);
        EXPECT_EQ(current_names, original_names);
    }

    // case: partial matching (some atoms matched, some not)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::Model& model = st.models[0];
        gemmi::Chain& chain = model.chains[0];
        gemmi::Residue& res_ala = chain.residues[9];

        std::vector<std::string> original_names = get_atom_names(res_ala);

        // Modify only some atom names
        std::vector<std::string> modified_names;
        for (size_t i = 0; i < res_ala.atoms.size(); i++) {
            if (i < 4) {
                std::string new_name = res_ala.atoms[i].name + "_tmp";
                res_ala.atoms[i].name = new_name;
                modified_names.push_back(new_name);
            } else {
                modified_names.push_back(res_ala.atoms[i].name);
            }
        }

        // Create matching pairs only for modified atoms
        std::vector<std::pair<std::pair<std::string, std::string>,
                              std::pair<std::string, std::string>>> matching_names;
        for (int i = 0; i < 4; i++) {
            matching_names.push_back(std::make_pair(
                std::make_pair(modified_names[i], ""),
                std::make_pair(original_names[i], "")
            ));
        }

        coot::graph_match_info_t gmi;
        gmi.success = true;
        gmi.matching_atom_names = matching_names;

        gmi.match_names_gemmi(res_ala);

        // Verify - check if modified names were restored
        std::vector<std::string> restored_names = get_atom_names(res_ala);
        bool all_correct = true;
        for (size_t i = 0; i < original_names.size(); i++) {
            if (restored_names[i] != original_names[i]) {
                all_correct = false;
                break;
            }
        }
        EXPECT_TRUE(all_correct) << "Partial matching not handled correctly";
    }

    // case: collision scenario (reference name already exists in residue)
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
        gemmi::Model& model = st.models[0];
        gemmi::Chain& chain = model.chains[0];
        gemmi::Residue& res_ala = chain.residues[9];

        std::vector<std::string> original_names = get_atom_names(res_ala);

        // Modify all atom names
        std::vector<std::string> modified_names;
        for (auto& atom : res_ala.atoms) {
            std::string new_name = atom.name + "_swap";
            atom.name = new_name;
            modified_names.push_back(new_name);
        }

        // Create matching pairs where we swap names (CA->CBX, CB->CAX)
        std::vector<std::pair<std::pair<std::string, std::string>,
                              std::pair<std::string, std::string>>> matching_names;
        matching_names.push_back(std::make_pair(
            std::make_pair(modified_names[1], ""),
            std::make_pair(original_names[2], "")
        ));
        matching_names.push_back(std::make_pair(
            std::make_pair(modified_names[2], ""),
            std::make_pair(original_names[1], "")
        ));

        coot::graph_match_info_t gmi;
        gmi.success = true;
        gmi.matching_atom_names = matching_names;

        gmi.match_names_gemmi(res_ala);

        // At least the matched atoms should have correct names
        std::vector<std::string> restored_names = get_atom_names(res_ala);
        bool has_ca = false;
        bool has_cb = false;
        for (const auto& name : restored_names) {
            if (name == original_names[1]) has_ca = true;
            if (name == original_names[2]) has_cb = true;
        }
        EXPECT_TRUE(has_ca) << "CA atom name not found after collision handling";
        EXPECT_TRUE(has_cb) << "CB atom name not found after collision handling";
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}