#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "geometry/protein-geometry.hh"

TEST(OracleTest, charged_atoms) {
    // Minimal setup matching oracle exactly
    coot::protein_geometry geom;
    geom.init_standard();

    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    gemmi::Model& model = st.models[0];

    // Helper to find residue by CID (chain/res_seq_id)
    auto find_residue_by_cid = [&model](const std::string& cid) -> gemmi::Residue* {
        // CID format: //chain/res_seq_id, e.g. "//A/1" or "//A/50"
        size_t last_slash = cid.find_last_of('/');
        std::string chain_id = cid.substr(2, last_slash - 2);
        int res_seq = std::stoi(cid.substr(last_slash + 1));

        for (auto& chain : model.chains) {
            if (chain.name == chain_id) {
                for (auto& res : chain.residues) {
                    if (res.seqid.num.value == res_seq) {
                        return &res;
                    }
                }
            }
        }
        return nullptr;
    };

    // Case 1: MET at //A/1
    {
        std::string cid = "//A/1";
        gemmi::Residue* res = find_residue_by_cid(cid);
        ASSERT_NE(res, nullptr);

        std::string res_name = res->name;
        std::pair<bool, coot::dictionary_residue_restraints_t> rest_pair =
            geom.get_monomer_restraints(res_name, 0);

        std::vector<std::pair<coot::dictionary_residue_restraints_t, gemmi::Residue*>> dict_res_pairs;
        dict_res_pairs.push_back(std::make_pair(rest_pair.second, res));

        std::vector<std::pair<gemmi::Atom*, float>> charged_ats =
            coot::dipole::charged_atoms_gemmi(dict_res_pairs);
        EXPECT_EQ(charged_ats.size(), 0u);
    }

    // Case 2: ARG at //A/50
    {
        std::string cid = "//A/50";
        gemmi::Residue* res = find_residue_by_cid(cid);
        ASSERT_NE(res, nullptr);

        std::string res_name = res->name;
        std::pair<bool, coot::dictionary_residue_restraints_t> rest_pair =
            geom.get_monomer_restraints(res_name, 0);

        std::vector<std::pair<coot::dictionary_residue_restraints_t, gemmi::Residue*>> dict_res_pairs;
        dict_res_pairs.push_back(std::make_pair(rest_pair.second, res));

        std::vector<std::pair<gemmi::Atom*, float>> charged_ats =
            coot::dipole::charged_atoms_gemmi(dict_res_pairs);
        EXPECT_EQ(charged_ats.size(), 0u);
    }

    // Case 3: empty input vector
    {
        std::vector<std::pair<coot::dictionary_residue_restraints_t, gemmi::Residue*>> empty_pairs;
        std::vector<std::pair<gemmi::Atom*, float>> charged_ats =
            coot::dipole::charged_atoms_gemmi(empty_pairs);
        EXPECT_EQ(charged_ats.size(), 0u);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}