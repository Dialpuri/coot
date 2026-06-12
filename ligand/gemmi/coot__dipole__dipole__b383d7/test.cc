#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "geometry/protein-geometry.hh"
#include "function.hh"

TEST(OracleTest, dipole) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_GT(st.models.size(), 0u);

    coot::protein_geometry geom;
    geom.init_standard();

    // Find residue by chain and seqnum
    auto find_residue = [&](const std::string& chain_id, int seq_num) -> gemmi::Residue* {
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                if (chain.name != chain_id) continue;
                for (auto& res : chain.residues) {
                    if (res.seqid.num.value == seq_num) {
                        return &res;
                    }
                }
            }
        }
        return nullptr;
    };

    // Case 1: valid residue — MET at //A/1, restraints found but no partial charges → exception
    {
        gemmi::Residue* res = find_residue("A", 1);
        EXPECT_NE(res, nullptr);

        EXPECT_EQ(res->name, "MET");

        std::pair<bool, coot::dictionary_residue_restraints_t> rest_pair =
            geom.get_monomer_restraints("MET", 0);

        EXPECT_TRUE(rest_pair.first);

        coot::dictionary_residue_restraints_t &dict_res = rest_pair.second;

        // Build CRA for this residue
        gemmi::CRA cra{nullptr, res, nullptr};
        EXPECT_THROW({
            coot::dipole::dipole_gemmi(dict_res, cra);
        }, std::exception);
    }

    // Case 2: empty restraints — should throw
    {
        coot::dictionary_residue_restraints_t empty_dict;
        gemmi::Residue* res2 = find_residue("A", 50);
        EXPECT_NE(res2, nullptr);

        gemmi::CRA cra{nullptr, res2, nullptr};
        EXPECT_THROW({
            coot::dipole::dipole_gemmi(empty_dict, cra);
        }, std::exception);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}