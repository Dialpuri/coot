#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, nomenclature_list) {
    // case: main test
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_GE(st.models.size(), 1u);

    gemmi::Model &model = st.models[0];

    // Call list_gemmi (equivalent to nomen.list(&geom) in MMDB)
    // Note: geom is unused in the original because apply_swap=false in fix_and_swap_maybe
    std::vector<gemmi::CRA> residues = coot::nomenclature::list_gemmi(model, nullptr);

    EXPECT_EQ(residues.size(), 8u);

    // Check each residue
    struct ExpectedResidue {
        std::string chain;
        int seqnum;
        std::string resname;
    };

    ExpectedResidue expected[] = {
        {"A", 4, "PHE"},
        {"A", 8, "GLU"},
        {"A", 19, "TYR"},
        {"A", 28, "GLU"},
        {"A", 82, "PHE"},
        {"A", 208, "GLU"},
        {"A", 236, "TYR"},
        {"A", 247, "ASP"},
    };

    for (size_t i = 0; i < residues.size(); i++) {
        EXPECT_NE(residues[i].residue, nullptr);
        if (residues[i].residue) {
            gemmi::CRA cra = residues[i];
            std::string chain_id = cra.chain ? cra.chain->name : "";
            int seqnum = cra.residue ? cra.residue->seqid.num.value : 0;
            std::string resname = cra.residue ? cra.residue->name : "";

            EXPECT_EQ(chain_id, expected[i].chain);
            EXPECT_EQ(seqnum, expected[i].seqnum);
            EXPECT_EQ(resname, expected[i].resname);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}