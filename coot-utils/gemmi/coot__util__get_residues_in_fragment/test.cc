#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "geometry/residue-and-atom-specs.hh"

TEST(OracleTest, get_residues_in_fragment) {
    // case: first residue
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_GE(st.models.size(), 1) << "no models";

    gemmi::Model& model = st.models[0];
    ASSERT_GE(model.chains.size(), 1) << "no chains";

    gemmi::Chain& chain = model.chains[0];
    
    ASSERT_GE(chain.residues.size(), 1) << "no residues in chain";

    gemmi::Residue& first_res = chain.residues[0];
    coot::residue_spec_t spec(std::string(chain.name), 
                              first_res.seqid.num.value, 
                              std::string(1, first_res.seqid.icode));

    std::vector<gemmi::Residue*> result = coot::util::get_residues_in_fragment_gemmi(&chain, spec);

    EXPECT_EQ(result.size(), chain.residues.size());
    EXPECT_EQ(result.size(), 267u);

    int count = std::min(3, static_cast<int>(result.size()));
    for (int i = 0; i < count; i++) {
        gemmi::Residue* r = result[i];
        if (i == 0) {
            EXPECT_EQ(chain.name, "A");
            EXPECT_EQ(r->seqid.num.value, 1);
            EXPECT_EQ(std::string(1, r->seqid.icode), " ");
            EXPECT_EQ(r->name, "MET");
        } else if (i == 1) {
            EXPECT_EQ(chain.name, "A");
            EXPECT_EQ(r->seqid.num.value, 2);
            EXPECT_EQ(std::string(1, r->seqid.icode), " ");
            EXPECT_EQ(r->name, "GLU");
        } else if (i == 2) {
            EXPECT_EQ(chain.name, "A");
            EXPECT_EQ(r->seqid.num.value, 3);
            EXPECT_EQ(std::string(1, r->seqid.icode), " ");
            EXPECT_EQ(r->name, "ASN");
        }
    }

    // case: middle residue
    gemmi::Residue* mid_res = &chain.residues[chain.residues.size()/2];
    coot::residue_spec_t mid_spec(std::string(chain.name),
                                  mid_res->seqid.num.value,
                                  std::string(1, mid_res->seqid.icode));
    
    std::vector<gemmi::Residue*> result2 = coot::util::get_residues_in_fragment_gemmi(&chain, mid_spec);

    EXPECT_EQ(result2.size(), 267u);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}