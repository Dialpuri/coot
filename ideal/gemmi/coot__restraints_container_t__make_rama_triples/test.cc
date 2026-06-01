#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

TEST(OracleTest, make_rama_triples_gemmi) {
    // Load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_FALSE(st.models.empty());

    gemmi::Model& model = st.models[0];

    // Build CRA vector for chain A only (mirrors Select("//A/"))
    std::vector<gemmi::CRA> chain_a;
    for (gemmi::Chain& chain : model.chains) {
        if (chain.name == "A") {
            for (gemmi::Residue& res : chain.residues) {
                chain_a.push_back(gemmi::CRA{&chain, &res, nullptr});
            }
        }
    }

    EXPECT_EQ(chain_a.size(), 267u);

    coot::protein_geometry pg;
    pg.init_standard();

    // Call make_rama_triples_gemmi
    std::vector<coot::rama_triple_t_gemmi> triples =
        coot::make_rama_triples_gemmi(chain_a, pg);

    EXPECT_EQ(triples.size(), 265u);

    // Verify first few triples
    EXPECT_NE(triples[0].r_1, nullptr);
    EXPECT_NE(triples[0].r_2, nullptr);
    EXPECT_NE(triples[0].r_3, nullptr);

    // triple_0: A/1/MET A/2/GLU A/3/ASN link_type=TRANS
    EXPECT_EQ(std::string(triples[0].chain_1->name), "A");
    EXPECT_EQ(triples[0].r_1->seqid.num.value, 1);
    EXPECT_EQ(std::string(triples[0].r_1->name), "MET");
    EXPECT_EQ(triples[0].r_2->seqid.num.value, 2);
    EXPECT_EQ(std::string(triples[0].r_2->name), "GLU");
    EXPECT_EQ(triples[0].r_3->seqid.num.value, 3);
    EXPECT_EQ(std::string(triples[0].r_3->name), "ASN");
    EXPECT_EQ(triples[0].link_type, "TRANS");

    // triple_1: A/2/GLU A/3/ASN A/4/PHE link_type=TRANS
    EXPECT_EQ(triples[1].r_1->seqid.num.value, 2);
    EXPECT_EQ(std::string(triples[1].r_1->name), "GLU");
    EXPECT_EQ(triples[1].r_2->seqid.num.value, 3);
    EXPECT_EQ(std::string(triples[1].r_2->name), "ASN");
    EXPECT_EQ(triples[1].r_3->seqid.num.value, 4);
    EXPECT_EQ(std::string(triples[1].r_3->name), "PHE");
    EXPECT_EQ(triples[1].link_type, "TRANS");

    // triple_2: A/3/ASN A/4/PHE A/5/GLN link_type=TRANS
    EXPECT_EQ(triples[2].r_1->seqid.num.value, 3);
    EXPECT_EQ(std::string(triples[2].r_1->name), "ASN");
    EXPECT_EQ(triples[2].r_2->seqid.num.value, 4);
    EXPECT_EQ(std::string(triples[2].r_2->name), "PHE");
    EXPECT_EQ(triples[2].r_3->seqid.num.value, 5);
    EXPECT_EQ(std::string(triples[2].r_3->name), "GLN");
    EXPECT_EQ(triples[2].link_type, "TRANS");

    // Verify some gap-jumping triples (e.g., triple_34 jumps from 36 to 46)
    // triple_34: A/35/ILE A/36/ARG A/46/SER link_type=TRANS
    EXPECT_EQ(triples[34].r_1->seqid.num.value, 35);
    EXPECT_EQ(std::string(triples[34].r_1->name), "ILE");
    EXPECT_EQ(triples[34].r_2->seqid.num.value, 36);
    EXPECT_EQ(std::string(triples[34].r_2->name), "ARG");
    EXPECT_EQ(triples[34].r_3->seqid.num.value, 46);
    EXPECT_EQ(std::string(triples[34].r_3->name), "SER");
    EXPECT_EQ(triples[34].link_type, "TRANS");

    // triple_61: A/71/HIS A/72/THR A/76/LEU link_type=TRANS
    EXPECT_EQ(triples[61].r_1->seqid.num.value, 71);
    EXPECT_EQ(std::string(triples[61].r_1->name), "HIS");
    EXPECT_EQ(triples[61].r_2->seqid.num.value, 72);
    EXPECT_EQ(std::string(triples[61].r_2->name), "THR");
    EXPECT_EQ(triples[61].r_3->seqid.num.value, 76);
    EXPECT_EQ(std::string(triples[61].r_3->name), "LEU");
    EXPECT_EQ(triples[61].link_type, "TRANS");

    // triple_132: A/145/ASP A/146/PHE A/166/LEU link_type=TRANS
    EXPECT_EQ(triples[132].r_1->seqid.num.value, 145);
    EXPECT_EQ(std::string(triples[132].r_1->name), "ASP");
    EXPECT_EQ(triples[132].r_2->seqid.num.value, 146);
    EXPECT_EQ(std::string(triples[132].r_2->name), "PHE");
    EXPECT_EQ(triples[132].r_3->seqid.num.value, 166);
    EXPECT_EQ(std::string(triples[132].r_3->name), "LEU");
    EXPECT_EQ(triples[132].link_type, "TRANS");

    // Last triple: triple_264: A/296/LEU A/297/ARG A/298/LEU link_type=TRANS
    size_t last_idx = triples.size() - 1;
    EXPECT_EQ(triples[last_idx].r_1->seqid.num.value, 296);
    EXPECT_EQ(std::string(triples[last_idx].r_1->name), "LEU");
    EXPECT_EQ(triples[last_idx].r_2->seqid.num.value, 297);
    EXPECT_EQ(std::string(triples[last_idx].r_2->name), "ARG");
    EXPECT_EQ(triples[last_idx].r_3->seqid.num.value, 298);
    EXPECT_EQ(std::string(triples[last_idx].r_3->name), "LEU");
    EXPECT_EQ(triples[last_idx].link_type, "TRANS");

    // All link types should be TRANS
    for (size_t i = 0; i < triples.size(); i++) {
        EXPECT_EQ(triples[i].link_type, "TRANS");
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}