#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

// Find the CA atom in a residue
static gemmi::CRA find_ca(gemmi::Chain& chain, gemmi::Residue& res) {
    for (gemmi::Atom& atom : res.atoms) {
        if (atom.name == "CA") {
            return gemmi::CRA{&chain, &res, &atom};
        }
    }
    return gemmi::CRA{&chain, &res, nullptr};
}

TEST(OracleTest, delete_upstream) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    EXPECT_GE(st.models.size(), 1);
    gemmi::Model& model = st.models[0];

    EXPECT_GE(model.chains.size(), 1);
    gemmi::Chain& chain = model.chains[0];

    std::string chain_id = chain.name;

    // Get residues by index (gemmi is 0-indexed, like C++)
    EXPECT_GE(chain.residues.size(), 10);

    gemmi::Residue& res_0 = chain.residues[0];
    gemmi::Residue& res_2 = chain.residues[2];
    gemmi::Residue& res_4 = chain.residues[4];
    gemmi::Residue& res_5 = chain.residues[5];
    gemmi::Residue& res_9 = chain.residues[9];

    // Build sel_atoms1: atoms in residues 0..4
    int res0_seq = res_0.seqid.num.value;
    int res4_seq = res_4.seqid.num.value;

    std::vector<gemmi::CRA> sel1_atoms;
    for (gemmi::Residue& res : chain.residues) {
        int seq = res.seqid.num.value;
        if (seq >= res0_seq && seq <= res4_seq) {
            for (gemmi::Atom& atom : res.atoms) {
                sel1_atoms.push_back(gemmi::CRA{&chain, &res, &atom});
            }
        }
    }

    // Build sel_atoms2: atoms in residues 5..9
    int res5_seq = res_5.seqid.num.value;
    int res9_seq = res_9.seqid.num.value;

    std::vector<gemmi::CRA> sel2_atoms;
    for (gemmi::Residue& res : chain.residues) {
        int seq = res.seqid.num.value;
        if (seq >= res5_seq && seq <= res9_seq) {
            for (gemmi::Atom& atom : res.atoms) {
                sel2_atoms.push_back(gemmi::CRA{&chain, &res, &atom});
            }
        }
    }

    // Find CA atoms in res_2 and res_5 for atom_pairs
    gemmi::CRA ca_res2 = find_ca(chain, res_2);
    gemmi::CRA ca_res5 = find_ca(chain, res_5);

    EXPECT_NE(ca_res2.atom, nullptr);
    EXPECT_NE(ca_res5.atom, nullptr);

    // === CASE 1: delete_upstream with from_first=true ===
    int residues_before = static_cast<int>(chain.residues.size());

    coot::match_container_for_residues_t container;
    container.add(ca_res2, ca_res5);

    // INPUT assertions
    EXPECT_EQ(res_2.seqid.num.value, 3);
    EXPECT_EQ(static_cast<int>(sel1_atoms.size()), 45);
    EXPECT_EQ(static_cast<int>(sel2_atoms.size()), 42);
    EXPECT_EQ(residues_before, 267);

    container.delete_upstream_gemmi(st, true, sel1_atoms, sel2_atoms);

    int residues_after_case1 = static_cast<int>(chain.residues.size());
    EXPECT_EQ(residues_after_case1, 265);
    EXPECT_EQ(residues_before - residues_after_case1, 2);

    // Check what the first residue is now
    EXPECT_GE(chain.residues.size(), 1);
    gemmi::Residue& new_first = chain.residues[0];
    EXPECT_EQ(new_first.seqid.num.value, 3);
    EXPECT_EQ(new_first.name, "ASN");

    // === CASE 2: empty selection (guard: n_atoms == 0 returns early) ===
    {
        std::vector<gemmi::CRA> sel_empty;
        std::vector<gemmi::CRA> sel_empty2;

        coot::match_container_for_residues_t container2;
        gemmi::CRA null_cra = gemmi::CRA{&chain, &new_first, nullptr};
        container2.add(null_cra, null_cra);

        // INPUT assertions
        EXPECT_EQ(static_cast<int>(sel_empty.size()), 0);

        int residues_before_empty = static_cast<int>(chain.residues.size());
        container2.delete_upstream_gemmi(st, true, sel_empty, sel_empty2);
        int residues_after_empty = static_cast<int>(chain.residues.size());

        // OUTPUT assertions
        EXPECT_TRUE(residues_before_empty == residues_after_empty);
        EXPECT_EQ(residues_after_empty, 265);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}