#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

namespace {

const gemmi::Atom* find_atom(const gemmi::Residue& res, const std::string& name) {
    for (const auto& a : res.atoms) {
        std::string trimmed_name = a.name;
        // Trim whitespace from both ends
        size_t start = trimmed_name.find_first_not_of(' ');
        if (start == std::string::npos) continue;
        size_t end = trimmed_name.find_last_not_of(' ');
        trimmed_name = trimmed_name.substr(start, end - start + 1);
        std::string target = name;
        size_t ts = target.find_first_not_of(' ');
        if (ts == std::string::npos) continue;
        size_t te = target.find_last_not_of(' ');
        target = target.substr(ts, te - ts + 1);
        if (trimmed_name == target) {
            return &a;
        }
    }
    return nullptr;
}

} // anonymous namespace

TEST(OracleTest, add_amino_hydrogens) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-hydrogen.pdb");

    gemmi::Model *model = &st.models[0];
    ASSERT_NE(model, nullptr);

    // Find chain "A" by iterating
    gemmi::Chain *chain = nullptr;
    for (auto &c : model->chains) {
        if (c.name == "A") {
            chain = &c;
            break;
        }
    }
    ASSERT_NE(chain, nullptr);

    // --- Case 1: ASN residue 3 ---
    {
        gemmi::Residue *asn_res = nullptr;
        for (auto &r : chain->residues) {
            if (r.seqid.num.value == 3 && r.name == "ASN") {
                asn_res = &r;
                break;
            }
        }

        if (asn_res) {
            EXPECT_STREQ(asn_res->name.c_str(), "ASN");
            EXPECT_STREQ(chain->name.c_str(), "A");
            EXPECT_EQ(asn_res->seqid.num.value, 3);

            const gemmi::Atom *nd2 = find_atom(*asn_res, "ND2");
            const gemmi::Atom *cg  = find_atom(*asn_res, "CG");
            const gemmi::Atom *od1 = find_atom(*asn_res, "OD1");
            EXPECT_NE(nd2, nullptr);
            EXPECT_NE(cg, nullptr);
            EXPECT_NE(od1, nullptr);

            int atoms_before = static_cast<int>(asn_res->atoms.size());
            EXPECT_EQ(atoms_before, 14);

            coot::reduce::add_amino_hydrogens_gemmi("HD21", "HD22", "ND2", "CG", "OD1", 1.01, *asn_res, st, chain);

            int atoms_after = static_cast<int>(asn_res->atoms.size());
            EXPECT_EQ(atoms_after, 14);

            const gemmi::Atom *hd21 = find_atom(*asn_res, "HD21");
            const gemmi::Atom *hd22 = find_atom(*asn_res, "HD22");
            EXPECT_NE(hd21, nullptr);
            EXPECT_NE(hd22, nullptr);
        }
    }

    // --- Case 2: GLN residue 29 (not present in PDB — oracle skips this case) ---
    // The PDB has VAL at residue 29, not GLN, so this case is omitted
    // to match the oracle's observed output.

    // --- Case 3: Edge case — residue without required atoms (should skip) ---
    {
        gemmi::Residue *gly_res = nullptr;
        for (auto &r : chain->residues) {
            if (r.name == "GLY") {
                gly_res = &r;
                break;
            }
        }

        if (gly_res) {
            EXPECT_STREQ(gly_res->name.c_str(), "GLY");
            EXPECT_EQ(gly_res->seqid.num.value, 11);

            const gemmi::Atom *fake1 = find_atom(*gly_res, "ND2");
            EXPECT_EQ(fake1, nullptr);

            int atoms_before = static_cast<int>(gly_res->atoms.size());
            EXPECT_EQ(atoms_before, 7);

            coot::reduce::add_amino_hydrogens_gemmi("HD21", "HD22", "ND2", "CG", "OD1", 1.01, *gly_res, st, chain);

            int atoms_after = static_cast<int>(gly_res->atoms.size());
            EXPECT_EQ(atoms_after, 7);
            EXPECT_EQ(atoms_after - atoms_before, 0);
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}