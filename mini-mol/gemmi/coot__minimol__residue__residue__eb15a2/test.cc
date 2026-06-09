#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/input.hpp>
#include "function.hh"

// --- Case 1: residue() default constructor ---
TEST(ResidueGemmiTest, default_constructor) {
    coot::minimol::residue r;
    EXPECT_EQ(r.seqnum, -2147483647);
    EXPECT_EQ(r.ins_code, "");
    EXPECT_EQ(r.name, "");
    EXPECT_EQ(r.n_atoms(), 0u);
}

// --- Case 2: residue(int i) ---
TEST(ResidueGemmiTest, int_constructor) {
    int i = 42;
    coot::minimol::residue r(i);
    EXPECT_EQ(r.seqnum, 42);
    EXPECT_EQ(r.ins_code, "");
    EXPECT_EQ(r.name, "");
    EXPECT_EQ(r.n_atoms(), 0u);
}

// --- Case 3: residue(int i, const std::string & resname) ---
TEST(ResidueGemmiTest, int_string_constructor) {
    int i = 100;
    std::string resname = "ALA";
    coot::minimol::residue r(i, resname);
    EXPECT_EQ(r.seqnum, 100);
    EXPECT_EQ(r.ins_code, "");
    EXPECT_EQ(r.name, "ALA");
    EXPECT_EQ(r.n_atoms(), 0u);
}

// --- Case 4: residue(mmdb::Residue * residue_p) ---
TEST(ResidueGemmiTest, gemmi_residue_constructor) {
    gemmi::Structure st = gemmi::read_structure_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
    ASSERT_FALSE(st.models.empty());
    gemmi::Model &model = st.models[0];

    // Find chain "A"
    gemmi::Chain *chain_ptr = nullptr;
    for (auto &ch : model.chains) {
        if (ch.name == "A") {
            chain_ptr = &ch;
            break;
        }
    }
    ASSERT_NE(chain_ptr, nullptr);

    // Find residue with seqnum 1299 in chain A
    gemmi::Residue *res_ptr = nullptr;
    for (auto &res : chain_ptr->residues) {
        if (res.seqid.num.value == 1299) {
            res_ptr = &res;
            break;
        }
    }
    ASSERT_NE(res_ptr, nullptr);

    coot::minimol::residue r(res_ptr);
    EXPECT_EQ(r.seqnum, 1299);
    EXPECT_EQ(r.ins_code, "");
    EXPECT_EQ(r.name, "LZA");
    EXPECT_EQ(r.n_atoms(), 25u);
}

// --- Case 5: residue(mmdb::Residue *, keep_only_these_atoms) ---
TEST(ResidueGemmiTest, gemmi_residue_constructor_keep_atoms) {
    gemmi::Structure st = gemmi::read_structure_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example-protein-ligand.cif");
    ASSERT_FALSE(st.models.empty());
    gemmi::Model &model = st.models[0];

    // Find chain "A"
    gemmi::Chain *chain_ptr = nullptr;
    for (auto &ch : model.chains) {
        if (ch.name == "A") {
            chain_ptr = &ch;
            break;
        }
    }
    ASSERT_NE(chain_ptr, nullptr);

    // Find residue with seqnum 1299 in chain A
    gemmi::Residue *res_ptr = nullptr;
    for (auto &res : chain_ptr->residues) {
        if (res.seqid.num.value == 1299) {
            res_ptr = &res;
            break;
        }
    }
    ASSERT_NE(res_ptr, nullptr);

    // Build keep list with first 3 atom names (padded to PDB format)
    auto pad_atom_name = [](const std::string &name) {
        if (name.size() >= 4) return name;
        return std::string(1, ' ') + name + std::string(4 - 1 - name.size(), ' ');
    };
    std::vector<std::string> keep;
    for (const auto &g_atom : res_ptr->atoms) {
        keep.push_back(pad_atom_name(g_atom.name));
        if (keep.size() >= 3) break;
    }

    coot::minimol::residue r(res_ptr, keep);
    EXPECT_EQ(r.seqnum, 1299);
    EXPECT_EQ(r.name, "LZA");
    EXPECT_EQ(r.n_atoms(), 3u);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
