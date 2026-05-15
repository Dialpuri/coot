#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include <map>
#include <string>

// Pad a gemmi atom name to MMDB's 4-character width
inline std::string pad_atom_name(const std::string &name) {
    if (name.size() == 1) return " " + name + "  ";
    if (name.size() == 2) return " " + name + " ";
    return name;  // already 4 chars
}

TEST(OracleTest, get_c_beta_deviations) {
    {
        gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

        std::map<const gemmi::Residue*, std::map<std::string, coot::c_beta_deviation_t>> dev_map = coot::get_c_beta_deviations_gemmi(st);

        EXPECT_EQ(dev_map.size(), 253u);

        // Build lookup maps from residue pointer to chain_id, seqnum, resname
        std::map<const gemmi::Residue*, std::string> res_chain;
        std::map<const gemmi::Residue*, int> res_seqnum;
        std::map<const gemmi::Residue*, std::string> res_resname;
        for (auto& model : st.models) {
            for (auto& chain : model.chains) {
                for (auto& res : chain.residues) {
                    res_chain[&res] = chain.name;
                    res_seqnum[&res] = res.seqid.num.value;
                    res_resname[&res] = res.name;
                }
            }
        }

        for (auto const &entry : dev_map) {
            const gemmi::Residue *res = entry.first;
            std::string chain_id = res_chain[res];
            int seqnum = res_seqnum[res];
            std::string resname = res_resname[res];

            for (auto const &kv : entry.second) {
                coot::c_beta_deviation_t dev = kv.second;
                std::string atom_name = dev.at ? pad_atom_name(dev.at->name) : "nullptr";

                if (seqnum == 6 && resname == "LYS") {
                    EXPECT_NEAR(kv.second.dist, 0.0145626, 1e-7);
                    EXPECT_EQ(atom_name, " CB ");
                } else if (seqnum == 5 && resname == "GLN") {
                    EXPECT_NEAR(kv.second.dist, 0.0755736, 1e-7);
                    EXPECT_EQ(atom_name, " CB ");
                } else if (seqnum == 4 && resname == "PHE") {
                    EXPECT_NEAR(kv.second.dist, 0.0855131, 1e-7);
                    EXPECT_EQ(atom_name, " CB ");
                } else if (seqnum == 10 && resname == "ILE") {
                    EXPECT_NEAR(kv.second.dist, 0.20076, 1e-5);
                    EXPECT_EQ(atom_name, " CB ");
                } else if (seqnum == 132 && resname == "ASN") {
                    EXPECT_NEAR(kv.second.dist, 0.300884, 1e-5);
                    EXPECT_EQ(atom_name, " CB ");
                } else if (seqnum == 295 && resname == "HIS") {
                    EXPECT_NEAR(kv.second.dist, 0.341183, 1e-5);
                    EXPECT_EQ(atom_name, " CB ");
                } else if (seqnum == 195 && resname == "GLU") {
                    EXPECT_NEAR(kv.second.dist, 0.347569, 1e-5);
                    EXPECT_EQ(atom_name, " CB ");
                }
            }
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}