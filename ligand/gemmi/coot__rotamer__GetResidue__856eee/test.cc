#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"
#include "geometry/protein-geometry.hh"

static FILE* info_old = nullptr;

TEST(OracleTest, richardson_rotamer) {
    // Suppress INFO:: messages to stderr for clarity
    info_old = freopen("/dev/null", "w", stderr);

    coot::protein_geometry pg;
    pg.init_standard();

    // Read the PDB with gemmi
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    EXPECT_EQ(st.models.size(), 1u);

    gemmi::Model& model = st.models[0];

    // Restore stderr so we see debug output
    fflush(stderr);
    freopen("/proc/self/fd/2", "w", stderr);

    // Helper to find a residue by chain and sequence number
    auto find_residue = [&model](const std::string& chain_id, int seq_num) -> gemmi::CRA {
        for (gemmi::Chain& chain : model.chains) {
            if (chain.name != chain_id) continue;
            for (gemmi::Residue& res : chain.residues) {
                if (res.seqid.num.value == seq_num) {
                    return gemmi::CRA{&chain, &res, nullptr};
                }
            }
        }
        return gemmi::CRA{nullptr, nullptr, nullptr};
    };

    // case 1: valid residue with rotamers (MET at //A/1, index 0)
    {
        gemmi::CRA cra = find_residue("A", 1);
        ASSERT_TRUE(cra.residue != nullptr);

        std::string res_type = cra.residue->name;
        std::cerr << "DEBUG case1: res_type=[" << res_type << "]" << std::endl;
        EXPECT_EQ(res_type, "MET");

        std::pair<short int, coot::dictionary_residue_restraints_t> p =
            pg.get_monomer_restraints(res_type, 0);

        std::cerr << "DEBUG case1: p.first=" << p.first << std::endl;
        EXPECT_TRUE(p.first);

        auto moving_res = coot::GetResidue_gemmi(cra, p.second, 0);
        std::cerr << "DEBUG case1: moving_res=" << (moving_res ? "non-null" : "null") << std::endl;
        EXPECT_NE(moving_res, nullptr);
        if (moving_res) {
            int n_atoms = static_cast<int>(moving_res->atoms.size());
            std::cerr << "DEBUG case1: n_atoms=" << n_atoms << std::endl;
            EXPECT_EQ(n_atoms, 8);
            if (n_atoms > 0) {
                std::string atom_name = moving_res->atoms[0].name;
                std::cerr << "DEBUG case1: first_atom_name=[" << atom_name << "] len=" << atom_name.size() << std::endl;
            }
        }
    }

    // case 2: rotamer index out of range
    {
        gemmi::CRA cra = find_residue("A", 1);
        ASSERT_TRUE(cra.residue != nullptr);

        std::string res_type = cra.residue->name;
        EXPECT_EQ(res_type, "MET");

        std::pair<short int, coot::dictionary_residue_restraints_t> p =
            pg.get_monomer_restraints(res_type, 0);

        if (p.first) {
            auto moving_res = coot::GetResidue_gemmi(cra, p.second, 9999);
            std::cerr << "DEBUG case2: moving_res=" << (moving_res ? "non-null" : "null") << std::endl;
            EXPECT_NE(moving_res, nullptr);
            if (moving_res) {
                int n_atoms = static_cast<int>(moving_res->atoms.size());
                std::cerr << "DEBUG case2: n_atoms=" << n_atoms << std::endl;
                EXPECT_EQ(n_atoms, 8);
            }
        }
    }

    // case 3: GLY at //A/11
    {
        gemmi::CRA cra = find_residue("A", 11);
        ASSERT_TRUE(cra.residue != nullptr);

        std::string res_type = cra.residue->name;
        std::cerr << "DEBUG case3: res_type=[" << res_type << "]" << std::endl;
        EXPECT_EQ(res_type, "GLY");

        std::pair<short int, coot::dictionary_residue_restraints_t> p =
            pg.get_monomer_restraints(res_type, 0);

        std::cerr << "DEBUG case3: p.first=" << p.first << std::endl;
        if (p.first) {
            auto moving_res = coot::GetResidue_gemmi(cra, p.second, 0);
            std::cerr << "DEBUG case3: moving_res=" << (moving_res ? "non-null" : "null") << std::endl;
            EXPECT_NE(moving_res, nullptr);
            if (moving_res) {
                int n_atoms = static_cast<int>(moving_res->atoms.size());
                std::cerr << "DEBUG case3: n_atoms=" << n_atoms << std::endl;
                EXPECT_EQ(n_atoms, 4);
            }
        }
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
