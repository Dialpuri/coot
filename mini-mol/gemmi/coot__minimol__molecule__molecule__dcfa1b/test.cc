#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <gemmi/calculate.hpp>
#include "function.hh"

#include <vector>
#include <cstdlib>
#include <iostream>

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    int r = RUN_ALL_TESTS();
    std::fflush(nullptr);
    _exit(r);
}

TEST(OracleTest, diagnostic) {
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

    std::cerr << "Models: " << st.models.size() << std::endl;
    std::cerr << "count_atom_sites: " << gemmi::count_atom_sites(st) << std::endl;
    for (size_t mi = 0; mi < st.models.size(); mi++) {
        auto& mod = st.models[mi];
        std::cerr << "  Model " << mi << ": " << mod.chains.size() << " chains, " << gemmi::count_atom_sites(mod) << " sites" << std::endl;
        for (auto& chain : mod.chains) {
            int res_count = 0;
            for (auto& res : chain.residues) res_count++;
            int lig_count = 0;
            for (auto& res : chain.get_ligands()) lig_count++;
            int wat_count = 0;
            for (auto& res : chain.get_waters()) wat_count++;
            std::cerr << "    Chain " << chain.name << ": " << res_count << " res, " << lig_count << " lig, " << wat_count << " wat" << std::endl;
        }
    }
}
