#include <gtest/gtest.h>
#include <fstream>
#include <string>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include "function.hh"

static void create_test_res_file(const std::string &filename) {
    std::ofstream f(filename.c_str());
    f << "TITL Test structure\n";
    f << "CELL 0.0 10.0 10.0 10.0 90.0 90.0 90.0\n";
    f << "SFAC C H O\n";
    f << "UNIT 6 12 6\n";
    f << "LATT 1\n";
    f << "SYMM X, Y, Z\n";
    f << "RESI 1\n";
    f << "C1  1  0.5000  0.5000  0.5000  0.50  2.00  C\n";
    f << "O1  1  0.6000  0.5000  0.5000  0.50  2.00  O\n";
    f << "N1  1  0.4000  0.5000  0.5000  0.50  2.00  N\n";
    f << "END\n";
    f.close();
}

TEST(OracleTest, ShelxReadFile) {
    {
        std::string res_file = "/tmp/test_shelx_valid.res";
        create_test_res_file(res_file);
        
        coot::ShelxIns s;
        coot::shelx_read_file_info_t_gemmi p = coot::read_file_gemmi(s, res_file);
        
        EXPECT_EQ(p.status, 1);
        EXPECT_NE(p.udd_afix_handle, -1);
        EXPECT_NE(p.structure, nullptr);
        
        if (p.structure) {
            const gemmi::Structure& st = *p.structure;
            EXPECT_EQ(st.models.size(), 1);
            
            const gemmi::Model& model = st.models[0];
            EXPECT_EQ(model.chains.size(), 1);
            
            const gemmi::Chain& chain = model.chains[0];
            EXPECT_EQ(chain.name, "A");
            
            EXPECT_EQ(chain.residues.size(), 1);
        }
    }
    
    {
        std::string nonexistent = "/tmp/nonexistent_file.res";
        
        coot::ShelxIns s2;
        coot::shelx_read_file_info_t_gemmi p2 = coot::read_file_gemmi(s2, nonexistent);
        
        EXPECT_EQ(p2.status, 0);
        EXPECT_EQ(p2.udd_afix_handle, -1);
        EXPECT_EQ(p2.structure, nullptr);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}