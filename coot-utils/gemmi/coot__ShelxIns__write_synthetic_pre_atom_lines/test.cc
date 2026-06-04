#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <fstream>
#include <sstream>
#include <string>
#include "function.hh"

TEST(OracleTest, write_synthetic_pre_atom_lines) {
    // Load PDB
    gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");
    ASSERT_EQ(st.models.size(), 1u) << "Expected at least one model";

    // Set up cell and space group properly on the gemmi structure
    st.cell = gemmi::UnitCell(54.0, 70.0, 80.0, 90.0, 90.0, 90.0);
    st.spacegroup_hm = "P 21 21 21";

    EXPECT_EQ(st.spacegroup_hm, "P 21 21 21");

    // --- Case 1: with cell info (directly set have_cell_flag + cell) ---
    {
        coot::ShelxIns shelx_ins;

        // Copy cell into ShelxIns
        clipper::Cell shelx_cell(clipper::Cell_descr(54.0, 70.0, 80.0, 90.0, 90.0, 90.0));
        shelx_ins.cell = shelx_cell;
        shelx_ins.have_cell_flag = true;

        EXPECT_TRUE(shelx_ins.have_cell_flag);
        EXPECT_NEAR(shelx_cell.descr().a(), 54.0, 1e-4);

        int sfac_size_before = static_cast<int>(shelx_ins.sfac.size());
        EXPECT_EQ(sfac_size_before, 0);

        const char* tmpfile = "/tmp/oracle_shelx_test1.ins";
        std::ofstream f(tmpfile);
        ASSERT_TRUE(f.is_open()) << "Failed to open temp file";

        shelx_ins.write_synthetic_pre_atom_lines_gemmi(st, f);
        f.close();

        int sfac_size_after = static_cast<int>(shelx_ins.sfac.size());
        EXPECT_EQ(sfac_size_after, 4);

        std::string sfac_contents;
        for (size_t i = 0; i < shelx_ins.sfac.size(); i++) {
            if (i > 0) sfac_contents += " ";
            sfac_contents += shelx_ins.sfac[i];
        }
        EXPECT_EQ(sfac_contents, "C N O S");

        std::ifstream fin(tmpfile);
        std::string written((std::istreambuf_iterator<char>(fin)),
                            std::istreambuf_iterator<char>());

        EXPECT_FALSE(written.empty());
        EXPECT_GT(written.size(), 100); // Large output string
        EXPECT_NE(written.find("TITL"), std::string::npos);
        EXPECT_NE(written.find("CELL"), std::string::npos);
        EXPECT_NE(written.find("ZERR"), std::string::npos);
        EXPECT_NE(written.find("LATT"), std::string::npos);
        EXPECT_NE(written.find("SYMM"), std::string::npos);
        EXPECT_NE(written.find("SFAC"), std::string::npos);
        EXPECT_NE(written.find("UNIT"), std::string::npos);
        EXPECT_NE(written.find("CGLS"), std::string::npos);
        EXPECT_NE(written.find("SHEL"), std::string::npos);
        EXPECT_NE(written.find("FMAP"), std::string::npos);
        EXPECT_NE(written.find("PLAN"), std::string::npos);
        EXPECT_NE(written.find("LIST"), std::string::npos);
        EXPECT_NE(written.find("WPDB"), std::string::npos);
        EXPECT_NE(written.find("CONN"), std::string::npos);
        EXPECT_NE(written.find("BUMP"), std::string::npos);
        EXPECT_NE(written.find("54.0000"), std::string::npos);
        EXPECT_NE(written.find("70.0000"), std::string::npos);
        EXPECT_NE(written.find("80.0000"), std::string::npos);
    }

    // --- Case 2: Without cell flag ---
    {
        coot::ShelxIns shelx_ins2;
        EXPECT_FALSE(shelx_ins2.have_cell_flag);

        const char* tmpfile2 = "/tmp/oracle_shelx_test2.ins";
        std::ofstream f2(tmpfile2);
        ASSERT_TRUE(f2.is_open()) << "Failed to open temp file 2";
        shelx_ins2.write_synthetic_pre_atom_lines_gemmi(st, f2);
        f2.close();

        std::ifstream fin2(tmpfile2);
        std::string written2((std::istreambuf_iterator<char>(fin2)),
                             std::istreambuf_iterator<char>());

        EXPECT_FALSE(written2.empty());
        EXPECT_NE(written2.find("TITL"), std::string::npos);
        EXPECT_NE(written2.find("CGLS"), std::string::npos);
        EXPECT_NE(written2.find("SHEL"), std::string::npos);
        EXPECT_NE(written2.find("FMAP"), std::string::npos);
        EXPECT_NE(written2.find("PLAN"), std::string::npos);
        EXPECT_NE(written2.find("LIST"), std::string::npos);
        EXPECT_NE(written2.find("WPDB"), std::string::npos);
        EXPECT_NE(written2.find("CONN"), std::string::npos);
        EXPECT_NE(written2.find("BUMP"), std::string::npos);

        // Should NOT have CELL-related lines when have_cell_flag is false
        EXPECT_EQ(written2.find("CELL"), std::string::npos);
        EXPECT_EQ(written2.find("ZERR"), std::string::npos);
        EXPECT_EQ(written2.find("LATT"), std::string::npos);
        EXPECT_EQ(written2.find("SYMM"), std::string::npos);
        EXPECT_EQ(written2.find("SFAC"), std::string::npos);
        EXPECT_EQ(written2.find("UNIT"), std::string::npos);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}