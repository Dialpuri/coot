#include <gtest/gtest.h>
#include <gemmi/pdb.hpp>
#include <gemmi/model.hpp>
#include <filesystem>
#include <fstream>
#include "function.hh"

TEST(OracleTest, MakeBackup) {
   // Load the PDB file with gemmi
   gemmi::Structure st = gemmi::read_pdb_file("/lmb/home/jdialpuri/Development/coot-tooling/test-data/example.pdb");

   // Create our standalone modification_info_t (equivalent to coot::molecule_t::modification_info_t)
   coot::modification_info_t info;
   info.mol_name = "test";
   info.backup_dir = "coot-backup";

   // Clean up any previous backup directory
   std::filesystem::remove_all("coot-backup");

   // Test 1: Normal backup with valid molecule
   {
      std::cout << "=== Test 1: Normal backup with valid molecule ===" << std::endl;
      std::cout << "INPUT mol: valid pointer" << std::endl;
      std::cout << "INPUT modification_info_string: \"test backup\"" << std::endl;
      std::cout << "INPUT save_info.size() before: " << info.save_info.size() << std::endl;

      std::string result = info.make_backup_gemmi(st, "test backup");
      std::cout << "OUTPUT result: \"" << result << "\"" << std::endl;
      std::cout << "OUTPUT save_info.size() after: " << info.save_info.size() << std::endl;
      std::cout << "OUTPUT modification_index: " << info.modification_index << std::endl;

      EXPECT_EQ(result, "");
      EXPECT_EQ(info.save_info.size(), 1);
      EXPECT_EQ(info.modification_index, 1);
   }

   // Test 2: Null molecule (empty structure)
   {
      std::cout << "\n=== Test 2: Null molecule ===" << std::endl;
      std::cout << "INPUT mol: nullptr (empty structure)" << std::endl;
      std::cout << "INPUT modification_info_string: \"null test\"" << std::endl;

      gemmi::Structure empty_st;  // empty = no models

      std::string result = info.make_backup_gemmi(empty_st, "null test");
      std::cout << "OUTPUT result: \"" << result << "\"" << std::endl;

      EXPECT_EQ(result, "null molecule");
      // size should remain unchanged
      EXPECT_EQ(info.save_info.size(), 1);
   }

   // Test 3: Multiple backups
   {
      std::cout << "\n=== Test 3: Multiple backups ===" << std::endl;
      std::cout << "INPUT save_info.size() before: " << info.save_info.size() << std::endl;

      info.make_backup_gemmi(st, "backup 2");
      info.make_backup_gemmi(st, "backup 3");
      std::cout << "OUTPUT save_info.size() after 2 more: " << info.save_info.size() << std::endl;
      std::cout << "OUTPUT modification_index: " << info.modification_index << std::endl;

      EXPECT_EQ(info.save_info.size(), 3);
      EXPECT_EQ(info.modification_index, 3);
   }

   // Test 4: Verify backup files were created
   {
      std::cout << "\n=== Test 4: Backup files ===" << std::endl;
      std::cout << "OUTPUT backup_dir: " << info.backup_dir << std::endl;
      std::cout << "OUTPUT save_info[0].file_name: " << info.save_info[0].file_name << std::endl;
      std::cout << "OUTPUT save_info[0].modification_info_string: " << info.save_info[0].modification_info_string << std::endl;

      EXPECT_EQ(info.backup_dir, "coot-backup");
      EXPECT_NE(info.save_info[0].file_name.find("coot-backup/test-0.pdb"), std::string::npos);
      EXPECT_EQ(info.save_info[0].modification_info_string, "test backup");
   }

   // Verify backup files actually exist on disk
   {
      std::cout << "\n=== Test 5: Backup files exist on disk ===" << std::endl;
      for (const auto& entry : info.save_info) {
         std::cout << "Checking file: " << entry.file_name << std::endl;
         EXPECT_TRUE(std::filesystem::exists(entry.file_name))
            << "Backup file does not exist: " << entry.file_name;
      }
   }

   // Clean up
   std::filesystem::remove_all("coot-backup");
}

int main(int argc, char** argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}