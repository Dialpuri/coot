#include "function.hh"
#include <iostream>

namespace coot {

void add_chem_links_gemmi(
    gemmi::cif::Loop& loop,
    std::map<unsigned int, std::vector<chem_link>>& chem_link_map) 
{
    for (size_t j = 0; j < loop.length(); j++) {
        std::string chem_link_id;
        std::string chem_link_comp_id_1;
        std::string chem_link_mod_id_1;
        std::string chem_link_group_comp_1;
        std::string chem_link_comp_id_2;
        std::string chem_link_mod_id_2;
        std::string chem_link_group_comp_2;
        std::string chem_link_name;

        int ierr_tot = 0;

        // Get column indices
        int col_id = loop.find_tag("id");
        int col_comp_id_1 = loop.find_tag("comp_id_1");
        int col_mod_id_1 = loop.find_tag("mod_id_1");
        int col_group_comp_1 = loop.find_tag("group_comp_1");
        int col_comp_id_2 = loop.find_tag("comp_id_2");
        int col_mod_id_2 = loop.find_tag("mod_id_2");
        int col_group_comp_2 = loop.find_tag("group_comp_2");
        int col_name = loop.find_tag("name");

        // Extract values - in gemmi, val() returns const std::string&
        // We consider a field missing if its column index is -1
        if (col_id >= 0) chem_link_id = loop.val(j, col_id);
        else {
            std::cout << "WARNING add_chem_links error getting field \"id\"" << std::endl;
            ierr_tot++;
        }

        if (col_comp_id_1 >= 0) chem_link_comp_id_1 = loop.val(j, col_comp_id_1);
        else {
            std::cout << "WARNING add_chem_links error getting field \"comp_id_1\"" << std::endl;
            ierr_tot++;
        }

        if (col_mod_id_1 >= 0) chem_link_mod_id_1 = loop.val(j, col_mod_id_1);
        else {
            std::cout << "WARNING add_chem_links error getting field \"mod_id_1\"" << std::endl;
            ierr_tot++;
        }

        if (col_group_comp_1 >= 0) chem_link_group_comp_1 = loop.val(j, col_group_comp_1);
        else {
            std::cout << "WARNING add_chem_links error getting field \"group_comp_1\"" << std::endl;
            ierr_tot++;
        }

        if (col_comp_id_2 >= 0) chem_link_comp_id_2 = loop.val(j, col_comp_id_2);
        else {
            std::cout << "WARNING add_chem_links error getting field \"comp_id_2\"" << std::endl;
            ierr_tot++;
        }

        if (col_mod_id_2 >= 0) chem_link_mod_id_2 = loop.val(j, col_mod_id_2);
        else {
            std::cout << "WARNING add_chem_links error getting field \"mod_id_2\"" << std::endl;
            ierr_tot++;
        }

        if (col_group_comp_2 >= 0) chem_link_group_comp_2 = loop.val(j, col_group_comp_2);
        else {
            std::cout << "WARNING add_chem_links error getting field \"group_comp_2\"" << std::endl;
            ierr_tot++;
        }

        if (col_name >= 0) chem_link_name = loop.val(j, col_name);
        else {
            std::cout << "WARNING add_chem_links error getting field \"name\"" << std::endl;
            ierr_tot++;
        }

        if (ierr_tot == 0) {
            chem_link clink(chem_link_id,
                           chem_link_comp_id_1, chem_link_mod_id_1, chem_link_group_comp_1,
                           chem_link_comp_id_2, chem_link_mod_id_2, chem_link_group_comp_2,
                           chem_link_name);

            chem_link_map[clink.get_hash_code()].push_back(clink);
        } else {
            std::cout << "WARNING:: an error occurred when trying to add link: "
                      << "\"" << chem_link_id << "\" "
                      << "\"" << chem_link_comp_id_1 << "\" "
                      << "\"" << chem_link_mod_id_1 << "\" "
                      << "\"" << chem_link_group_comp_1 << "\" "
                      << "\"" << chem_link_comp_id_2 << "\" "
                      << "\"" << chem_link_mod_id_2 << "\" "
                      << "\"" << chem_link_group_comp_2 << "\" "
                      << "\"" << chem_link_name << "\" "
                      << std::endl;
        }
    }
}

} // namespace coot
