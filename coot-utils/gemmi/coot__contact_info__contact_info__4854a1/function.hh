#pragma once

#include <vector>
#include <utility>

namespace coot {

struct contacts_pair {
    int id1;
    int id2;
    contacts_pair(int i1, int i2) : id1(i1), id2(i2) {}
};

class contact_info {
public:
    std::vector<contacts_pair> contacts;

    // Gemmi-port of the original constructor that took mmdb::Contact*
    static contact_info contact_info_gemmi(std::vector<std::pair<int, int>> const& con_in) {
        contact_info ci;
        for (const auto& con : con_in) {
            ci.contacts.push_back(contacts_pair(con.first, con.second));
        }
        return ci;
    }

    int n_contacts() const {
        return static_cast<int>(contacts.size());
    }
};

} // namespace coot