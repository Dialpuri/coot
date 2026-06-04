#pragma once

#include <ostream>
#include <gemmi/model.hpp>

namespace {
    std::string pad_atom_name(const std::string& name) {
        // PDB-style: 1 leading space + name + trailing spaces to fill 4 chars
        std::string result = " " + name;
        while (result.size() < 4) result += " ";
        return result.substr(0, 4);
    }
}

std::ostream& stream_cra_gemmi(std::ostream& s, gemmi::CRA cra, int model_num) {
    if (cra.atom) {
        s << model_num << "/" << cra.chain->name << "/"
          << cra.residue->seqid.num.value;
        if (cra.residue->seqid.icode != ' ') {
            s << cra.residue->seqid.icode;
        }
        s << " {" << cra.residue->name << "}/"
          << pad_atom_name(cra.atom->name) << " altLoc :";
        if (cra.atom->altloc != '\0') {
            s << cra.atom->altloc;
        }
        s << ": segid :: pos: ("
          << cra.atom->pos.x << "," << cra.atom->pos.y << "," << cra.atom->pos.z
          << ") B-factor: " << cra.atom->b_iso;
    } else {
        s << "NULL";
    }
    return s;
}