#pragma once
#include <ostream>
#include <string>
#include <gemmi/model.hpp>

// Pad an atom name to 4 chars the PDB way: leading space, then name, then trailing spaces
// e.g. "N"  -> " N  "
//      "CA" -> " CA "
inline std::string pad_atom_name(const std::string& name) {
    std::string padded = " " + name;
    while (padded.size() < 4)
        padded += ' ';
    return padded.substr(0, 4);
}

/**
 * gemmi port of the MMDB operator<<(std::ostream&, mmdb::Atom&).
 *
 * In MMDB the atom object carries parent context (chain, residue, model).
 * In gemmi the atom is standalone, so the caller must provide a CRA that
 * carries the chain pointer, residue pointer, and atom pointer. We also
 * pass model_num (0-based index; the original printed 1-based).
 */
inline std::ostream& stream_atom_gemmi(
    std::ostream& s,
    gemmi::CRA cra,
    int model_num)   // 0-based model index; original printed model_num+1
{
    // Guard against nulls
    if (!cra.atom || !cra.residue) {
        s << "(null atom)";
        return s;
    }

    const gemmi::Atom&    atom     = *cra.atom;
    const gemmi::Residue& residue  = *cra.residue;

    // Model number: original used 1-based (GetModelNum returned 1 for the first model)
    int printed_model = model_num + 1;

    // Chain ID
    const char* chain_id = "(null)";
    if (cra.chain)
        chain_id = cra.chain->name.c_str();

    // Insertion code — MMDB returned "" when none, gemmi stores ' '
    // Only print icode if it's not the "empty" space character
    char icode = residue.seqid.icode;

    // Atom name: MMDB used a 4-char padded name; gemmi's atom.name is unpadded.
    // Pad to match PDB convention (space + name + trailing spaces).
    std::string padded_name = pad_atom_name(atom.name);

    // altloc: gemmi stores ' ' for no alternate location; MMDB used '\0'
    // so we print nothing when gemmi's altloc is ' ' to match MMDB output
    std::string altloc_str = (atom.altloc != ' ' && atom.altloc != '\0')
                              ? std::string(1, atom.altloc) : "";

    s << printed_model << "/" << chain_id << "/"
      << residue.seqid.num.value
      << (icode != ' ' ? std::string(1, icode) : "")
      << "/" << residue.name
      << "/" << padded_name << " altLoc :" << altloc_str << ": pos: ("
      << atom.pos.x << "," << atom.pos.y << "," << atom.pos.z
      << ") B-factor: " << atom.b_iso;

    return s;
}