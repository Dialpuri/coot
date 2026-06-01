#pragma once

#include <coot/geometry/energy-lib.hh>
#include <coot/geometry/hb-types.hh>
#include <gemmi/cifdoc.hpp>
#include <string>
#include <cstdlib>

namespace coot {

inline void add_energy_lib_atoms_gemmi(energy_lib_t& lib, gemmi::cif::Loop& loop) {
    for (size_t j = 0; j < loop.length(); j++) {
        std::string type;
        double weight = -1;
        hb_t hb_type = HB_UNASSIGNED;
        double vdw_radius = -1;
        double vdwh_radius = -1;
        double ion_radius = -1;
        std::string element;
        int valency = -1;
        int sp_hybridisation = -1;

        int ierr = 0;
        int ierr_tot = 0;

        // type (required — contributes to ierr_tot)
        if (loop.has_tag("type")) {
            int col = loop.find_tag("type");
            const std::string& v = loop.val(j, col);
            if (!gemmi::cif::is_null(v)) {
                type = gemmi::cif::as_string(v);
            } else {
                ierr = 1;
            }
        } else {
            ierr = 1;
        }
        ierr_tot += ierr;
        ierr = 0;

        // weight (optional — failure does NOT contribute to ierr_tot)
        if (loop.has_tag("weight")) {
            int col = loop.find_tag("weight");
            const std::string& v = loop.val(j, col);
            if (!gemmi::cif::is_null(v)) {
                try {
                    weight = std::stod(v);
                } catch (...) {
                    weight = -1;
                }
            }
        }

        // hb_type (required — contributes to ierr_tot)
        if (loop.has_tag("hb_type")) {
            int col = loop.find_tag("hb_type");
            const std::string& v = loop.val(j, col);
            if (!gemmi::cif::is_null(v)) {
                std::string ss = gemmi::cif::as_string(v);
                if (ss == "D")
                    hb_type = HB_DONOR;
                if (ss == "A")
                    hb_type = HB_ACCEPTOR;
                if (ss == "B")
                    hb_type = HB_BOTH;
                if (ss == "N")
                    hb_type = HB_NEITHER;
                if (ss == "H")
                    hb_type = HB_HYDROGEN;
            } else {
                ierr = 1;
            }
        } else {
            ierr = 1;
        }
        ierr_tot += ierr;
        ierr = 0;

        // vdw_radius (optional)
        if (loop.has_tag("vdw_radius")) {
            int col = loop.find_tag("vdw_radius");
            const std::string& v = loop.val(j, col);
            if (!gemmi::cif::is_null(v)) {
                try {
                    vdw_radius = std::stod(v);
                } catch (...) {
                    vdw_radius = -1;
                }
            }
        }

        // vdwh_radius (optional)
        if (loop.has_tag("vdwh_radius")) {
            int col = loop.find_tag("vdwh_radius");
            const std::string& v = loop.val(j, col);
            if (!gemmi::cif::is_null(v)) {
                try {
                    vdwh_radius = std::stod(v);
                } catch (...) {
                    vdwh_radius = -1;
                }
            }
        }

        // ion_radius (optional)
        if (loop.has_tag("ion_radius")) {
            int col = loop.find_tag("ion_radius");
            const std::string& v = loop.val(j, col);
            if (!gemmi::cif::is_null(v)) {
                try {
                    ion_radius = std::stod(v);
                } catch (...) {
                    ion_radius = -1;
                }
            }
        }

        // element (required — contributes to ierr_tot)
        if (loop.has_tag("element")) {
            int col = loop.find_tag("element");
            const std::string& v = loop.val(j, col);
            if (!gemmi::cif::is_null(v)) {
                element = gemmi::cif::as_string(v);
            } else {
                ierr = 1;
            }
        } else {
            ierr = 1;
        }
        ierr_tot += ierr;
        ierr = 0;

        // valency (required — contributes to ierr_tot)
        if (loop.has_tag("valency")) {
            int col = loop.find_tag("valency");
            const std::string& v = loop.val(j, col);
            if (!gemmi::cif::is_null(v)) {
                valency = gemmi::cif::as_int(v);
            } else {
                ierr = 1;
            }
        } else {
            ierr = 1;
        }
        ierr_tot += ierr;
        ierr = 0;

        // sp_hybridisation (optional)
        if (loop.has_tag("sp")) {
            int col = loop.find_tag("sp");
            const std::string& v = loop.val(j, col);
            if (!gemmi::cif::is_null(v)) {
                sp_hybridisation = gemmi::cif::as_int(v);
            }
        }

        if (ierr_tot == 0) {
            energy_lib_atom at(type, hb_type, weight, vdw_radius, vdwh_radius,
                              ion_radius, element, valency, sp_hybridisation);
            lib.add_energy_lib_atom(at);
        }
    }
}

} // namespace coot