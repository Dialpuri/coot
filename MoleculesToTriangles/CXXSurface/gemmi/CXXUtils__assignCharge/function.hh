#pragma once

#include <gemmi/model.hpp>
#include <map>
#include <string>
#include <vector>

typedef std::map<std::string, double> stringDoubleMap;
typedef std::map<std::string, stringDoubleMap> stringStringDoubleMap;

class CXXChargeTable {
private:
    stringStringDoubleMap atomNameToChargeMap;
    int addResidueMap(std::string residueName, stringDoubleMap theAtomCharges);
    int addChargeToResidueMap(stringDoubleMap& theMap, std::string theAtomName, double theCharge);
public:
    CXXChargeTable();
    double getCharge(std::string residueName, std::string atomName);
};

namespace CXXUtils {
    int assignCharge_gemmi(std::vector<gemmi::CRA> selected, CXXChargeTable *theChargeTable, std::map<const gemmi::Atom*, double>& charges_out);
}