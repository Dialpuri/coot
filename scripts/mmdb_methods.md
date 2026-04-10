# MMDB2 API Reference (AST-driven)

---

## API Index

Okay, here's a concise API index for the MMDB2 library based on the documentation you provided, incorporating the requested format and focusing on essential information for Coot developers.

**MMDB2 API Index**

**Classes**

*   `Alignment` (header.h) - Represents a sequence alignment, facilitating the alignment of two sequences based on scoring parameters and potentially affine gap penalties.
    *   `Alignment::Alignment()` - Constructs an empty alignment object.
    *   `Alignment::Alignment(io::RPStream Object)` - Constructs an alignment object from an existing `io::RPStream` object.
    *   `Alignment::~Alignment()` - Destroys the alignment object, releasing allocated memory.
    *   `Alignment::SetAffineModel(realtype WGap, realtype WSpace)` - Sets the affine gap penalty parameters.
    *   `Alignment::SetScores(realtype SEqual, realtype SNEqual)` - Sets the scoring parameters for match and mismatch.
    *   `Alignment::Align(cpstr S, cpstr T, ALIGN_METHOD Method=ALIGN_GLOBAL)` - Aligns the sequences `S` and `T`.
    *   `Alignment::GetAlignedS() -> pstr` - Returns a pointer to the aligned sequence `S`.
    *   `Alignment::GetAlignedT() -> pstr` - Returns a pointer to the aligned sequence `T`.
    *   `Alignment::GetScore() -> realtype` - Returns the achieved alignment score.
    *   `Alignment::GetSpace() -> char` - Returns the space character used in the alignment.
    *   `Alignment::GetSimilarity() -> realtype` - Calculates a similarity score.
    *   `Alignment::GetSeqId() -> realtype` - Returns the sequence ID.
    *   `Alignment::GetNAlign() -> int` - Returns the number of aligned residues.
    *   `Alignment::OutputResults(io::RFile f, cpstr S, cpstr T)` - Outputs the alignment results to a file.
    *   `Alignment::read(io::RFile f)` - Reads alignment data from a file.
    *   `Alignment::write(io::RFile f)` - Writes alignment data to a file.
    *   `Alignment::InitAlignment()` - Initializes the alignment data structures.
    *   `Alignment::FreeMemory()` - Frees any dynamically allocated memory used by the alignment object.
    *   `Alignment::Score(char A, char B) -> realtype` - Calculates the score for a match or mismatch between two characters.
    *   `Alignment::BuildGATable(cpstr S, cpstr T, bool FreeSEnd, bool FreeTEnd)` - Builds a gap opening/extension table.
    *   `Alignment::BuildLATable(cpstr S, cpstr T)` - Builds a local alignment table.
    *   `Alignment::Backtrace(cpstr S, cpstr T, int J, int I, bool StopAtZero)` - Performs the traceback to reconstruct the aligned sequences.
    *   `Alignment::AdjustEnds(cpstr S, cpstr T, int J, int I)` - Adjusts the ends of the aligned sequences.
    *   `Alignment::PrintVT(cpstr S, cpstr T)` - Prints the VT table for debugging.

*   `Atom` (header.h) - Represents an atom in a structural biology file, storing its coordinates, bond information, and related data.
    *   `Atom::AddBond(PAtom bond_atom, int bond_order, int nAdd_bonds=1)` - Adds a bond to the atom.
    *   `Atom::Copy()` - Creates a deep copy of the `Atom` object.
    *   `Atom::GetAtomID() -> pstr` - Generates a unique string identifier for the atom.
    *   `Atom::GetAtomName() -> pstr` - Returns the atom's name.
    *   `Atom::GetBonds() -> void` - Returns a pointer to the internal structure representing the atom's bonds.
    *   `Atom::GetChain() -> pstr` - Returns the name of the chain to which the atom belongs.
    *   `Atom::GetChainID() -> pstr` - Returns the chain ID within the chain.
    *   `Atom::GetElementName() -> pstr` - Returns the element symbol associated with the atom.
    *   `Atom::GetIndex() -> int` - Returns the atom's index within the residue.
    *   `Atom::GetInsCode() -> pstr` - Returns the insertion code for the atom.
    *   `Atom::GetModel() -> pstr` - Returns the model number to which the atom belongs.
    *   `Atom::GetNBonds() -> int` - Returns the number of bonds associated with the atom.
    *   `Atom::GetResName() -> pstr` - Returns the name of the residue to which the atom belongs.
    *   `Atom::GetResidue() -> PResidue` - Returns a pointer to the `Residue` object to which the atom belongs.
    *   `Atom::GetSeqNum() -> int` - Returns the sequence number of the residue to which the atom belongs.
    *   `Atom::GetUDData() -> int` - Returns the user-defined data mask indicating which fields of the atom are set.
    *   `Atom::MakeTer() -> void` - Marks the atom as a chain terminator (ter).
    *   `Atom::PutUDData(int UDDhandle, int iudd)` - Writes user-defined data to a specific handle.
    *   `Atom::PutUDData(int UDDhandle, realtype rudd)` - Writes a real-valued user-defined data to a specific handle.
    *   `Atom::PutUDData(int UDDhandle, cpstr sudd)` - Writes a string user-defined data to a specific handle.
    *   `Atom::PutUDData(int UDDhandle, int& iudd)` - Reads user-defined data from a specific handle.
    *   `Atom::PutUDData(int UDDhandle, realtype& rudd)` - Reads a real-valued user-defined data from a specific handle.
    *   `Atom::PutUDData(int UDDhandle, pstr& sudd)` - Reads a string user-defined data from a specific handle.
    *   `Atom::InitAlignment()` - Initializes the alignment data structures.
    *   `Atom::FreeMemory()` - Frees any dynamically allocated memory used by the atom.
    *   `Atom::Score(char A, char B) -> realtype` - Calculates the score for a match or mismatch between two characters.
    *   `Atom::BuildGATable(cpstr S, cpstr T, bool FreeSEnd, bool FreeTEnd)` - Builds a gap opening/extension table.
    *   `Atom::BuildLATable(cpstr S, cpstr T)` - Builds a local alignment table.
    *   `Atom::Backtrace(cpstr S, cpstr T, int J, int I, bool StopAtZero)` - Performs the traceback to reconstruct the aligned sequences.
    *   `Atom::AdjustEnds(cpstr S, cpstr T, int J, int I)` - Adjusts the ends of the aligned sequences after the traceback.
    *   `Atom::PrintVT(cpstr S, cpstr T)` - Prints the VT table to the console for debugging.

*   `Category` (header.h) - Manages tags and their associated data within a mmCIF file.
    *   `Category::GetCategoryName() -> pstr` - Returns the category name.
    *   `Category::GetNofTags() -> int` - Returns the total number of tags in the category.
    *   `Category::GetTag(int tagNo) -> pstr` - Returns the tag at the specified serial number.
    *   `Category::AddTag(cpstr ttag) -> int` - Adds a tag to the category.
    *   `Category::GetTagNo(cpstr ttag) -> int` - Returns the serial number of a tag.
    *   `Category::PrintTags()` - Prints the list of tags to standard output.
    *   `Category::CheckTags(cpstr * tagList)` - Checks if all tags from the list are found in the category.
    *   `Category::Copy(PCategory Category)` - Performs a deep copy of the category.
    *   `Category::WriteMMCIF(io::RFile)` - Writes category's content into mmCIF file.
    *   `Category::ReadMMCIF(io::RFile)` - Reads category's content from mmCIF file.
    *   `Category::Optimize()` - Optimizes the data structures within the category.
    *   `Category::Sort()` - Sorts the category's data.
    *   `Category::SetCategoryName(cpstr N)` - Sets the category name.
    *   `Category::GetCategoryID() -> MMCIF_ITEM` - Returns the category type.

*   `Chain` (header.h) - Represents a chain of atoms in a structural biology model.
    *   `Chain::AddResidue(PResidue res) -> int` - Adds a residue to the chain.
    *   `Chain::Copy()` - Creates a complete copy of the chain.
    *   `Chain::DeleteResidue(int resNo) -> int` - Deletes a residue from the chain.
    *   `Chain::GetChainID() -> cpstr` - Returns the chain ID as a string.
    *   `Chain::GetCoordHierarchy() -> PManager` - Returns a pointer to the coordinate hierarchy for this chain.
    *   `Chain::ConvertPDBASCII(cpstr S)` - Converts a PDB ASCII string to a `Chain` object.
    *   `Chain::GetCIF(mmcif::PData CIF, int & n)` - Retrieves the CIF data for the chain.
    *   `Chain::GetClassID() -> CLASS_ID` - Returns the class ID for the `Chain` object.
    *   `Chain::MakeCIF(mmcif::PData CIF)` - Creates a mmCIF data structure for the chain.
    *   `Chain::Copy(PChain Chain)` - Copies the data from another `Chain` object.
    *   `Chain::write(io::RFile f)` - Writes the chain data to an `io::RFile` object.
    *   `Chain::read(io::RFile f)` - Reads the chain data from an `io::RFile` object.
    *   `Chain::InitChain()` - Initializes the chain object.
    *   `Chain::FreeMemory()` - Frees any dynamically allocated memory used by the chain.
    *   `Chain::SortResidues() -> void` - Sorts the residues in the chain.
    *   `Chain::TrimResidueTable() -> void` - Trims the residue table.
    *   `Chain::isAminoacidChain() -> bool` - Checks if the chain contains only amino acid residues.
    *   `Chain::isInSelection(int selHnd) -> bool` - Checks if the chain is present in a selection.
    *   `Chain::isNucleotideChain() -> bool` - Checks if the chain contains only nucleotide residues.
    *   `Chain::isSolventChain() -> bool` - Checks if the chain contains only solvent residues.
    *   `Chain::isNTerminus() -> bool` - Checks if the chain is the N-terminus.
    *   `Chain::isCTerminus() -> bool` - Checks if the chain is the C-terminus.

*   `CoorManager` (header.h) - Manages structural biology data, coordinates, and related structures.
    *   `CoorManager::MakeBonds(bool calc_only) -> void` - Calculates and creates bonds based on the current structure.
    *   `CoorManager::RemoveBonds() -> void` - Removes all currently existing bonds from the structure.

*   `Data` (header.h) - Represents a data block in an MMDB2 file.
    *   `Data::Data()` - Basic constructor.
    *   `Data::Data(cpstr N)` - Constructor that assigns data block name.
    *   `Data::Data(io::RPStream Object)` - Constructor for MMDB data streaming functions.
    *   `Data::~Data()` - Destructor.
    *   `Data::SetPrintWarnings(bool SPW)` - Sets flag to print warnings on reading mmCIF files.
    *   `Data::SetStopOnWarning(bool SOW)` - Sets flag to stop on warning when reading an mmCIF file.
    *   `Data::SetFlag(word Flag)` - Sets a flag in the structure.
    *   `Data::DeleteFlag(word Flag)` - Deletes the flag in the structure.
    *   `Data::GetFlag(word Flag)` - Gets the flag in the structure.
    *   `Data::PutData(int NoDataType, cpstr TName, cpstr S)` - Sets string value for the specified tag name.
    *   `Data::PutData(int NoDataType, realtype R, cpstr TName)` - Sets floating-point value for the specified tag name.
    *   `Data::PutData(int NoDataType, int I, cpstr TName)` - Sets integer value for the specified tag name.
    *   `Data::GetField(int tagNo, int FType, cpstr TName)` - Returns the value of the field corresponding to the specified tag number.
    *   `Data::GetField(int tagNo, realtype& R, cpstr TName)` - Returns the value of the field associated with the given tag number.
    *   `Data::GetField(int tagNo, int &I, cpstr TName)` - Returns the value of the field associated with the given tag number.
    *   `Data::GetField(int tagNo, cpstr &S, cpstr TName)` - Returns the value of the field associated with the given tag number.
    *   `Data::WriteMMCIFData(cpstr FName)` - Writes the data in mmCIF format into the given file.
    *   `Data::ReadMMCIFData(cpstr FName)` - Reads the data from the given file.
    *   `Data::InitData()` - Initializes the data block.
    *   `Data::FreeMemory()` - Frees any memory allocated by the data block.
    *   `Data::Sort()` - Sorts the data.
    *   `Data::GetNumberOfCategories() -> int` - Returns the number of categories in the data block.
    *   `Data::GetCategory(int categoryNo)` - Returns a pointer to the specified category.
    *   `Data::GetCategoryName(cpstr & categoryName)` - Returns the name of the specified category.
    *   `Data::GetSequenceNum(cpstr TName)` - Returns the sequence number for the specified tag name.
    *   `Data::GetDataName(pstr & dname)` - Returns the data block name.
    *   `Data::CheckData(cpstr C, cpstr T)` - Checks if a field exists in the data block.
    *   `Data::AddCategory(cpstr cname)` - Adds a category to the data block.
    *   `Data::DeleteCategory(cpstr cname)` - Deletes a category from the data block.
    *   `Data::PutCategoryName(cpstr n)` - Sets the category name.
    *   `Data::WriteMMCIF(io::RFile f)` - Writes the data to the mmCIF file.
    *   `Data::ReadMMCIF(io::RFile f)` - Reads the data from the mmCIF file.
    *   `Data::MemRead(cpstr S, int & count)` - Reads a line of input from a memory string.
    *   `Data::MemWrite(pstr S, int & count)` - Writes a string to a memory string.

*   `Graph` (header.h) - Represents a graph of atoms in a structural biology file.
    *   `Graph::AddVertex(PVertex V)` - Adds a vertex to the graph.
    *   `Graph::AddEdge(PEdge G)` - Adds an edge to the graph.
    *   `Graph::SetVertex(int index, PVertex V)` - Sets the vertex at the given index.
    *   `Graph::DeleteVertex(int index)` - Deletes the vertex at the given index.
    *   `Graph::SetEdge(int index, PEdge G)` - Sets the edge at the given index.
    *   `Graph::DeleteEdge(int index)` - Deletes the edge at the given index.
    *   `Graph::MakeGraph(PResidue R, cpstr altLoc=NULL)` - Creates a graph representing the residue.
    *   `Graph::MakeGraph(PPAtom a, int n)` - Creates a graph from a list of atoms.
    *   `Graph::Copy(PGraph G)` - Copies the contents of another graph.
    *   `Graph::write(io::RFile f)` - Writes the graph data to an RFile.
    *   `Graph::read(io::RFile f)` - Reads the graph data from an RFile.
    *   `Graph::InitGraph()` - Initializes the graph object.
    *   `Graph::FreeMemory()` - Frees the memory allocated for the graph.
    *   `Graph::Sort()` - Sorts the graph's data.

*   `GraphMatch` (header.h) - This class implements an algorithm for finding maximal common subgraphs between two graphs.
    *   `GraphMatch::SetTimeLimit(int maxTimeToRun=0)` - Sets a maximum time limit for the matching process.
    *   `GraphMatch::GetMatch(int matchNo, ivector & FV1, ivector & FV2, int & nv, realtype & p1, realtype & p2)` - Returns the vertices of a specific matching subgraph.
    *   `GraphMatch::MatchGraphs(PGraph Gh1, PGraph Gh2, int minMatch, bool vertexType, VERTEX_EXT_TYPE vertexExt)` - Finds maximal common subgraphs of a minimum size between two graphs.
    *   `GraphMatch::GetNofMatches() -> int` - Returns the number of found matches.
    *   `GraphMatch::GetStopSignal() -> bool` - Returns the stop signal status.
    *   `GraphMatch::Print()` - Prints the GraphMatch object's information to the standard output stream.

*   `Helix` (header.h) - Represents a helix in a protein structure.
    *   `Helix::Helix()` - Constructs a new `Helix` object.
    *   `Helix::Helix(io::RPStream Object)` - Constructs a new `Helix` object from an `io::RPStream` object.
    *   `Helix::~Helix()` - Destroys the `Helix` object, releasing any allocated memory.
    *   `Helix::PDBASCIIDump(io::RFile f)` - Writes the helix data to PDB ASCII.
    *   `Helix::MakeCIF(mmcif::PData CIF)` - Generates a mmCIF data structure.
    *   `Helix::ConvertPDBASCII(cpstr S)` - Converts a PDB ASCII string to a `Helix` object.
    *   `Helix::GetCIF(mmcif::PData CIF, int & n)` - Retrieves the CIF data for the helix.
    *   `Helix::GetClassID() -> CLASS_ID` - Returns the class ID for the `Helix` object.
    *   `Helix::Copy(PSheet Sheet)` - Copies the data from another `Helix` object.
    *   `Helix::write(io::RFile f)` - Writes the helix data to an `io::RFile` object.
    *   `Helix::read(io::RFile f)` - Reads the helix data from an `io::RFile` object.
    *   `Helix::InitHelix()` - Initializes the `Helix` object.
    *   `Helix::FreeMemory()` - Frees memory allocated by the `Helix` object.
    *   `Helix::OrderSheet()` - Orders the strands within the `Helix` object.
    *   `Helix::GetStrand(int strand_no)` - Returns the strand at the given index.

*   `Link` (header.h) - Represents a single bond or link between two atoms in a structure.
    *   `Link::Link()` - Constructs a new `Link` object.
    *   `Link::Link(io::RPStream Object)` - Constructs a `Link` object from an `RPStream` object.
    *   `Link::~Link()` - Destroys the `Link` object, releasing allocated memory.
    *   `Link::PDBASCIIDump(pstr S, int N)` - Writes the `Link` data to a PDB ASCII string.
    *   `Link::MakeCIF(mmcif::PData CIF, int N)` - Creates a mmCIF data structure for the `Link`.
    *   `Link::ConvertPDBASCII(cpstr S)` - Converts a PDB ASCII string to a `Link` object.
    *   `Link::GetCIF(mmcif::PData CIF, int & n)` - Retrieves the CIF data for the `Link`.
    *   `Link::GetClassID() -> CLASS_ID` - Returns the class ID for the `Link` object.
    *   `Link::Copy(PLink Link)` - Copies the data from another `Link` object.
    *   `Link::write(io::RFile f)` - Writes the `Link` data to an `io::RFile` object.
    *   `Link::read(io::RFile f)` - Reads the `Link` data from an `io::RFile` object.
    *   `Link::InitLink()` - Initializes the `Link` object.
    *   `Link::FreeMemory()` - Frees any dynamically allocated memory used by the `Link` object.

*   `LinkContainer` (header.h) - This class represents a collection of linked structural elements.
    *   `LinkContainer::LinkContainer()` - Constructs a new `LinkContainer` object.
    *   `LinkContainer::LinkContainer(io::RPStream Object)` - Constructs a new `LinkContainer` object from an `io::RPStream` object.
    *   `LinkContainer::~LinkContainer()` - Destroys the `LinkContainer` object, releasing any associated resources.
    *   `LinkContainer::LinkContainer(PContainerClass LinkContainer)` - Copy the `LinkContainer` object.
    *   `LinkContainer::read(io::RFile f)` - Reads data from an `io::RFile` object.
    *   `LinkContainer::write(io::RFile f)` - Writes data to an `io::RFile` object.

*   `Root` (header.h) - The `Root` class provides a central interface for accessing and manipulating structural biology data in the MMDB2 format.
    *   `Root::CrystReady() -> bool` - Checks if crystallographic information is complete, noting any imprecise or translation-containing details.
    *   `Root::DeleteAllModels() -> int` - Deletes all models from the structure, resetting the structure to an empty state.
    *   `Root::FinishStructEdit() -> int` - Completes structural editing operations.
    *   `Root::GetAtomI(int index) -> PAtom` - Returns the atom at the given serial number index in the structure’s atom array.
    *   `Root::GetCell() -> mat44` - Returns the cell transformation matrix for the current structure.
    *   `Root::GetFileType() -> int` - Returns the file type of the currently loaded structure.
    *   `Root::GetInputBuffer(pstr Line, int & count) -> void` - Retrieves a line of input from the input buffer and stores it in the specified string.
    *   `Root::GetNumberOfSymOps() -> int` - Returns the number of symmetry operations defined in the structure.
    *   `Root::GetSpaceGroup() -> pstr` - Returns the name of the space group for the current structure.
    *   `Root::GetStructureTitle(pstr & L) -> void` - Retrieves the structure title and stores it in the specified string.
    *   `Root::GetSymOp(int Nop) -> pstr` - Returns the XYZ symmetry operation name corresponding to the given symmetry operation number.
    *   `Root::GetTMatrix(mat44 & TMatrix, int Nop) -> int` - Calculates and returns the coordinate transformation matrix for the given symmetry operation number.
    *   `Root::GetUDDHandle(UDR_TYPE udr_type, cpstr UDDataID) -> int` - Retrieves the handle for a user-defined data record.
    *   `Root::PDBCleanup(word CleanKey) -> word` - Performs PDB cleanup operations on the structure, ensuring compliance with PDB standards.
    *   `Root::PutAtom(int index, int serNum, const AtomName atomName, const ResName resName, const ChainID chainID, int seqNum, const InsCode insCode, const AltLoc altLoc, const SegID segID, const Element element) -> int` - Adds an atom to the structure at the specified index, setting its properties.
    *   `Root::PutPDBString(cpstr PDBString) -> ERROR_CODE` - Adds a PDB-keyworded string to the existing structure.
    *   `Root::ReadCoorFile(cpstr LFName) -> ERROR_CODE` - Reads coordinate data from a file in a coordinate format.
    *   `Root::RegisterUDInteger(UDR_TYPE udr_type, cpstr UDDataID) -> int` - Registers a user-defined integer data record.
    *   `Root::SetFlag(word Flag) -> void` - Sets a flag in the structure.
    *   `Root::SetSpaceGroup(cpstr spGroup) -> int` - Sets the space group for the current structure.
    *   `Root::WriteCIFASCII(cpstr CIFFileName) -> ERROR_CODE` - Writes CIF data to an ASCII file.
    *   `Root::WritePDBASCII(cpstr PDBFileName) -> ERROR_CODE` - Writes PDB data to an ASCII file.
    *   `Root::DeleteAllModels() -> int` - Deletes all models from the structure.
    *   `Root::FinishStructEdit() -> int` - Completes structural editing operations.

*   `Remark` (header.h) - Represents a remark within a PDB file, providing methods for dumping, converting, and accessing its associated data.
    *   `Remark::Remark()` - Constructs a new `Remark` object with default values.
    *   `Remark::Remark(cpstr S)` - Constructs a `Remark` object from a C-style string.
    *   `Remark::Remark(io::RPStream Object)` - Constructs a `Remark` object from an `io::RPStream` object.
    *   `Remark::~Remark()` - Destroys the `Remark` object, releasing any allocated memory.
    *   `Remark::PDBASCIIDump(pstr S, int N)` - Writes the `Remark` data to a given string `S` in PDB ASCII format.
    *   `Remark::MakeCIF(mmcif::PData CIF, int N)` - Generates a mmCIF data structure `CIF` from the current `Remark` data.
    *   `Remark::ConvertPDBASCII(cpstr S)` - Converts a PDB ASCII string `S` to a `Remark` object.
    *   `Remark::GetCIF(mmcif::PData CIF, int & n)` - Retrieves the mmCIF data object `CIF` from the `Remark` object.
    *   `Remark::GetClassID() -> CLASS_ID` - Returns the class ID for the `Remark` object.
    *   `Remark::Copy(PContainerClass RemarkClass)` - Copies the contents of another `Remark` object into this `Remark` object.
    *   `Remark::write(io::RFile f)` - Writes the `Remark` data to an `io::RFile` object.
    *   `Remark::read(io::RFile f)` - Reads the `Remark` data from an `io::RFile` object.
    *   `Remark::InitRemark()` - Initializes the `Remark` object.
    *   `Remark::FreeMemory()` - Frees any dynamically allocated memory used by the `Remark` object.

*   `SelManager` (header.h) - Manages selection masks and atom selection within a structural biology file.
    *   `SelManager::NewSelection() -> int` - Creates a new selection mask and returns its handle.
    *   `SelManager::GetSelType(int selHnd)` - Returns the selection type associated with a given selection handle.
    *   `SelManager::DeleteSelection(int selHnd)` - Deletes the specified selection mask and removes the corresponding selection attributes from all atoms selected with this mask.
    *   `SelManager::DeleteAllSelections()` - Deletes all selection masks and unselects all atoms in the file.
    *   `SelManager::SelectAtoms(int selHnd, int iSer1, int iSer2, SELECTION_KEY selKey)` - Selects atoms in the serial number range of `iSer1` to `iSer2` by adding them to the set of atoms marked by the given mask.
    *   `SelManager::SelectAtoms(int selHnd, ivector asn, int nsn, SELECTION_KEY selKey)` - Selects atoms based on the provided sequence number vector `asn` of length `nsn`.
    *   `SelManager::UnselectAtoms(int selHnd, int iSer1, int iSer2)` - Clears the specified mask for atoms in the serial number range of `iSer1` to `iSer2`.
    *   `SelManager::SelectAtom(int selHnd, PAtom A, SELECTION_KEY selKey, bool makeIndex)` - Selects a single atom according to the value of selection key.
    *   `SelManager::SelectResidue(int selHnd, PResidue Res, SELECTION_TYPE sType, SELECTION_KEY sKey, bool makeIndex)` - Selects a single residue, chain or model, or all their hierarchical descendants depending on the value of `sType`.
    *   `SelManager::SelectChain(int selHnd, PChain chain, SELECTION_TYPE sType, SELECTION_KEY sKey, bool makeIndex)` - Selects a single chain, or all their hierarchical descendants depending on the value of `sType`.
    *   `SelManager::SelectModel(int selHnd, PModel mdl, SELECTION_TYPE sType, SELECTION_KEY sKey, bool makeIndex)` - Selects a single model, or all their hierarchical descendants depending on the value of `sType`.
    *   `SelManager::MakeSelIndex(int selHnd)` - Calculates the selection index for the selection handled by `selHnd`.
    *   `SelManager::GetSelIndex(PPAtom & SelAtom, int & nSelAtoms)` - Returns an array of atom handles of atoms selected by a particular selection mask.
    *   `SelManager::DeleteAllSelections()` - Deletes all selection masks and unselects all atoms in the file.

*   `Sheet` (header.h) - Represents a sheet in a protein structure.
    *   `Sheet::Sheet()` - Constructs a new `Sheet` object.
    *   `Sheet::Sheet(io::RPStream Object)` - Constructs a new `Sheet` object from an `io::RPStream` object.
    *   `Sheet::~Sheet()` - Destroys the `Sheet` object, releasing any allocated memory.
    *   `Sheet::PDBASCIIDump(io::RFile f)` - Writes the sheet data to an RFile in PDB ASCII format.
    *   `Sheet::MakeCIF(mmcif::PData CIF)` - Creates a mmCIF data structure `CIF` from the current `Sheet` data.
    *   `Sheet::ConvertPDBASCII(cpstr S)` - Converts a PDB ASCII string `S` to a `Sheet` object.
    *   `Sheet::GetCIF(mmcif::PData CIF, int & n)` - Retrieves the CIF data for the sheet.
    *   `Sheet::GetClassID() -> CLASS_ID` - Returns the class ID for the `Sheet` object.
    *   `Sheet::Copy(PSheet Sheet)` - Copies the data from another `Sheet` object.
    *   `Sheet::write(io::RFile f)` - Writes the sheet data to an `io::RFile` object.
    *   `Sheet::read(io::RFile f)` - Reads the sheet data from an `io::RFile` object.
    *   `Sheet::InitSheet()` - Initializes the `Sheet` object.
    *   `Sheet::FreeMemory()` - Frees any dynamically allocated memory used by the `Sheet` object.
    *   `Sheet::OrderSheet()` - Orders the strands within the `Sheet` object.
    *   `Sheet::GetStrand(int strand_no)` - Returns the strand at the given index.

*   `Sheets` (header.h) - This class represents a collection of protein sheets within a PDB or CIF file.
    *   `Sheets::Sheets()` - Constructs a new `Sheets` object.
    *   `Sheets::Sheets(io::RPStream Object)` - Constructs a new `Sheets` object from an `io::RPStream` object.
    *   `Sheets::~Sheets()` - Destroys the `Sheets` object, releasing any associated resources.
    *   `Sheets::FreeMemory()` - Frees any memory allocated by the `Sheets` object.
    *   `Sheets::PDBASCIIDump(io::RFile f)` - Writes the sheet data to an RFile in PDB ASCII format.
    *   `Sheets::MakeCIF(mmcif::PData CIF)` - Creates a mmCIF data structure `CIF` from the current `Sheets` data.
    *   `Sheets::ConvertPDBASCII(cpstr S)` - Converts a PDB ASCII string `S` to a `Sheets` object.
    *   `Sheets::GetCIF(mmcif::PData CIF, cpstr sheet_id)` - Retrieves the CIF data for the sheets.
    *   `Sheets::GetClassID() -> CLASS_ID` - Returns the class ID for the `Sheets` object.
    *   `Sheets::Copy(PSheets Sheets)` - Copies the contents of another `Sheets` object.
    *   `Sheets::write(io::RFile f)` - Writes the sheet data to an `io::RFile` object.
    *   `Sheets::read(io::RFile f)` - Reads the sheet data from an `io::RFile` object.
    *   `Sheets::InitSheets()` - Initializes the `Sheets` object.
    *   `Sheets::FreeMemory()` - Frees any dynamically allocated memory used by the `Sheets` object.

*   `SymOps` (header.h) - This class provides functions for manipulating symmetry operations within an MMDB2 structure.
    *   `SymOps::SymOps()` - Constructs a new `SymOps` object.
    *   `SymOps::SymOps(io::RPStream Object)` - Constructs a new `SymOps` object from an `io::RPStream` object.
    *   `SymOps::~SymOps()` - Destroys the `SymOps` object, releasing any allocated memory.
    *   `SymOps::FreeMemory()` - Frees the memory allocated by the `SymOps` object.
    *   `SymOps::SetGroupSymopLib(cpstr SpaceGroup, cpstr symop_lib=NULL) -> int` - Sets the SpaceGroup name from a symop.lib file.
    *   `SymOps::SetGroup(cpstr SpaceGroup, cpstr syminfo_lib=NULL) -> int` - Sets the SpaceGroup name from a syminfo.lib file.
    *   `SymOps::Reset()` - Removes all

---

## Per-Class Method Documentation

## `Alignment`
This class represents a sequence alignment, facilitating the alignment of two sequences based on scoring parameters and potentially affine gap penalties.

#### `Alignment::Alignment()`
Constructs an empty alignment object.
Params: None
Returns: None
Side effects: Initializes member variables.

#### `Alignment::Alignment(io::RPStream Object)`
Constructs an alignment object from an existing `io::RPStream` object.
Params: `Object` — An `io::RPStream` object.
Returns: None
Side effects: Initializes member variables using the data from the `io::RPStream` object.

#### `Alignment::~Alignment()`
Destroys the alignment object, releasing allocated memory.
Params: None
Returns: None
Side effects: Releases all dynamically allocated memory.

#### `Alignment::SetAffineModel(realtype WGap, realtype WSpace)`
Sets the affine gap penalty parameters for gap opening (`WGap`) and extension (`WSpace`).
Params: `WGap` — The penalty for opening a gap. `WSpace` — The penalty for extending a gap.
Returns: None
Side effects: Updates internal gap penalty parameters.

#### `Alignment::SetScores(realtype SEqual, realtype SNEqual)`
Sets the scoring parameters for match (`SEqual`) and mismatch (`SNEqual`).
Params: `SEqual` — The score for a matching residue. `SNEqual` — The score for a non-matching residue.
Returns: None
Side effects: Updates internal scoring parameters.

#### `Alignment::Align(cpstr S, cpstr T, ALIGN_METHOD Method=ALIGN_GLOBAL)`
Aligns the sequences `S` and `T` based on specified scoring parameters and alignment method.
Params: `S` — The first sequence. `T` — The second sequence. `Method` — The alignment method to use (default: `ALIGN_GLOBAL`).
Returns: None
Side effects: Performs the sequence alignment, populating the `AlgnS` and `AlgnT` member variables.

#### `Alignment::GetAlignedS() -> pstr`
Returns a pointer to the aligned sequence `S`.
Params: None
Returns: A pointer to the aligned sequence `AlgnS`.

#### `Alignment::GetAlignedT() -> pstr`
Returns a pointer to the aligned sequence `T`.
Params: None
Returns: A pointer to the aligned sequence `AlgnT`.

#### `Alignment::GetScore() -> realtype`
Returns the achieved alignment score.
Params: None
Returns: The total alignment score.

#### `Alignment::GetSpace() -> char`
Returns the space character used in the alignment.
Params: None
Returns: The space character `Space`.

#### `Alignment::GetSimilarity() -> realtype`
Calculates a similarity score based on the alignment, taking into account sequence lengths and scoring.
Params: None
Returns: The calculated similarity score.

#### `Alignment::GetSeqId() -> realtype`
Returns the primitive sequence ID.
Params: None
Returns: The ID of the sequence.

#### `Alignment::GetNAlign() -> int`
Returns the number of aligned residues.
Params: None
Returns: The number of aligned residues.

#### `Alignment::OutputResults(io::RFile f, cpstr S, cpstr T)`
Outputs the alignment results to a file.
Params: `f` — An `io::RFile` object to write to. `S` — The first sequence. `T` — The second sequence.
Returns: None
Side effects: Writes the aligned sequences and their score to the specified file.

#### `Alignment::read(io::RFile f)`
Reads alignment data from a file.
Params: `f` — An `io::RFile` object to read from.
Returns: None
Side effects: Populates the alignment object with data read from the file.

#### `Alignment::write(io::RFile f)`
Writes alignment data to a file.
Params: `f` — An `io::RFile` object to write to.
Returns: None
Side effects: Writes the alignment data to the specified file.

#### `Alignment::InitAlignment()`
Initializes the alignment data structures.
Params: None
Returns: None
Side effects: Sets up the alignment data, including sequence lengths, scoring parameters, and gap penalty matrices.

#### `Alignment::FreeMemory()`
Frees any dynamically allocated memory used by the alignment object.
Params: None
Returns: None
Side effects: Releases memory allocated to the alignment object.

#### `Alignment::Score(char A, char B) -> realtype`
Calculates the score for a match or mismatch between two characters.
Params: `A` — The first character. `B` — The second character.
Returns: The score for the match/mismatch.

#### `Alignment::BuildGATable(cpstr S, cpstr T, bool FreeSEnd, bool FreeTEnd)`
Builds a gap opening/extension table for use in the alignment algorithm.
Params: `S` — The first sequence. `T` — The second sequence. `FreeSEnd` — Flag to indicate free end of sequence S. `FreeTEnd` — Flag to indicate free end of sequence T.
Returns: None
Side effects: Populates internal gap penalty matrices.

#### `Alignment::BuildLATable(cpstr S, cpstr T)`
Builds a local alignment table for use in the alignment algorithm.
Params: `S` — The first sequence. `T` — The second sequence.
Returns: None
Side effects: Populates internal scoring and traceback matrices.

#### `Alignment::Backtrace(cpstr S, cpstr T, int J, int I, bool StopAtZero)`
Performs the traceback to reconstruct the aligned sequences.
Params: `S` — The first sequence. `T` — The second sequence. `J` — The column index. `I` — The row index. `StopAtZero` — Flag to stop at zero score.
Returns: None
Side effects: Constructs the aligned sequences based on the traceback path.

#### `Alignment::AdjustEnds(cpstr S, cpstr T, int J, int I)`
Adjusts the ends of the aligned sequences after the traceback.
Params: `S` — The first sequence. `T` — The second sequence. `J` — The column index. `I` — The row index.
Returns: None
Side effects: Fine-tunes the alignment endpoints.

#### `Alignment::PrintVT(cpstr S, cpstr T)`
Prints the VT table to the console for debugging.
Params: `S` — The first sequence. `T` — The second sequence.
Returns: None
Side effects: Prints the VT table, useful for debugging.


---

## `Atom`
Represents an atom in a structural biology file, storing its coordinates, bond information, and related data.

#### `Atom::AddBond(PAtom bond_atom, int bond_order, int nAdd_bonds=1)`
Adds a bond to the atom, connecting it to another atom with a specified bond order.  
Params: `bond_atom` — Pointer to the bonded atom. `bond_order` — The type of bond (e.g., single, double, triple). `nAdd_bonds` — Number of bonds to add, defaults to 1.
Returns: `<=0` - Error (e.g., atom already bonded). `>0` - Success, returns the current number of bonds.
Side effects: Updates the internal bond list of the atom.
Errors: Returns an error code if the target atom is already bonded to this atom.

#### `Atom::Copy()`
Creates a deep copy of the `Atom` object, duplicating all its attributes.
Params: None.
Returns: A new `Atom` object.
Side effects: Creates a new `Atom` object with the same data as the original.

#### `Atom::GetAtomID() -> pstr`
Generates a unique string identifier for the atom in the format `/m/c/r(rn).i/n[e]:a`.
Params: None.
Returns: A string containing the atom ID.
Side effects: None.

#### `Atom::GetAtomName() -> pstr`
Returns the atom's name as stored in the file.
Params: None.
Returns: A pointer to the atom's name.
Side effects: None.

#### `Atom::GetBonds() -> void`
Returns a pointer to the internal structure representing the atom's bonds.
Params: None.
Returns: None.
Side effects: None.

#### `Atom::GetChain() -> pstr`
Returns the name of the chain to which the atom belongs.
Params: None.
Returns: A pointer to the chain name.
Side effects: None.

#### `Atom::GetChainID() -> pstr`
Returns the chain ID (index) within the chain.
Params: None.
Returns: A pointer to the chain ID.
Side effects: None.

#### `Atom::GetElementName() -> pstr`
Returns the element symbol associated with the atom (e.g., "C", "O", "N").
Params: None.
Returns: A pointer to the element name.
Side effects: None.

#### `Atom::GetIndex() -> int`
Returns the atom's index within the residue.
Params: None.
Returns: An integer representing the atom's index.
Side effects: None.

#### `Atom::GetInsCode() -> pstr`
Returns the insertion code (residue sequence number) for the atom.
Params: None.
Returns: A pointer to the insertion code string.
Side effects: None.

#### `Atom::GetModel() -> pstr`
Returns the model number to which the atom belongs.
Params: None.
Returns: A pointer to the model name.
Side effects: None.

#### `Atom::GetModelNum() -> int`
Returns the model number to which the atom belongs.
Params: None.
Returns: An integer representing the model number.
Side effects: None.

#### `Atom::GetNBonds() -> int`
Returns the number of bonds associated with the atom.
Params: None.
Returns: An integer representing the number of bonds.
Side effects: None.

#### `Atom::GetResName() -> pstr`
Returns the name of the residue to which the atom belongs.
Params: None.
Returns: A pointer to the residue name.
Side effects: None.

#### `Atom::GetResidue() -> PResidue`
Returns a pointer to the residue object to which the atom belongs.
Params: None.
Returns: A pointer to the `Residue` object.
Side effects: None.

#### `Atom::GetSeqNum() -> int`
Returns the sequence number of the residue to which the atom belongs.
Params: None.
Returns: An integer representing the sequence number.
Side effects: None.

#### `Atom::GetUDData() -> int`
Returns the user-defined data mask indicating which fields of the atom are set.
Params: None.
Returns: An integer representing the data mask.
Side effects: None.

#### `Atom::MakeTer() -> void`
Marks the atom as a chain terminator (ter).
Params: None.
Returns: None.
Side effects: Sets the `Ter` flag to `true`.

#### `Atom::PutUDData(int UDDhandle, int iudd)`
Writes user defined data to a specific handle.
Params: `UDDhandle` - Handle to the user defined data. `iudd` - Integer data.
Returns: None.
Side effects: Updates the user defined data.

#### `Atom::PutUDData(int UDDhandle, realtype rudd)`
Writes a real-valued user defined data to a specific handle.
Params: `UDDhandle` - Handle to the user defined data. `rudd` - Real value.
Returns: None.
Side effects: Updates the user defined data.

#### `Atom::PutUDData(int UDDhandle, cpstr sudd)`
Writes a string user defined data to a specific handle.
Params: `UDDhandle` - Handle to the user defined data. `suddt` - String.
Returns: None.
Side effects: Updates the user defined data.

#### `Atom::PutUDData(int UDDhandle, int& iudd)`
Reads user defined data from a specific handle.
Params: `UDDhandle` - Handle to the user defined data. `iudd` - Integer data.
Returns: None.
Side effects: Updates the user defined data.

#### `Atom::PutUDData(int UDDhandle, realtype& rudd)`
Reads a real-valued user defined data from a specific handle.
Params: `UDDhandle` - Handle to the user defined data. `rudd` - Real value.
Returns: None.
Side effects: Updates the user defined data.

#### `Atom::PutUDData(int UDDhandle, pstr & sudd)`
Reads a string user defined data from a specific handle.
Params: `UDDhandle` - Handle to the user defined data. `suddt` - String.
Returns: None.
Side effects: Updates the user defined data.

#### `Atom::SetAtomName(int ix, int sN, const AtomName aName, const AltLoc aLoc, const SegID sID, const Element eName)`
Sets the atom's name (serial number, name, altLoc, segID, element).
Params: `ix` - Atom index. `sN` - Serial number. `aName` - Atom name. `aLoc` - Alternative location. `sID` - Segment ID. `eName` - Element name.
Returns: None.
Side effects: Updates the atom's name and associated data.

#### `Atom::SetAtomName(const AtomName atomName)`
Sets the atom's name to the provided name.
Params: `atomName` - The name to set.
Returns: None.
Side effects: Updates the atom's name.

#### `Atom::SetElementName(const Element elName)`
Sets the atom's element to the provided element.
Params: `elName` - The element to set.
Returns: None.
Side effects: Updates the atom's element.

#### `Atom::SetResidue(PResidue res)`
Sets the atom's residue to the provided residue.
Params: `res` - Pointer to the residue object.
Returns: None.
Side effects: Updates the atom's residue.

#### `Atom::SetCoordinates(realtype xx, realtype yy, realtype zz, realtype occ, realtype tFac)`
Sets the atom's coordinates and occupancy/temperature factor.
Params: `xx`, `yy`, `zz` - X, Y, and Z coordinates. `occ` - Occupancy. `tFac` - Temperature factor.
Returns: None.
Side effects: Updates the atom's coordinates, occupancy, and temperature factor.

#### `Atom::Transform(const mat33 & tm, vect3 & v)`
Transforms the atom's coordinates using the given 3x3 transformation matrix.
Params: `tm` - The 3x3 transformation matrix. `v` - The vector to transform.
Returns: None.
Side effects: Updates the atom's coordinates based on the transformation.

#### `Atom::Transform(const mat44 & tm)`
Transforms the atom's coordinates using the given 4x4 transformation matrix.
Params: `tm` - The 4x4 transformation matrix.
Returns: None.
Side effects: Updates the atom's coordinates based on the transformation.

#### `Atom::TransformCopy(const mat44 & tm, realtype & xx, realtype & yy, realtype & zz)`
Transforms the atom's coordinates using the given 4x4 transformation matrix and stores result in `xx`, `yy`, `zz`.
Params: `tm` - The 4x4 transformation matrix. `xx`, `yy`, `zz` - Coordinate X, Y, and Z values.
Returns: None.
Side effects: Updates the atom's coordinates based on the transformation.

#### `Atom::TransformCopy(const mat44 & tm, vect3 & xyz)`
Transforms the atom's coordinates using the given 4x4 transformation matrix and stores result in `xyz`.
Params: `tm` - The 4x4 transformation matrix. `xyz` - Vector to store result.
Returns: None.
Side effects: Updates the atom's coordinates based on the transformation.

#### `Atom::isInSelection(int selHnd)`
Checks if the atom is currently selected in the selection set.
Params: `selHnd` - The handle of the selection set.
Returns: `true` if the atom is selected, `false` otherwise.
Side effects: None.

#### `Atom::isMetal() -> bool`
Checks if the atom is a metal atom.
Params: None.
Returns: `true` if the atom is a metal, `false` otherwise.
Side effects: None.

#### `Atom::isTer() -> bool`
Checks if the atom is a chain terminator.
Params: None.
Returns: `true` if the atom is a chain terminator, `false` otherwise.
Side effects: None.

#### `Atom::operator=`
Assigns another `Atom` object to this object.
Params: None.
Returns: `*this` (self-assignment).
Side effects: Copies the data from the other `Atom` object to this object.


---

## `Author`
This class represents an author associated with a PDB file.

#### `Author::Author()`
Constructs a new `Author` object.
Params: None
Returns: None
Side effects: Initializes the `Author` object.

#### `Author::Author(cpstr S)`
Constructs a new `Author` object from a C-style string.
Params: `cpstr S` — The PDB author string.
Returns: None
Side effects: Initializes the `Author` object with the provided string.

#### `Author::Author(io::RPStream Object)`
Constructs a new `Author` object from an `io::RPStream`.
Params: `io::RPStream Object` — The input stream.
Returns: None
Side effects: Initializes the `Author` object from the stream.

#### `Author::~Author()`
Destroys the `Author` object, releasing any resources.
Params: None
Returns: None
Side effects: Cleans up the `Author` object's internal state.

#### `Author::PDBASCIIDump(pstr S, int N)`
Writes the `Author`'s PDB ASCII string to a specified string.
Params: `pstr S` — The destination string buffer, `int N` — The maximum length of the buffer.
Returns: None
Side effects: Writes the `Author`'s PDB ASCII string to the string `S`.

#### `Author::ConvertPDBASCII(cpstr S)`
Converts a PDB ASCII string to an `Author` object.
Params: `cpstr S` — The PDB ASCII string.
Returns: `ERROR_CODE` — The error code if the conversion fails, otherwise `NoError`.
Side effects: Converts the input string to an `Author` object.

#### `Author::GetClassID()`
Returns the class identifier for an `Author` object.
Params: None
Returns: `CLASS_ID` — The `ClassID_Author` constant.
Side effects: None

#### `Author::write(io::RFile f)`
Writes the `Author` object's data to an `io::RFile`.
Params: `io::RFile f` — The output file stream.
Returns: None
Side effects: Writes the `Author` object's data to the file stream.

#### `Author::read(io::RFile f)`
Reads the `Author` object's data from an `io::RFile`.
Params: `io::RFile f` — The input file stream.
Returns: None
Side effects: Reads the `Author` object's data from the file stream.


---

## `BondManager`
This class manages bonds within an MMDB structure.

#### `BondManager::MakeBonds(bool calc_only) -> void`
Calculates and creates bonds based on the current structure and atom selections, or only calculates the bonds if `calc_only` is true.
Params: `calc_only` — A boolean indicating whether to calculate only or create bonds.
Returns: `void`.
Side effects: Creates or updates bonds in the MMDB file.

#### `BondManager::RemoveBonds() -> void`
Removes all currently existing bonds from the MMDB structure.
Params: None.
Returns: `void`.
Side effects: Removes bonds from the MMDB file.


---

## `Category`
Represents a category of mmCIF data, managing tags and their associated data.

#### `Category::GetCategoryName() -> pstr`
Returns the category name if it was set, otherwise returns `NULL`.

#### `Category::GetNofTags() -> int`
Returns the total number of tags in the category.

#### `Category::GetTag(int tagNo) -> pstr`
Returns the tag at the specified serial number (0..nTags-1), or `NULL` if tagNo is out of range.

#### `Category::AddTag(cpstr ttag) -> int`
Adds a tag to the category and returns the tag's serial number (0..nTags-1) or a negative value if the tag was not added.

#### `Category::GetTagNo(cpstr ttag) -> int`
Returns the serial number of a tag in the category, or a negative value if the tag was not found.

#### `Category::PrintTags()`
Prints the list of tags to standard output.

#### `Category::CheckTags(cpstr * tagList)`
Checks if all tags from the list are found in the category, returns true or false.

#### `Category::Copy(PCategory Category)`
Performs a deep copy of the category.

#### `Category::WriteMMCIF(io::RFile)`
Writes category's content into mmCIF file. Default implementation does nothing.

#### `Category::ReadMMCIF(io::RFile)`
Reads category's content from mmCIF file. Default implementation does nothing.

#### `Category::Optimize()`
Optimizes the data structures within the category for faster access.

#### `Category::Sort()`
Sorts the category's data for quicker access, preserving tag order as added.

#### `Category::SetCategoryName(cpstr N)`
Sets the category name.

#### `Category::GetCategoryID() -> MMCIF_ITEM`
Returns the category type (MMCIF_Category, MMCIF_Struct, or MMCIF_Loop).

#### `Category::FreeMemory()`
Frees all memory allocated by the category.

#### `Category::ExpandTags(int nTagsNew)`
Expands the tag vector to accommodate new tags.

#### `Category::PutCategoryName(cpstr newName)`
Internal method to put category name.


---

## `Chain`
Represents a chain of atoms in a structural biology model.

#### `Chain::AddResidue(PResidue res) -> int`
Adds a residue to the chain.
Params: `res` — A pointer to the residue to be added.
Returns: The index of the added residue in the chain's residue table.
Side effects: Adds the residue to the chain's residue table and updates internal counters.
Errors: None.

#### `Chain::Copy() -> void`
Creates a complete copy of this chain.
Params: None.
Returns: None.
Side effects: Allocates new memory for the copy and copies all data from the original chain.
Errors: None.

#### `Chain::DeleteResidue(int resNo) -> int`
Deletes a residue from the chain.
Params: `resNo` — The index of the residue to be deleted.
Returns: The index of the deleted residue.
Side effects: Removes the specified residue from the chain's residue table and adjusts internal counters.
Errors: None.

#### `Chain::GetChainID() -> cpstr`
Returns the chain ID as a string.
Params: None.
Returns: A C-string representing the chain ID.
Side effects: None.
Errors: None.

#### `Chain::GetCoordHierarchy() -> PManager`
Returns a pointer to the coordinate hierarchy for this chain.
Params: None.
Returns: A pointer to a `PManager` object representing the coordinate hierarchy.
Side effects: None.
Errors: None.

#### `Chain::ConvertDBREF(cpstr PDBString) -> ERROR_CODE`
Converts a database reference string to the internal format.
Params: `PDBString` — A C-style string containing the database reference.
Returns: An `ERROR_CODE` indicating success or failure.
Side effects: Updates the `DBRef` container with the converted database reference.
Errors: Returns an error code if the string is invalid or cannot be parsed.

#### `Chain::ConvertSEQADV(cpstr PDBString) -> ERROR_CODE`
Converts a sequence advancement record string to the internal format.
Params: `PDBString` — A C-style string containing the sequence advancement record.
Returns: An `ERROR_CODE` indicating success or failure.
Side effects: Updates the `seqAdv` container with the converted sequence advancement record.
Errors: Returns an error code if the string is invalid or cannot be parsed.

#### `Chain::ConvertSEQRES(cpstr PDBString) -> ERROR_CODE`
Converts a sequence residue record string to the internal format.
Params: `PDBString` — A C-style string containing the sequence residue record.
Returns: An `ERROR_CODE` indicating success or failure.
Side effects: Updates the `seqRes` container with the converted sequence residue record.
Errors: Returns an error code if the string is invalid or cannot be parsed.

#### `Chain::ConvertMODRES(cpstr PDBString) -> ERROR_CODE`
Converts a modification record string to the internal format.
Params: `PDBString` — A C-style string containing the modification record.
Returns: An `ERROR_CODE` indicating success or failure.
Side effects: Updates the `modRes` container with the converted modification record.
Errors: Returns an error code if the string is invalid or cannot be parsed.

#### `Chain::ConvertHET(cpstr PDBString) -> ERROR_CODE`
Converts a heterocyclic record string to the internal format.
Params: `PDBString` — A C-style string containing the heterocyclic record.
Returns: An `ERROR_CODE` indicating success or failure.
Side effects: Updates the `Het` container with the converted heterocyclic record.
Errors: Returns an error code if the string is invalid or cannot be parsed.

#### `Chain::PDBASCIIDump(io::RFile f) -> void`
Dumps the chain's data into a PDB ASCII file.
Params: `f` — A file stream object to write the PDB ASCII data to.
Returns: None.
Side effects: Writes the chain's data to the specified file in PDB ASCII format.
Errors: None.

#### `Chain::PDBASCIIAtomDump(io::RFile f) -> void`
Dumps the chain's atoms into a PDB ASCII file.
Params: `f` — A file stream object to write the PDB ASCII data to.
Returns: None.
Side effects: Writes the chain's atoms to the specified file in PDB ASCII format.
Errors: None.

#### `Chain::MakeAtomCIF(mmcif::PData CIF) -> void`
Creates a MMCIF data structure from the chain's data.
Params: `CIF` — A pointer to the MMCIF data structure to be populated.
Returns: None.
Side effects: Populates the MMCIF data structure with the chain's data.
Errors: None.

#### `Chain::GetNumberOfResidues() -> int`
Returns the number of residues in the chain.
Params: None.
Returns: The number of residues in the chain.
Side effects: None.
Errors: None.

#### `Chain::GetResidue(int seqNum, const InsCode insCode) -> PResidue`
Returns a pointer to a residue in the chain.
Params: `seqNum` — The sequence number of the residue. `insCode` — The insertion code.
Returns: A pointer to a `PResidue` object representing the specified residue.
Side effects: None.
Errors: Returns `NULL` if the residue is not found in the chain.

#### `Chain::GetResidue(int resNo) -> PResidue`
Returns a pointer to a residue in the chain by residue number.
Params: `resNo` — The index of the residue in the chain's residue table.
Returns: A pointer to a `PResidue` object representing the specified residue.
Side effects: None.
Errors: Returns `NULL` if the residue is not found in the chain.

#### `Chain::GetResidueNo(int seqNum, const InsCode insCode) -> int`
Returns the residue number in the chain's residue table.
Params: `seqNum` — The sequence number of the residue. `insCode` — The insertion code.
Returns: The residue number in the chain's residue table.
Side effects: None.
Errors: Returns -1 if the residue is not found in the chain.

#### `Chain::GetResidueTable(PPResidue & resTable, int & NumberOfResidues) -> void`
Returns the chain's residue table.
Params: `resTable` — A pointer to an array of `PResidue` objects to be populated. `NumberOfResidues` — A reference to an integer to store the number of residues in the table.
Returns: None.
Side effects: Populates the `resTable` array with pointers to the chain's residues and sets `NumberOfResidues` to the number of residues.
Errors: None.

#### `Chain::GetResidueCreate(const ResName resName, int seqNum, const InsCode insCode, bool Enforce) -> PResidue`
Creates a residue and adds it to the chain.
Params: `resName` — The name of the residue. `seqNum` — The sequence number of the residue. `insCode` — The insertion code. `Enforce` — A flag indicating whether to enforce existing residue names.
Returns: A pointer to the newly created `PResidue` object.
Side effects: Adds the new residue to the chain's residue table.
Errors: Returns `NULL` if the residue is already present in the chain or if `Enforce` is set to `false`.

#### `Chain::GetCoordHierarchy() -> PManager`
Returns a pointer to the coordinate hierarchy for this chain.
Params: None.
Returns: A pointer to a `PManager` object representing the coordinate hierarchy.
Side effects: None.
Errors: None.

#### `Chain::SetChain(const ChainID chID) -> void`
Sets the chain ID for this chain.
Params: `chID` — The chain ID to set.
Returns: None.
Side effects: Updates the chain's ID and adjusts internal references accordingly.
Errors: None.

#### `Chain::SortResidues() -> void`
Sorts the residues in the chain based on their sequence numbers.
Params: None.
Returns: None.
Side effects: Sorts the chain's residue table, updating the residue pointers.
Errors: None.

#### `Chain::TrimResidueTable() -> void`
Removes residues from the end of the chain's residue table.
Params: None.
Returns: None.
Side effects: Removes all residues from the chain's residue table.
Errors: None.

#### `Chain::isAminoacidChain() -> bool`
Checks if the chain contains only amino acid residues.
Params: None.
Returns: `true` if the chain contains only amino acid residues, `false` otherwise.
Side effects: None.
Errors: None.

#### `Chain::isInSelection(int selHnd) -> bool`
Checks if the chain is present in a selection.
Params: `selHnd` — The handle of the selection.
Returns: `true` if the chain is present in the selection, `false` otherwise.
Side effects: None.
Errors: None.

#### `Chain::isNucleotideChain() -> bool`
Checks if the chain contains only nucleotide residues.
Params: None.
Returns: `true` if the chain contains only nucleotide residues, `false` otherwise.
Side effects: None.
Errors: None.

#### `Chain::isSolventChain() -> bool`
Checks if the chain contains only solvent residues.
Params: None.
Returns: `true` if the chain contains only solvent residues, `false` otherwise.
Side effects: None.
Errors: None.

#### `Chain::AddResidue(PResidue res) -> int`
Adds a residue to the chain.
Params: `res` — A pointer to the residue to be added.
Returns: The index of the added residue in the chain's residue table.
Side effects: Adds the residue to the chain's residue table and updates internal counters.
Errors: None.

#### `Chain::InsResidue(PResidue res, int pos) -> int`
Inserts a residue at a specified position in the chain.
Params: `res` — A pointer to the residue to be inserted. `pos` — The index of the position to insert the residue.
Returns: The index of the inserted residue in the chain's residue table.
Side effects: Inserts the residue at the specified position in the chain's residue table and adjusts internal counters.
Errors: None.

#### `Chain::InsResidue(PResidue res, int seqNum, const InsCode insCode) -> int`
Inserts a residue at a specified sequence number in the chain.
Params: `res` — A pointer to the residue to be inserted. `seqNum` — The sequence number of the position to insert the residue. `insCode` — The insertion code.
Returns: The index of the inserted residue in the chain's residue table.
Side effects: Inserts the residue at the specified sequence number in the chain's residue table and adjusts internal counters.
Errors: None.

#### `Chain::DeleteResidue(int resNo) -> int`
Deletes a residue from the chain.
Params: `resNo` — The index of the residue to be deleted.
Returns: The index of the deleted residue.
Side effects: Removes the specified residue from the chain's residue table and adjusts internal counters.
Errors: None.

#### `Chain::DeleteResidue(int seqNum, const InsCode insCode) -> int`
Deletes a residue from the chain by sequence number and insertion code.
Params: `seqNum` — The sequence number of the residue to be deleted. `insCode` — The insertion code.
Returns: The index of the deleted residue.
Side effects: Removes the specified residue from the chain's residue table and adjusts internal counters.
Errors: None.

#### `Chain::DeleteAllResidues() -> int`
Deletes all residues from the chain.
Params: None.
Returns: The number of residues deleted.
Side effects: Removes all residues from the chain's residue table and resets internal counters.
Errors: None.

#### `Chain::DeleteSolvent() -> int`
Deletes all solvent residues from the chain.
Params: None.
Returns: The number of solvent residues deleted.
Side effects: Removes all solvent residues from the chain's residue table and resets internal counters.
Errors: None.

#### `Chain::TrimResidueTable() -> void`
Trims the residue table, removing any residues that are no longer present in the chain.
Params: None.
Returns: None.
Side effects: Removes any residues from the chain's residue table that are no longer present in the chain, adjusting internal counters.
Errors: None.

#### `Chain::GetNumberOfResidues() -> int`
Returns the number of residues in the chain.
Params: None.
Returns: The number of residues in the chain.
Side effects: None.
Errors: None.

#### `Chain::GetResidue(int resNo) -> PResidue`
Returns a pointer to a residue in the chain by residue number.
Params: `resNo` — The index of the residue in the chain's residue table.
Returns: A pointer to a `PResidue` object representing the specified residue.
Side effects: None.
Errors: Returns `NULL` if the residue is not found in the chain.

#### `Chain::GetResidue(int seqNum, const InsCode insCode) -> PResidue`
Returns a pointer to a residue in the chain by sequence number and insertion code.
Params: `seqNum` — The sequence number of the residue. `insCode` — The insertion code.
Returns: A pointer to a `PResidue` object representing the specified residue.
Side effects: None.
Errors: Returns `NULL` if the residue is not found in the chain.

#### `Chain::GetResidueNo(int seqNum, const InsCode insCode) -> int`
Returns the residue number in the chain's residue table.
Params: `seqNum` — The sequence number of the residue. `insCode` — The insertion code.
Returns: The residue number in the chain's residue table.
Side effects: None.
Errors: Returns -1 if the residue is not found in the chain.

#### `Chain::GetResidueTable(PPResidue & resTable, int & NumberOfResidues) -> void`
Returns the chain's residue table.
Params: `resTable` — A pointer to an array of `PResidue` objects to be populated. `NumberOfResidues` — A reference to an integer to store the number of residues in the table.
Returns: None.
Side effects: Populates the `resTable` array with pointers to the chain's residues and sets `NumberOfResidues` to the number of residues.
Errors: None.

#### `Chain::GetResidueCreate(const ResName resName, int seqNum, const InsCode insCode, bool Enforce) -> PResidue`
Creates a residue and adds it to the chain.
Params: `resName` — The name of the residue. `seqNum` — The sequence number of the residue. `insCode` — The insertion code. `Enforce` — A flag indicating whether to enforce existing residue names.
Returns: A pointer to the newly created `PResidue` object.
Side effects: Adds the new residue to the chain's residue table.
Errors: Returns `NULL` if the residue is already present in the chain or if `Enforce` is set to `false`.

#### `Chain::GetCoordHierarchy() -> PManager`
Returns a pointer to the coordinate hierarchy for this chain.
Params: None.
Returns: A pointer to a `PManager` object representing the coordinate hierarchy.
Side effects: None.
Errors: None.

#### `Chain::SetChain(const ChainID chID) -> void`
Sets the chain ID for this chain.
Params: `chID` — The chain ID to set.
Returns: None.
Side effects: Updates the chain's ID and adjusts internal references accordingly.
Errors: None.

#### `Chain::SortResidues() -> void`
Sorts the residues in the chain based on their sequence numbers.
Params: None.
Returns: None.
Side effects: Sorts the chain's residue table, updating the residue pointers.
Errors: None.

#### `Chain::TrimResidueTable() -> void`
Trims the residue table, removing any residues that are no longer present in the chain.
Params: None.
Returns: None.
Side effects: Removes any residues from the chain's residue table that are no longer present in the chain, adjusting internal counters.
Errors: None.

#### `Chain::isAminoacidChain() -> bool`
Checks if the chain contains only amino acid residues.
Params: None.
Returns: `true` if the chain contains only amino acid residues, `false` otherwise.
Side effects: None.
Errors: None.

#### `Chain::isInSelection(int selHnd) -> bool`
Checks if the chain is present in a selection.
Params: `selHnd` — The handle of the selection.
Returns: `true` if the chain is present in the selection, `false` otherwise.
Side effects: None.
Errors: None.

#### `Chain::isNucleotideChain() -> bool`
Checks if the chain contains only nucleotide residues.
Params: None.
Returns: `true` if the chain contains only nucleotide residues, `false` otherwise.
Side effects: None.
Errors: None.

#### `Chain::isSolventChain() -> bool`
Checks if the chain contains only solvent residues.
Params: None.
Returns: `true` if the chain contains only solvent residues, `false` otherwise.
Side effects: None.
Errors: None.

#### `Chain::AddResidue(PResidue res) -> int`
Adds a residue to the chain.
Params: `res` — A pointer to the residue to be added.
Returns: The index of the added residue in the chain's residue table.
Side effects: Adds the residue to the chain's residue table and updates internal counters.
Errors: None.

#### `Chain::InsResidue(PResidue res, int pos) -> int`
Inserts a residue at a specified position in the chain.
Params: `res` — A pointer to the residue to be inserted. `pos` — The index of the position to insert the residue.
Returns: The index of the inserted residue in the chain's residue table.
Side effects: Inserts the residue at the specified position in the chain's residue table and adjusts internal counters.
Errors: None.

#### `Chain::InsResidue(PResidue res, int seqNum, const InsCode insCode) -> int`
Inserts a residue at a specified sequence number in the chain.
Params: `res` — A pointer to the residue to be inserted. `seqNum` — The sequence number of the residue. `insCode` — The insertion code.
Returns: The index of the inserted residue in the chain's residue table.
Side effects: Inserts the residue at the specified sequence number in the chain's residue table and adjusts internal counters.
Errors: None.

#### `Chain::DeleteResidue(int resNo) -> int`
Deletes a residue from the chain.
Params: `resNo` — The index of the residue to be deleted.
Returns: The index of the deleted residue.
Side effects: Removes the specified residue from the chain's residue table and adjusts internal counters.
Errors: None.

#### `Chain::DeleteResidue(int seqNum, const InsCode insCode) -> int`
Deletes a residue from the chain by sequence number and insertion code.
Params: `seqNum` — The sequence number of the residue to be deleted. `insCode` — The insertion code.
Returns: The index of the deleted residue.
Side effects: Removes the specified residue from the chain's residue table and adjusts internal counters.
Errors: None.

#### `Chain::DeleteAllResidues() -> int`
Deletes all residues from the chain.
Params: None.
Returns: The number of residues deleted.
Side effects: Removes all residues from the chain's residue table and resets internal counters.
Errors: None.

#### `Chain::DeleteSolvent() -> int`
Deletes all solvent residues from the chain.
Params: None.
Returns: The number of solvent residues deleted.
Side effects: Removes all solvent residues from the chain's residue table and resets internal counters.
Errors: None.

#### `Chain::TrimResidueTable() -> void`
Trims the residue table, removing any residues that are no longer present in the chain.
Params: None.
Returns: None.
Side effects: Removes any residues from the chain's residue table that are no longer present in the chain, adjusting internal counters.
Errors: None.

#### `Chain::GetNumberOfResidues() -> int`
Returns the number of residues in the chain.
Params: None.
Returns: The number of residues in the chain.
Side effects: None.
Errors: None.

#### `Chain::GetResidue(int seqNum, const InsCode insCode) -> PResidue`
Returns a pointer to a residue in the chain by sequence number and insertion code.
Params: `seqNum` — The sequence number of the residue. `insCode` — The insertion code.
Returns: A pointer to a `PResidue` object representing the specified residue.
Side effects: None.
Errors: Returns `NULL` if the residue is not found in the chain.

#### `Chain::GetResidue(int resNo) -> PResidue`
Returns a pointer to a residue in the chain by residue number.
Params: `resNo` — The index of the residue in the chain's residue table.
Returns: A pointer to a `PResidue` object representing the specified residue.
Side effects: None.
Errors: Returns `NULL` if the residue is not found in the chain.

#### `Chain::GetResidueNo(int seqNum, const InsCode insCode) -> int`
Returns the residue number in the chain's residue table.
Params: `seqNum` — The sequence number of the residue. `insCode` — The insertion code.
Returns: The residue number in the chain's residue table.
Side effects: None.
Errors: Returns -1 if the residue is not found in the chain.

#### `Chain::GetResidueTable(PPResidue & resTable, int & NumberOfResidues) -> void`
Returns the chain's residue table.
Params: `resTable` — A pointer to an array of `PResidue` objects to be populated. `NumberOfResidues` — A reference to an integer to store the number of residues in the table.
Returns: None.
Side effects: Populates the `resTable` array with pointers to the chain's residues and sets `NumberOfResidues` to the number of residues.
Errors: None.

#### `Chain::GetResidueCreate(const ResName resName, int seqNum, const InsCode insCode, bool Enforce) -> PResidue`
Creates a residue and adds it to the chain.
Params: `resName` — The name of the residue. `seqNum` — The sequence number of the residue. `insCode` — The insertion code. `Enforce` — A flag indicating whether to enforce existing residue names.
Returns: A pointer to the newly created `PResidue` object.
Side effects: Adds the new residue to the chain's residue table.
Errors: Returns `NULL` if the residue is already present in the chain or if `Enforce` is set to `false`.

#### `Chain::GetCoordHierarchy() -> PManager`
Returns a pointer to the coordinate hierarchy for this chain.
Params: None.
Returns: A pointer to a `PManager` object representing the coordinate hierarchy.
Side effects: None.
Errors: None.

#### `Chain::SetChain(const ChainID chID) -> void`
Sets the chain ID for this chain.
Params: `chID` — The chain ID to set.
Returns: None.
Side effects: Updates the chain's ID and adjusts internal references accordingly.
Errors: None.

#### `Chain::SortResidues() -> void`
Sorts the residues in the chain based on their sequence numbers.
Params: None.
Returns: None.
Side effects: Sorts the chain's residue table, updating the residue pointers.
Errors: None.

#### `Chain::TrimResidueTable() -> void`
Trims the residue table, removing any residues that are no longer present in the chain.
Params: None.
Returns: None.
Side effects: Removes any residues from the chain's residue table that are no longer present in the chain, adjusting internal counters.
Errors: None.

#### `Chain::isAminoacidChain() -> bool`
Checks if the chain contains only amino acid residues.
Params: None.
Returns: `true` if the chain contains only amino acid residues, `false` otherwise.
Side effects: None.
Errors: None.

#### `Chain::isInSelection(int selHnd) -> bool`
Checks if the chain is present in a selection.
Params: `selHnd` — The handle of the selection.
Returns: `true` if the chain is present in the selection, `false` otherwise.
Side effects: None.
Errors: None.

#### `Chain::isNucleotideChain() -> bool`
Checks if the chain contains only nucleotide residues.
Params: None.
Returns: `true` if the chain contains only nucleotide residues, `false` otherwise.
Side effects: None.
Errors: None.

#### `Chain::isSolventChain() -> bool`
Checks if the chain contains only solvent residues.
Params: None.
Returns: `true` if the chain contains only solvent residues, `false` otherwise.
Side effects: None.
Errors: None.

#### `Chain::AddResidue(PResidue res) -> int`
Adds a residue to the chain.
Params: `res` — A pointer to the residue to be added.
Returns: The index of the added residue in the chain's residue table.
Side effects: Adds the residue to the chain's residue table and updates internal counters.
Errors: None.

#### `Chain::InsResidue(PResidue res, int pos) -> int`
Inserts a residue at a specified position in the chain.
Params: `res` — A pointer to the residue to be inserted. `pos` — The index of the position to insert the residue.
Returns: The index of the inserted residue in the chain's residue table.
Side effects: Inserts the residue at the specified position in the chain's residue table and adjusts internal counters.
Errors: None.

#### `Chain::InsResidue(PResidue res, int seqNum, const InsCode insCode) -> int`
Inserts a residue at a specified sequence number in the chain.
Params: `res` — A pointer to the residue to be inserted. `seqNum` — The sequence number of the residue. `insCode` — The insertion code.
Returns: The index of the inserted residue in the chain's residue table.
Side effects: Inserts the residue at the specified sequence number in the chain's residue table and adjusts internal counters.
Errors: None.

#### `Chain::DeleteResidue(int resNo) -> int`
Deletes a residue from the chain.
Params: `resNo` — The index of the residue to be deleted.
Returns: The index of the deleted residue.
Side effects: Removes the specified residue from the chain's residue table and adjusts internal counters.
Errors: None.

#### `Chain::DeleteResidue(int seqNum, const InsCode insCode) -> int`
Deletes a residue from the chain by sequence number and insertion code.
Params: `seqNum` — The sequence number of the residue to be deleted. `insCode` — The insertion code.
Returns: The index of the deleted residue.
Side effects: Removes the specified residue from the chain's residue table and adjusts internal counters.
Errors: None.

#### `Chain::DeleteAllResidues() -> int`
Deletes all residues from the chain.
Params: None.
Returns: The number of residues deleted.
Side effects: Removes all residues from the chain's residue table and resets internal counters.
Errors: None.

#### `Chain::DeleteSolvent() -> int`
Deletes all solvent residues from the chain.
Params: None.
Returns: The number of solvent residues deleted.
Side effects: Removes all solvent residues from the chain's residue table and resets internal counters.
Errors: None.

#### `Chain::TrimResidueTable() -> void`
Trims the residue table, removing any residues that are no longer present in the chain.
Params: None.
Returns: None.
Side effects: Removes any residues from the chain's residue table that are no longer present in the chain, adjusting internal counters.
Errors: None.

#### `Chain::GetNumberOfResidues() -> int`
Returns the number of residues in the chain.
Params: None.
Returns: The number of residues in the chain.
Side effects: None.
Errors: None.

#### `Chain::GetResidue(int seqNum, const InsCode insCode) -> PResidue`
Returns a pointer to a residue in the chain by sequence number and insertion code.
Params: `seqNum` — The sequence number of the residue. `insCode` — The insertion code.
Returns: A pointer to a `PResidue` object representing the specified residue.
Side effects: None.
Errors: Returns `NULL` if the residue is not found in the chain.

#### `Chain::GetResidue(int resNo) -> PResidue`
Returns a pointer to a residue in the chain by residue number.
Params: `resNo` — The index of the residue in the chain's residue table.
Returns: A pointer to a `PResidue` object representing the specified residue.
Side effects: None.
Errors: Returns `NULL` if the residue is not found in the chain.

#### `Chain::GetResidueNo(int seqNum, const InsCode insCode) -> int`
Returns the residue number in the chain's residue table.
Params: `seqNum` — The sequence number of the residue. `insCode` — The insertion code.
Returns: The residue number in the chain's residue table.
Side effects: None.
Errors: Returns -1 if the residue is not found in the chain.

#### `Chain::GetResidueTable(PPResidue & resTable, int & NumberOfResidues) -> void`
Returns the chain's residue table.
Params: `resTable` — A pointer to an array of `PResidue` objects to be populated. `NumberOfResidues` — A reference to an integer to store the number of residues in the table.
Returns: None.
Side effects: Populates the `resTable` array with pointers to the chain's residues and sets `NumberOfResidues` to the number of residues.
Errors: None.

#### `Chain::GetResidueCreate(const ResName resName, int seqNum, const InsCode insCode, bool Enforce) -> PResidue`
Creates a residue and adds it to the chain.
Params: `resName` — The name of the residue. `seqNum` — The sequence number of the residue. `insCode` — The insertion code. `Enforce` — A flag indicating whether to enforce existing residue names.
Returns: A pointer to the newly created `PResidue` object.
Side effects: Adds the new residue to the chain's residue table.
Errors: Returns `NULL` if the residue is already present in the chain or if `Enforce` is set to `false`.

#### `Chain::GetCoordHierarchy() -> PManager`
Returns a pointer to the coordinate hierarchy for this chain.
Params: None.
Returns: A pointer to a `PManager` object representing the coordinate hierarchy.
Side effects: None.
Errors: None.

#### `Chain::SetChain(const ChainID chID) -> void`
Sets the chain ID for this chain.
Params: `chID` — The chain ID to set.
Returns: None.
Side effects: Updates the chain's ID and adjusts internal references accordingly.
Errors: None.

#### `Chain::SortResidues() -> void`
Sorts the residues in the chain based on their sequence numbers.
Params: None.
Returns: None.
Side effects: Sorts the chain's residue table, updating the residue pointers.
Errors: None.

#### `Chain::TrimResidueTable() -> void`
Trims the residue table, removing any residues that are no longer present in the chain.
Params: None.
Returns: None.
Side effects: Removes any residues from the chain's residue table that are no longer present in the chain, adjusting internal counters.
Errors: None.

#### `Chain::isAminoacidChain() -> bool`
Checks if the chain contains only amino acid residues.
Params: None.
Returns: `true` if the chain contains only amino acid residues, `false` otherwise.
Side effects: None.
Errors: None.

#### `Chain::isInSelection(int selHnd) -> bool`
Checks if the chain is present in a selection.
Params: `selHnd` — The handle of the selection.
Returns: `true` if the chain is present in the selection, `false` otherwise.
Side effects: None.
Errors: None.

#### `Chain::isNucleotideChain() -> bool`
Checks if the chain contains only nucleotide residues.
Params: None.
Returns: `true` if the chain contains only nucleotide residues, `false` otherwise.
Side effects: None.
Errors: None.

#### `Chain::isSolventChain() -> bool`
Checks if the chain contains only solvent residues.
Params: None.
Returns: `true` if the chain contains only solvent residues, `false` otherwise.
Side effects: None.
Errors: None.

#### `Chain::AddResidue(PResidue res) -> int`
Adds a residue to the chain.
Params: `res` — A pointer to the residue to be added.
Returns: The index of the added residue in the chain's residue table.
Side effects: Adds the residue to the chain's residue table and updates internal counters.
Errors: None.

#### `Chain::InsResidue(PResidue res, int pos) -> int`
Inserts a residue at a specified position in the chain.
Params: `res` — A pointer to the residue to be inserted. `pos` — The index of the position to insert the residue.
Returns: The index of the inserted residue in the chain's residue table.
Side effects: Inserts the residue at the specified position in the chain's residue table and adjusts internal counters.
Errors: None.

#### `Chain::InsResidue(PResidue res, int seqNum, const InsCode insCode) -> int`
Inserts a residue at a specified sequence number in the chain.
Params: `res` — A pointer to the residue to be inserted. `seqNum` — The sequence number of the residue. `insCode` — The insertion code.
Returns: The index of the inserted residue in the chain's residue table.
Side effects: Inserts the residue at the specified sequence number in the chain's residue table and adjusts internal counters.
Errors: None.

#### `Chain::DeleteResidue(int resNo) -> int`
Deletes a residue from the chain.
Params: `resNo` — The index of the residue to be deleted.
Returns: The index of the deleted residue.
Side effects: Removes the specified residue from the chain's residue table and adjusts internal counters.
Errors: None.

#### `Chain::DeleteResidue(int seqNum, const InsCode insCode) -> int`
Deletes a residue from the chain by sequence number and insertion code

---

## `CisPep`
Represents a cis-peptide structure with associated measurement data.

#### `CisPep::CisPep()`
Constructs a new `CisPep` object with all member variables initialized to their default values.
Params: None
Returns: `void`
Side effects: Initializes `CisPep` object.

#### `CisPep::CisPep(cpstr S)`
Constructs a `CisPep` object from a C-style string representation.
Params: `S` — A pointer to the C-style string.
Returns: `void`
Side effects: Parses the string and populates the `CisPep` object.

#### `CisPep::CisPep(io::RPStream Object)`
Constructs a `CisPep` object from an `io::RPStream` object.
Params: `Object` — An `io::RPStream` object.
Returns: `void`
Side effects: Parses data from the `RPStream` object and populates the `CisPep` object.

#### `CisPep::~CisPep()`
Destructs a `CisPep` object, releasing any resources.
Params: None
Returns: `void`
Side effects: Destroys the `CisPep` object.

#### `CisPep::Copy(PContainerClass cisPep)`
Copies the contents of one `CisPep` object to another.
Params: `cisPep` — A pointer to the destination `CisPep` object.
Returns: `void`
Side effects:  Updates the destination `CisPep` object with the values of the source.

#### `CisPep::ConvertPDBASCII(cpstr S)`
Converts a PDB ASCII string to the internal `CisPep` representation.
Params: `S` — A pointer to the PDB ASCII string.
Returns: `CLASS_ID` — Returns `ClassID_CisPep` if successful, or an `ERROR_CODE` on failure.
Side effects: Parses the PDB ASCII string and populates the `CisPep` object.

#### `CisPep::GetClassID() -> CLASS_ID`
Returns the class identifier for the `CisPep` object.
Params: None
Returns: `CLASS_ID` — Returns `ClassID_CisPep`.
Side effects: None.

#### `CisPep::InitCisPep()`
Initializes the `CisPep` object to a default state.
Params: None
Returns: `void`
Side effects: Initializes the `CisPep` object’s internal state.

#### `CisPep::PDBASCIIDump(pstr S, int N)`
Writes the `CisPep` data to a PDB ASCII string.
Params: `S` — A pointer to the destination PDB ASCII string. `N` — The maximum length of the string.
Returns: `void`
Side effects: Formats the `CisPep` data into a PDB ASCII string and writes it to `S`.

#### `CisPep::read(io::RFile f)`
Reads data from an `io::RFile` object into the `CisPep` object.
Params: `f` — An `io::RFile` object.
Returns: `void`
Side effects: Reads data from the input file and populates the `CisPep` object.

#### `CisPep::write(io::RFile f)`
Writes the `CisPep` data to an `io::RFile` object.
Params: `f` — An `io::RFile` object.
Returns: `void`
Side effects: Writes the `CisPep` data to the output file.

#### `CisPep::operator=`
Assigns the values from another `CisPep` object to this object.
Params: None
Returns: `CisPep*` — Returns a pointer to the current `CisPep` object.
Side effects: Copies the contents of the source object into the destination.


---

## `ClassContainer`
This class manages a collection of container classes, likely representing structural classes within a molecular data file.

#### `ClassContainer::ClassContainer()`
Constructs a new `ClassContainer` object.

#### `ClassContainer::ClassContainer(io::RPStream Object)`
Constructs a new `ClassContainer` object from an `io::RPStream` object.

#### `ClassContainer::FreeContainer()`
Releases the memory allocated for the container data.

#### `ClassContainer::AddData(PContainerClass Data)`
Adds a `PContainerClass` object to the container.

#### `ClassContainer::GetContainerClass(int ContClassNo)` // *** USED IN COOT ***
Returns a pointer to the `PContainerClass` object at the specified index.

#### `ClassContainer::Length() // *** USED IN COOT ***`
Returns the number of `PContainerClass` objects currently stored in the container.

#### `ClassContainer::MakeCIF(mmcif::PData CIF)`
Creates a `mmcif::PData` object representing the CIF format.

#### `ClassContainer::GetCIF(mmcif::PData CIF, int ClassID)` // *** USED IN COOT ***
Creates a CIF data structure, failing if the ClassID is invalid.

#### `ClassContainer::MakeContainerClass(int ClassID)`
Creates a `PContainerClass` object with the specified `ClassID`.

#### `ClassContainer::Copy(PClassContainer CContainer)`
Copies the contents of another `ClassContainer` object, clearing the container on `NULL` parameter.

#### `ClassContainer::write(io::RFile f)`
Writes the internal state of the container to the given `io::RFile`.

#### `ClassContainer::read(io::RFile f)`
Reads the internal state of the container from the given `io::RFile`.

#### `ClassContainer::Init()`
Initializes the internal state of the `ClassContainer` object.


---

## `Compound`
Represents a compound object within an MMDB2 file.

#### `Compound::Compound()`
Constructs a new compound object.
Params: None
Returns: None
Side effects: Initializes a new compound object.
Errors: None

#### `Compound::Compound(cpstr S)`
Constructs a compound object from a PDBASCII string.
Params: `cpstr S` — A PDBASCII string representation of the compound.
Returns: None
Side effects: Initializes a compound object with the provided string.
Errors: None

#### `Compound::Compound(io::RPStream Object)`
Constructs a compound object from an RPStream.
Params: `io::RPStream Object` — An RPStream object containing the compound data.
Returns: None
Side effects: Initializes a compound object from the stream.
Errors: None

#### `Compound::~Compound()`
Destroys the compound object.
Params: None
Returns: None
Side effects: Releases resources associated with the compound object.
Errors: None

#### `Compound::PDBASCIIDump(pstr S, int N)`
Generates a PDBASCII string representation of the compound and writes it to the provided string.
Params: `pstr S` — The string to write the PDBASCII representation to, `int N` — The maximum number of characters to write.
Returns: None
Side effects: Writes the PDBASCII representation of the compound to the string `S`.
Errors: None

#### `Compound::ConvertPDBASCII(cpstr S)`
Converts a PDBASCII string to a Compound object.
Params: `cpstr S` — The PDBASCII string to convert.
Returns: `ERROR_CODE` — The error code if the conversion fails, otherwise `CLASS_ID_Compound`.
Side effects: Converts the PDBASCII string into a Compound object.
Errors: Returns `CLASS_ID_Compound` on success, otherwise returns an error code.

#### `Compound::GetClassID()`
Returns the class ID for a Compound object.
Params: None
Returns: `CLASS_ID` — The class ID for a Compound object.
Side effects: None
Errors: None

#### `Compound::write(io::RFile f)`
Writes the compound data to an RPStream.
Params: `io::RFile f` — The RPStream object to write to.
Returns: None
Side effects: Writes the compound data to the RPStream.
Errors: None

#### `Compound::read(io::RFile f)`
Reads the compound data from an RPStream.
Params: `io::RFile f` — The RPStream object to read from.
Returns: None
Side effects: Reads the compound data from the RPStream.
Errors: None


---

## `CoorManager`
Manages structural biology data, specifically models, chains, residues, and atoms, within an MMDB2 structure.

#### `CoorManager::CoorIDCode`
Integer representing the last returned atom extraction ID.
Params: None
Returns: None
Side effects: None
Errors: None

#### `CoorManager::CoorIDCode()`
Constructor.
Params: None
Returns: None
Side effects: Initializes `CoorIDCode` to 0.
Errors: None

#### `CoorManager::CoorIDCode(io::RPStream Object)`
Constructor from an `RPStream` object.
Params: `Object` — The `RPStream` object.
Returns: None
Side effects: Initializes `CoorIDCode` to 0.
Errors: None

#### `CoorManager::~CoorIDCode()`
Destructor.
Params: None
Returns: None
Side effects: No resource deallocation is performed.
Errors: None

#### `CoorManager::SetDefaultCoorID(cpstr CID)`
Sets the default atom extraction ID to the given `CID`.
Params: `CID` — The `CID` string.
Returns: Integer representing the new atom extraction ID.
Side effects: Updates the `CoorIDCode`.
Errors: None

#### `CoorManager::RemoveBricks()`
Removes all brick data.
Params: None
Returns: None
Side effects: Clears the `brick` and associated variables.
Errors: None

#### `CoorManager::areBricks() -> bool`
Checks if brick data is active.
Params: None
Returns: `true` if brick data is active, `false` otherwise.
Side effects: None
Errors: None

#### `CoorManager::MakeBricks(PPAtom atmvec, int avlen, realtype Margin, realtype BrickSize)`
Creates brick data based on given parameters.
Params: `atmvec` — Pointer to an array of atoms. `alen` — Length of the array. `Margin` — Margin for brick calculation. `BrickSize` — Size of the brick.
Returns: None
Side effects: Allocates and populates the `brick` array.
Errors: None

#### `CoorManager::GetBrickDimension(int & nxmax, int & nymax, int & nzmax)`
Gets the dimensions of the brick.
Params: `nxmax` — Pointer to the maximum x-dimension. `nymax` — Pointer to the maximum y-dimension. `nzmax` — Pointer to the maximum z-dimension.
Returns: None
Side effects: Populates `nxmax`, `nymax`, and `nzmax` with the dimensions of the brick.
Errors: None

#### `CoorManager::GetBrickCoor(PAtom A, int & nx, int & ny, int & nz)`
Gets the coordinates of the brick corners.
Params: `A` — A pointer to an atom. `nx` — Pointer to the x-coordinate. `ny` — Pointer to the y-coordinate. `nz` — Pointer to the z-coordinate.
Returns: None
Side effects: Populates `nx`, `ny`, and `nz` with the coordinates of the brick corners.
Errors: None

#### `CoorManager::GetBrickCoor(realtype x, realtype y, realtype z, int & nx, int & ny, int & nz)`
Gets the coordinates of the brick corners using decimal coordinates.
Params: `x` — X-coordinate. `y` — Y-coordinate. `z` — Z-coordinate. `nx` — Pointer to the x-coordinate. `ny` — Pointer to the y-coordinate. `nz` — Pointer to the z-coordinate.
Returns: None
Side effects: Populates `nx`, `ny`, and `nz` with the coordinates of the brick corners.
Errors: None

#### `CoorManager::GetBrickCoor(vect3 & xyz, int & nx, int & ny, int & nz)`
Gets the coordinates of the brick corners using a `vect3` object.
Params: `xyz` — A `vect3` object to populate with the coordinates. `nx` — Pointer to the x-coordinate. `ny` — Pointer to the y-coordinate. `nz` — Pointer to the z-coordinate.
Returns: None
Side effects: Populates `nx`, `ny`, and `nz` with the coordinates of the brick corners.
Errors: None

#### `CoorManager::GetBrick(int nx, int ny, int nz)`
Gets a brick data structure.
Params: `nx` — X-coordinate. `ny` — Y-coordinate. `nz` — Z-coordinate.
Returns: A pointer to a `Brick` object.
Side effects: None
Errors: None

#### `CoorManager::RemoveMBricks()`
Removes all MBrick data.
Params: None
Returns: None
Side effects: Clears the `mbrick` and associated variables.
Errors: None

#### `CoorManager::areMBricks() -> bool`
Checks if MBrick data is active.
Params: None
Returns: `true` if MBrick data is active, `false` otherwise.
Side effects: None
Errors: None

#### `CoorManager::MakeMBricks(PPAtom * atmvec, ivector avlen, int nStructures, realtype Margin, realtype BrickSize)`
Creates MBrick data based on given parameters.
Params: `atmvec` — Pointer to an array of atom pointers. `alen` — Length of the array. `nStructures` — Number of structures. `Margin` — Margin for brick calculation. `BrickSize` — Size of the brick.
Returns: None
Side effects: Allocates and populates the `mbrick` array.
Errors: None

#### `CoorManager::GetMBrickDimension(int & nxmax, int & nymax, int & nzmax)`
Gets the dimensions of the MBrick.
Params: `nxmax` — Pointer to the maximum x-dimension. `nymax` — Pointer to the maximum y-dimension. `nzmax` — Pointer to the maximum z-dimension.
Returns: None
Side effects: Populates `nxmax`, `nymax`, and `nzmax` with the dimensions of the MBrick.
Errors: None

#### `CoorManager::GetMBrickCoor(PAtom A, int & nx, int & ny, int & nz)`
Gets the coordinates of the MBrick corners.
Params: `A` — A pointer to an atom. `nx` — Pointer to the x-coordinate. `ny` — Pointer to the y-coordinate. `nz` — Pointer to the z-coordinate.
Returns: None
Side effects: Populates `nx`, `ny`, and `nz` with the coordinates of the MBrick corners.
Errors: None

#### `CoorManager::GetMBrickCoor(realtype x, realtype y, realtype z, int & nx, int & ny, int & nz)`
Gets the coordinates of the MBrick corners using decimal coordinates.
Params: `x` — X-coordinate. `y` — Y-coordinate. `z` — Z-coordinate. `nx` — Pointer to the x-coordinate. `ny` — Pointer to the y-coordinate. `nz` — Pointer to the z-coordinate.
Returns: None
Side effects: Populates `nx`, `ny`, and `nz` with the coordinates of the MBrick corners.
Errors: None

#### `CoorManager::GetMBrickCoor(vect3 & xyz, int & nx, int & ny, int & nz)`
Gets the coordinates of the MBrick corners using a `vect3` object.
Params: `xyz` — A `vect3` object to populate with the coordinates. `nx` — Pointer to the x-coordinate. `ny` — Pointer to the y-coordinate. `nz` — Pointer to the z-coordinate.
Returns: None
Side effects: Populates `nx`, `ny`, and `nz` with the coordinates of the MBrick corners.
Errors: None

#### `CoorManager::GetMBrick(int nx, int ny, int nz)`
Gets an MBrick data structure.
Params: `nx` — X-coordinate. `ny` — Y-coordinate. `nz` — Z-coordinate.
Returns: A pointer to a `MBrick` object.
Side effects: None
Errors: None

#### `CoorManager::GetNumberOfModels() -> int`
Gets the number of models.
Params: None
Returns: The number of models.
Side effects: None
Errors: None

#### `CoorManager::GetFirstModelNum() -> int`
Gets the number of the first model.
Params: None
Returns: The number of the first model.
Side effects: None
Errors: None

#### `CoorManager::GetFirstDefinedModel() -> PModel`
Gets the first defined model.
Params: None
Returns: A pointer to the first defined model.
Side effects: None
Errors: None

#### `CoorManager::GetModel(int modelNo) -> PModel`
Gets a model.
Params: `modelNo` — The model number.
Returns: A pointer to the model.
Side effects: None
Errors: None

#### `CoorManager::GetModel(cpstr CID) -> PModel`
Gets a model from a CID.
Params: `CID` — The CID string.
Returns: A pointer to the model.
Side effects: None
Errors: None

#### `CoorManager::GetModelTable(PPModel & modTable, int & NumberOfModels)`
Gets a table of models.
Params: `modTable` — Pointer to a pointer to a `Model` object. `NumberOfModels` — Pointer to the number of models.
Returns: None
Side effects: Populates `modTable` and `NumberOfModels` with model data.
Errors: None

#### `CoorManager::DeleteModel(cpstr CID)`
Deletes a model.
Params: `CID` — The CID string.
Returns: The model number.
Side effects: Deletes the specified model.
Errors: None

#### `CoorManager::DeleteModel(int modelNo)`
Deletes a model by its model number.
Params: `modelNo` — The model number.
Returns: The model number.
Side effects: Deletes the specified model.
Errors: None

#### `CoorManager::AddModel(PModel mdl)`
Adds a model.
Params: `mdl` — A pointer to the model to be added.
Returns: The model number.
Side effects: Adds the specified model.
Errors: None

#### `CoorManager::InsModel(PModel mdl, int modelNo)`
Inserts a model.
Params: `mdl` — A pointer to the model to be inserted. `modelNo` — The model number to insert it into.
Returns: The model number.
Side effects: Inserts the specified model into the model list.
Errors: None

#### `CoorManager::RotateModels(int modelNo1, int modelNo2, int rotdir)`
Rotates two models.
Params: `modelNo1` — The first model number. `modelNo2` — The second model number. `rotdir` — The rotation direction.
Returns: None
Side effects: Rotates the specified models.
Errors: None

#### `CoorManager::SwapModels(int modelNo1, int modelNo2)`
Swaps two models.
Params: `modelNo1` — The first model number. `modelNo2` — The second model number.
Returns: None
Side effects: Swaps the specified models.
Errors: None

#### `CoorManager::GetNumberOfChains(int modelNo) -> int`
Gets the number of chains in a model.
Params: `modelNo` — The model number.
Returns: The number of chains.
Side effects: None
Errors: None

#### `CoorManager::GetNumberOfChains(cpstr CID) -> int`
Gets the number of chains in a model from a CID.
Params: `CID` — The CID string.
Returns: The number of chains.
Side effects: None
Errors: None

#### `CoorManager::GetChain(int modelNo, const ChainID chainID) -> PChain`
Gets a chain.
Params: `modelNo` — The model number. `chainID` — The chain ID.
Returns: A pointer to the chain.
Side effects: None
Errors: None

#### `CoorManager::GetChain(int modelNo, int chainNo) -> PChain`
Gets a chain by its chain number.
Params: `modelNo` — The model number. `chainNo` — The chain number.
Returns: A pointer to the chain.
Side effects: None
Errors: None

#### `CoorManager::GetChain(cpstr CID) -> PChain`
Gets a chain from a CID.
Params: `CID` — The CID string.
Returns: A pointer to the chain.
Side effects: None
Errors: None

#### `CoorManager::GetChainTable(int modelNo, PPChain & chainTable, int & NumberOfChains)`
Gets a table of chains.
Params: `modelNo` — The model number. `chainTable` — Pointer to a pointer to a `Chain` object. `NumberOfChains` — Pointer to the number of chains.
Returns: None
Side effects: Populates `chainTable` and `NumberOfChains` with chain data.
Errors: None

#### `CoorManager::GetChainTable(cpstr CID, PPChain & chainTable, int & NumberOfChains)`
Gets a table of chains from a CID.
Params: `CID` — The CID string. `chainTable` — Pointer to a pointer to a `Chain` object. `NumberOfChains` — Pointer to the number of chains.
Returns: None
Side effects: Populates `chainTable` and `NumberOfChains` with chain data.
Errors: None

#### `CoorManager::DeleteChain(int modelNo, const ChainID chID)`
Deletes a chain.
Params: `modelNo` — The model number. `chID` — The chain ID.
Returns: The chain number.
Side effects: Deletes the specified chain.
Errors: None

#### `CoorManager::DeleteChain(int modelNo, int chainNo)`
Deletes a chain by its chain number.
Params: `modelNo` — The model number. `chainNo` — The chain number.
Returns: The chain number.
Side effects: Deletes the specified chain.
Errors: None

#### `CoorManager::DeleteAllChains(int modelNo)`
Deletes all chains in a model.
Params: `modelNo` — The model number.
Returns: The number of chains deleted.
Side effects: Deletes all chains in the specified model.
Errors: None

#### `CoorManager::DeleteAllChains()`
Deletes all chains in all models.
Params: None
Returns: The number of chains deleted.
Side effects: Deletes all chains in all models.
Errors: None

#### `CoorManager::AddChain(int modelNo, PChain chain)`
Adds a chain.
Params: `modelNo` — The model number. `chain` — A pointer to the chain to be added.
Returns: None
Side effects: Adds the specified chain.
Errors: None

#### `CoorManager::GetNumberOfResidues(int modelNo, const ChainID chainID) -> int`
Gets the number of residues in a chain.
Params: `modelNo` — The model number. `chainID` — The chain ID.
Returns: The number of residues.
Side effects: None
Errors: None

#### `CoorManager::GetNumberOfResidues(int modelNo, int chainNo) -> int`
Gets the number of residues in a chain by its chain number.
Params: `modelNo` — The model number. `chainNo` — The chain number.
Returns: The number of residues.
Side effects: None
Errors: None

#### `CoorManager::GetNumberOfResidues(cpstr CID) -> int`
Gets the number of residues in a chain from a CID.
Params: `CID` — The CID string.
Returns: The number of residues.
Side effects: None
Errors: None

#### `CoorManager::GetResidue(int modelNo, const ChainID chainID, int seqNo, const InsCode insCode) -> PResidue`
Gets a residue.
Params: `modelNo` — The model number. `chainID` — The chain ID. `seqNo` — The sequence number. `insCode` — The insertion code.
Returns: A pointer to the residue.
Side effects: None
Errors: None

#### `CoorManager::GetResidue(int modelNo, int chainNo, int seqNo, const InsCode insCode) -> PResidue`
Gets a residue by its chain number.
Params: `modelNo` — The model number. `chainNo` — The chain number. `seqNo` — The sequence number. `insCode` — The insertion code.
Returns: A pointer to the residue.
Side effects: None
Errors: None

#### `CoorManager::GetResidue(int modelNo, const ChainID chainID, int resNo, const InsCode insCode) -> PResidue`
Gets a residue by its residue number.
Params: `modelNo` — The model number. `chainID` — The chain ID. `resNo` — The residue number. `insCode` — The insertion code.
Returns: A pointer to the residue.
Side effects: None
Errors: None

#### `CoorManager::GetResidue(int modelNo, int chainNo, int resNo) -> PResidue`
Gets a residue by its chain number and residue number.
Params: `modelNo` — The model number. `chainNo` — The chain number. `resNo` — The residue number.
Returns: A pointer to the residue.
Side effects: None
Errors: None

#### `CoorManager::GetResidue(cpstr CID) -> PResidue`
Gets a residue from a CID.
Params: `CID` — The CID string.
Returns: A pointer to the residue.
Side effects: None
Errors: None

#### `CoorManager::GetResidueNo(int modelNo, const ChainID chainID, int seqNo, const InsCode insCode) -> int`
Gets the residue number.
Params: `modelNo` — The model number. `chainID` — The chain ID. `seqNo` — The sequence number. `insCode` — The insertion code.
Returns: The residue number.
Side effects: None
Errors: None

#### `CoorManager::GetResidueNo(int modelNo, int chainNo, int seqNo, const InsCode insCode) -> int`
Gets the residue number by its chain number.
Params: `modelNo` — The model number. `chainNo` — The chain number. `seqNo` — The sequence number. `insCode` — The insertion code.
Returns: The residue number.
Side effects: None
Errors: None

#### `CoorManager::GetResidueTable(PPResidue & resTable, int & NumberOfResidues)`
Gets a table of residues.
Params: `resTable` — Pointer to a pointer to a `Residue` object. `NumberOfResidues` — Pointer to the number of residues.
Returns: None
Side effects: Populates `resTable` and `NumberOfResidues` with residue data.
Errors: None

#### `CoorManager::GetResidueTable(int modelNo, const ChainID chainID, int seqNo, const InsCode insCode, PPResidue & resTable, int & NumberOfResidues)`
Gets a table of residues for a given chain.
Params: `modelNo` — The model number. `chainID` — The chain ID. `seqNo` — The sequence number. `insCode` — The insertion code. `resTable` — Pointer to a pointer to a `Residue` object. `NumberOfResidues` — Pointer to the number of residues.
Returns: None
Side effects: Populates `resTable` and `NumberOfResidues` with residue data.
Errors: None

#### `CoorManager::GetResidueTable(int modelNo, int chainNo, int seqNo, const InsCode insCode, PPResidue & resTable, int & NumberOfResidues)`
Gets a table of residues for a given chain by its chain number.
Params: `modelNo` — The model number. `chainNo` — The chain number. `seqNo` — The sequence number. `insCode` — The insertion code. `resTable` — Pointer to a pointer to a `Residue` object. `NumberOfResidues` — Pointer to the number of residues.
Returns: None
Side effects: Populates `resTable` and `NumberOfResidues` with residue data.
Errors: None

#### `CoorManager::GetResidueTable(int modelNo, const ChainID chainID, int resNo, int seqNo, const InsCode insCode, PPResidue & resTable, int & NumberOfResidues)`
Gets a table of residues for a given chain by its residue number.
Params: `modelNo` — The model number. `chainID` — The chain ID. `resNo` — The residue number. `seqNo` — The sequence number. `insCode` — The insertion code. `resTable` — Pointer to a pointer to a `Residue` object. `NumberOfResidues` — Pointer to the number of residues.
Returns: None
Side effects: Populates `resTable` and `NumberOfResidues` with residue data.
Errors: None

#### `CoorManager::GetResidueTable(cpstr CID, PPResidue & resTable, int & NumberOfResidues)`
Gets a table of residues from a CID.
Params: `CID` — The CID string. `resTable` — Pointer to a pointer to a `Residue` object. `NumberOfResidues` — Pointer to the number of residues.
Returns: None
Side effects: Populates `resTable` and `NumberOfResidues` with residue data.
Errors: None

#### `CoorManager::DeleteResidue(int modelNo, const ChainID chainID, int seqNo, const InsCode insCode)`
Deletes a residue.
Params: `modelNo` — The model number. `chainID` — The chain ID. `seqNo` — The sequence number. `insCode` — The insertion code.
Returns: The residue number.
Side effects: Deletes the specified residue.
Errors: None

#### `CoorManager::DeleteResidue(int modelNo, int chainNo, int seqNo, const InsCode insCode)`
Deletes a residue by its chain number.
Params: `modelNo` — The model number. `chainNo` — The chain number. `seqNo` — The sequence number. `insCode` — The insertion code.
Returns: The residue number.
Side effects: Deletes the specified residue.
Errors: None

#### `CoorManager::DeleteAllResidues(int modelNo, const ChainID chainID)`
Deletes all residues in a chain.
Params: `modelNo` — The model number. `chainID` — The chain ID.
Returns: The number of residues deleted.
Side effects: Deletes all residues in the specified chain.
Errors: None

#### `CoorManager::DeleteAllResidues(int modelNo)`
Deletes all residues in a model.
Params: `modelNo` — The model number.
Returns: The number of residues deleted.
Side effects: Deletes all residues in the specified model.
Errors: None

#### `CoorManager::DeleteAllResidues()`
Deletes all residues in all models.
Params: None
Returns: The number of residues deleted.
Side effects: Deletes all residues in all models.
Errors: None

#### `CoorManager::DeleteSolvent() `
Deletes solvent atoms.
Params: None
Returns: None
Side effects: Deletes all solvent atoms.
Errors: None

#### `CoorManager::AddAtom(int modelNo, const ChainID chainID, int seqNo, const InsCode insCode, PAtom atom)`
Adds an atom.
Params: `modelNo` — The model number. `chainID` — The chain ID. `seqNo` — The sequence number. `insCode` — The insertion code. `atom` — A pointer to the atom to be added.
Returns: None
Side effects: Adds the specified atom.
Errors: None

#### `CoorManager::AddAtom(int modelNo, int chainNo, int seqNo, const InsCode insCode, PAtom atom)`
Adds an atom by its chain number.
Params: `modelNo` — The model number. `chainNo` — The chain number. `seqNo` — The sequence number. `insCode` — The insertion code. `atom` — A pointer to the atom to be added.
Returns: None
Side effects: Adds the specified atom.
Errors: None

#### `CoorManager::GenerateSymMates(PGenSym genSym)`
Generates symmetry mates for all models.
Params: `genSym` — Pointer to a `GenSym` object.
Returns: None
Side effects: Generates symmetry mates.
Errors: None

#### `CoorManager::ApplyTransform(const mat44 & TMatrix)`
Applies the transformation matrix to all coordinates.
Params: `TMatrix` — The transformation matrix.
Returns: None
Side effects: Transforms all coordinates using the provided matrix.
Errors: None

#### `CoorManager::BringToUnitCell() -> int`
Brings all chains into the 0th unit cell.
Params: None
Returns: Number of chains that were moved to the first unit cell.
Side effects: Moves all chains into the 0th unit cell.
Errors: None

#### `CoorManager::Orth2Frac(mat44 & F, mat44 & T)`
Transforms from orthogonal to fractional coordinates.
Params: `F` — Pointer to a `mat44` object (fractional coordinates). `T` — Pointer to a `mat44` object (orthogonal coordinates).
Returns: True if the transformation was successful, False otherwise.
Side effects: Transforms coordinates.
Errors: None

#### `CoorManager::Frac2Orth(realtype & xfrac, realtype & yfrac, realtype & zfrac)`
Transforms from fractional to orthogonal coordinates.
Params: `xfrac` — Pointer to the x-coordinate. `yfrac` — Pointer to the y-coordinate. `zfrac` — Pointer to the z-coordinate.
Returns: True if the transformation was successful, False otherwise.
Side effects: Transforms coordinates.
Errors: None

#### `CoorManager::SeekContacts(PPAtom AIndex, int ilen, int atomNum, realtype dist1, realtype dist2, int seqDist, RPContact contact, int bricking, bool doSqrt)`
Searches for contacts between atoms.
Params: `AIndex` — Pointer to an array of atom indices. `ilen` — Length of the array. `atomNum` — The number of the 1st atom in the contact. `dist1` — Minimal contact distance. `dist2` — Maximal contact distance. `seqDist` — Sequence distance to neglect. `contact` — Pointer to a `RPContact` object. `bricking` - bricking control flag. `doSqrt` - enable square root calculation.
Returns: None
Side effects: Finds contacts between atoms.
Errors: None

#### `CoorManager::SeekContacts(PAtom A, PPAtom AIndex, int ilen, realtype dist1, realtype dist2, int seqDist, RPContact contact, int bricking, bool doSqrt)`
Searches for contacts between atoms using a single atom.
Params: `A` — A pointer to an atom. `AIndex` — Pointer to an array of atom indices. `ilen` — Length of the array. `dist1` — Minimal contact distance. `dist2` — Maximal contact distance. `seqDist` — Sequence distance to neglect. `contact` — Pointer to a `RPContact` object. `bricking` - bricking control flag. `doSqrt` - enable square root calculation.
Returns: None
Side effects: Finds contacts between atoms.
Errors: None

#### `CoorManager::SeekContacts(PPAtom AIndex1, int ilen1, PPAtom AIndex2, int ilen2, realtype dist1, realtype dist2, int seqDist, RPContact contact, int bricking, bool doSqrt)`
Searches for contacts between two sets of atoms.
Params: `AIndex1` — Pointer to the first atom index array. `AIndex2` — Pointer to the second atom index array. `ilen1` — Length of the first array. `ilen2` — Length of the second array. `dist1` — Minimal contact distance. `dist2` — Maximal contact distance. `seqDist` — Sequence distance to neglect. `contact` — Pointer to a `RPContact` object. `bricking` - bricking control flag. `doSqrt` - enable square root calculation.
Returns: None
Side effects: Finds contacts between two sets of atoms.
Errors: None

#### `CoorManager::SeekContacts(PPAtom AIndex1, int ilen1, PPAtom AIndex2, int ilen2, realtype contDist, PContact contact, int bricking)`
Searches for contacts between two sets of atoms optimized for speed.
Params: `AIndex1` — Pointer to the first atom index array. `AIndex2` — Pointer to the second atom index array. `ilen1` — Length of the first array. `ilen2` — Length of the second array. `contDist` — Maximal contact distance. `contact` — Pointer to a `PContact` object. `bricking` - bricking control flag.
Returns: None
Side effects: Finds contacts between two sets of atoms.
Errors: None

#### `CoorManager::SeekContacts(PPAtom AIndex1, int ilen1, PPAtom AIndex2, int ilen2, realtype contDist, PContact contact, int bricking)`
Searches for contacts between two sets of atoms optimized for speed and convenience.
Params: `AIndex1` — Pointer to the first atom index array. `AIndex2` — Pointer to the second atom index array. `ilen1` — Length of the first array. `ilen2` — Length of the second array. `contDist` — Maximal contact distance. `contact` — Pointer to a `PContact` object. `bricking` - bricking control flag.
Returns: None
Side effects: Finds contacts between two sets of atoms.
Errors: None

#### `CoorManager::SeekContacts(PPAtom AIndex1, int ilen1, PPAtom AIndex2, int ilen2, realtype contDist, PContact contact, int bricking)`
Searches for contacts between two sets of atoms optimized for speed and convenience.
Params: `AIndex1` — Pointer to the first atom index array. `AIndex2` — Pointer to the second atom index array. `ilen1` — Length of the first array. `ilen2` — Length of the second array. `contDist` — Maximal contact distance. `contact` — Pointer to a `PContact` object. `bricking` - bricking control flag.
Returns: None
Side effects: Finds contacts between two sets of atoms.
Errors: None

#### `CoorManager::InitMMDBCoorManager()`
Initializes the MMDB2 CoorManager.
Params: None
Returns: None
Side effects: Initializes the CoorManager object.
Errors: None

#### `CoorManager::ApplySymTransform(int SymMatrixNo, PGenSym genSym)`
Applies symmetry transformation.
Params: `SymMatrixNo` — Number of symmetry matrix. `genSym` — Pointer to `GenSym` object.
Returns: None
Side effects: Applies the symmetry matrix.
Errors: None

#### `CoorManager::ResetManager()`
Resets all data in the CoorManager.
Params: None
Returns: None
Side effects: Clears all data.
Errors: None

#### `CoorManager::FindSeqSection(PAtom atom, int seqDist, int & seq1, int & seq2)`
Finds the sequence section for a given atom.
Params: `atom` — A pointer to an atom. `seqDist` — Sequence distance. `seq1` — Pointer to the beginning sequence number. `seq2` — Pointer to the end sequence number.
Returns: None
Side effects: Sets the sequence numbers.
Errors: None

#### `CoorManager::areBricks() -> bool`
Checks if brick data is active.
Params: None
Returns: `true` if brick data is active, `false` otherwise.
Side effects: None
Errors: None

#### `CoorManager::areMBricks() -> bool`
Checks if MBrick data is active.
Params: None
Returns: `true` if MBrick data is active, `false` otherwise.
Side effects: None
Errors: None


---

## `Cryst`
Represents a crystal structure, including its cell parameters, symmetry, and atomic coordinates.

#### `Cryst::GetNumberOfSymOps()`
Returns the number of symmetry operations defined in the crystal.
Params: None
Returns: The number of symmetry operations (int).
Side effects: None.
Errors: None.

#### `Cryst::GetTMatrix(mat44 &TMatrix, int Nop, int cellshift_a, int cellshift_b, int cellshift_c, PSymOps symOpers=NULL)` // *** USED IN COOT ***
Calculates and returns the coordinate transformation matrix, which converts orthogonal coordinates according to the symmetry operation number Nop and places them into unit cell shifted by cellshift_a a's, cellshift_b b's and cellshift_c c's.
Params: `TMatrix` — The output transformation matrix (mat44 &). `Nop` — The symmetry operation number (int). `cellshift_a`, `cellshift_b`, `cellshift_c` — Cell shift parameters (int). `symOpers` — Symmetry operation set (PSymOps).
Returns: 0 if successful, 1 if Nop is not defined, 2 if transformation matrices were not calculated, 3 if cell parameters were not set up.
Side effects: Sets the `TMatrix` field.
Errors: None.

#### `Cryst::ConvertPDBString(pstr PDBString)`
Interprets an ASCII PDB line and fills the corresponding data fields.
Params: `PDBString` — The PDB string to interpret (pstr).
Returns: 0 if successful, a non-negative error code otherwise.
Side effects: Updates crystal data fields.
Errors: `PDBString` must be at least 81 characters long.

#### `Cryst::RWBROOKReadPrintout()`
May be invoked after reading PDB file for simulating the old RWBROOK messages and warnings.
Params: None
Returns: None.
Side effects: Prints RWBROOK messages.
Errors: None.

#### `Cryst::SetCell(realtype cell_a, realtype cell_b, realtype cell_c, realtype cell_alpha, realtype cell_beta, realtype cell_gamma, int OrthCode)`
Sets the cell parameters and orthogonalization code.
Params: `cell_a`, `cell_b`, `cell_c`, `cell_alpha`, `cell_beta`, `cell_gamma` — Cell parameters (realtype). `OrthCode` — Orthogonalization code (int).
Returns: None.
Side effects: Updates crystal data fields.
Errors: None.

#### `Cryst::PutCell(realtype cell_a, realtype cell_b, realtype cell_c, realtype cell_alpha, realtype cell_beta, realtype cell_gamma, int OrthCode)`
Sets the cell parameters and orthogonalization code.
Params: `cell_a`, `cell_b`, `cell_c`, `cell_alpha`, `cell_beta`, `cell_gamma` — Cell parameters (realtype). `OrthCode` — Orthogonalization code (int).
Returns: None.
Side effects: Updates crystal data fields.
Errors: None.

#### `Cryst::GetCell(realtype & cell_a, realtype & cell_b, realtype & cell_c, realtype & cell_alpha, realtype & cell_beta, realtype & cell_gamma, realtype & vol)`
Retrieves the cell parameters and volume.
Params: `cell_a`, `cell_b`, `cell_c`, `cell_alpha`, `cell_beta`, `cell_gamma` — Cell parameters (realtype &). `vol` — Cell volume (realtype &).
Returns: None.
Side effects: None.
Errors: None.

#### `Cryst::GetRCell(realtype & cell_as, realtype & cell_bs, realtype & cell_cs, realtype & cell_alphas, realtype & cell_betas, realtype & cell_gammas, realtype & vols)`
Retrieves the cell parameters and volume in reciprocal space.
Params: `cell_as`, `cell_bs`, `cell_cs`, `cell_alphas`, `cell_betas`, `cell_gammas`, `vols` — Cell parameters in reciprocal space (realtype &).
Returns: None.
Side effects: None.
Errors: None.

#### `Cryst::SetSpaceGroup(cpstr spGroup)`
Sets the space group.
Params: `spGroup` — The space group string (cpstr).
Returns: None.
Side effects: Updates the `spaceGroup` field.
Errors: None.

#### `Cryst::GetSpaceGroup(pstr spGroup)`
Retrieves the space group string.
Params: `spGroup` — The space group string (pstr).
Returns: None.
Side effects: Populates the `spGroup` field.
Errors: None.

#### `Cryst::GetSpaceGroupFix(pstr spGroup)`
Retrieves the space group fix.
Params: `spGroup` — The space group fix string (pstr).
Returns: None.
Side effects: Populates the `spaceGroupFix` field.
Errors: None.

#### `Cryst::CalcCoordTransforms()`
Calculates the coordinate transformation matrices after reading and converting data.
Params: None
Returns: None.
Side effects: Updates `CellCheck` field based on errors.
Errors: None.

#### `Cryst::PDBASCIIDump(io::RFile f)`
Dumps the crystal data into a PDB ASCII file.
Params: `f` — The output file stream (io::RFile).
Returns: None.
Side effects: Writes crystal data to the file.
Errors: None.

#### `Cryst::GetCIF(mmcif::PData CIF)`
Retrieves the crystal data in mmCIF format.
Params: `CIF` — The output CIF data structure (mmcif::PData).
Returns: None.
Side effects: Populates the CIF data structure.
Errors: None.

#### `Cryst::MakeCIF(mmcif::PData CIF)`
Creates the crystal data in mmCIF format.
Params: `CIF` — The output CIF data structure (mmcif::PData).
Returns: None.
Side effects: Populates the CIF data structure.
Errors: None.

#### `Cryst::areMatrices()`
Returns True if the orthogonal-to-fractional and fractional-to-orthogonal matrices are defined.
Params: None
Returns: True or False.
Side effects: None.
Errors: None.

#### `Cryst::Frac2Orth(realtype x, realtype y, realtype z, realtype & xx, realtype & yy, realtype & zz)`
Transforms fractional coordinates to orthogonal coordinates.
Params: `x`, `y`, `z` — Fractional coordinates (realtype). `xx`, `yy`, `zz` — Output orthogonal coordinates (realtype &).
Returns: True if transformation was done, False if matrices were not calculated.
Side effects: Updates the `xx`, `yy`, `zz` fields.
Errors: None.

#### `Cryst::Orth2Frac(realtype x, realtype y, realtype z, realtype & xx, realtype & yy, realtype & zz)`
Transforms orthogonal coordinates to fractional coordinates.
Params: `x`, `y`, `z` — Orthogonal coordinates (realtype). `xx`, `yy`, `zz` — Output fractional coordinates (realtype &).
Returns: True if transformation was done, False if matrices were not calculated.
Side effects: Updates the `xx`, `yy`, `zz` fields.
Errors: None.

#### `Cryst::Frac2Orth(mat44 & F, mat44 & T)`
Transforms fractional coordinates to orthogonal coordinates using matrices.
Params: `F` — Output orthogonal transformation matrix (mat44 &). `T` — Output orthogonal transformation matrix (mat44 &).
Returns: True if transformation was done, False if matrices were not calculated.
Side effects: Updates `F` and `T`.
Errors: None.

#### `Cryst::Orth2Frac(mat44 & T, mat44 & F)`
Transforms orthogonal coordinates to fractional coordinates using matrices.
Params: `T` — Output orthogonal transformation matrix (mat44 &). `F` — Output fractional transformation matrix (mat44 &).
Returns: True if transformation was done, False if matrices were not calculated.
Side effects: Updates `T` and `F`.
Errors: None.

#### `Cryst::CalcOrthMatrices()`
Calculates the RR, AC, cella's and Vol.
Params: None
Returns: None.
Side effects: Updates RR, AC, cella's and Vol.
Errors: None.

#### `Cryst::isNCSMatrix()`
Returns True if the non-crystallographic symmetry matrices are defined.
Params: None
Returns: True or False.
Side effects: None.
Errors: None.

#### `Cryst::isScaleMatrix()`
Returns True if the scale matrix is defined.
Params: None
Returns: True or False.
Side effects: None.
Errors: None.

#### `Cryst::isCellParameters()`
Returns True if the cell parameters are defined.
Params: None
Returns: True or False.
Side effects: None.
Errors: None.

#### `Cryst::GetNumberOfSymOps()` // *** USED IN COOT ***
Returns the number of symmetry operations defined in the crystal.
Params: None
Returns: The number of symmetry operations (int).
Side effects: None.
Errors: None.

#### `Cryst::GetSymOp(int Nop)` // *** USED IN COOT ***
Retrieves the symmetry operation data for a given symmetry operation number.
Params: `Nop` — The symmetry operation number (int).
Returns: None.
Side effects: None.
Errors: None.

#### `Cryst::GetNumberOfNCSMatrices()` // Returns the number of NCS mates not given in the file (iGiven==0)
Returns the number of non-crystallographic symmetry matrices.
Params: None
Returns: The number of NCS matrices (int).
Side effects: None.
Errors: None.

#### `Cryst::GetNCSMatrix(int NCSMatrixNo, mat33 & ncs_m, vect3 & ncs_v)`
Retrieves the non-crystallographic symmetry matrix data.
Params: `NCSMatrixNo` — The NCS matrix number (int). `ncs_m` — The output NCS matrix (mat33 &). `ncs_v` — The output NCS vector (vect3 &).
Returns: None.
Side effects: Updates `ncs_m` and `ncs_v`.
Errors: None.

#### `Cryst::GetNCSMatrix(int NCSMatrixNo, mat44 & ncs_m, int & iGiven)` // no=0..N-1
Retrieves the non-crystallographic symmetry matrix data.
Params: `NCSMatrixNo` — The NCS matrix number (int). `ncs_m` — The output NCS matrix (mat44 &). `iGiven` — Indicator of whether this matrix was given in the file (int &).
Returns: None.
Side effects: Updates `ncs_m` and `iGiven`.
Errors: None.

#### `Cryst::AddNCSMatrix(mat33 & ncs_m, vect3 & ncs_v, int iGiven)`
Adds a non-crystallographic symmetry matrix to the list.
Params: `ncs_m` — The NCS matrix (mat33 &). `ncs_v` — The NCS vector (vect3 &). `iGiven` — Indicator of whether this matrix was given in the file (int).
Returns: None.
Side effects: Updates `ncs_m` and `ncs_v`.
Errors: None.

#### `Cryst::GetTMatrix(mat44 &TMatrix, int Nop, int cellshift_a, int cellshift_b, int cellshift_c, PSymOps symOpers=NULL)` // *** USED IN COOT ***
Calculates and returns the coordinate transformation matrix.
Params: `TMatrix` — The output transformation matrix (mat44 &). `Nop` — The symmetry operation number (int). `cellshift_a`, `cellshift_b`, `cellshift_c` — Cell shift parameters (int). `symOpers` — Symmetry operation set (PSymOps).
Returns: 0 if successful, 1 if Nop is not defined, 2 if transformation matrices were not calculated, 3 if cell parameters were not set up.
Side effects: Sets the `TMatrix` field.
Errors: None.

#### `Cryst::GetUCTMatrix(mat44 &TMatrix, int Nop, realtype x, realtype y, realtype z, int cellshift_a, int cellshift_b, int cellshift_c, PSymOps symOpers=NULL)`
Calculates and returns the unit-cell transformation matrix.
Params: `TMatrix` — The output transformation matrix (mat44 &). `Nop` — The symmetry operation number (int). `x`, `y`, `z` — Fractional coordinates (realtype). `cellshift_a`, `cellshift_b`, `cellshift_c` — Cell shift parameters (int). `symOpers` — Symmetry operation set (PSymOps).
Returns: 0 if successful, 1 if Nop is not defined, 2 if transformation matrices were not calculated, 3 if cell parameters were not set up.
Side effects: Sets the `TMatrix` field.
Errors: None.

#### `Cryst::GetFractMatrix(mat44 &TMatrix, int Nop, int cellshift_a, int cellshift_b, int cellshift_c, PSymOps symOpers=NULL)`
Calculates and returns the fractional coordinate transformation matrix.
Params: `TMatrix` — The output transformation matrix (mat44 &). `Nop` — The symmetry operation number (int). `cellshift_a`, `cellshift_b`, `cellshift_c` — Cell shift parameters (int). `symOpers` — Symmetry operation set (PSymOps).
Returns: 0 if successful, 1 if Nop is not defined, 2 if transformation matrices were not calculated, 3 if cell parameters were not set up.
Side effects: Sets the `TMatrix` field.
Errors: None.

#### `Cryst::GetSymOpMatrix(mat44 &TMatrix, int Nop)`
Calculates and returns the transformation matrix for the given symmetry operation number.
Params: `TMatrix` — The output transformation matrix (mat44 &). `Nop` — The symmetry operation number (int).
Returns: 0 if successful, 1 if Nop is not defined, 2 if transformation matrices were not calculated, 3 if cell parameters were not set up.
Side effects: Sets the `TMatrix` field.
Errors: None.

#### `Cryst::Copy(PCryst Cryst)`
Copies crystal data from another `Cryst` object.
Params: `Cryst` — The source `Cryst` object (PCryst).
Returns: None.
Side effects: Copies data to this object.
Errors: None.

#### `Cryst::write(io::RFile f)`
Writes header to PDB binary file.
Params: `f` — The output file stream (io::RFile).
Returns: None.
Side effects: Writes header to the file.
Errors: None.

#### `Cryst::read(io::RFile f)`
Reads header from PDB binary file.
Params: `f` — The input file stream (io::RFile).
Returns: None.
Side effects: Reads header from the file.
Errors: None.

#### `Cryst::Init(bool fullInit)`
Initializes the crystal data.
Params: `fullInit` — Boolean flag indicating if full initialization is needed.
Returns: None.
Side effects: Initializes various fields of the `Cryst` object.
Errors: None.


---

## `DBReference`
Represents a reference to a sequence database entry, typically a PDB file.

#### `DBReference::DBReference()`
Constructs a DBReference object with default values.
Params: None
Returns: None
Side effects: Initializes the DBReference object.
Errors: None

#### `DBReference::DBReference(PChain Chain_Owner)`
Constructs a DBReference object from a PChain.
Params: `Chain_Owner` — The PChain object to which this DBReference belongs.
Returns: None
Side effects: Initializes the DBReference object, linking it to the PChain.
Errors: None

#### `DBReference::DBReference(PChain Chain_Owner, cpstr S)`
Constructs a DBReference object from a PChain and a string.
Params: `Chain_Owner` — The PChain object to which this DBReference belongs, `S` — A string containing PDB ASCII data.
Returns: None
Side effects: Initializes the DBReference object, linking it to the PChain and parsing the PDB ASCII data.
Errors: None

#### `DBReference::DBReference(io::RPStream Object)`
Constructs a DBReference object from an RPStream.
Params: `Object` — An io::RPStream object.
Returns: None
Side effects: Initializes the DBReference object, reading data from the RPStream.
Errors: None

#### `DBReference::~DBReference()`
Destroys a DBReference object, releasing any associated resources.
Params: None
Returns: None
Side effects: Deallocates memory and cleans up resources used by the DBReference object.
Errors: None

#### `DBReference::PDBASCIIDump(pstr S, int N)`
Writes the PDB ASCII data of the DBReference to a string.
Params: `S` — A pointer to a string where the PDB ASCII data will be written, `N` — The maximum number of characters to write to the string.
Returns: None
Side effects: Writes the PDB ASCII data to the string `S`.
Errors: None

#### `DBReference::MakeCIF(mmcif::PData CIF, int N)`
Generates a CIF (Collaborative Indexing Facility) data object from the DBReference.
Params: `CIF` — A mmcif::PData object to populate, `N` — A dummy parameter.
Returns: None
Side effects: Constructs a CIF data object using information from the DBReference.
Errors: None

#### `DBReference::ConvertPDBASCII(cpstr S)`
Converts a PDB ASCII string to the internal representation of the DBReference.
Params: `S` — A pointer to the PDB ASCII string.
Returns: `ERROR_CODE` — The conversion result, typically 0 on success.
Side effects: Parses the PDB ASCII data and populates the DBReference object.
Errors: Returns an error code if the PDB ASCII string is invalid.

#### `DBReference::GetCIF(mmcif::PData CIF, int & n)`
Retrieves a CIF data object from the DBReference.
Params: `CIF` — A mmcif::PData object to populate, `n` — A reference to an integer that will store the number of atoms in the CIF.
Returns: `ERROR_CODE` — The retrieval result, typically 0 on success.
Side effects: Constructs a CIF data object using information from the DBReference and populates `n` with the number of atoms.
Errors: Returns an error code if the CIF creation fails.

#### `DBReference::GetClassID() -> ClassID`
Returns the class identifier for the DBReference.
Params: None
Returns: `ClassID_DBReference` — The class identifier.
Side effects: None
Errors: None

#### `DBReference::Copy(PContainerClass DBRef)`
Copies the contents of one DBReference to another.
Params: `DBRef` — A pointer to another DBReference object to copy to.
Returns: None
Side effects: Copies the data from the current DBReference to the destination DBReference.
Errors: None

#### `DBReference::write(io::RFile f)`
Writes the DBReference data to an RFile.
Params: `f` — An io::RFile object to write to.
Returns: None
Side effects: Serializes the DBReference object to the RFile.
Errors: None

#### `DBReference::read(io::RFile f)`
Reads the DBReference data from an RFile.
Params: `f` — An io::RFile object to read from.
Returns: None
Side effects: Deserializes the DBReference object from the RFile.
Errors: None

#### `DBReference::InitDBReference()`
Initializes the DBReference object's member variables.
Params: None
Returns: None
Side effects: Sets the initial values for the DBReference object's member variables.
Errors: None


---

## `Data`
The `Data` class provides a flexible interface for reading and writing data from MMDB2 structure files.

#### `Data::Data()`
Basic constructor.

#### `Data::Data(cpstr N)`
Constructor that assigns data block name.

#### `Data::Data(io::RPStream Object)`
Constructor for MMDB data streaming functions.

#### `Data::~Data()`
Destructor.

#### `Data::SetPrintWarnings(bool SPW)`
Sets flag to print warnings on reading mmCIF files.

#### `Data::SetStopOnWarning(bool SOW)`
Sets flag to stop on warning when reading an mmCIF file.

#### `Data::SetFlag(CIF_FLAG F)`
Sets optional flag(s) for reading mmCIF files.

#### `Data::RemoveFlag(CIF_FLAG F)`
Removes optional flag(s) for reading mmCIF files.

#### `Data::GetWarnings() -> int`
Returns bit-encoded warnings issued at last file read.

#### `Data::SetWrongFields(cpstr *cats, cpstr *tags)`
Sets category names and tags that are to be ignored on file read.

#### `Data::ReadMMCIFData(cpstr FName, io::GZ_MODE gzipMode)`
Reads mmCIF data block from file.

#### `Data::ReadMMCIFData(io::RFile f, pstr S, int & lcount)`
Reads sequential mmCIF data blocks from file.

#### `Data::WriteMMCIFData(cpstr FName, io::GZ_MODE gzipMode)`
Writes mmCIF data block into file.

#### `Data::WriteMMCIF(io::RFile f)`
Writes (next) mmCIF data block into file.

#### `Data::GetNumberOfCategories() -> int`
Retrieves the number of categories (structures and loops) in data block.

#### `Data::GetCategory(int categoryNo)`
Retrieves pointer to category (a structure or a loop) by category number.

#### `Data::GetStructure(cpstr CName)`
Retrieves mmCIF structure with given name.

#### `Data::GetLoop(cpstr CName)`
Retrieves mmCIF loop with given name.

#### `Data::FindLoop(cpstr * tagList)`
Finds loop containing all tags from the tag list.

#### `Data::GetDataName(pstr & dname, bool Remove)`
Retrieves the data block name into dynamically-allocated string.

#### `Data::GetDataName() -> pstr`
Returns the data block name.

#### `Data::CheckData(cpstr CName, cpstr TName)`
Checks whether a field exists in the CIF file.

#### `Data::DeleteCategory(cpstr CName)`
Deletes a category.

#### `Data::Copy(PData Data)`
Copies data from another Data object.

#### `Data::CopyCategory(PData Data, cpstr CName, cpstr newCName)`
Copies a category from another Data object.

#### `Data::PrintCategories()`
Prints category names for debugging purposes.

#### `Data::write(io::RFile f)`
Writes the current data block to a file.

#### `Data::read(io::RFile f)`
Reads data from a file.

#### `Data::InitData()`
Initializes the data block.

#### `Data::FreeWrongFields()`
Frees memory allocated for incorrect fields.

#### `Data::CheckWrongField(cpstr C, cpstr T)`
Checks if a particular field is wrong.

#### `Data::Sort()`
Sorts the data (implementation detail).

#### `Data::GetCategoryNo(cpstr cname)`
Gets the index of the category in the category list.

#### `Data::AddCategory(cpstr cname)`
Adds a category to the data block.

#### `Data::DeleteCategory(int CatNo)`
Deletes a category from the data block.

#### `Data::GetDataItem(io::RFile f, pstr S, pstr & L, pstr & p, int & lcount, int & llen)`
Reads data from a file, using internal storage.

#### `Data::GetLoop(io::RFile f, pstr S, pstr & L, pstr & p, int & lcount, int & llen)`
Reads data from a loop.

#### `Data::GetField(io::RFile f, pstr S, pstr & L, pstr & p, int & lcount, int & llen)`
Reads a field from the file, handling potential errors.

#### `Data::PutNoData(int NoDataType, cpstr CName, cpstr TName, int nrow)`
Puts data into a category without using the specific field name.

#### `Data::PutString(cpstr S, cpstr CName, cpstr TName, bool Concatenate)`
Puts a string value into a category and tag.

#### `Data::PutDate(cpstr CName, cpstr TName)`
Puts a date value into a category and tag.

#### `Data::PutReal(realtype R, cpstr CName, cpstr TName, int prec)`
Puts a real value into a category and tag.

#### `Data::PutInteger(int I, cpstr CName, cpstr TName)`
Puts an integer value into a category and tag.

#### `Data::PutLoopNoData(int NoDataType, cpstr CName, cpstr TName, int nrow)`
Puts data into a loop category without using the specific field name.

#### `Data::PutLoopString(cpstr S, cpstr CName, cpstr TName, int nrow)`
Puts a string value into a loop category and tag.

#### `Data::PutLoopReal(realtype R, cpstr CName, cpstr TName, int nrow, int prec)`
Puts a real value into a loop category and tag.

#### `Data::PutLoopInteger(int I, cpstr CName, cpstr TName, int nrow)`
Puts an integer value into a loop category and tag.

#### `Data::PutLoopSVector(psvector S, cpstr CName, cpstr TName, int i1, int i2)`
Puts strings from a vector into a loop.

#### `Data::PutLoopRVector(rvector R, cpstr CName, cpstr TName, int i1, int i2, int prec)`
Puts reals from a vector into a loop.

#### `Data::PutLoopIVector(ivector I, cpstr CName, cpstr TName, int i1, int i2)`
Puts integers from a vector into a loop.

#### `Data::RenameCategory(cpstr CName, cpstr newCName)`
Renames a category.

#### `Data::InitData()`
Initializes the data block.

#### `Data::FreeWrongFields()`
Frees memory allocated for incorrect fields.

#### `Data::CheckWrongField(cpstr C, cpstr T)`
Checks if a particular field is wrong.

#### `Data::Sort()`
Sorts the data (implementation detail).


---

## `Edge`
Represents an edge in a graph structure.

#### `Edge::GetType() -> int`
Returns the type of the edge.
Params: None
Returns: The edge type as an integer.
Side effects: None
Errors: None

#### `Edge::GetVertex1() -> int`
Returns the index of the first vertex connected by the edge (numbered starting from 1).
Params: None
Returns: The index of the first vertex.
Side effects: None
Errors: None

#### `Edge::GetVertex2() -> int`
Returns the index of the second vertex connected by the edge (numbered starting from 1).
Params: None
Returns: The index of the second vertex.
Side effects: None
Errors: None


---

## `File`
This class represents a data file, typically a PDB or MMCIF file, for structural biology data.

#### `File::File()`
Constructs a new `File` object.
Params: None
Returns: None
Side effects: Initializes the `File` object.

#### `File::File(cpstr FName, io::GZ_MODE gzipMode=io::GZM_CHECK)`
Constructs a new `File` object from a file name.
Params: `cpstr FName` — The name of the file to load. `io::GZ_MODE gzipMode` — The gzip mode to use.
Returns: None
Side effects: Loads the data from the specified file.

#### `File::File(io::RPStream Object)`
Constructs a new `File` object from an `RPStream` object.
Params: `io::RPStream Object` — The `RPStream` object to use.
Returns: None
Side effects: Loads the data from the `RPStream` object.

#### `File::~File()`
Destroys the `File` object, releasing associated resources.
Params: None
Returns: None
Side effects: Frees memory allocated by the `File` object.

#### `File::AddCIFData(cpstr DName)`
Adds a data entry to the file, using its name.
Params: `cpstr DName` — The name of the data entry to add.
Returns: `int` — The index of the newly added data entry, or -1 on error.
Side effects: Adds the specified data entry to the file’s internal data structure.

#### `File::DeleteCIFData(int dataNo)`
Deletes a data entry from the file, using its index.
Params: `int dataNo` — The index of the data entry to delete.
Returns: `int` — 0 on success, -1 on error.
Side effects: Removes the data entry specified by index from the file's internal data structure.

#### `File::DeleteCIFData(cpstr DName)`
Deletes a data entry from the file, using its name.
Params: `cpstr DName` — The name of the data entry to delete.
Returns: `int` — 0 on success, -1 on error.
Side effects: Removes the data entry specified by name from the file's internal data structure.

#### `File::DeleteCIFData(int dataNo)`
Deletes a data entry from the file, using its index.
Params: `int dataNo` — The index of the data entry to delete.
Returns: `int` — 0 on success, -1 on error.
Side effects: Removes the data entry specified by index from the file's internal data structure.

#### `File::GetCIFData(int dataNo)`
Returns a pointer to a data entry in the file, given its index.
Params: `int dataNo` — The index of the data entry to retrieve.
Returns: `PData` — A pointer to the data entry, or `nullptr` if invalid index.
Side effects: Returns the data entry specified by index.

#### `File::GetCIFData(cpstr DName)`
Returns a pointer to a data entry in the file, given its name.
Params: `cpstr DName` — The name of the data entry to retrieve.
Returns: `PData` — A pointer to the data entry, or `nullptr` if not found.
Side effects: Returns the data entry specified by name.

#### `File::GetCIFDataNo(cpstr DName)`
Returns the index of a data entry in the file, given its name.
Params: `cpstr DName` — The name of the data entry to retrieve.
Returns: `int` — The index of the data entry, or -1 if not found.
Side effects: Returns the index of the data entry specified by name.

#### `File::GetNofData() -> int`
Returns the total number of data entries in the file.
Params: None
Returns: `int` — The number of data entries.
Side effects: None.

#### `File::ReadMMCIFFile(cpstr FName, io::GZ_MODE gzipMode=io::GZM_CHECK)`
Reads an MMCIF file into the `File` object.
Params: `cpstr FName` — The name of the MMCIF file to read. `io::GZ_MODE gzipMode` — The gzip mode to use.
Returns: `int` — 0 on success, -1 on error.
Side effects: Loads the data from the MMCIF file.

#### `File::Rewrite()`
Rewrites the file to the state of the current `File` object.
Params: None
Returns: None
Side effects: Writes the data to the specified file in the current state.

#### `File::Shut()`
Shuts down the `File` object, releasing any associated resources.
Params: None
Returns: None
Side effects: Releases resources and prepares the `File` object for destruction.

#### `File::WriteMMCIF(io::RFile f)`
Writes the current file contents to the specified `RFile` object.
Params: `io::RFile f` — The `RFile` object to write to.
Returns: None
Side effects: Writes the data to the `RFile` object.

#### `File::WriteMMCIFFile(cpstr FName, io::GZ_MODE gzipMode=io::GZM_CHECK)`
Writes the current file contents to a new MMCIF file.
Params: `cpstr FName` — The name of the MMCIF file to write to. `io::GZ_MODE gzipMode` — The gzip mode to use.
Returns: `int` — 0 on success, -1 on error.
Side effects: Writes the data to the specified MMCIF file.

#### `File::assign(PFile File)`
Assigns the contents of another `File` object to the current object.
Params: `PFile File` — The `File` object to copy from.
Returns: None
Side effects: Copies the data from the source `File` object to the current object.

#### `File::copy(PFile File)`
Copies the contents of another `File` object to the current object.
Params: `PFile File` — The `File` object to copy from.
Returns: None
Side effects: Copies the data from the source `File` object to the current object.

#### `File::FreeMemory()`
Frees the memory allocated by the `File` object.
Params: None
Returns: None
Side effects: Releases the memory used by the `File` object.

#### `File::InitFile()`
Initializes the `File` object, setting up internal data structures.
Params: None
Returns: None
Side effects: Sets up the initial state of the `File` object.

#### `File::Read(io::RFile f)`
Reads the data from the input stream `RFile` object.
Params: `io::RFile f` — The input stream to read from.
Returns: None
Side effects: Reads the data from the input stream, updating the internal state of the `File` object.

#### `File::Sort()`
Sorts the internal list of data entries.
Params: None
Returns: None
Side effects: Sorts the internal data structure.

#### `File::TakeFilePool()`
Takes ownership of the file pool.
Params: None
Returns: None
Side effects: Sets the file pool to the current `File` object.

#### `File::Write(io::RFile f)`
Writes the data from the `File` object to the output stream `RFile` object.
Params: `io::RFile f` — The output stream to write to.
Returns: None
Side effects: Writes the data to the output stream, updating the internal state of the `File` object.


---

## `Graph`
Represents a graph of atoms, useful for structural analysis and symmetry relief.

#### `Graph::SetName(cpstr gname)`
Sets the name of the graph.
Params: `gname` — The name to assign to the graph.
Returns: `void`.
Side effects: Sets the internal `name` member.
Errors: None.

#### `Graph::AddVertex(PVertex V)`
Adds a vertex to the graph.
Params: `V` — A pointer to the vertex to be added.
Returns: `void`.
Side effects: Adds the vertex to the graph’s internal vertex list.
Errors: None.

#### `Graph::AddEdge(PEdge G)`
Adds an edge to the graph.
Params: `G` — A pointer to the edge to be added.
Returns: `void`.
Side effects: Adds the edge to the graph’s internal edge list.
Errors: None.

#### `Graph::SetVertices(PPVertex V, int vlen)`
Sets the list of vertices in the graph.
Params: `V` — A pointer to an array of vertex pointers. `vlen` — The number of vertices in the array.
Returns: `void`.
Side effects: Updates the internal vertex list with the provided vertices.
Errors: None.

#### `Graph::SetEdges(PEdge G, int glen)`
Sets the list of edges in the graph.
Params: `G` — A pointer to an array of edge pointers. `glen` — The number of edges in the array.
Returns: `void`.
Side effects: Updates the internal edge list with the provided edges.
Errors: None.

#### `Graph::RemoveChirality()`
Removes chirality information from the graph.
Params: None
Returns: `void`.
Side effects: Removes chirality data, potentially simplifying the graph.
Errors: None.

#### `Graph::LeaveChirality(int eltype)`
Specifies the type of chirality to leave for the graph.
Params: `eltype` — The type of chirality to keep.
Returns: `void`.
Side effects: Sets the chirality type that will be used.
Errors: None.

#### `Graph::MakeGraph(PResidue R, cpstr altLoc=NULL)`
Creates a graph representing the specified residue.
Params: `R` — A pointer to the residue. `altLoc` — An optional alternate location specification.
Returns: `int` — Returns `MKGRAPH_Ok` on success, or an error code if the graph cannot be constructed (e.g., no atoms, incorrect altcode).
Side effects: Populates the graph with vertices and edges based on the residue’s atoms.
Errors:  `MKGRAPH_NoAtoms`, `MKGRAPH_ChangedAltLoc`, `MKGRAPH_MaxOccupancy`.

#### `Graph::MakeGraph(PPAtom atom, int nAtoms)`
Creates a graph from a list of atoms.
Params: `atom` — An array of atom pointers. `nAtoms` — The number of atoms in the array.
Returns: `int` — Returns `MKGRAPH_Ok` on success, or an error code if the graph cannot be constructed.
Side effects: Populates the graph with vertices and edges based on the provided atoms.
Errors: None.

#### `Graph::HideType(int bond_vx_type)`
Hides a specific bond type.
Params: `bond_vx_type` — The type of bond to hide.
Returns: `void`.
Side effects: Removes the specified bond type from the graph, preventing it from being considered during matching.
Errors: None.

#### `Graph::ExcludeType(int type)`
Excludes a specific bond type.
Params: `type` — The type of bond to exclude.
Returns: `void`.
Side effects: Removes the specified bond type from the graph, preventing it from being considered during matching.
Errors: None.

#### `Graph::MakeSymmetryRelief(bool noCO2)`
Sets up symmetry relief parameters.
Params: `noCO2` — A boolean indicating whether to include CO2 in symmetry calculations.
Returns: `void`.
Side effects: Modifies the graph to consider symmetry relationships.
Errors: None.

#### `Graph::IdentifyRings()`
Identifies and labels rings in the graph.
Params: None
Returns: `void`.
Side effects: Analyzes the graph to identify and label rings, potentially affecting subsequent graph matching.
Errors: None.

#### `Graph::IdentifyConnectedComponents()`
Identifies and counts connected components in the graph.
Params: None
Returns: `int` — The number of connected components (greater than or equal to 1).
Side effects: Analyzes the graph to find all connected components.
Errors: None.

#### `Graph::Build(bool bondOrder)`
Builds the graph structure, establishing edges based on atom distances.
Params: `bondOrder` — A boolean indicating whether to build bonds based on distance.
Returns: `int` — Returns `0` on success, or an error code if the graph building fails.
Side effects: Calculates and establishes chemical bonds between atoms based on specified cut-off distances.
Errors: None.

#### `Graph::MakeVertexIDs()`
Assigns unique integer IDs to all vertices in the graph.
Params: None
Returns: `void`.
Side effects: Numerates the vertices sequentially starting from 1.
Errors: None.

#### `Graph::GetVertexID(int vertexNo)`
Gets the integer ID of a vertex.
Params: `vertexNo` — The integer ID of the vertex.
Returns: `int` — The integer ID of the vertex.
Side effects: Returns the vertex ID.
Errors: `vertexNo` must be between 1 and the number of vertices in the graph.

#### `Graph::GetVertexNo(cpstr vname)`
Gets the name of a vertex.
Params: `vname` — The name of the vertex.
Returns: `int` — The integer ID of the vertex.
Side effects: Returns the vertex ID.
Errors: None.

#### `Graph::GetNBondedVertices(int vertexNo)`
Gets the number of bonded vertices to a given vertex.
Params: `vertexNo` — The vertex ID.
Returns: `int` — The number of bonded vertices.
Side effects: Returns the number of vertices that are connected to the specified vertex through bonds.
Errors: `vertexNo` must be between 1 and the number of vertices in the graph.

#### `Graph::GetBondedVertexID(int vertexNo, int bond_vx_type, int bondNo)`
Gets the ID of a bonded vertex given the vertex number, bond type and bond number.
Params: `vertexNo` — The vertex ID. `bond_vx_type` — The bond type. `bondNo` — The bond number.
Returns: `int` — The ID of the bonded vertex.
Side effects: Returns the ID of the vertex bonded to the specified vertex through the specified bond type and bond number.
Errors: None.

#### `Graph::GetVertex(int vertexNo)`
Retrieves a vertex object by its integer ID.
Params: `vertexNo` — The integer ID of the vertex.
Returns: `PVertex` — A pointer to the vertex object.
Side effects: Returns the vertex object.
Errors: `vertexNo` must be between 1 and the number of vertices in the graph.

#### `Graph::GetNofVertices()`
Returns the number of vertices in the graph.
Params: None
Returns: `int` — The number of vertices.
Side effects: Returns the number of vertices in the graph.
Errors: None.

#### `Graph::GetEdge(int edgeNo)`
Retrieves an edge object by its integer ID.
Params: `edgeNo` — The integer ID of the edge.
Returns: `PEdge` — A pointer to the edge object.
Side effects: Returns the edge object.
Errors: `edgeNo` must be between 1 and the number of edges in the graph.

#### `Graph::GetNofEdges()`
Returns the number of edges in the graph.
Params: None
Returns: `int` — The number of edges.
Side effects: Returns the number of edges in the graph.
Errors: None.

#### `Graph::GetVertices(PPVertex & V, int & nV)`
Retrieves all vertices in the graph.
Params: `V` — A pointer to an array of vertex pointers. `nV` — A reference to the number of vertices to be stored.
Returns: `void`.
Side effects: Populates the `V` array with pointers to all vertices in the graph and sets `nV` to the number of vertices.
Errors: None.

#### `Graph::GetEdges(PEdge & E, int & nE)`
Retrieves all edges in the graph.
Params: `E` — A pointer to an array of edge pointers. `nE` — A reference to the number of edges to be stored.
Returns: `void`.
Side effects: Populates the `E` array with pointers to all edges in the graph and sets `nE` to the number of edges.
Errors: None.

#### `Graph::Print()`
Prints the contents of the graph to the console.
Params: None
Returns: `void`.
Side effects: Prints a textual representation of the graph’s structure.
Errors: None.

#### `Graph::Print1()`
Prints one vertex of the graph to the console.
Params: None
Returns: `void`.
Side effects: Prints the information on a single vertex of the graph.
Errors: None.

#### `Graph::Copy(PGraph G)`
Copies the graph to another graph.
Params: `G` — A pointer to the destination graph.
Returns: `void`.
Side effects: Creates a deep copy of the current graph, populating the destination graph with the same vertices and edges.
Errors: None.

#### `Graph::read(io::RFile f)`
Reads data from an input file.
Params: `f` — An input file stream.
Returns: `void`.
Side effects: Reads data from the input file and updates the graph’s internal state.
Errors: None.

#### `Graph::write(io::RFile f)`
Writes data to an output file.
Params: `f` — An output file stream.
Returns: `void`.
Side effects: Writes data from the graph’s internal state to the output file.
Errors: None.

#### `Graph::mem_read(cpstr S, int & l)`
Reads a string from memory.
Params: `S` — A pointer to the memory buffer. `l` — A reference to the length of the string.
Returns: `void`.
Side effects: Reads a string from the memory buffer and updates the length.
Errors: None.

#### `Graph::mem_write(pstr S, int & l)`
Writes a string to memory.
Params: `S` — A pointer to the memory buffer. `l` — A reference to the length of the string.
Returns: `void`.
Side effects: Writes a string to the memory buffer and updates the length.
Errors: None.

#### `Graph::InitGraph()`
Initializes the graph data structures.
Params: None
Returns: `void`.
Side effects: Initializes the internal data structures of the graph, such as vertex and edge lists.
Errors: None.

#### `Graph::FreeMemory()`
Frees the memory allocated for the graph.
Params: None
Returns: `void`.
Side effects: Deallocates all memory allocated for the graph’s internal data structures.
Errors: None.

#### `Graph::markConnected(int vno, int cno)`
Marks vertices in a connected component.
Params: `vno` — The ID of the vertex. `cno` — The component number.
Returns: `void`.
Side effects: Marks all vertices belonging to the same connected component as `vno`.
Errors: None.


---

## `GraphMatch`
This class implements an algorithm for finding maximal common subgraphs between two graphs.

#### `GraphMatch::SetTimeLimit(int maxTimeToRun=0) -> void`
Sets a maximum time limit for the matching process, preventing excessive runtime.
Params: `maxTimeToRun` — The maximum execution time in seconds; defaults to 0 (no limit).
Returns: `void`.
Side effects: Sets the `timeLimit` member variable.

#### `GraphMatch::GetMatch(int MatchNo, ivector & FV1, ivector & FV2, int & nv, realtype & p1, realtype & p2) -> void`
Returns the vertices of a specific matching subgraph, along with associated edge weights, and the number of vertices.
Params: `MatchNo` — The index of the matching subgraph to retrieve.
FV1: `ivector &` — A reference to a vector containing the vertex IDs of the first subgraph.
FV2: `ivector &` — A reference to a vector containing the vertex IDs of the second subgraph.
nv: `int &` — A reference to an integer storing the number of vertices in the subgraph.
p1: `realtype &` — A reference to a floating-point value representing the edge weight between vertices in the first subgraph.
p2: `realtype &` — A reference to a floating-point value representing the edge weight between vertices in the second subgraph.
Returns: `void`.
Side effects: Populates the provided vectors with vertex IDs and edge weights.

#### `GraphMatch::MatchGraphs(PGraph Gh1, PGraph Gh2, int minMatch, bool vertexType, VERTEX_EXT_TYPE vertexExt) -> void`
Finds maximal common subgraphs of a minimum size between two graphs, considering vertex types.
Params: `Gh1` — Pointer to the first graph.
`Gh2` — Pointer to the second graph.
`minMatch` — The minimum size of the common subgraph to find.
`vertexType` — Whether vertex types should be taken into account.
`vertexExt` — The type of extended vertex type to use.
Returns: `void`.
Side effects:  Internal call to the matching algorithm, potentially modifying internal data structures.

#### `GraphMatch::GetNofMatches() -> int`
Returns the number of matched subgraphs found.
Params: None.
Returns: `int` — The number of matched subgraphs.
Side effects: None.

#### `GraphMatch::GetNofMatches() -> int`
Returns the number of found matching subgraphs.
Params: None.
Returns: `int` — The number of matching subgraphs.
Side effects: None.

#### `GraphMatch::GetMaxMatchSize() -> int`
Returns the size of the largest matching subgraph found.
Params: None.
Returns: `int` — The size of the largest matching subgraph.
Side effects: None.

#### `GraphMatch::GetStopSignal() -> bool`
Returns the stop signal status, indicating whether the matching process has been interrupted.
Params: None.
Returns: `bool` — The current stop signal status.
Side effects: None.

#### `GraphMatch::PrintMatches() -> void`
Prints the details of the found matching subgraphs to the console.
Params: None.
Returns: `void`.
Side effects: Performs output to the console.

#### `GraphMatch::SetFlag(word flag) -> void`
Sets a specified flag for internal use by the matching algorithm.
Params: `flag` — The flag to set.
Returns: `void`.
Side effects: Modifies the `flags` member variable.

#### `GraphMatch::RemoveFlag(word flag) -> void`
Removes a specified flag from the internal flags.
Params: `flag` — The flag to remove.
Returns: `void`.
Side effects: Modifies the `flags` member variable.

#### `GraphMatch::SetMaxNofMatches(int maxNofMatches, bool stopOnMaxN) -> void`
Sets the maximum number of matching subgraphs to find and whether to stop on reaching that limit.
Params: `maxNofMatches` — The maximum number of matches to find.
`stopOnMaxN` — If true, stops the search when `maxNofMatches` matches are found.
Returns: `void`.
Side effects: Updates the `maxNofMatches` and `stopOnMaxNMathches` member variables.

#### `GraphMatch::SetMaxNofMatches(int maxNofMatches, bool stopOnMaxN) -> void`
Sets the maximum number of matching subgraphs to find and whether to stop on reaching that limit.
Params: `maxNofMatches` — The maximum number of matches to find.
`stopOnMaxN` — If true, stops the search when `maxNofMatches` matches are found.
Returns: `void`.
Side effects: Updates the `maxNofMatches` and `stopOnMaxNMathches` member variables.

#### `GraphMatch::MatchGraphs(PGraph Gh1, PGraph Gh2, int minMatch,  bool vertexType=true, VERTEX_EXT_TYPE vertexExt=EXTTYPE_Ignore ) -> void`
Finds maximal common subgraphs of a minimum size between two graphs, considering vertex types.
Params: `Gh1` — Pointer to the first graph.
`Gh2` — Pointer to the second graph.
`minMatch` — The minimum size of the common subgraph to find.
`vertexType` — Whether vertex types should be taken into account.
`vertexExt` — The type of extended vertex type to use.
Returns: `void`.
Side effects: Internal call to the matching algorithm, potentially modifying internal data structures.

#### `GraphMatch::InitGraphMatch() -> void`
Initializes the internal data structures for the GraphMatch object.
Params: None.
Returns: `void`.
Side effects:  Initializes member variables, sets initial states.

#### `GraphMatch::FreeMemory() -> void`
Frees the allocated memory associated with the GraphMatch object.
Params: None.
Returns: `void`.
Side effects: Releases memory, including those of internal structures.

#### `GraphMatch::FreeRecHeap() -> void`
Frees the memory allocated for the recursive heap.
Params: None.
Returns: `void`.
Side effects: Releases memory used by the recursive matching heap.

#### `GraphMatch::GetMemory() -> void`
Allocates memory needed by the GraphMatch object.
Params: None.
Returns: `void`.
Side effects: Allocates memory for internal data structures.

#### `GraphMatch::GetRecHeap() -> void`
Allocates memory for the recursive heap used during the matching process.
Params: None.
Returns: `void`.
Side effects: Allocates memory for the matching heap.

#### `GraphMatch::Initialize(bool vertexType, int vertexExt) -> int`
Initializes the GraphMatch object, including setting default values and performing initial setup.
Params: `vertexType` — Whether to consider vertex types.
`vertexExt` — The type of extended vertex type.
Returns: `int` — The number of matches found.
Side effects: Initializes internal variables, sets up the matching process.

#### `GraphMatch::Backtrack(int i) -> void`
Performs backtracking in the graph matching algorithm.
Params: `i` — The index of the node to backtrack from.
Returns: `void`.
Side effects: Modifies the internal matching state during backtracking.

#### `GraphMatch::Backtrack1(int i, int k0) -> void`
Performs a more detailed backtracking step in the graph matching algorithm.
Params: `i` — The index of the node to backtrack from.
`k0` — A value used in the backtracking process.
Returns: `void`.
Side effects: Modifies the internal matching state during backtracking.

#### `GraphMatch::CollectMatch(int nm) -> void`
Collects the vertices of a matching subgraph with a specific number of matches.
Params: `nm` — The number of matches to collect.
Returns: `void`.
Side effects: Populates the match data structures.

#### `GraphMatch::mem_read(cpstr S, int & l) -> void`
Reads data from a character pointer into the GraphMatch object.
Params: `S` — Pointer to the data to read.
`l` — Length of the data to read.
Returns: `void`.
Side effects: Reads data into the internal data structures.

#### `GraphMatch::mem_write(pstr  S, int & l) -> void`
Writes data from the GraphMatch object to a character pointer.
Params: `S` — Pointer to the buffer to write to.
`l` — Length of the data to write.
Returns: `void`.
Side effects: Writes data from the internal data structures to a character pointer.

#### `GraphMatch::Reset() -> void`
Resets the GraphMatch object to its initial state.
Params: None.
Returns: `void`.
Side effects: Clears internal data structures and resets state variables.


---

## `Helix`
Represents a helix structure in a molecular data file.

#### `Helix::Helix()`
Constructs a `Helix` object with default values.
Params: None
Returns: `void`
Side effects: Initializes member variables.
Errors: None

#### `Helix::Helix(cpstr S)`
Constructs a `Helix` object from a C-style string.
Params: `S` — The C-style string representation of the helix.
Returns: `void`
Side effects: Parses the string and initializes member variables.
Errors: `ConvertPDBASCII` must succeed for parsing to be valid.

#### `Helix::Helix(io::RPStream Object)`
Constructs a `Helix` object from an `io::RPStream`.
Params: `Object` — The input stream.
Returns: `void`
Side effects: Reads data from the stream and initializes member variables.
Errors: Dependent on the stream's read operation.

#### `Helix::~Helix()`
Destructs a `Helix` object, releasing any resources.
Params: None
Returns: `void`
Side effects: Releases allocated memory and resets member variables to default values.
Errors: None

#### `Helix::PDBASCIIDump(pstr S, int N)`
Writes the helix data to a PDB ASCII string.
Params: `S` — The output PDB ASCII string. `N` — The maximum length of the string.
Returns: `void`
Side effects: Formats and writes data to the string `S`.
Errors: `N` must be large enough to hold the resulting string.

#### `Helix::MakeCIF(mmcif::PData CIF, int N)`
Generates a CIF (Crystallographic Information File) data structure.
Params: `CIF` — The CIF data structure to populate. `N` — The maximum length of the string.
Returns: `void`
Side effects: Populates the CIF data structure with helix information.
Errors: `N` must be large enough to hold the resulting string.

#### `Helix::ConvertPDBASCII(cpstr S)`
Converts a PDB ASCII string to a `Helix` object.
Params: `S` — The PDB ASCII string.
Returns: `ERROR_CODE` — The error code, or `NoError` on success.
Side effects: Parses the string and initializes member variables.
Errors: Returns `ErrorPDBASCII` on invalid PDB ASCII format.

#### `Helix::GetCIF(mmcif::PData CIF, int & n)`
Generates a CIF (Crystallographic Information File) data structure.
Params: `CIF` — The CIF data structure to populate. `n` — The number of residues in the helix.
Returns: `ERROR_CODE` — The error code, or `NoError` on success.
Side effects: Populates the CIF data structure with helix information.
Errors: Returns `ErrorCIF` on issues writing to CIF.

#### `Helix::GetClassID() -> CLASS_ID`
Returns the class ID for the `Helix` class.
Params: None
Returns: `CLASS_ID` — `ClassID_Helix`
Side effects: None
Errors: None

#### `Helix::Copy(PContainerClass Helix)`
Copies the data from one `Helix` object to another.
Params: `Helix` — The source `Helix` object.
Returns: `void`
Side effects: Copies member variables to the destination object.
Errors: None

#### `Helix::write(io::RFile f)`
Writes the `Helix` data to an `io::RFile`.
Params: `f` — The output file.
Returns: `void`
Side effects: Writes the `Helix` data to the file.
Errors: Dependent on the file's write operation.

#### `Helix::read(io::RFile f)`
Reads the `Helix` data from an `io::RFile`.
Params: `f` — The input file.
Returns: `void`
Side effects: Reads the `Helix` data from the file and initializes member variables.
Errors: Dependent on the file's read operation.

#### `Helix::InitHelix()`
Initializes the `Helix` object with default values.
Params: None
Returns: `void`
Side effects: Performs initialization of member variables.
Errors: None


---

## `Journal`
The `Journal` class represents a PDBASCII entry within an MMDB file.

#### `Journal::Journal()`
Constructs a new `Journal` object with no initial data.
Params: None
Returns: None
Side effects: Initializes the `Journal` object.

#### `Journal::Journal(cpstr S)`
Constructs a new `Journal` object from a PDBASCII string.
Params: `cpstr S` — The PDBASCII string to parse.
Returns: None
Side effects: Parses the PDBASCII string and populates the `Journal` object.

#### `Journal::Journal(io::RPStream Object)`
Constructs a new `Journal` object from an `RPStream` object.
Params: `io::RPStream Object` — The `RPStream` object containing the PDBASCII data.
Returns: None
Side effects: Parses the PDBASCII data from the `RPStream` object and populates the `Journal` object.

#### `Journal::~Journal()`
Destroys the `Journal` object, releasing any associated resources.
Params: None
Returns: None
Side effects: Releases memory and performs necessary cleanup.

#### `Journal::PDBASCIIDump(pstr S, int N)`
Writes the PDBASCII representation of the journal to a given string.
Params: `pstr S` — The string to write the PDBASCII data to, `N` — The maximum number of characters to write.
Returns: None
Side effects: Writes the journal's PDBASCII data to the string `S`.

#### `Journal::ConvertPDBASCII(cpstr S)`
Converts a PDBASCII string to a `Journal` object.
Params: `cpstr S` — The PDBASCII string to convert.
Returns: `ERROR_CODE` — The error code if conversion fails, otherwise 0.
Side effects: Parses the PDBASCII string and populates the `Journal` object.

#### `Journal::GetClassID()`
Returns the class ID for the `Journal` object.
Params: None
Returns: `CLASS_ID` — `ClassID_Journal`.
Side effects: None

#### `Journal::write(io::RFile f)`
Writes the journal's data to an `RFile` object.
Params: `io::RFile f` — The `RFile` object to write to.
Returns: None
Side effects: Serializes the journal data to the `RFile`.

#### `Journal::read(io::RFile f)`
Reads the journal's data from an `RFile` object.
Params: `io::RFile f` — The `RFile` object to read from.
Returns: None
Side effects: Deserializes the journal data from the `RFile`.


---

## `Link`
Represents a single bond or link between two atoms in a molecular structure.

#### `Link::Copy(PContainerClass link)`  // *** USED IN COOT ***
Copies the contents of another `Link` object into the current object.
Params: `link` — Another `Link` object to copy from.
Returns: `void`.
Side effects: The current `Link` object is updated with the values from `link`.
Errors: None.

#### `Link::operator=(Link other)`  // *** USED IN COOT ***
Assigns the contents of another `Link` object to the current object.
Params: `other` — Another `Link` object to copy from.
Returns: `*this`.
Side effects: The current `Link` object is updated with the values from `other`.
Errors: None.

#### `Link::ConvertPDBASCII(cpstr S)`  // *** USED IN COOT ***
Converts a PDB ASCII string into a `Link` object.
Params: `S` — The PDB ASCII string to convert.
Returns: `ERROR_CODE` — The error code if the conversion fails.
Side effects: The `Link` object is populated with data from the PDB ASCII string.
Errors: Returns an error code if the input string is invalid.

#### `Link::GetCIF(mmcif::PData CIF, int & n)`  // *** USED IN COOT ***
Generates a CIF data structure representing the `Link`.
Params: `CIF` — The CIF data structure to populate. `n` — The number of atoms in the CIF.
Returns: `ERROR_CODE` — The error code if the generation fails.
Side effects: The `CIF` data structure is populated.

#### `Link::GetClassID()`
Returns the class ID for `Link`.
Params: None
Returns: `CLASS_ID` — `ClassID_Link`.
Side effects: None.
Errors: None.

#### `Link::MakeCIF(mmcif::PData CIF, int N)`  // *** USED IN COOT ***
Generates a CIF data structure representing the `Link`.
Params: `CIF` — The CIF data structure to populate. `N` — The number of atoms in the CIF.
Returns: `void`.
Side effects: The `CIF` data structure is populated.
Errors: None.

#### `Link::PDBASCIIDump(pstr S, int N)`  // *** USED IN COOT ***
Generates a PDB ASCII string representation of the `Link`.
Params: `S` — The PDB ASCII string to populate. `N` — The length of the string.
Returns: `void`.
Side effects: The `S` string is populated with the PDB ASCII representation.
Errors: None.

#### `Link::read(io::RFile f)`
Reads a `Link` object from an input stream.
Params: `f` — An input file stream.
Returns: `void`.
Side effects: The `Link` object is populated with data from the input stream.
Errors: None.

#### `Link::write(io::RFile f)`
Writes a `Link` object to an output stream.
Params: `f` — An output file stream.
Returns: `void`.
Side effects: The `Link` object is serialized to the output stream.
Errors: None.

#### `Link::InitLink()`
Initializes the internal state of the `Link` object.
Params: None
Returns: `void`.
Side effects: The internal state of the `Link` object is initialized.
Errors: None.


---

## `LinkContainer`
This class represents a container for linked structural elements.

#### `LinkContainer::LinkContainer()`
Constructs a `LinkContainer` object.
Params: None
Returns: `void`
Side effects: Initializes the `ClassContainer` base class.

#### `LinkContainer::LinkContainer(io::RPStream Object)`
Constructs a `LinkContainer` object from an `RPStream` object.
Params: `io::RPStream Object` — The input stream to read the container data from.
Returns: `void`
Side effects: Initializes the `ClassContainer` base class, reading data from the stream.

#### `LinkContainer::~LinkContainer()`
Destroys the `LinkContainer` object.
Params: None
Returns: `void`
Side effects: Destroys the `ClassContainer` base class.


---

## `LinkR`
Represents a single inter-chain or inter-protein link in a structural biology file.

#### `LinkR::LinkR()`
Constructs a `LinkR` object with default values.
Params: None
Returns: `void`
Side effects: Initializes member variables to their default values.
Errors: None

#### `LinkR::LinkR(cpstr S)`
Constructs a `LinkR` object from a C-style string representing the link data.
Params: `cpstr S` — The input string.
Returns: `void`
Side effects: Parses the input string and initializes the `LinkR` object.
Errors: None

#### `LinkR::LinkR(io::RPStream Object)`
Constructs a `LinkR` object from an `io::RPStream` object.
Params: `io::RPStream Object` — The input stream.
Returns: `void`
Side effects: Reads the link data from the stream and initializes the `LinkR` object.
Errors: None

#### `LinkR::~LinkR()`
Destructs the `LinkR` object, releasing any resources.
Params: None
Returns: `void`
Side effects: Cleans up the `LinkR` object by deallocating any memory.
Errors: None

#### `LinkR::PDBASCIIDump(pstr S, int N)`
Formats the `LinkR` data and writes it to a C-style string.
Params: `pstr S` — The output string buffer. `int N` — The buffer size.
Returns: `void`
Side effects: Writes the link data to the specified string buffer.
Errors: None

#### `LinkR::MakeCIF(mmcif::PData CIF, int N)`
Creates a mmCIF data object representing the `LinkR` information.
Params: `mmcif::PData CIF` — The mmCIF data object. `int N` — The size of the CIF data.
Returns: `void`
Side effects: Populates the mmCIF data object with the `LinkR` data.
Errors: None

#### `LinkR::ConvertPDBASCII(cpstr S)`
Converts a PDB ASCII string into a `LinkR` object.
Params: `cpstr S` — The input PDB ASCII string.
Returns: `ERROR_CODE` —  Indicates success or failure.
Side effects: Parses the input string and initializes the `LinkR` object.
Errors: Returns an error code if parsing fails.

#### `LinkR::GetCIF(mmcif::PData CIF, int & n)`
Retrieves the `LinkR` data and populates an mmCIF data object.
Params: `mmcif::PData CIF` — The mmCIF data object. `int & n` — The size of the CIF data.
Returns: `ERROR_CODE` — Indicates success or failure.
Side effects: Populates the mmCIF data object with the `LinkR` data and sets the size of the CIF data.
Errors: Returns an error code if the CIF object is not properly populated.

#### `LinkR::GetClassID() -> CLASS_ID`
Returns the class ID for the `LinkR` object.
Params: None
Returns: `CLASS_ID` —  `ClassID_LinkR`.
Side effects: None
Errors: None

#### `LinkR::Copy(PContainerClass LinkR)`
Creates a copy of the `LinkR` object.
Params: `PContainerClass LinkR` — The `LinkR` object to copy.
Returns: `void`
Side effects: Duplicates the member variables of the source `LinkR` object.
Errors: None

#### `LinkR::read(io::RFile f)`
Reads `LinkR` data from an input file stream.
Params: `io::RFile f` — The input file stream.
Returns: `void`
Side effects: Reads the data from the input stream and initializes the `LinkR` object.
Errors: None

#### `LinkR::write(io::RFile f)`
Writes `LinkR` data to an output file stream.
Params: `io::RFile f` — The output file stream.
Returns: `void`
Side effects: Writes the data to the output stream.
Errors: None

#### `LinkR::InitLinkR()`
Initializes the member variables of the `LinkR` object.
Params: None
Returns: `void`
Side effects: Sets the member variables to their default values.
Errors: None


---

## `Loop`
Represents a loop in a mmCIF structure, containing data fields and tags.

#### `Loop::Loop()`
Basic constructor that initializes the loop.
Params: None
Returns: None
Side effects: Initializes the internal state of the loop.
Errors: None

#### `Loop::Loop(cpstr N)`
Constructor that assigns the loop name.
Params: `N` — The loop name as a C-style string.
Returns: None
Side effects: Initializes the internal state of the loop with the given name.
Errors: None

#### `Loop::Loop(io::RPStream Object)`
Constructor for MMDB data streaming functions.
Params: `Object` — An input/output stream object.
Returns: None
Side effects: Initializes the internal state of the loop using the stream object.
Errors: None

#### `Loop::~Loop()`
Destructor that releases all resources allocated by the loop.
Params: None
Returns: None
Side effects: Frees memory and resets internal state.
Errors: None

#### `Loop::AddLoopTag(cpstr T, bool Remove=true)`
Adds a tag to the loop.
Params: `T` — The tag name as a C-style string. `Remove` — Flag to remove all fields in the loop.
Returns: None
Side effects: Appends the tag to the list of tags;  order cannot be changed.
Errors: None

#### `Loop::AddString(cpstr S, bool NonBlankOnly=false)`
Sets string value at current loop position.
Params: `S` — Character string with value to be set. `NonBlankOnly` — Flag to treat white-space-only strings.
Returns: None
Side effects: Sets the string value at the current loop position, handling "data not given" or "data not available" cases.
Errors: None

#### `Loop::AddNoData(int NoDataType)`
Sets \"data not given\" or \"data not available\" at current loop position.
Params: `NoDataType` — Integer key specifying which type of data absence should be set as a value.
Returns: None
Side effects: Sets the appropriate "data not" value at the current loop position.
Errors: None

#### `Loop::AddReal(realtype R, int prec=8)`
Sets float-point value at current loop position.
Params: `R` — Real number with value to be set. `prec` — Float-point precision.
Returns: None
Side effects: Sets the float-point value at the current loop position, using the specified precision.
Errors: None

#### `Loop::AddReal(realtype R, cpstr format)`
Sets float-point value at current loop position in given format.
Params: `R` — Real number with value to be set. `format` — Format string to convert R.
Returns: None
Side effects: Sets the float-point value at the current loop position, using the specified format.
Errors: None

#### `Loop::AddInteger(int I)`
Sets integer value at current loop position.
Params: `I` — Integer number with value to be set.
Returns: None
Side effects: Sets the integer value at the current loop position.
Errors: None

#### `Loop::GetLoopLength()`
Returns the current length of the loop.
Params: None
Returns: `int` — The number of rows in the loop.
Side effects: Returns the length of the loop.
Errors: None

#### `Loop::GetField(int rowNo, int tagNo)`
Returns the value at the specified tag and row.
Params: `rowNo` — Row number. `tagNo` — Tag number.
Returns: `pstr` — Pointer to string value, or `NULL` if tag or row do not exist.
Side effects: Retrieves the value from the loop based on the row and tag numbers.
Errors: None

#### `Loop::GetString(pstr & S, cpstr TName, int nrow, bool Remove=false)`
Fetches string value from the loop.
Params: `S` — Pointer to string to store the value. `TName` — Tag name. `nrow` — Row number. `Remove` — Flag to remove the field.
Returns: `int` — Return code.
Side effects: Retrieves the string value from the loop, handling "data not given" or "data not available" cases.
Errors: None

#### `Loop::CopyString(pstr buf, int maxlength, cpstr TName, int nrow, int & RC)`
Copies value for given tag, row, into a buffer.
Params: `buf` — Buffer to copy data to. `maxlength` — Maximum number of bytes to copy. `TName` — Tag name. `nrow` — Row number. `RC` — Return code.
Returns: None
Side effects: Copies string data from the loop into the specified buffer.
Errors: None

#### `Loop::DeleteField(cpstr TName, int nrow)`
Deletes the field associated with the given tag and row.
Params: `TName` — Tag name. `nrow` — Row number.
Returns: `int` — Return code.
Side effects: Deletes the field from the loop.
Errors: None

#### `Loop::DeleteRow(int nrow)`
Deletes the row corresponding to the given row number.
Params: `nrow` — Row number.
Returns: `int` — Return code.
Side effects: Deletes the specified row from the loop.
Errors: None

#### `Loop::Optimize()`
Optimizes loop for RAM and data access speed.
Params: None
Returns: None
Side effects: Sorts loop structure to reduce memory use and improve access speeds.
Errors: None

#### `Loop::WriteMMCIFLoop(cpstr FName, io::GZ_MODE gzipMode=io::GZM_CHECK)`
Writes loop data in mmCIF format into file.
Params: `FName` — File name. `gzipMode` — Compression mode.
Returns: `bool` — `true` on success, `false` on failure.
Side effects: Writes the loop data to the specified file in mmCIF format.
Errors: None

#### `Loop::WriteMMCIF(io::RFile f)`
Writes loop data into given file.
Params: `f` — MMDB's file class.
Returns: None
Side effects: Writes loop data into the specified file.
Errors: None

#### `Loop::Copy(PCategory Loop)`
Deep copy of loops.
Params: `Loop` — A pointer to mmcif::Loop, whose content is copied into 'this' loop.
Returns: None
Side effects: Performs a deep copy of the loop data.
Errors: None

#### `Loop::write(io::RFile f)`
MMDB stream writer.
Params: `f` — MMDB's file class.
Returns: None
Side effects: Writes loop data to the specified file stream.
Errors: None

#### `Loop::read(io::RFile f)`
MMDB stream reader.
Params: `f` — MMDB's file class.
Returns: None
Side effects: Reads loop data from the specified file stream.
Errors: None


---

## `Manager`

This class manages MMDB structural biology data.

#### `Manager::Copy(PManager MMDB, COPY_MASK CopyMask)`
Transfers information between two MMDBs based on the specified copy mask, replacing existing data.
Params: `MMDB` — The source MMDB to copy from. `CopyMask` — A mask specifying which data to copy.
Returns: `void`.
Side effects: The current MMDB is updated with data from the source MMDB.
Errors: None.

#### `Manager::Delete(word DelMask)`
Deletes data from the MMDB according to the specified delete mask.
Params: `DelMask` — A mask specifying which data to delete.
Returns: `void`.
Side effects: The current MMDB is modified by removing the specified data.
Errors: None.

#### `Manager::GetRemarks() -> PTitleContainer`
Returns a container of titles, providing remarks associated with the MMDB.
Returns: `PTitleContainer` — A pointer to a title container containing remarks.
Side effects: None.
Errors: None.


---

## `Model`
Represents a structural biology model loaded from an MMDB2 file.

#### `Model::AddChain(PChain chn) -> int`
Adds a chain to the model.
Params: `chn` — Pointer to the chain to be added.
Returns: The serial number of the newly added chain.
Side effects: Adds the chain to the model's chain table.

#### `Model::AddCisPep(cpstr cisPepID) -> int`
Adds a cis peptide to the model.
Params: `cisPepID` —  The identifier of the cis peptide.
Returns: The serial number of the newly added cis peptide.
Side effects: Adds the cis peptide to the model's cispep table.

#### `Model::AddLink(PLink link) -> int`
Adds a link (a covalent bond) to the model.
Params: `link` — Pointer to the link to be added.
Returns: The serial number of the newly added link.
Side effects: Adds the link to the model's link table.

#### `Model::CalcSecStructure(bool flagBulge=true, int aminoSelHnd=-1) -> int`
Calculates secondary structure elements (helices, strands, sheets) for the model.
Params: `flagBulge` —  Boolean indicating whether to calculate bulge elements.
`aminoSelHnd` —  Handle to amino acid selection (unused).
Returns: The number of helices and sheets calculated.
Side effects: Updates the model's secondary structure information.

#### `Model::Copy(PModel model) -> void`
Copies the entire model structure into another `Model` object.
Params: `model` — Pointer to the destination `Model` object.
Returns: `void`.
Side effects: Duplicates all chains, residues, atoms, links, and other data in the model.

#### `Model::DeleteChain(const ChainID chID) -> int`
Deletes a chain from the model.
Params: `chID` — The identifier of the chain to be deleted.
Returns: The serial number of the deleted chain.
Side effects: Removes the chain from the model's chain table.

#### `Model::GetChain(const ChainID chID) -> PChain`
Retrieves a chain from the model based on its identifier.
Params: `chID` — The identifier of the chain.
Returns: A pointer to the chain if found, otherwise `NULL`.
Side effects: None.

#### `Model::GetCisPep(int CisPepNum) -> cpstr`
Retrieves a cis peptide from the model based on its serial number.
Params: `CisPepNum` — The serial number of the cis peptide.
Returns: The identifier of the cis peptide.
Side effects: None.

#### `Model::GetCoordHierarchy() -> PManager`
Returns the `Manager` object, which is used for coordinate access.
Params: None
Returns: A pointer to the `Manager` class.
Side effects: None.

#### `Model::GetChainCreate(const ChainID chID, bool enforceUniqueChainID) -> PChain`
Creates a new chain if one with the given identifier does not exist.
Params: `chID` — The identifier of the chain to create.
`enforceUniqueChainID` — Boolean controlling uniqueness of chain IDs.
Returns: A pointer to the newly created chain.
Side effects: Creates the chain and adds it to the model's chain table.

#### `Model::CreateChain(const ChainID chID) -> PChain`
Creates a new chain with the given identifier, regardless of existing chains.
Params: `chID` — The identifier of the chain to create.
Returns: A pointer to the newly created chain.
Side effects: Creates the chain and adds it to the model's chain table.

#### `Model::GetEntryID() -> cpstr`
Returns the entry ID for the model.
Params: None
Returns: A C-style string representing the entry ID.
Side effects: None.

#### `Model::SetEntryID(const IDCode idCode) -> void`
Sets the entry ID for the model.
Params: `idCode` — The IDCode to set.
Returns: `void`.
Side effects: Updates the model's entry ID.

#### `Model::GetSerNum() -> int`
Returns the model's serial number.
Params: None
Returns: The serial number of the model.
Side effects: None.

#### `Model::GetModelID(pstr modelID) -> cpstr`
Returns the model's ID.
Params: `modelID` — A C-style string to store the model ID.
Returns: A C-style string representing the model ID ("/mdl").
Side effects: None.

#### `Model::GetNumberOfModels() -> int`
Returns the total number of models in the MMDB2 file.
Params: None
Returns: The number of models.
Side effects: None.

#### `Model::GetNumberOfAtoms(bool countTers) -> int`
Returns the number of atoms in the model.
Params: `countTers` — Boolean indicating whether to count terminal atoms.
Returns: The number of atoms.
Side effects: None.

#### `Model::GetNumberOfResidues() -> int`
Returns the number of residues in the model.
Params: None
Returns: The number of residues.
Side effects: None.

#### `Model::GetNumberOfChains() -> int`
Returns the number of chains in the model.
Params: None
Returns: The number of chains.
Side effects: None.

#### `Model::GetNewChainID(ChainID chID, int length=1) -> bool`
Generates a new unique chain ID if one does not exist.
Params: `chID` — The requested chain ID.
`length` — The desired length of the new chain ID (default 1).
Returns: `true` if a new ID was generated, `false` otherwise.
Side effects: Updates the model's chain table with a unique ID.

#### `Model::GetChain(const ChainID chID) -> PChain`
Returns a pointer to the chain identified by chID.
Params: `chID` — The chain identifier.
Returns: Pointer to the chain if found, NULL if not.
Side effects: None.

#### `Model::GetChain(int chainNo) -> PChain`
Returns the chain at the specified index.
Params: `chainNo` — The index of the chain (0-based).
Returns: Pointer to the chain.
Side effects: None.

#### `Model::GetChainTable(PPChain & chainTable, int & NumberOfChains)`
Returns a pointer to the chain table, and the number of chains.
Params: `chainTable` — Pointer to an array of `PChain` objects.
`NumberOfChains` —  Integer to hold the number of chains.
Returns: `void`.
Side effects: Populates `chainTable` with pointers to the chains and sets `NumberOfChains`.

#### `Model::DeleteChain(const ChainID chID) -> int`
Deletes a chain from the model.
Params: `chID` — The identifier of the chain to be deleted.
Returns: The serial number of the deleted chain.
Side effects: Removes the chain from the model's chain table.

#### `Model::DeleteChain(int chainNo) -> int`
Deletes a chain from the model based on its index.
Params: `chainNo` — The index of the chain to be deleted (0-based).
Returns: The serial number of the deleted chain.
Side effects: Removes the chain from the model's chain table.

#### `Model::DeleteAllChains() -> int`
Deletes all chains from the model.
Params: None
Returns: The number of chains that were deleted.
Side effects: Clears the model's chain table.

#### `Model::DeleteSolventChains() -> int`
Deletes all solvent chains from the model.
Params: None
Returns: The number of solvent chains that were deleted.
Side effects: Removes solvent chains from the model's chain table.

#### `Model::TrimChainTable() -> void`
Reduces the size of the chain table to match the actual number of chains.
Params: None
Returns: `void`.
Side effects: Deallocates unused space in the chain table.

#### `Model::SortChains(int sortKey) -> void`
Sorts the chains in the model based on the specified key.
Params: `sortKey` — The sorting key (e.g., SORT_CHAIN_NAME).
Returns: `void`.
Side effects: Reorders the chains in the model's chain table.

#### `Model::GetNumberOfResidues(const ChainID chainID) -> int`
Returns the number of residues in a specific chain.
Params: `chainID` — The identifier of the chain.
Returns: The number of residues in the chain.
Side effects: None.

#### `Model::GetNumberOfResidues(int chainNo) -> int`
Returns the number of residues in a specific chain based on its index.
Params: `chainNo` — The index of the chain (0-based).
Returns: The number of residues in the chain.
Side effects: None.

#### `Model::GetResidue(const ChainID chainID, int seqNo, const InsCode insCode) -> PResidue`
Retrieves a residue from the model based on its chain identifier, sequence number, and insertion code.
Params: `chainID` — The identifier of the chain.
`seqNo` — The sequence number of the residue.
`insCode` — The insertion code of the residue.
Returns: A pointer to the residue.
Side effects: None.

#### `Model::GetResidue(const ChainID chainID, int resNo) -> PResidue`
Retrieves a residue from the model based on its chain identifier and residue number.
Params: `chainID` — The identifier of the chain.
`resNo` — The residue number.
Returns: A pointer to the residue.
Side effects: None.

#### `Model::GetResidue(int chainNo, int seqNo, const InsCode insCode) -> PResidue`
Retrieves a residue from the model based on its chain index, sequence number, and insertion code.
Params: `chainNo` — The index of the chain (0-based).
`seqNo` — The sequence number of the residue.
`insCode` — The insertion code of the residue.
Returns: A pointer to the residue.
Side effects: None.

#### `Model::GetResidue(int chainNo, int resNo) -> PResidue`
Retrieves a residue from the model based on its chain index and residue number.
Params: `chainNo` — The index of the chain (0-based).
`resNo` — The residue number.
Returns: A pointer to the residue.
Side effects: None.

#### `Model::GetResidueNo(const ChainID chainID, int seqNo, const InsCode insCode) -> int`
Returns the residue number for a residue given its chain ID, sequence number, and insertion code.
Params: `chainID` — The identifier of the chain.
`seqNo` — The sequence number of the residue.
`insCode` — The insertion code of the residue.
Returns: The residue number.
Side effects: None.

#### `Model::GetResidueNo(int chainNo, int seqNo, const InsCode insCode) -> int`
Returns the residue number for a residue given its chain index, sequence number, and insertion code.
Params: `chainNo` — The index of the chain (0-based).
`seqNo` — The sequence number of the residue.
`insCode` — The insertion code of the residue.
Returns: The residue number.
Side effects: None.

#### `Model::GetResidueTable(PPChain & chainTable, int & NumberOfChains) -> void`
Returns a pointer to the residue table, and the number of chains.
Params: `chainTable` — Pointer to an array of `PChain` objects.
`NumberOfChains` —  Integer to hold the number of chains.
Returns: `void`.
Side effects: Populates `chainTable` with pointers to the chains and sets `NumberOfChains`.

#### `Model::GetResidueTable(const ChainID chainID, PPChain & chainTable, int & NumberOfChains) -> void`
Returns a pointer to the residue table, and the number of chains.
Params: `chainID` — The identifier of the chain.
`chainTable` — Pointer to an array of `PChain` objects.
`NumberOfChains` —  Integer to hold the number of chains.
Returns: `void`.
Side effects: Populates `chainTable` with pointers to the chains and sets `NumberOfChains`.

#### `Model::GetResidueTable(int chainNo, PPChain & chainTable, int & NumberOfChains) -> void`
Returns a pointer to the residue table, and the number of chains.
Params: `chainNo` — The index of the chain (0-based).
`chainTable` — Pointer to an array of `PChain` objects.
`NumberOfChains` —  Integer to hold the number of chains.
Returns: `void`.
Side effects: Populates `chainTable` with pointers to the chains and sets `NumberOfChains`.

#### `Model::DeleteAtom(const ChainID chID, int seqNo, const InsCode insCode, const AtomName aname, const Element elmnt, const AltLoc aloc) -> int`
Deletes an atom from the model.
Params: `chID` — The identifier of the chain.
`seqNo` — The sequence number of the residue.
`insCode` — The insertion code of the residue.
`aname` — The name of the atom.
`elmnt` — The element of the atom.
`aloc` — The alternate location of the atom.
Returns: The serial number of the deleted atom.
Side effects: Removes the atom from the model's atom table.

#### `Model::DeleteAtom(const ChainID chID, int seqNo, const InsCode insCode, int atomNo) -> int`
Deletes an atom from the model based on its chain identifier, sequence number, insertion code, and atom number.
Params: `chID` — The identifier of the chain.
`seqNo` — The sequence number of the residue.
`insCode` — The insertion code of the residue.
`atomNo` — The number of the atom.
Returns: The serial number of the deleted atom.
Side effects: Removes the atom from the model's atom table.

#### `Model::DeleteAtom(const ChainID chID, int resNo, const AtomName aname, const Element elmnt, const AltLoc aloc) -> int`
Deletes an atom from the model based on its chain identifier, residue number, atom name, element, and alternate location.
Params: `chID` — The identifier of the chain.
`resNo` — The residue number.
`aname` — The name of the atom.
`elmnt` — The element of the atom.
`aloc` — The alternate location of the atom.
Returns: The serial number of the deleted atom.
Side effects: Removes the atom from the model's atom table.

#### `Model::DeleteAtom(const ChainID chID, int resNo, int atomNo) -> int`
Deletes an atom from the model based on its chain identifier, residue number, and atom number.
Params: `chID` — The identifier of the chain.
`resNo` — The residue number.
`atomNo` — The number of the atom.
Returns: The serial number of the deleted atom.
Side effects: Removes the atom from the model's atom table.

#### `Model::DeleteAllAtoms(const ChainID chID, int seqNo, const InsCode insCode) -> int`
Deletes all atoms from a specific chain.
Params: `chID` — The identifier of the chain.
`seqNo` — The sequence number of the residue.
`insCode` — The insertion code of the residue.
Returns: The number of atoms deleted.
Side effects: Clears the model's atom table for the given chain.

#### `Model::DeleteAllAtoms(const ChainID chID, int resNo) -> int`
Deletes all atoms from a specific chain based on its residue number.
Params: `chID` — The identifier of the chain.
`resNo` — The residue number.
Returns: The number of atoms deleted.
Side effects: Clears the model's atom table for the given chain.

#### `Model::DeleteAllAtoms(const ChainID chID) -> int`
Deletes all atoms from a specific chain based on its identifier.
Params: `chID` — The identifier of the chain.
Returns: The number of atoms deleted.
Side effects: Clears the model's atom table for the given chain.

#### `Model::DeleteAllAtoms() -> int`
Deletes all atoms from the model.
Params: None
Returns: The number of atoms deleted.
Side effects: Clears the model's atom table.

#### `Model::DeleteAltLocs() -> int`
Deletes all alternate locations.
Params: None
Returns: The number of deleted alternative locations.
Side effects: Removes all alternate location information from the model.

#### `Model::AddAtom(const ChainID chID, int seqNo, const InsCode insCode, PAtom atom) -> int`
Adds an atom to the model.
Params: `chID` — The identifier of the chain.
`seqNo` — The sequence number of the residue.
`insCode` — The insertion code of the residue.
`atom` — Pointer to the atom to be added.
Returns: The serial number of the newly added atom.
Side effects: Adds the atom to the model's atom table.

#### `Model::AddResidue(const ChainID chainID, PResidue res) -> int`
Adds a residue to the model.
Params: `chainID` — The identifier of the chain.
`res` — Pointer to the residue to be added.
Returns: The serial number of the newly added residue.
Side effects: Adds the residue to the model's residue table.

#### `Model::AddResidue(int chainNo, PResidue res) -> int`
Adds a residue to the model based on its chain index.
Params: `chainNo` — The index of the chain (0-based).
`res` — Pointer to the residue to be added.
Returns: The serial number of the newly added residue.
Side effects: Adds the residue to the model's residue table.

#### `Model::ConvertPDBString(pstr PDBString) -> ERROR_CODE`
Converts a PDB string to a model representation.
Params: `PDBString` — The PDB string to convert.
Returns: An `ERROR_CODE` indicating the success or failure of the conversion.
Side effects: Updates the model's data based on the PDB string.

#### `Model::PDBASCIIDumpPS(io::RFile f) -> void`
Dumps the primary structure records from a PDB file to a file.
Params: `f` — Pointer to the output file.
Returns: `void`.
Side effects: Writes the PDB records to the file.

#### `Model::PDBASCIIDumpCP(io::RFile f) -> void`
Dumps the cispep records from a PDB file to a file.
Params: `f` — Pointer to the output file.
Returns: `void`.
Side effects: Writes the cispep records to the file.

#### `Model::PDBASCIIDump(io::RFile f) -> void`
Dumps the coordinate (ATOM etc.) records from a PDB file to a file.
Params: `f` — Pointer to the output file.
Returns: `void`.
Side effects: Writes the PDB records to the file.

#### `Model::MakeAtomCIF(mmcif::PData CIF) -> void`
Generates a CIF (Crystallographic Information File) from the model.
Params: `CIF` — Pointer to the CIF data structure.
Returns: `void`.
Side effects: Creates a CIF file containing the model's data.

#### `Model::MakePSCIF(mmcif::PData CIF) -> void`
Generates a PSCIF (Protein Structure Crystallographic File) from the model.
Params: `CIF` — Pointer to the CIF data structure.
Returns: `void`.
Side effects: Creates a PSCIF file containing the model's data.

#### `Model::GetCIF(mmcif::PData CIF) -> ERROR_CODE`
Gets CIF from the model.
Params: `CIF` — Pointer to the CIF data structure.
Returns: An `ERROR_CODE` indicating the success or failure of the conversion.
Side effects: Populates the CIF with the model's data.

#### `Model::MoveChain(PChain & m_chain, PPAtom m_atom, PPAtom atom, int & atom_index, int chain_ext) -> void`
Moves a chain to the top chain array.
Params: `m_chain` — Pointer to the chain to be moved.
`m_atom` — Pointer to the atom array.
`atom` — Pointer to the atom array.
`atom_index` — Index of the first free position in atom.
`chain_ext` — Chain extension number.
Returns: `void`.
Side effects: Moves the chain and its corresponding atoms to the top of the chain table.

#### `Model::GetAIndexRange(int & i1, int & i2)`
Returns the index range of the atoms.
Params: `i1` — Index of the first atom.
`i2` — Index of the last atom.
Returns: `void`.
Side effects: Populates `i1` and `i2`.

#### `Model::MaskAtoms(PMask mask) -> void`
Masks atoms.
Params: `mask` — Pointer to the mask.
Returns: `void`.
Side effects: Updates the model's atom table based on the mask.

#### `Model::MaskResidues(PMask mask) -> void`
Masks residues.
Params: `mask` — Pointer to the mask.
Returns: `void`.
Side effects: Updates the model's residue table based on the mask.

#### `Model::MaskChains(PMask mask) -> void`
Masks chains.
Params: `mask` — Pointer to the mask.
Returns: `void`.
Side effects: Updates the model's chain table based on the mask.

#### `Model::UnmaskAtoms(PMask mask) -> void`
Unmasks atoms.
Params: `mask` — Pointer to the mask.
Returns: `void`.
Side effects: Updates the model's atom table based on the mask.

#### `Model::UnmaskResidues(PMask mask) -> void`
Unmasks residues.
Params: `mask` — Pointer to the mask.
Returns: `void`.
Side effects: Updates the model's residue table based on the mask.

#### `Model::UnmaskChains(PMask mask) -> void`
Unmasks chains.
Params: `mask` — Pointer to the mask.
Returns: `void`.
Side effects: Updates the model's chain table based on the mask.

#### `Model::GetAtomStatistics(RAtomStat AS)`
Returns atom statistics.
Params: `AS` — Pointer to RAtomStat.
Returns: `void`.
Side effects: Populates the RAtomStat struct.

#### `Model::CalAtomStatistics(RAtomStat AS) -> void`
Calculates atom statistics.
Params: `AS` — Pointer to RAtomStat.
Returns: `void`.
Side effects: Updates the RAtomStat struct with calculated values.

#### `Model::DeleteAltLocs() -> int`
Deletes all alternate locations.
Params: None
Returns: The number of deleted alternative locations.
Side effects: Removes all alternate location information from the model.

#### `Model::AddAtom(const ChainID chID, int seqNo, const InsCode insCode, PAtom atom) -> int`
Adds an atom to the model.
Params: `chID` — The identifier of the chain.
`seqNo` — The sequence number of the residue.
`insCode` — The insertion code of the residue.
`atom` — Pointer to the atom to be added.
Returns: The serial number of the newly added atom.
Side effects: Adds the atom to the model's atom table.

#### `Model::AddResidue(int chainNo, PResidue res) -> int`
Adds a residue to the model based on its chain index.
Params: `chainNo` — The index of the chain (0-based).
`res` — Pointer to the residue to be added.
Returns: The serial number of the newly added residue.
Side effects: Adds the residue to the model's residue table.


---

## `Remark`
The `Remark` class represents a remark within a PDB or CIF file, providing methods for dumping, converting, and accessing its associated data.

#### `Remark::Remark()`
Constructs a `Remark` object with default values.
Params: None
Returns: None
Side effects: Initializes the `Remark` object.
Errors: None

#### `Remark::Remark(cpstr S)`
Constructs a `Remark` object from a PDB ASCII string.
Params: `S` — A pointer to a character string representing the PDB ASCII data.
Returns: None
Side effects: Parses the PDB ASCII string and populates the `Remark` object.
Errors: Returns `ERROR_CODE` if parsing fails.

#### `Remark::Remark(io::RPStream Object)`
Constructs a `Remark` object from an `io::RPStream` object.
Params: `Object` — An `io::RPStream` object for reading the remark data.
Returns: None
Side effects: Reads the remark data from the `io::RPStream` object and populates the `Remark` object.
Errors: None

#### `Remark::~Remark()`
Destroys the `Remark` object, releasing any allocated resources.
Params: None
Returns: None
Side effects: Releases resources used by the `Remark` object.
Errors: None

#### `Remark::PDBASCIIDump(pstr S, int N)`
// *** USED IN COOT ***
Dumps the `Remark` data into a PDB ASCII string.
Params: `S` — A pointer to the output PDB ASCII string, `N` — The maximum length of the string.
Returns: None
Side effects: Writes the `Remark` data to the specified string.
Errors: None

#### `Remark::MakeCIF(mmcif::PData CIF, int N)`
// *** USED IN COOT ***
Converts the `Remark` data into a CIF data structure.
Params: `CIF` — A pointer to the CIF data structure, `N` — The maximum length of the CIF.
Returns: None
Side effects: Populates the CIF data structure with the `Remark` data.
Errors: None

#### `Remark::ConvertPDBASCII(cpstr S)`
// *** USED IN COOT ***
Converts a PDB ASCII string to a `Remark` object.
Params: `S` — A pointer to a character string representing the PDB ASCII data.
Returns: `ERROR_CODE` — The error code, or `ERROR_CODE_OK` if conversion succeeds.
Side effects: Parses the PDB ASCII string and populates the `Remark` object.
Errors: Returns `ERROR_CODE_OK` on success, other `ERROR_CODE` values on failure.

#### `Remark::GetCIF(mmcif::PData CIF, int & n)`
// *** USED IN COOT ***
Retrieves the CIF data structure associated with the `Remark`.
Params: `CIF` — A pointer to the CIF data structure, `n` — A reference to an integer that will store the CIF ID.
Returns: `ERROR_CODE` — The error code, or `ERROR_CODE_OK` if retrieval succeeds.
Side effects: Populates the CIF data structure and sets the CIF ID in the `n` parameter.
Errors: Returns `ERROR_CODE_OK` on success, other `ERROR_CODE` values on failure.

#### `Remark::GetClassID() -> ClassID`
Returns the class ID for the `Remark` object.
Params: None
Returns: `ClassID_Remark` — The class ID.
Side effects: None
Errors: None

#### `Remark::Copy(PContainerClass RemarkClass)`
// *** USED IN COOT ***
Copies the contents of one `Remark` object to another.
Params: `RemarkClass` — A pointer to the `Remark` object to copy to.
Returns: None
Side effects: Copies the data from one `Remark` object to another.
Errors: None

#### `Remark::write(io::RFile f)`
// *** USED IN COOT ***
Writes the `Remark` data to an `io::RFile` object.
Params: `f` — An `io::RFile` object for writing the data.
Returns: None
Side effects: Writes the `Remark` data to the specified file.
Errors: None

#### `Remark::read(io::RFile f)`
// *** USED IN COOT ***
Reads the `Remark` data from an `io::RFile` object.
Params: `f` — An `io::RFile` object for reading the data.
Returns: None
Side effects: Reads the `Remark` data from the specified file.
Errors: None


---

## `Residue`
Represents a single residue in a protein structure, storing atom coordinates, sequence information, and other relevant data.

#### `Residue::AddAtom(PAtom atm)`
Adds an atom to the residue's atom table.
Params: `atm` — A pointer to the atom to be added.
Returns: The index of the newly added atom in the atom table.
Side effects: Adds the atom to the `atom` array and updates the `nAtoms` counter.
Errors: None.

#### `Residue::DeleteAtom(int atomNo)`
Deletes an atom from the residue's atom table by its index.
Params: `atomNo` — The index of the atom to delete.
Returns: 1 if the atom was successfully deleted, 0 otherwise.
Side effects: Removes the atom from the `atom` array and updates the `nAtoms` counter.
Errors: Returns 0 if the atom index is invalid.

#### `Residue::GetAtom(const AtomName aname, const Element elname=NULL, const AltLoc aloc=cpstr(""))`
Retrieves an atom from the residue's atom table by its name and optional coordinates.
Params: `aname` — The atom name. `elname` —  The element name (optional). `aloc` — The alternative location (optional).
Returns: A pointer to the atom if found, or `NULL` if not found.
Side effects: None.
Errors: Returns `NULL` if the atom with the given name does not exist.

#### `Residue::GetAtomTable(PPAtom & atomTable, int & NumberOfAtoms)`
Retrieves the residue's atom table.
Params: `atomTable` — A pointer to a dynamically allocated array to store the atom pointers. `NumberOfAtoms` — A reference to an integer to store the number of atoms in the table.
Returns: None.
Side effects: Allocates memory for the atom table and populates it with pointers to the atoms in the residue.
Errors: None.

#### `Residue::GetChain()`
Returns a pointer to the residue's chain.
Params: None.
Returns: A pointer to the `Chain` object.
Side effects: None.
Errors: Returns `NULL` if the residue has no chain assigned.

#### `Residue::GetModel()`
Returns a pointer to the residue's model.
Params: None.
Returns: A pointer to the `Model` object.
Side effects: None.
Errors: Returns `NULL` if the residue has no model assigned.

#### `Residue::GetModelNum()`
Returns the model number of the residue.
Params: None.
Returns: The model number as an integer.
Side effects: None.
Errors: None.

#### `Residue::GetChainID()`
Returns the chain ID of the residue.
Params: None.
Returns: The chain ID as a string.
Side effects: None.
Errors: None.

#### `Residue::GetLabelAsymID()`
Returns the asymmetric ID of the chain.
Params: None.
Returns: The asymmetric ID as a string.
Side effects: None.
Errors: None.

#### `Residue::GetResName()`
Returns the residue name.
Params: None.
Returns: The residue name as a string.
Side effects: None.
Errors: None.

#### `Residue::SetChain(PChain Chain_Owner)`
Sets the chain for the residue.
Params: `Chain_Owner` — A pointer to the `Chain` object to associate with the residue.
Returns: None.
Side effects: Sets the `chain` member to point to the given `Chain` object.
Errors: None.

#### `Residue::SetChainID(const ChainID chID)`
Sets the chain ID of the residue.
Params: `chID` — The chain ID to set.
Returns: None.
Side effects: Sets the `label_asym_id` member to the given `ChainID`.
Errors: None.

#### `Residue::PDBASCIIAtomDump(io::RFile f)`
Writes the atom coordinates in PDB ASCII format to a file.
Params: `f` — A pointer to the output file.
Returns: None.
Side effects: Writes the atom coordinates to the specified file in PDB ASCII format.
Errors: None.

#### `Residue::MakeAtomCIF(mmcif::PData CIF)`
Creates a mmCIF data structure for the residue.
Params: `CIF` — A pointer to the mmCIF data structure to populate.
Returns: None.
Side effects: Populates the `CIF` data structure with the residue's data.
Errors: None.

#### `Residue::GetSeqNum()`
Returns the sequence number of the residue.
Params: None.
Returns: The sequence number as an integer.
Side effects: None.
Errors: None.

#### `Residue::GetLabelSeqID()`
Returns the label sequence ID of the residue.
Params: None.
Returns: The label sequence ID as an integer.
Side effects: None.

#### `Residue::GetLabelEntityID()`
Returns the label entity ID of the residue.
Params: None.
Returns: The label entity ID as an integer.
Side effects: None.

#### `Residue::GetInsCode()`
Returns the insertion code of the residue.
Params: None.
Returns: The insertion code as an integer.
Side effects: None.

#### `Residue::GetResidueNo()`
Returns the residue number.
Params: None.
Returns: The residue number as an integer.
Side effects: None.

#### `Residue::GetCenter(realtype & x, realtype & y, realtype & z)`
Calculates the center coordinates of the residue.
Params: `x`, `y`, `z` — References to realtype variables to store the coordinates.
Returns: None.
Side effects: Calculates the weighted average coordinates of the atoms in the residue and stores them in the specified variables.
Errors: None.

#### `Residue::GetCoordHierarchy()`
Returns the coordinate hierarchy associated with the residue.
Params: None.
Returns: A pointer to the coordinate hierarchy.
Side effects: None.
Errors: Returns `NULL` if no coordinate hierarchy is associated with the residue.

#### `Residue::GetAtomStatistics(RAtomStat AS)`
Populates the `RAtomStat` structure with atom statistics.
Params: `AS` — A pointer to the `RAtomStat` structure to populate.
Returns: None.
Side effects: Populates the `RAtomStat` structure with information about the atoms in the residue.
Errors: None.

#### `Residue::CalAtomStatistics(RAtomStat AS)`
Calculates the atom statistics.
Params: `AS` — A pointer to the `RAtomStat` structure to populate.
Returns: None.
Side effects: Calculates the atom statistics and populates the `RAtomStat` structure.
Errors: None.

#### `Residue::GetResidueID(pstr ResidueID)`
Returns the residue ID.
Params: `ResidueID` - A pointer to a string where the residue ID will be written.
Returns: 1 if the residue ID was successfully created, 0 otherwise.
Side effects: Creates a string representing the residue ID and writes it to the `ResidueID` parameter.
Errors: Returns 0 if the residue ID cannot be created.

#### `Residue::GetAltLocations(int & nAltLocs, PAltLoc & aLoc, rvector & occupancy, int & alflag)`
Retrieves alternative location data for the residue.
Params: `nAltLocs` — A reference to an integer to store the number of alternative locations. `aLoc` — A pointer to the `PAltLoc` structure to store the alternative location data. `occupancy` — A reference to the rvector to store the occupancy values. `alflag` — A reference to an integer to store the alternative location flag.
Returns: None.
Side effects: Populates the specified data structures with information about the alternative locations of the atoms in the residue.
Errors: None.

#### `Residue::GetNofAltLocations()`
Returns the number of alternative locations in the residue.
Params: None.
Returns: The number of alternative locations as an integer.
Side effects: None.
Errors: None.

#### `Residue::isAminoacid()`
Returns true if the residue is an amino acid.
Params: None.
Returns: True if the residue is an amino acid, false otherwise.
Side effects: None.
Errors: None.

#### `Residue::isNucleotide()`
Returns true if the residue is a nucleotide.
Params: None.
Returns: True if the residue is a nucleotide, false otherwise.
Side effects: None.

#### `Residue::isDNARNA()`
Returns the type of nucleotide (DNA, RNA or neither).
Params: None.
Returns: 0 for neither, 1 for DNA, 2 for RNA.
Side effects: None.
Errors: None.

#### `Residue::isSugar()`
Returns true if the residue is a sugar.
Params: None.
Returns: True if the residue is a sugar, false otherwise.
Side effects: None.
Errors: None.

#### `Residue::isSolvent()`
Returns true if the residue is a solvent.
Params: None.
Returns: True if the residue is a solvent, false otherwise.
Side effects: None.
Errors: None.

#### `Residue::isModRes()`
Returns true if the residue is a modified residue.
Params: None.
Returns: True if the residue is a modified residue, false otherwise.
Side effects: None.
Errors: None.

#### `Residue::isInSelection(int selHnd)`
Returns true if the residue is in the selection.
Params: `selHnd` — The handle of the selection.
Returns: True if the residue is in the selection, false otherwise.
Side effects: None.
Errors: None.

#### `Residue::isNTerminus()`
Returns true if the residue is the N-terminus of the chain.
Params: None.
Returns: True if the residue is the N-terminus, false otherwise.
Side effects: None.
Errors: None.

#### `Residue::isCTerminus()`
Returns true if the residue is the C-terminus of the chain.
Params: None.
Returns: True if the residue is the C-terminus, false otherwise.
Side effects: None.

#### `Residue::CheckID(int * snum, const InsCode inscode=pstr(""), const ResName resname=NULL)`
Checks if a residue is identified by the given parameters.
Params: `snum` — A pointer to the sequence number. `inscode` — The insertion code. `resname` — The residue name.
Returns: 1 if the residue is identified, 0 otherwise.
Side effects: Compares the provided parameters with the residue's data.
Errors: Returns 0 if the parameters do not match.

#### `Residue::CheckIDS(cpstr ID)`
Checks if a residue is identified by the given ID string.
Params: `ID` — The residue ID string.
Returns: 1 if the residue is identified, 0 otherwise.
Side effects: Compares the given ID string with the residue's data.
Errors: Returns 0 if the ID does not match.

#### `Residue::AddAtom(PAtom atm)`
Adds an atom to the residue’s atom table.
Params: `atm` — A pointer to the atom to add.
Returns: The index of the newly added atom.
Side effects: Adds the atom to the atom array and updates the atom count.
Errors: None

#### `Residue::DeleteAtom(int atomNo)`
Deletes an atom from the residue’s atom table.
Params: `atomNo` — The index of the atom to delete.
Returns: 1 if successful, 0 otherwise.
Side effects: Removes the atom from the atom array and updates the atom count.
Errors: Returns 0 if the atom index is invalid.

#### `Residue::GetAtom(const AtomName aname, const Element elname=NULL, const AltLoc aloc=cpstr(""))`
Retrieves an atom from the residue’s atom table.
Params: `aname` — The atom name. `elname` — The element name (optional). `aloc` — The alternative location (optional).
Returns: A pointer to the atom if found, or NULL if not found.
Side effects: None.
Errors: Returns NULL if the atom with the given name does not exist.

#### `Residue::GetAtomTable(PPAtom & atomTable, int & NumberOfAtoms)`
Retrieves the residue’s atom table.
Params: `atomTable` — A pointer to a dynamically allocated array to store the atom pointers. `NumberOfAtoms` — A reference to an integer to store the number of atoms in the table.
Returns: None.
Side effects: Allocates memory for the atom table and populates it with pointers to the atoms in the residue.
Errors: None.

#### `Residue::TrimAtomTable()`
Removes empty atom entries from the residue’s atom table.
Params: None.
Returns: None.
Side effects: Reduces the size of the atom table by removing any atoms without coordinates.
Errors: None.

#### `Residue::isAminoacid()`
Checks if the residue is an amino acid.
Params: None.
Returns: True if the residue is an amino acid, false otherwise.
Side effects: None.
Errors: None.

#### `Residue::isNucleotide()`
Checks if the residue is a nucleotide.
Params: None.
Returns: True if the residue is a nucleotide, false otherwise.
Side effects: None.

#### `Residue::isDNARNA()`
Returns the type of nucleotide (DNA, RNA or neither).
Params: None.
Returns: 0 for neither, 1 for DNA, 2 for RNA.
Side effects: None.
Errors: None.

#### `Residue::isSugar()`
Checks if the residue is a sugar.
Params: None.
Returns: True if the residue is a sugar, false otherwise.
Side effects: None.
Errors: None.

#### `Residue::isSolvent()`
Checks if the residue is a solvent.
Params: None.
Returns: True if the residue is a solvent, false otherwise.
Side effects: None.
Errors: None.

#### `Residue::isModRes()`
Checks if the residue is a modified residue.
Params: None.
Returns: True if the residue is a modified residue, false otherwise.
Side effects: None.
Errors: None.

#### `Residue::isInSelection(int selHnd)`
Checks if the residue is present in the specified selection.
Params: `selHnd` — The handle of the selection.
Returns: True if the residue is in the selection, false otherwise.
Side effects: None.
Errors: None.

#### `Residue::isNTerminus()`
Checks if the residue is the N-terminus of the chain.
Params: None.
Returns: True if the residue is the N-terminus, false otherwise.
Side effects: None.
Errors: None.

#### `Residue::isCTerminus()`
Checks if the residue is the C-terminus of the chain.
Params: None.
Returns: True if the residue is the C-terminus, false otherwise.
Side effects: None.

#### `Residue::SetChain(PChain Chain_Owner)`
Sets the chain for the residue.
Params: `Chain_Owner` — A pointer to the `Chain` object to associate with the residue.
Returns: None.
Side effects: Sets the `chain` member to point to the given `Chain` object.
Errors: None.

#### `Residue::SetChainID(const ChainID chID)`
Sets the chain ID of the residue.
Params: `chID` — The chain ID to set.
Returns: None.
Side effects: Sets the `label_asym_id` member to the given `ChainID`.
Errors: None.

#### `Residue::PDBASCIIAtomDump(io::RFile f)`
Writes the atom coordinates in PDB ASCII format to a file.
Params: `f` — A pointer to the output file.
Returns: None.
Side effects: Writes the atom coordinates to the specified file in PDB ASCII format.
Errors: None.

#### `Residue::MakeAtomCIF(mmcif::PData CIF)`
Creates a mmCIF data structure for the residue.
Params: `CIF` — A pointer to the mmCIF data structure to populate.
Returns: None.
Side effects: Populates the `CIF` data structure with the residue's data.
Errors: None.

#### `Residue::GetModel()`
Returns a pointer to the residue's model.
Params: None.
Returns: A pointer to the `Model` object.
Side effects: None.
Errors: Returns NULL if the residue has no model assigned.

#### `Residue::GetModelNum()`
Returns the model number of the residue.
Params: None.
Returns: The model number as an integer.
Side effects: None.
Errors: None.

#### `Residue::GetNumberOfAtoms()`
Returns the number of atoms in the residue.
Params: None.
Returns: The number of atoms in the residue.
Side effects: None.
Errors: None.

#### `Residue::GetNumberOfAtoms(bool countTers)`
Returns the number of atoms in the residue, including terminal atoms.
Params: `countTers` — A boolean flag indicating whether to include terminal atoms.
Returns: The number of atoms in the residue.
Side effects: None.
Errors: None.

#### `Residue::GetAtom(const AtomName aname, const Element elname=NULL, const AltLoc aloc=cpstr(""))`
Retrieves an atom from the residue’s atom table by its name and optional coordinates.
Params: `aname` — The atom name. `elname` — The element name (optional). `aloc` — The alternative location (optional).
Returns: A pointer to the atom if found, or NULL if not found.
Side effects: None.
Errors: Returns NULL if the atom with the given name does not exist.

#### `Residue::GetAtomTable(PPAtom & atomTable, int & NumberOfAtoms)`
Retrieves the residue’s atom table.
Params: `atomTable` — A pointer to a dynamically allocated array to store the atom pointers. `NumberOfAtoms` — A reference to an integer to store the number of atoms in the table.
Returns: None.
Side effects: Allocates memory for the atom table and populates it with pointers to the atoms in the residue.
Errors: None.

#### `Residue::InsertAtom(PAtom atm, int position)`
Inserts an atom into the residue’s atom table at the specified position.
Params: `atm` — A pointer to the atom to insert. `position` — The index at which to insert the atom.
Returns: The index of the newly added atom.
Side effects: Inserts the atom into the atom table.
Errors: None.

#### `Residue::InsertAtom(PAtom atm, const AtomName aname)`
Inserts an atom into the residue’s atom table before the atom with the given name.
Params: `atm` — A pointer to the atom to insert. `aname` — The atom name of the atom before which to insert the new atom.
Returns: The index of the newly added atom.
Side effects: Inserts the atom into the atom table before the atom with the specified name.
Errors: None.

#### `Residue::DeleteAltLocs()`
Deletes all alternative locations from the residue's atom table.
Params: None.
Returns: The number of deleted alternative locations.
Side effects: Removes all alternative locations from the atom table.
Errors: None.

#### `Residue::TrimAtomTable()`
Removes empty atom entries from the residue’s atom table.
Params: None.
Returns: None.
Side effects: Reduces the size of the atom table by removing any atoms without coordinates.
Errors: None.

#### `Residue::ApplyTransform(const mat44 & TMatrix)`
Applies a transformation matrix to the residue’s atoms.
Params: `TMatrix` — A pointer to the transformation matrix.
Returns: None.
Side effects: Transforms the coordinates of all atoms in the residue by multiplying them with the transformation matrix.
Errors: None.

#### `Residue::MaskAtoms(PMask Mask)`
Masks atoms in the residue according to the provided mask.
Params: `Mask` — A pointer to the mask.
Returns: None.
Side effects: Masks the atoms in the residue according to the provided mask.
Errors: None.

#### `Residue::UnmaskAtoms(PMask Mask)`
Unmasks atoms in the residue according to the provided mask.
Params: `Mask` — A pointer to the mask.
Returns: None.
Side effects: Unmasks the atoms in the residue according to the provided mask.
Errors: None.

#### `Residue::PutUDData(int UDDhandle, int iudd)`
Stores user-defined data in the residue’s UDData structure.
Params: `UDDhandle` — The handle of the UDData structure. `iudd` — The value to store.
Returns: None.
Side effects: Stores the user-defined data in the specified UDData structure.
Errors: None.

#### `Residue::PutUDData(int UDDhandle, realtype rudd)`
Stores user-defined data in the residue’s UDData structure.
Params: `UDDhandle` — The handle of the UDData structure. `rudd` — The value to store.
Returns: None.
Side effects: Stores the user-defined data in the specified UDData structure.
Errors: None.

#### `Residue::PutUDData(int UDDhandle, cpstr sudd)`
Stores user-defined data in the residue’s UDData structure.
Params: `UDDhandle` — The handle of the UDData structure. `suddt` — The string value to store.
Returns: None.
Side effects: Stores the user-defined data in the specified UDData structure.
Errors: None.

#### `Residue::GetUDData(int UDDhandle, int & iudd)`
Retrieves user-defined data from the residue’s UDData structure.
Params: `UDDhandle` — The handle of the UDData structure. `iudd` — A reference to an integer to store the retrieved value.
Returns: None.
Side effects: Retrieves the user-defined data from the specified UDData structure and stores it in the specified integer variable.
Errors: None.

#### `Residue::GetUDData(int UDDhandle, realtype & rudd)`
Retrieves user-defined data from the residue’s UDData structure.
Params: `UDDhandle` — The handle of the UDData structure. `rudd` — A reference to a realtype variable to store the retrieved value.
Returns: None.
Side effects: Retrieves the user-defined data from the specified UDData structure and stores it in the specified realtype variable.
Errors: None.

#### `Residue::GetUDData(int UDDhandle, pstr sudd, int maxLen)`
Retrieves user-defined data from the residue’s UDData structure.
Params: `UDDhandle` — The handle of the UDData structure. `suddt` — A pointer to a string to store the retrieved value. `maxLen` — The maximum length of the string to store.
Returns: None.
Side effects: Retrieves the user-defined data from the specified UDData structure and stores it in the specified string variable.
Errors: None.

#### `Residue::GetUDData(int UDDhandle, pstr & sudd)`
Retrieves user-defined data from the residue’s UDData structure.
Params: `UDDhandle` — The handle of the UDData structure. `suddt` — A pointer to a string to store the retrieved value.
Returns: None.
Side effects: Retrieves the user-defined data from the specified UDData structure and stores it in the specified string variable.
Errors: None.

#### `Residue::Copy(PResidue res)`
Copies the residue data.
Params: `res` — A pointer to the residue to copy.
Returns: None.
Side effects: Copies the residue data from the source to the destination.
Errors: None.

#### `Residue::write(io::RFile f)`
Writes the residue data to the output file.
Params: `f` — A pointer to the output file.
Returns: None.
Side effects: Writes the residue data to the specified file.
Errors: None.

#### `Residue::read(io::RFile f)`
Reads the residue data from the input file.
Params: `f` — A pointer to the input file.
Returns: None.
Side effects: Reads the residue data from the specified file.
Errors: None.

#### `Residue::SetResID(const ResName resName, int sqNum, const InsCode ins)`
Sets the residue ID.
Params: `resName` — The residue name. `sqNum` — The sequence number. `ins` — The insertion code.
Returns: None.
Side effects: Sets the residue ID, sequence number and insertion code.
Errors: None.


---

## `Root`
The `Root` class provides a central interface for accessing and manipulating structural biology data in the MMDB2 format.

#### `Root::CrystReady()` -> `bool`
Checks if crystallographic information is complete, noting any imprecise or translation-containing details.
Params: None
Returns: `true` if crystallographic information is complete, `false` otherwise.
Side effects: None.
Errors: None.

#### `Root::DeleteAllModels()` -> `int`
Deletes all models from the structure, resetting the structure to an empty state.
Params: None
Returns: The number of models deleted.
Side effects: Destroys all model objects.
Errors: None.

#### `Root::FinishStructEdit()` -> `int`
Completes structural editing operations, finalizing changes to the structure.
Params: None
Returns: 0 on success, error code on failure.
Side effects: Finalizes the current structure, potentially updating internal data structures.
Errors: None.

#### `Root::GetAtomI(int index)` -> `PAtom`
Returns the atom at the given serial number index in the structure's atom array.
Params: `index` — The index of the atom in the atom array (1-based).
Returns: A pointer to the atom at the specified index.
Side effects: None.
Errors: Returns `nullptr` if the index is out of bounds.

#### `Root::GetCell() -> mat44`
Returns the cell transformation matrix for the current structure.
Params: None
Returns: A 4x4 matrix representing the cell transformation.
Side effects: None.
Errors: None.

#### `Root::GetFileType() -> int`
Returns the file type of the currently loaded structure.
Params: None
Returns: An integer code representing the file type (e.g., 0 for PDB, 1 for CIF).
Side effects: None.
Errors: None.

#### `Root::GetInputBuffer(pstr Line, int & count)` -> `void`
Retrieves a line of input from the input buffer and stores it in the specified string.
Params: `Line` — A pointer to a string where the line will be stored. `count` — A reference to an integer that will be set to the number of characters in the line.
Returns: None
Side effects: Populates `Line` with a line from the input buffer.
Errors: None.

#### `Root::GetNumberOfSymOps() -> int`
Returns the number of symmetry operations defined in the structure.
Params: None
Returns: The number of symmetry operations.
Side effects: None.
Errors: None.

#### `Root::GetSpaceGroup() -> pstr`
Returns the name of the space group for the current structure.
Params: None
Returns: A pointer to a string containing the space group name.
Side effects: None.
Errors: Returns `nullptr` if the space group is not defined.

#### `Root::GetStructureTitle(pstr & L)` -> `void`
Retrieves the structure title and stores it in the specified string.
Params: `L` — A pointer to a string where the title will be stored.
Returns: None
Side effects: Populates `L` with the structure title.
Errors: None.

#### `Root::GetSymOp(int Nop) -> pstr`
Returns the XYZ symmetry operation name corresponding to the given symmetry operation number.
Params: `Nop` — The number of the symmetry operation.
Returns: A pointer to a string containing the name of the symmetry operation.
Side effects: None.
Errors: Returns `nullptr` if the operation number is invalid.

#### `Root::GetTMatrix(mat44 & TMatrix, int Nop, int cellshift_a, int cellshift_b, int cellshift_c) -> int`
Calculates and returns the coordinate transformation matrix for the given symmetry operation number.
Params: `TMatrix` — A 4x4 matrix to store the transformation matrix. `Nop` — The symmetry operation number. `cellshift_a`, `cellshift_b`, `cellshift_c` —  Cell shift parameters.
Returns: 0 on success, error code on failure.
Side effects: Populates `TMatrix` with the transformation matrix.
Errors: Returns a non-zero error code if the calculation fails.

#### `Root::GetUDDHandle(UDR_TYPE udr_type, cpstr UDDataID) -> int`
Retrieves the handle for a user-defined data record.
Params: `udr_type` — The type of user-defined data. `UDDataID` — The ID of the user-defined data record.
Returns: The handle for the user-defined data record.
Side effects: None.
Errors: Returns a negative error code if the record is not found.

#### `Root::PDBCleanup(word CleanKey) -> word`
Performs PDB cleanup operations on the structure, ensuring compliance with PDB standards.
Params: `CleanKey` —  A key used to determine which cleanup operations to perform.
Returns: 0 on success, error code on failure.
Side effects: Modifies the structure's atom and residue data to comply with PDB standards.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::PutAtom(int index, int serNum, const AtomName atomName, const ResName resName, const ChainID chainID, int seqNum, const InsCode insCode, const AltLoc altLoc, const SegID segID, const Element element) -> int`
Adds an atom to the structure at the specified index, setting its properties.
Params: `index` — The index at which to add the atom. `serNum` — The serial number of the atom. `atomName` — The name of the atom. `resName` — The name of the residue. `chainID` — The chain ID of the residue. `seqNum` — The sequence number of the atom. `insCode` — The insertion code. `altLoc` — The alternative location code. `segID` — The segment ID. `element` — The element symbol.
Returns: 0 on success, error code on failure.
Side effects: Adds an atom to the structure, updating related data structures.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::PutPDBString(cpstr PDBString) -> ERROR_CODE`
Adds a PDB-keyworded string to the existing structure.
Params: `PDBString` — The string to add.
Returns: 0 on success, error code on failure.
Side effects: Adds the string as a remark to the current structure.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::ReadCoorFile(cpstr LFName) -> ERROR_CODE`
Reads coordinate data from a file in a coordinate format.
Params: `LFName` — The name of the coordinate file.
Returns: 0 on success, error code on failure.
Side effects: Reads coordinate data from the file and updates the structure's atom coordinates.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::ReadCIFASCII(cpstr CIFFileName) -> ERROR_CODE`
Reads CIF data from an ASCII file.
Params: `CIFFileName` — The name of the CIF file.
Returns: 0 on success, error code on failure.
Side effects: Reads CIF data from the file and updates the structure's data.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::RegisterUDInteger(UDR_TYPE udr_type, cpstr UDDataID) -> int`
Registers a user-defined integer data record.
Params: `udr_type` — The type of user-defined data. `UDDataID` — The ID of the user-defined data record.
Returns: The handle for the user-defined data record.
Side effects: Adds a user-defined integer data record to the structure's data table.
Errors: Returns a negative error code if the record cannot be added.

#### `Root::SetFlag(word Flag)` -> `void`
Sets a flag in the structure.
Params: `Flag` — The flag to set.
Returns: None
Side effects: Modifies the structure's internal flags.
Errors: None.

#### `Root::SetSpaceGroup(cpstr spGroup) -> int`
Sets the space group for the current structure.
Params: `spGroup` — The name of the space group.
Returns: 0 on success, error code on failure.
Side effects: Updates the structure's space group information.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::WriteCIFASCII(cpstr CIFFileName) -> ERROR_CODE`
Writes CIF data to an ASCII file.
Params: `CIFFileName` — The name of the CIF file.
Returns: 0 on success, error code on failure.
Side effects: Writes CIF data to the file.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::WritePDBASCII(cpstr PDBFileName) -> ERROR_CODE`
Writes PDB data to an ASCII file.
Params: `PDBFileName` — The name of the PDB file.
Returns: 0 on success, error code on failure.
Side effects: Writes PDB data to the file.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::DeleteAllModels() -> int`
Deletes all models from the structure.
Params: None
Returns: The number of models deleted.
Side effects: Destroys all model objects.
Errors: None.

#### `Root::FinishStructEdit() -> int`
Completes structural editing operations.
Params: None
Returns: 0 on success, error code on failure.
Side effects: Finalizes the current structure, potentially updating internal data structures.
Errors: None.

#### `Root::GetAtomI(int index) -> PAtom`
Returns the atom at the given serial number index in the structure's atom array.
Params: `index` — The index of the atom in the atom array (1-based).
Returns: A pointer to the atom at the specified index.
Side effects: None.
Errors: Returns `nullptr` if the index is out of bounds.

#### `Root::GetCell() -> mat44`
Returns the cell transformation matrix for the current structure.
Params: None
Returns: A 4x4 matrix representing the cell transformation.
Side effects: None.
Errors: None.

#### `Root::GetFileType() -> int`
Returns the file type of the currently loaded structure.
Params: None
Returns: An integer code representing the file type (e.g., 0 for PDB, 1 for CIF).
Side effects: None.
Errors: None.

#### `Root::GetInputBuffer(pstr Line, int & count) -> void`
Retrieves a line of input from the input buffer and stores it in the specified string.
Params: `Line` — A pointer to a string where the line will be stored. `count` — A reference to an integer that will be set to the number of characters in the line.
Returns: None
Side effects: Populates `Line` with a line from the input buffer.
Errors: None.

#### `Root::GetNumberOfSymOps() -> int`
Returns the number of symmetry operations defined in the structure.
Params: None
Returns: The number of symmetry operations.
Side effects: None.
Errors: None.

#### `Root::GetSpaceGroup() -> pstr`
Returns the name of the space group for the current structure.
Params: None
Returns: A pointer to a string containing the space group name.
Side effects: None.
Errors: Returns `nullptr` if the space group is not defined.

#### `Root::GetStructureTitle(pstr & L) -> void`
Retrieves the structure title and stores it in the specified string.
Params: `L` — A pointer to a string where the title will be stored.
Returns: None
Side effects: Populates `L` with the structure title.
Errors: None.

#### `Root::GetSymOp(int Nop) -> pstr`
Returns the XYZ symmetry operation name corresponding to the given symmetry operation number.
Params: `Nop` — The number of the symmetry operation.
Returns: A pointer to a string containing the name of the symmetry operation.
Side effects: None.
Errors: Returns `nullptr` if the operation number is invalid.

#### `Root::GetTMatrix(mat44 & TMatrix, int Nop, int cellshift_a, int cellshift_b, int cellshift_c) -> int`
Calculates and returns the coordinate transformation matrix for the given symmetry operation number.
Params: `TMatrix` — A 4x4 matrix to store the transformation matrix. `Nop` — The symmetry operation number. `cellshift_a`, `cellshift_b`, `cellshift_c` —  Cell shift parameters.
Returns: 0 on success, error code on failure.
Side effects: Populates `TMatrix` with the transformation matrix.
Errors: Returns a non-zero error code if the calculation fails.

#### `Root::GetUDDHandle(UDR_TYPE udr_type, cpstr UDDataID) -> int`
Retrieves the handle for a user-defined data record.
Params: `udr_type` — The type of user-defined data. `UDDataID` — The ID of the user-defined data record.
Returns: The handle for the user-defined data record.
Side effects: None.
Errors: Returns a negative error code if the record is not found.

#### `Root::PDBCleanup(word CleanKey) -> word`
Performs PDB cleanup operations on the structure, ensuring compliance with PDB standards.
Params: `CleanKey` —  A key used to determine which cleanup operations to perform.
Returns: 0 on success, error code on failure.
Side effects: Modifies the structure's atom and residue data to comply with PDB standards.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::PutAtom(int index, int serNum, const AtomName atomName, const ResName resName, const ChainID chainID, int seqNum, const InsCode insCode, const AltLoc altLoc, const SegID segID, const Element element) -> int`
Adds an atom to the structure at the specified index, setting its properties.
Params: `index` — The index at which to add the atom. `serNum` — The serial number of the atom. `atomName` — The name of the atom. `resName` — The name of the residue. `chainID` — The chain ID of the residue. `seqNum` — The sequence number of the atom. `insCode` — The insertion code. `altLoc` — The alternative location code. `segID` — The segment ID. `element` — The element symbol.
Returns: 0 on success, error code on failure.
Side effects: Adds an atom to the structure, updating related data structures.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::PutPDBString(cpstr PDBString) -> ERROR_CODE`
Adds a PDB-keyworded string to the existing structure.
Params: `PDBString` — The string to add.
Returns: 0 on success, error code on failure.
Side effects: Adds the string as a remark to the current structure.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::ReadCoorFile(cpstr LFName) -> ERROR_CODE`
Reads coordinate data from a file in a coordinate format.
Params: `LFName` — The name of the coordinate file.
Returns: 0 on success, error code on failure.
Side effects: Reads coordinate data from the file and updates the structure's atom coordinates.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::RegisterUDInteger(UDR_TYPE udr_type, cpstr UDDataID) -> int`
Registers a user-defined integer data record.
Params: `udr_type` — The type of user-defined data. `UDDataID` — The ID of the user-defined data record.
Returns: The handle for the user-defined data record.
Side effects: Adds a user-defined integer data record to the structure's data table.
Errors: Returns a negative error code if the record cannot be added.

#### `Root::SetFlag(word Flag) -> void`
Sets a flag in the structure.
Params: `Flag` — The flag to set.
Returns: None
Side effects: Modifies the structure's internal flags.
Errors: None.

#### `Root::SetSpaceGroup(cpstr spGroup) -> int`
Sets the space group for the current structure.
Params: `spGroup` — The name of the space group.
Returns: 0 on success, error code on failure.
Side effects: Updates the structure's space group information.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::WriteCIFASCII(cpstr CIFFileName) -> ERROR_CODE`
Writes CIF data to an ASCII file.
Params: `CIFFileName` — The name of the CIF file.
Returns: 0 on success, error code on failure.
Side effects: Writes CIF data to the file.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::WritePDBASCII(cpstr PDBFileName) -> ERROR_CODE`
Writes PDB data to an ASCII file.
Params: `PDBFileName` — The name of the PDB file.
Returns: 0 on success, error code on failure.
Side effects: Writes PDB data to the file.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::DeleteAllModels() -> int`
Deletes all models from the structure.
Params: None
Returns: The number of models deleted.
Side effects: Destroys all model objects.
Errors: None.

#### `Root::FinishStructEdit() -> int`
Completes structural editing operations.
Params: None
Returns: 0 on success, error code on failure.
Side effects: Finalizes the current structure, potentially updating internal data structures.
Errors: None.

#### `Root::GetAtomI(int index) -> PAtom`
Returns the atom at the given serial number index in the structure's atom array.
Params: `index` — The index of the atom in the atom array (1-based).
Returns: A pointer to the atom at the specified index.
Side effects: None.
Errors: Returns `nullptr` if the index is out of bounds.

#### `Root::GetCell() -> mat44`
Returns the cell transformation matrix for the current structure.
Params: None
Returns: A 4x4 matrix representing the cell transformation.
Side effects: None.
Errors: None.

#### `Root::GetFileType() -> int`
Returns the file type of the currently loaded structure.
Params: None
Returns: An integer code representing the file type (e.g., 0 for PDB, 1 for CIF).
Side effects: None.
Errors: None.

#### `Root::GetInputBuffer(pstr Line, int & count) -> void`
Retrieves a line of input from the input buffer and stores it in the specified string.
Params: `Line` — A pointer to a string where the line will be stored. `count` — A reference to an integer that will be set to the number of characters in the line.
Returns: None
Side effects: Populates `Line` with a line from the input buffer.
Errors: None.

#### `Root::GetNumberOfSymOps() -> int`
Returns the number of symmetry operations defined in the structure.
Params: None
Returns: The number of symmetry operations.
Side effects: None.
Errors: None.

#### `Root::GetSpaceGroup() -> pstr`
Returns the name of the space group for the current structure.
Params: None
Returns: A pointer to a string containing the space group name.
Side effects: None.
Errors: Returns `nullptr` if the space group is not defined.

#### `Root::GetStructureTitle(pstr & L) -> void`
Retrieves the structure title and stores it in the specified string.
Params: `L` — A pointer to a string where the title will be stored.
Returns: None
Side effects: Populates `L` with the structure title.
Errors: None.

#### `Root::GetSymOp(int Nop) -> pstr`
Returns the XYZ symmetry operation name corresponding to the given symmetry operation number.
Params: `Nop` — The number of the symmetry operation.
Returns: A pointer to a string containing the name of the symmetry operation.
Side effects: None.
Errors: Returns `nullptr` if the operation number is invalid.

#### `Root::GetTMatrix(mat44 & TMatrix, int Nop, int cellshift_a, int cellshift_b, int cellshift_c) -> int`
Calculates and returns the coordinate transformation matrix for the given symmetry operation number.
Params: `TMatrix` — A 4x4 matrix to store the transformation matrix. `Nop` — The symmetry operation number. `cellshift_a`, `cellshift_b`, `cellshift_c` —  Cell shift parameters.
Returns: 0 on success, error code on failure.
Side effects: Populates `TMatrix` with the transformation matrix.
Errors: Returns a non-zero error code if the calculation fails.

#### `Root::GetUDDHandle(UDR_TYPE udr_type, cpstr UDDataID) -> int`
Retrieves the handle for a user-defined data record.
Params: `udr_type` — The type of user-defined data. `UDDataID` — The ID of the user-defined data record.
Returns: The handle for the user-defined data record.
Side effects: None.
Errors: Returns a negative error code if the record is not found.

#### `Root::PDBCleanup(word CleanKey) -> word`
Performs PDB cleanup operations on the structure, ensuring compliance with PDB standards.
Params: `CleanKey` —  A key used to determine which cleanup operations to perform.
Returns: 0 on success, error code on failure.
Side effects: Modifies the structure's atom and residue data to comply with PDB standards.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::PutAtom(int index, int serNum, const AtomName atomName, const ResName resName, const ChainID chainID, int seqNum, const InsCode insCode, const AltLoc altLoc, const SegID segID, const Element element) -> int`
Adds an atom to the structure at the specified index, setting its properties.
Params: `index` — The index at which to add the atom. `serNum` — The serial number of the atom. `atomName` — The name of the atom. `resName` — The name of the residue. `chainID` — The chain ID of the residue. `seqNum` — The sequence number of the atom. `insCode` — The insertion code. `altLoc` — The alternative location code. `segID` — The segment ID. `element` — The element symbol.
Returns: 0 on success, error code on failure.
Side effects: Adds an atom to the structure, updating related data structures.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::PutPDBString(cpstr PDBString) -> ERROR_CODE`
Adds a PDB-keyworded string to the existing structure.
Params: `PDBString` — The string to add.
Returns: 0 on success, error code on failure.
Side effects: Adds the string as a remark to the current structure.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::ReadCoorFile(cpstr LFName) -> ERROR_CODE`
Reads coordinate data from a file in a coordinate format.
Params: `LFName` — The name of the coordinate file.
Returns: 0 on success, error code on failure.
Side effects: Reads coordinate data from the file and updates the structure's atom coordinates.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::RegisterUDInteger(UDR_TYPE udr_type, cpstr UDDataID) -> int`
Registers a user-defined integer data record.
Params: `udr_type` — The type of user-defined data. `UDDataID` — The ID of the user-defined data record.
Returns: The handle for the user-defined data record.
Side effects: Adds a user-defined integer data record to the structure's data table.
Errors: Returns a negative error code if the record cannot be added.

#### `Root::SetFlag(word Flag) -> void`
Sets a flag in the structure.
Params: `Flag` — The flag to set.
Returns: None
Side effects: Modifies the structure's internal flags.
Errors: None.

#### `Root::SetSpaceGroup(cpstr spGroup) -> int`
Sets the space group for the current structure.
Params: `spGroup` — The name of the space group.
Returns: 0 on success, error code on failure.
Side effects: Updates the structure's space group information.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::WriteCIFASCII(cpstr CIFFileName) -> ERROR_CODE`
Writes CIF data to an ASCII file.
Params: `CIFFileName` — The name of the CIF file.
Returns: 0 on success, error code on failure.
Side effects: Writes CIF data to the file.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::WritePDBASCII(cpstr PDBFileName) -> ERROR_CODE`
Writes PDB data to an ASCII file.
Params: `PDBFileName` — The name of the PDB file.
Returns: 0 on success, error code on failure.
Side effects: Writes PDB data to the file.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::DeleteAllModels() -> int`
Deletes all models from the structure.
Params: None
Returns: The number of models deleted.
Side effects: Destroys all model objects.
Errors: None.

#### `Root::FinishStructEdit() -> int`
Completes structural editing operations.
Params: None
Returns: 0 on success, error code on failure.
Side effects: Finalizes the current structure, potentially updating internal data structures.
Errors: None.

#### `Root::GetAtomI(int index) -> PAtom`
Returns the atom at the given serial number index in the structure's atom array.
Params: `index` — The index of the atom in the atom array (1-based).
Returns: A pointer to the atom at the specified index.
Side effects: None.
Errors: Returns `nullptr` if the index is out of bounds.

#### `Root::GetCell() -> mat44`
Returns the cell transformation matrix for the current structure.
Params: None
Returns: A 4x4 matrix representing the cell transformation.
Side effects: None.
Errors: None.

#### `Root::GetFileType() -> int`
Returns the file type of the currently loaded structure.
Params: None
Returns: An integer code representing the file type (e.g., 0 for PDB, 1 for CIF).
Side effects: None.
Errors: None.

#### `Root::GetInputBuffer(pstr Line, int & count) -> void`
Retrieves a line of input from the input buffer and stores it in the specified string.
Params: `Line` — A pointer to a string where the line will be stored. `count` — A reference to an integer that will be set to the number of characters in the line.
Returns: None
Side effects: Populates `Line` with a line from the input buffer.
Errors: None.

#### `Root::GetNumberOfSymOps() -> int`
Returns the number of symmetry operations defined in the structure.
Params: None
Returns: The number of symmetry operations.
Side effects: None.
Errors: None.

#### `Root::GetSpaceGroup() -> pstr`
Returns the name of the space group for the current structure.
Params: None
Returns: A pointer to a string containing the space group name.
Side effects: None.
Errors: Returns `nullptr` if the space group is not defined.

#### `Root::GetStructureTitle(pstr & L) -> void`
Retrieves the structure title and stores it in the specified string.
Params: `L` — A pointer to a string where the title will be stored.
Returns: None
Side effects: Populates `L` with the structure title.
Errors: None.

#### `Root::GetSymOp(int Nop) -> pstr`
Returns the XYZ symmetry operation name corresponding to the given symmetry operation number.
Params: `Nop` — The number of the symmetry operation.
Returns: A pointer to a string containing the name of the symmetry operation.
Side effects: None.
Errors: Returns `nullptr` if the operation number is invalid.

#### `Root::GetTMatrix(mat44 & TMatrix, int Nop, int cellshift_a, int cellshift_b, int cellshift_c) -> int`
Calculates and returns the coordinate transformation matrix for the given symmetry operation number.
Params: `TMatrix` — A 4x4 matrix to store the transformation matrix. `Nop` — The symmetry operation number. `cellshift_a`, `cellshift_b`, `cellshift_c` —  Cell shift parameters.
Returns: 0 on success, error code on failure.
Side effects: Populates `TMatrix` with the transformation matrix.
Errors: Returns a non-zero error code if the calculation fails.

#### `Root::GetUDDHandle(UDR_TYPE udr_type, cpstr UDDataID) -> int`
Retrieves the handle for a user-defined data record.
Params: `udr_type` — The type of user-defined data. `UDDataID` — The ID of the user-defined data record.
Returns: The handle for the user-defined data record.
Side effects: None.
Errors: Returns a negative error code if the record is not found.

#### `Root::PDBCleanup(word CleanKey) -> word`
Performs PDB cleanup operations on the structure, ensuring compliance with PDB standards.
Params: `CleanKey` —  A key used to determine which cleanup operations to perform.
Returns: 0 on success, error code on failure.
Side effects: Modifies the structure's atom and residue data to comply with PDB standards.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::PutAtom(int index, int serNum, const AtomName atomName, const ResName resName, const ChainID chainID, int seqNum, const InsCode insCode, const AltLoc altLoc, const SegID segID, const Element element) -> int`
Adds an atom to the structure at the specified index, setting its properties.
Params: `index` — The index at which to add the atom. `serNum` — The serial number of the atom. `atomName` — The name of the atom. `resName` — The name of the residue. `chainID` — The chain ID of the residue. `seqNum` — The sequence number of the atom. `insCode` — The insertion code. `altLoc` — The alternative location code. `segID` — The segment ID. `element` — The element symbol.
Returns: 0 on success, error code on failure.
Side effects: Adds an atom to the structure, updating related data structures.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::PutPDBString(cpstr PDBString) -> ERROR_CODE`
Adds a PDB-keyworded string to the existing structure.
Params: `PDBString` — The string to add.
Returns: 0 on success, error code on failure.
Side effects: Adds the string as a remark to the current structure.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::ReadCoorFile(cpstr LFName) -> ERROR_CODE`
Reads coordinate data from a file in a coordinate format.
Params: `LFName` — The name of the coordinate file.
Returns: 0 on success, error code on failure.
Side effects: Reads coordinate data from the file and updates the structure's atom coordinates.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::RegisterUDInteger(UDR_TYPE udr_type, cpstr UDDataID) -> int`
Registers a user-defined integer data record.
Params: `udr_type` — The type of user-defined data. `UDDataID` — The ID of the user-defined data record.
Returns: The handle for the user-defined data record.
Side effects: Adds a user-defined integer data record to the structure's data table.
Errors: Returns a negative error code if the record cannot be added.

#### `Root::SetFlag(word Flag) -> void`
Sets a flag in the structure.
Params: `Flag` — The flag to set.
Returns: None
Side effects: Modifies the structure's internal flags.
Errors: None.

#### `Root::SetSpaceGroup(cpstr spGroup) -> int`
Sets the space group for the current structure.
Params: `spGroup` — The name of the space group.
Returns: 0 on success, error code on failure.
Side effects: Updates the structure's space group information.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::WriteCIFASCII(cpstr CIFFileName) -> ERROR_CODE`
Writes CIF data to an ASCII file.
Params: `CIFFileName` — The name of the CIF file.
Returns: 0 on success, error code on failure.
Side effects: Writes CIF data to the file.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::WritePDBASCII(cpstr PDBFileName) -> ERROR_CODE`
Writes PDB data to an ASCII file.
Params: `PDBFileName` — The name of the PDB file.
Returns: 0 on success, error code on failure.
Side effects: Writes PDB data to the file.
Errors: Returns a non-zero error code if an error occurs.

#### `Root::DeleteAllModels() -> int`
Deletes all models from the structure.
Params: None
Returns: The number of models deleted.
Side effects: Destroys all model objects.
Errors: None.

#### `Root::FinishStructEdit() -> int`
Completes structural editing operations.
Params: None
Returns: 0 on success, error code on failure.
Side effects: Finalizes the current structure, potentially updating internal data structures.
Errors: None.

#### `Root::GetAtomI(int index) -> PAtom`
Returns the atom at the given serial number index in the structure's atom array.
Params: `index` — The index of the atom in the atom array (1-based).
Returns: A pointer to the atom at the specified index.
Side effects: None.
Errors: Returns `nullptr` if the index is out of bounds.

#### `Root::GetCell() -> mat44`
Returns the cell transformation matrix for the current structure.
Params: None
Returns: A 4x4 matrix representing the cell transformation.
Side effects: None.
Errors: None.

#### `Root::GetFileType() -> int`
Returns the file type of the currently loaded structure.
Params: None
Returns: An integer code representing the file type (e.g., 0 for PDB, 1 for CIF).
Side effects: None.
Errors: None.

#### `Root::GetInputBuffer(pstr Line, int & count) -> void`
Retrieves a line of input from the input buffer and stores it in the specified string.
Params: `Line` — A pointer to a string where the line will be stored. `count` — A reference to an integer that will be set to the number of characters in the line.
Returns

---

## `SelManager`

Manages selection masks and atom selection within a structural biology file.

#### `SelManager::NewSelection()`
Creates a new selection mask and returns its handle. A handle is always a positive (non-zero) integer. Calling `NewSelection()` is the only way to create a new selection mask. Masks will be automatically copied from another MMDB if coordinates are copied; if this is the case, the mask handles will be inherited from the source MMDB as well. The masks will also be automatically deleted if coordinates are deleted. Returns: `int` - The handle of the newly created selection mask.

#### `SelManager::GetSelType(int selHnd)`
Returns the selection type associated with a given selection handle. Returns: `STYPE_XXXXX` - The selection type (e.g., STYPE_ATOM, STYPE_RESIDUE).

#### `SelManager::DeleteSelection(int selHnd)`
Deletes the specified selection mask and removes the corresponding selection attributes from all atoms selected with this mask. If an atom was selected also with other mask(s), the other selection(s) will remain, provided that the corresponding masks are valid. Returns: `void` - Deletes the mask and unselects atoms.

#### `SelManager::DeleteAllSelections()`
Deletes all selection masks and unselects all atoms in the file. All mask handles become invalid. Returns: `void` - Deletes all selection masks and resets all atoms to unselected.

#### `SelManager::SelectAtoms(int selHnd, int iSer1, int iSer2, SELECTION_KEY selKey)`
Selects atoms in the serial number range of `iSer1` to `iSer2` by adding them to the set of atoms marked by the given mask. If `iSer1=iSer2=0` then all atoms are selected. Each atom may be selected by a number of masks simultaneously. Returns: `void` - Adds atoms to the selection based on the serial number range.

#### `SelManager::SelectAtoms(int selHnd, ivector asn, int nsn, SELECTION_KEY selKey)`
Selects atoms based on the provided sequence number vector `asn` of length `nsn`. Returns: `void` - Adds atoms to the selection based on the sequence number vector.

#### `SelManager::UnselectAtoms(int selHnd, int iSer1, int iSer2)`
Clears the specified mask for atoms in the serial number range of `iSer1` to `iSer2`. If `iSer1=iSer2=0` then all atoms are cleared of the specified mask. If `selHnd` is set to 0, then the atoms are cleared of any mask. Returns: `void` - Unselects atoms based on the serial number range.

#### `SelManager::SelectAtom(int selHnd, PAtom A, SELECTION_KEY selKey, bool makeIndex)`
Selects a single atom according to the value of selection key. If `makeIndex` is false, then the routine does not update the selection index. In order to update the index after all single-atom selections are done, use `MakeSelIndex(selHnd)` found next. Returns: `void` - Adds a single atom to the selection based on the key.

#### `SelManager::SelectResidue(int selHnd, PResidue Res, SELECTION_TYPE sType, SELECTION_KEY sKey, bool makeIndex)`
Selects a single residue, chain or model, or all their hierarchical descendants depending on the value of `sType` (i.e. atoms, residues (in chain and model) and chains (in model only). Returns: `void` - Adds a single residue to the selection based on the key and type.

#### `SelManager::SelectChain(int selHnd, PChain chain, SELECTION_TYPE sType, SELECTION_KEY sKey, bool makeIndex)`
Selects a single chain, or all their hierarchical descendants depending on the value of `sType`. Returns: `void` - Adds a single chain to the selection based on the key and type.

#### `SelManager::SelectModel(int selHnd, PModel mdl, SELECTION_TYPE sType, SELECTION_KEY sKey, bool makeIndex)`
Selects a single model, or all their hierarchical descendants depending on the value of `sType`. Returns: `void` - Adds a single model to the selection based on the key and type.

#### `SelManager::MakeSelIndex(int selHnd)`
Calculates selection index for selection addressed by `selHnd`. All selection functions except the `SelectAtom(..)` above, update selection index automatically. This function is for use after a series of calls to `SelectAtom(..)` with `makeIndex` parameter set false. Returns: `int` - The number of selected objects.

#### `SelManager::MakeAllSelIndexes()`
Calculates the selection index for all selections in the object list. Returns: `void` - Updates all selection indices.

#### `SelManager::SelectAtoms(int selHnd, int iModel, cpstr Chains, int ResNo1, cpstr Ins1, int ResNo2, cpstr Ins2)`
Selects atoms based on a complex set of criteria, including model number, chain list, residue number range, and insertion code range. Returns: `void` - Adds atoms to the selection based on the specified criteria.

#### `SelManager::SelectAtoms(int selHnd, int iModel, cpstr Chains, int ResNo1, cpstr Ins1, int ResNo2, cpstr Ins2)`
Same as above, but using a `ivector` for the sequence numbers. Returns: `void` - Adds atoms to the selection based on the specified criteria.

#### `SelManager::SelectUDD(int selHnd, SELECTION_TYPE sType, int UDDhandle, int selMin, int selMax, SELECTION_KEY sKey)`
Selects atoms based on a user-defined data value, offering flexibility in selection criteria. Returns: `void` - Adds atoms to the selection based on the UDD value.

#### `SelManager::SelectUDD(int selHnd, SELECTION_TYPE sType, int UDDhandle, cpstr selStr, int cmpRule, SELECTION_KEY sKey)`
Same as above, but using a string for the selection criteria. Returns: `void` - Adds atoms to the selection based on the UDD value.

#### `SelManager::SelectSphere(int selHnd, SELECTION_TYPE sType, realtype x, realtype y, realtype z, realtype r, SELECTION_KEY sKey)`
Selects atoms within a sphere centered at the given point with a specified radius. Returns: `void` - Adds atoms to the selection based on the sphere parameters.

#### `SelManager::SelectCylinder(int selHnd, SELECTION_TYPE sType, realtype x1, realtype y1, realtype z1, realtype x2, realtype y2, realtype z2, realtype r, SELECTION_KEY sKey)`
Selects atoms within a cylinder defined by its axis and radius. Returns: `void` - Adds atoms to the selection based on the cylinder parameters.

#### `SelManager::SelectSlab(int selHnd, SELECTION_TYPE sType, realtype a, realtype b, realtype c, realtype d, realtype r, SELECTION_KEY sKey)`
Selects atoms located on one side of a plane. Returns: `void` - Adds atoms to the selection based on the plane parameters.

#### `SelManager::SelectNeighbours(int selHnd, SELECTION_TYPE sType, PPAtom sA, int alen, realtype d1, realtype d2, SELECTION_KEY sKey)`
Selects atoms within a specified distance from another selection of atoms. Returns: `void` - Adds atoms to the selection based on the neighbour distance.

#### `SelManager::DeleteSelObjects(int selHnd)`
Deletes the objects associated with a given selection handle. Returns: `void` - Removes the selected atoms from the MMDB.

#### `SelManager::GetSelLength(int selHnd)`
Returns the length of the selection mask. Returns: `int` - The number of atoms in the selection.

#### `SelManager::GetSelIndex(PPAtom & SelAtom, int & nSelAtoms)`
Returns an array of atom handles of atoms selected by a particular selection mask. Returns: `void` - Populates the provided array with the selection indexes.

#### `SelManager::DeleteAllSelections()`
Deletes all selection masks and unselects all atoms in the file. Returns: `void` - Clears all selections and resets all atoms.


---

## `Sheet`
Represents a sheet in a protein structure.

#### `Sheet::Sheet()`
Constructs a new `Sheet` object.
Params: None
Returns: `void`
Side effects: Initializes the `Sheet` object to a default state.

#### `Sheet::Sheet(io::RPStream Object)`
Constructs a new `Sheet` object from an `io::RPStream` object.
Params: `io::RPStream Object` — The input stream to read the sheet data from.
Returns: `void`
Side effects: Initializes the `Sheet` object with the data from the input stream.

#### `Sheet::~Sheet()`
Destroys the `Sheet` object, releasing any allocated memory.
Params: None
Returns: `void`
Side effects: Frees the memory associated with the `Sheet` object.

#### `Sheet::FreeMemory()`
Frees the memory allocated for the strands within the `Sheet` object.
Params: None
Returns: `void`
Side effects: Deallocates the `strand` array.

#### `Sheet::OrderSheet()`
Orders the strands within the `Sheet` object.
Params: None
Returns: `void`
Side effects:  Sorts the strands based on a defined order.

#### `Sheet::PDBASCIIDump(io::RFile f)`
Writes the sheet data in PDB ASCII format to the specified `io::RFile`.
Params: `io::RFile f` — The output file to write the PDB ASCII data to.
Returns: `void`
Side effects: Writes the sheet's structure data to a PDB ASCII file.

#### `Sheet::MakeCIF(mmcif::PData CIF)`
Creates a mmCIF data object (`mmcif::PData`) representing the sheet.
Params: `mmcif::PData CIF` — The output mmCIF data object.
Returns: `void`
Side effects: Populates the `mmcif::PData` object with information about the sheet.

#### `Sheet::ConvertPDBASCII(cpstr S)`
Converts a PDB ASCII string to a `Sheet` object.
Params: `cpstr S` — The input PDB ASCII string.
Returns: `int` — The number of strands in the sheet.
Errors: Returns `-1` if the conversion fails.

#### `Sheet::GetCIF(mmcif::PData CIF)`
Retrieves the mmCIF data object (`mmcif::PData`) associated with the sheet.
Params: `mmcif::PData CIF` — The output mmCIF data object.
Returns: `int` — The number of strands in the sheet.
Side effects: Populates the `mmcif::PData` object with information about the sheet.

#### `Sheet::Copy(PSheet sheet)`
Creates a copy of another `Sheet` object.
Params: `PSheet sheet` — The `Sheet` object to copy.
Returns: `void`
Side effects: Duplicates the data and structure of the input `Sheet` object.

#### `Sheet::write(io::RFile f)`
Writes the sheet data to the specified `io::RFile`.
Params: `io::RFile f` — The output file to write the sheet data to.
Returns: `void`
Side effects: Writes the sheet's internal representation to a file.

#### `Sheet::read(io::RFile f)`
Reads the sheet data from the specified `io::RFile`.
Params: `io::RFile f` — The input file to read the sheet data from.
Returns: `void`
Side effects: Populates the `Sheet` object with data read from the file.

#### `Sheet::InitSheet()`
Initializes the `Sheet` object, performing necessary setup steps.
Params: None
Returns: `void`
Side effects: Sets default values and initializes internal data structures.

#### `Sheet::CIFFindStrands(mmcif::PData CIF, cpstr Category)`
Finds the strands within a mmCIF data object based on a category.
Params: `mmcif::PData CIF` — The mmCIF data object to search within.
`cpstr Category` — The category to filter strands by.
Returns: `void`
Side effects: Populates the sheet with strands matching the category.

#### `Sheet::TryStrand(int strand_no)`
Attempts to set the given strand number to a valid state within the sheet.
Params: `int strand_no` — The index of the strand to process.
Returns: `void`
Side effects: Performs actions specific to the strand number.

#### `Sheet::GetStrand(int strand_no)`
Retrieves the `PPStrand` object for the specified strand number.
Params: `int strand_no` — The index of the strand to retrieve.
Returns: `PPStrand` — The `PPStrand` object for the specified strand.
Side effects: Returns a pointer to the corresponding strand.


---

## `Sheets`
This class represents a collection of protein sheets within a PDB or CIF file.

#### `Sheets::Sheets()`
Constructs a Sheets object.
Params: None
Returns: `void`
Side effects: Initializes the Sheets object.
Errors: None

#### `Sheets::Sheets(io::RPStream Object)`
Constructs a Sheets object from an RPStream.
Params: `Object` — An RPStream object.
Returns: `void`
Side effects: Initializes the Sheets object using the provided RPStream.
Errors: None

#### `Sheets::~Sheets()`
Destroys the Sheets object, freeing allocated memory.
Params: None
Returns: `void`
Side effects: Frees memory associated with the Sheets object.
Errors: None

#### `Sheets::FreeMemory()`
Frees the memory allocated by the Sheets object.
Params: None
Returns: `void`
Side effects: Releases all dynamically allocated memory.
Errors: None

#### `Sheets::PDBASCIIDump(io::RFile f)`
Writes the sheet data to an RFile in PDB ASCII format.
Params: `f` — An RFile object to write to.
Returns: `void`
Side effects: Writes the sheet data to the RFile.
Errors: None

#### `Sheets::MakeCIF(mmcif::PData CIF)`
Converts the sheet data into a CIF (Crystallographic Information File) format.
Params: `CIF` — An output mmcif::PData object to populate.
Returns: `void`
Side effects: Populates the `CIF` object with sheet data.
Errors: None

#### `Sheets::ConvertPDBASCII(cpstr S)`
Converts a PDB ASCII string into a Sheets object.
Params: `S` — A PDB ASCII string.
Returns: `int` — The number of sheets successfully parsed, or an error code if parsing fails.
Side effects: Parses the PDB ASCII string and populates the Sheets object.
Errors: Returns a non-zero error code on failure.

#### `Sheets::GetCIF(mmcif::PData CIF)`
Retrieves the CIF data for the sheets in the Sheets object.
Params: `CIF` — An input mmcif::PData object to populate.
Returns: `int` — The number of sheets successfully written to the CIF, or an error code if writing fails.
Side effects: Populates the `CIF` object with sheet data.
Errors: Returns a non-zero error code on failure.

#### `Sheets::Copy(PSheets Sheets)`
Copies the contents of one Sheets object to another.
Params: `Sheets` — Another Sheets object to copy data to.
Returns: `void`
Side effects: Copies sheet data from one Sheets object to another.
Errors: None

#### `Sheets::write(io::RFile f)`
Writes the sheet data to an RFile.
Params: `f` — An RFile object to write to.
Returns: `void`
Side effects: Writes sheet data to the RFile.
Errors: None

#### `Sheets::read(io::RFile f)`
Reads the sheet data from an RFile.
Params: `f` — An RFile object to read from.
Returns: `void`
Side effects: Reads sheet data from the RFile.
Errors: None

#### `Sheets::InitSheets()`
Initializes the internal representation of the sheets.
Params: None
Returns: `void`
Side effects: Initializes the sheet data.
Errors: None

#### `Sheets::CIFFindSheets(mmcif::PData CIF, cpstr Category)`
Searches the CIF data for sheets matching a given category.
Params: `CIF` — The CIF data to search. `Category` — The category to search for.
Returns: `void`
Side effects: Populates the `CIF` object with sheets matching the category.
Errors: None


---

## `Strand`
Represents a single strand within a protein structure, including its sequence, residues, and insertion codes.

#### `Strand::Strand()`
Constructs a new `Strand` object with default values.
Params: None
Returns: None
Side effects: Initializes the `Strand` object's member variables.
Errors: None

#### `Strand::Strand(io::RPStream Object)`
Constructs a new `Strand` object from an `io::RPStream` object.
Params: `io::RPStream Object` — an `RPStream` object containing the strand data.
Returns: None
Side effects: Initializes the `Strand` object's member variables based on the data in the `RPStream` object.
Errors: None

#### `Strand::~Strand()`
Destroys the `Strand` object, releasing any allocated memory.
Params: None
Returns: None
Side effects: Releases the `Strand` object's resources.
Errors: None

#### `Strand::PDBASCIIDump(pstr S)`
Writes the strand's data to a given string `S` in PDB ASCII format.
Params: `pstr S` — a pointer to a string to store the PDB ASCII representation.
Returns: None
Side effects: Modifies the string `S` with the strand's data in PDB ASCII format.
Errors: None

#### `Strand::MakeCIF(mmcif::PData CIF)`
Generates a mmCIF data structure `CIF` from the current `Strand` data.
Params: `mmcif::PData CIF` — a pointer to a mmCIF data structure object.
Returns: None
Side effects: Populates the `CIF` object with the strand's data.
Errors: None

#### `Strand::ConvertPDBASCII(cpstr S)`
Parses a PDB ASCII string `S` and constructs a `Strand` object.
Params: `cpstr S` — a pointer to a string containing PDB ASCII data.
Returns: `ERROR_CODE` — The error code if parsing fails, or zero on success.
Side effects: Initializes the `Strand` object with data parsed from the input string.
Errors: Returns an error code if the input string is invalid or parsing fails.

#### `Strand::GetCIF(mmcif::PData CIF, cpstr sheet_id)`
Writes the strand's data to the mmCIF data structure `CIF` with the specified sheet ID.
Params: `mmcif::PData CIF` — a pointer to a mmCIF data structure object.
Params: `cpstr sheet_id` — a string containing the sheet ID.
Returns: `int` — The number of bytes written to the mmCIF data structure.
Side effects: Populates the `CIF` object with the strand's data.
Errors: None

#### `Strand::Copy(PStrand Strand)`
Copies the data from another `Strand` object into this `Strand` object.
Params: `PStrand Strand` — a pointer to another `Strand` object.
Returns: None
Side effects: Copies the data from the input `Strand` object to this `Strand` object.
Errors: None

#### `Strand::write(io::RFile f)`
Writes the strand's data to an `io::RFile` object.
Params: `io::RFile f` — a pointer to an `RFile` object to write to.
Returns: None
Side effects: Writes the strand's data to the output stream.
Errors: None

#### `Strand::read(io::RFile f)`
Reads the strand's data from an `io::RFile` object.
Params: `io::RFile f` — a pointer to an `RFile` object to read from.
Returns: None
Side effects: Reads the strand's data from the input stream and initializes the object's member variables.
Errors: None

#### `Strand::InitStrand()`
Initializes the `Strand` object's member variables.
Params: None
Returns: None
Side effects: Sets the initial values of all member variables of the `Strand` object.
Errors: None


---

## `Struct`

This class represents a mmCIF structure, providing methods to manage its fields and optimize its data representation.

#### `Struct::AddField(cpstr F, cpstr T, bool Concatenate)`
Adds a field to the structure, allowing for concatenation of existing fields or replacement of the field's value.

Params: `F` — The field value (string). `T` — The tag name (string). `Concatenate` —  Boolean indicating whether to append to the existing field or replace it.

Returns: `void`

Side effects: Adds a new field to the structure's internal data structure.

Errors: None.

#### `Struct::GetCategoryID() -> MMCIF_ITEM`
Returns the category ID for this structure, which is `MMCIF_Struct`.

Params: None

Returns: `MMCIF_ITEM` — The category ID.

Side effects: None.

#### `Struct::Optimize()`
Optimizes the structure's data for RAM usage and data access speed by sorting the field data.

Params: None

Returns: `void`

Side effects: Sorts the internal field data structure, updates indexes, and improves data access.

Errors: None.

#### `Struct::GetField(int tagNo) -> pstr` // *** USED IN COOT ***
Returns the value of the field corresponding to the specified tag number (index).

Params: `tagNo` — The tag number (integer, 0-based indexing).

Returns: `pstr` — The string value of the field, or `NULL` if the tag doesn't exist.

Side effects: Reads the value from the structure's internal data.

Errors: `tagNo` out of range (returns `NULL`).

#### `Struct::GetString(pstr & S, cpstr TName, bool Remove)` // *** USED IN COOT ***
Fetches the value associated with the specified tag name and stores it in the provided string pointer.

Params: `S` — A pointer to the string where the retrieved value will be stored. `TName` — The tag name (string). `Remove` — Boolean flag to remove the tag and its value from the structure after reading.

Returns: `int` — Return code: `CIFRC_NoTag` if the tag doesn't exist, `CIFRC_NoField` if the field doesn't exist, or `CIFRC_Ok` on success.

Side effects: Reads the value from the structure's internal data and stores it in the `S` string pointer.  If `Remove` is true, the tag and corresponding data are deleted.

Errors: `S` is `NULL`.

#### `Struct::DeleteField(cpstr TName)` // *** USED IN COOT ***
Deletes the field associated with the specified tag name from the structure.

Params: `TName` — The tag name (string).

Returns: `int` — Return code: `0` if the field was deleted, negative value otherwise.

Side effects: Removes the field from the structure's internal data.

Errors: `TName` doesn't exist in the structure.

#### `Struct::GetReal(realtype & R, cpstr TName, bool Remove)` // *** USED IN COOT ***
Fetches the real number (floating-point) value associated with the specified tag name.

Params: `R` — A reference to the `realtype` variable to store the retrieved value. `TName` — The tag name (string). `Remove` — Boolean flag to remove the tag and its value from the structure after reading.

Returns: `int` — Return code: `CIFRC_NoTag` if the tag doesn't exist, `CIFRC_NoField` if the field doesn't exist, `CIFRC_WrongFormat` if the value isn't a real number, or `CIFRC_NoData` if the value is "data not available" or "data not given".

Side effects: Reads the value from the structure's internal data and stores it in the `R` real type. If `Remove` is true, the tag and corresponding data are deleted.

Errors: `TName` doesn't exist or the value is not a real number.

#### `Struct::GetInteger(int & I, cpstr TName, bool Remove)` // *** USED IN COOT ***
Fetches the integer value associated with the specified tag name.

Params: `I` — A reference to the `int` variable to store the retrieved value. `TName` — The tag name (string). `Remove` — Boolean flag to remove the tag and its value from the structure after reading.

Returns: `int` — Return code: `CIFRC_NoTag` if the tag doesn't exist, `CIFRC_NoField` if the field doesn't exist, `CIFRC_WrongFormat` if the value isn't an integer, or `CIFRC_NoData` if the value is "data not available" or "data not given".

Side effects: Reads the value from the structure's internal data and stores it in the `I` integer variable. If `Remove` is true, the tag and corresponding data are deleted.

Errors: `TName` doesn't exist or the value isn't an integer.

#### `Struct::PutString(cpstr S, cpstr TName, bool NonBlankOnly)` // *** USED IN COOT ***
Sets the string value for the specified tag name.

Params: `S` — The string value to be set. `TName` — The tag name (string). `NonBlankOnly` — Boolean flag controlling the value to be set when `S` is blank or white-space-only.

Returns: `void`

Side effects: Updates the structure's internal data with the provided string value. If `S` is `NULL`, the "data not given" value is set. If `S` is empty string, the "data not available" value is set.

Errors: `TName` doesn't exist in the structure.

#### `Struct::PutDate(cpstr T)` // *** USED IN COOT ***
Sets the date value for the specified tag name.

Params: `T` — The tag name (string).

Returns: `void`

Side effects: Updates the structure's internal data with the date string.

Errors: `TName` doesn't exist in the structure.

#### `Struct::PutNoData(int NoDataType, cpstr T)` // *** USED IN COOT ***
Sets either "data not given" or "data not available" value for the specified tag.

Params: `NoDataType` — The type of missing data to set (`CIF_NODATA_DOT` or `CIF_NODATA_QUESTION`). `T` — The tag name (string).

Returns: `void`

Side effects: Updates the structure's internal data with the specified missing data string.

Errors: `TName` doesn't exist in the structure.

#### `Struct::PutReal(realtype R, cpstr TName, int prec)` // *** USED IN COOT ***
Sets the floating-point value for the specified tag name.

Params: `R` — The floating-point value to be set. `TName` — The tag name (string). `prec` — The floating-point precision.

Returns: `void`

Side effects: Updates the structure's internal data with the provided real number.

Errors: `TName` doesn't exist in the structure.

#### `Struct::PutReal(realtype R, cpstr TName, cpstr format)` // *** USED IN COOT ***
Sets the floating-point value for the specified tag name.

Params: `R` — The floating-point value to be set. `TName` — The tag name (string). `format` — The format string for the conversion.

Returns: `void`

Side effects: Updates the structure's internal data with the provided real number.

Errors: `TName` doesn't exist in the structure.

#### `Struct::PutInteger(int I, cpstr TName)` // *** USED IN COOT ***
Sets the integer value for the specified tag name.

Params: `I` — The integer value to be set. `TName` — The tag name (string).

Returns: `void`

Side effects: Updates the structure's internal data with the provided integer.

Errors: `TName` doesn't exist in the structure.

#### `Struct::WriteMMCIFStruct(cpstr FName, io::GZ_MODE gzipMode)` // *** USED IN COOT ***
Writes the structure's data in mmCIF format into a file.

Params: `FName` — The file name (string). `gzipMode` — Controls file compression: `GZM_NONE`, `GZM_CHECK`, `GZM_ENFORCE_GZIP` or `GZM_ENFORCE_COMPRESS`.

Returns: `bool` — `true` on success, `false` on failure.

Side effects: Writes the structure's data to the specified file in mmCIF format.

Errors: Failure to open file for writing.

#### `Struct::Copy(PCategory Struct)` // *** USED IN COOT ***
Performs a deep copy of the structure.

Params: `Struct` — The structure to be copied from.

Returns: `void`

Side effects: Creates a completely new structure with the same data.

Errors: None.

#### `Struct::write(io::RFile f)` // *** USED IN COOT ***
Writes structure data to a stream.

Params: `f` — Reference to stream class.

Returns: `void`

Side effects: Writes data to a stream.

#### `Struct::read(io::RFile f)` // *** USED IN COOT ***
Reads structure data from a stream.

Params: `f` — Reference to stream class.

Returns: `void`

Side effects: Reads data from a stream.

#### `Struct::InitStruct()`
Initializes the structure's internal data.

Params: None

Returns: `void`

Side effects: Initializes the internal data structure of the `Struct` object.

Errors: None.

#### `Struct::FreeMemory()`
Frees memory allocated by the `Struct` object.

Params: None

Returns: `void`

Side effects: Releases all dynamically allocated memory held by the `Struct` object.

Errors: None.


---

## `SymOps`
This class provides functions for manipulating symmetry operations within an MMDB2 structure.

#### `SymOps::SymOps()`
Constructs a new SymOps object.
Params: None
Returns: None
Side effects: Initializes the SymOps object.

#### `SymOps::SymOps(io::RPStream Object)`
Constructs a new SymOps object from an RPStream object.
Params: `io::RPStream Object` — The RPStream object to initialize from.
Returns: None
Side effects: Initializes the SymOps object.

#### `SymOps::~SymOps()`
Destructs the SymOps object, freeing allocated memory.
Params: None
Returns: None
Side effects: Releases all allocated memory.

#### `SymOps::FreeMemory()`
// *** USED IN COOT ***
Releases any allocated memory associated with the SymOps object.
Params: None
Returns: None
Side effects: Frees memory allocated internally.

#### `SymOps::SetGroupSymopLib(cpstr SpaceGroup, cpstr symop_lib=NULL) -> int`
// *** USED IN COOT ***
Sets the SpaceGroup name from a symop.lib file.
Params: `cpstr SpaceGroup` — The name of the SpaceGroup. `cpstr symop_lib` — Optional symop.lib file (defaults to NULL).
Returns: `int` — Returns `SYMOP_Ok` on success.
Side effects: Updates the internal SpaceGroup representation.
Errors: Returns a non-zero value if the SpaceGroup cannot be found.

#### `SymOps::SetGroup(cpstr SpaceGroup, cpstr syminfo_lib=NULL) -> int`
// *** USED IN COOT ***
Sets the SpaceGroup name from a syminfo.lib file.
Params: `cpstr SpaceGroup` — The name of the SpaceGroup. `cpstr syminfo_lib` — Optional syminfo.lib file (defaults to NULL).
Returns: `int` — Returns `SYMOP_Ok` on success.
Side effects: Updates the internal SpaceGroup representation.
Errors: Returns a non-zero value if the SpaceGroup cannot be found.

#### `SymOps::Reset()`
// *** USED IN COOT ***
Removes all symmetry operations from the SymOps object.
Params: None
Returns: None
Side effects: Clears the internal symmetry operation list.

#### `SymOps::AddSymOp(cpstr XYZOperation) -> int`
// *** USED IN COOT ***
Adds a symmetry operation to the SymOps object.
Params: `cpstr XYZOperation` — The name of the symmetry operation.
Returns: `int` — Returns the index of the added operation, or a non-zero value if the operation is invalid.
Side effects: Adds the operation to the internal list.
Errors: Returns a non-zero value if the operation name is not valid.

#### `SymOps::PutGroupName(cpstr SpGroupName)`
// *** USED IN COOT ***
Sets the group name.
Params: `cpstr SpGroupName` — The group name.
Returns: None
Side effects: Updates the internal group name.

#### `SymOps::GetNofSymOps() -> int`
// *** USED IN COOT ***
Returns the number of symmetry operations.
Params: None
Returns: `int` — The number of symmetry operations.
Side effects: None

#### `SymOps::GetSymOp(int Nop) -> pstr`
// *** USED IN COOT ***
Retrieves the name of the symmetry operation with the given index.
Params: `int Nop` — The index of the symmetry operation.
Returns: `pstr` — The name of the symmetry operation.
Side effects: None
Errors: Returns NULL if the index is out of range.

#### `SymOps::Transform(realtype & x, realtype & y, realtype & z, int Nop) -> int`
// *** USED IN COOT ***
Transforms the coordinates according to the symmetry operation Nop.
Params: `realtype & x` — Output parameter for the x-coordinate. `realtype & y` — Output parameter for the y-coordinate. `realtype & z` — Output parameter for the z-coordinate. `int Nop` — The index of the symmetry operation.
Returns: `int` — Returns `SYMOP_Ok` on success, or a non-zero value if the operation is invalid.
Side effects: Transforms the coordinates.
Errors: Returns a non-zero value if the operation number is invalid.

#### `SymOps::GetTMatrix(mat44 & TMatrix, int Nop) -> int`
// *** USED IN COOT ***
Retrieves the transformation matrix for the symmetry operation Nop.
Params: `mat44 & TMatrix` — Output parameter for the transformation matrix. `int Nop` — The index of the symmetry operation.
Returns: `int` — Returns `SYMOP_Ok` on success, or a non-zero value if the operation is invalid.
Side effects: Populates the transformation matrix.
Errors: Returns a non-zero value if the operation number is invalid.

#### `SymOps::Print()`
// *** USED IN COOT ***
Prints the SymOps object's information to the standard output stream.
Params: None
Returns: None
Side effects: Outputs the object's state to the console.

#### `SymOps::Copy(PSymOps symOps)`
// *** USED IN COOT ***
Copies the SymOps object's state to another SymOps object.
Params: `PSymOps symOps` — The SymOps object to copy to.
Returns: None
Side effects: Copies the state of the current object to the destination object.

#### `SymOps::write(io::RFile f)`
// *** USED IN COOT ***
Writes the SymOps object's data to the given RFile object.
Params: `io::RFile f` — The RFile object to write to.
Returns: None
Side effects: Serializes the object's state to the file.

#### `SymOps::read(io::RFile f)`
// *** USED IN COOT ***
Reads the SymOps object's data from the given RFile object.
Params: `io::RFile f` — The RFile object to read from.
Returns: None
Side effects: Deserializes the object's state from the file.


---

## `Title`
Represents the title section of a PDB file, containing information about the structure.

#### `Title::GetJournal() -> pstr`
Returns the journal records associated with the structure.
Params: None
Returns: A pointer to a `TitleContainer` object containing the journal records.
Side effects: None.
Errors: None.


---

## `TitleContainer`
This class represents a title within an MMDB file.

#### `TitleContainer::TitleContainer()`
Initializes a `TitleContainer` object.

#### `TitleContainer::TitleContainer(io::RPStream Object)`
Initializes a `TitleContainer` object from an `RPStream`.

#### `TitleContainer::~TitleContainer()`
Destroys the `TitleContainer` object.

#### `TitleContainer::MakeContainerClass(int ClassID)` // *** USED IN COOT ***
Creates a `PContainerClass` object representing the title container, using the provided `ClassID`.


---

## `Vertex`
Represents a single vertex in the MMDB structure, storing information about atoms and bonds.

#### `Vertex::GetName() -> cpstr`
Returns the name of the vertex.
Params: None
Returns: The vertex name as a C-string.
Side effects: None
Errors: None

#### `Vertex::GetUserID() -> int`
Returns the user-defined ID of the vertex.
Params: None
Returns: The user-defined vertex ID.
Side effects: None
Errors: None

#### `Vertex::GetType() -> int`
Returns the type of the vertex.
Params: None
Returns: The vertex type.
Side effects: None
Errors: None

#### `Vertex::SetVertex(cpstr)`
Sets the vertex name.
Params: `cpstr` — The vertex name.
Returns: None
Side effects: Updates the `name` member.
Errors: None

#### `Vertex::SetVertex(int, cpstr)`
Sets the vertex name and type.
Params: `int` — The vertex type, `cpstr` — The vertex name.
Returns: None
Side effects: Updates the `name` and `type` members.
Errors: None

#### `Vertex::SetVertex(int)`
Sets the vertex type.
Params: `int` — The vertex type.
Returns: None
Side effects: Updates the `type` member.
Errors: None

#### `Vertex::SetType(int)`
Sets the vertex type.
Params: `int` — The vertex type.
Returns: None
Side effects: Updates the `type` member.
Errors: None

#### `Vertex::SetTypeExt(int)`
Sets the vertex type extension.
Params: `int` — The vertex type extension.
Returns: None
Side effects: Updates the `type_ext` member.
Errors: None

#### `Vertex::RemoveChirality()`
Removes chirality information from the vertex.
Params: None
Returns: None
Side effects: Clears chirality-related data.
Errors: None

#### `Vertex::LeaveChirality(int)`
Leaves chirality information for a specific vertex type.
Params: `int` — The vertex type.
Returns: None
Side effects: Sets chirality-related flags for the given vertex type.
Errors: None

#### `Vertex::SetName(cpstr)`
Sets the vertex name.
Params: `cpstr` — The vertex name.
Returns: None
Side effects: Updates the `name` member.
Errors: None

#### `Vertex::SetProperty(int)`
Sets a property flag for the vertex.
Params: `int` — The property flag.
Returns: None
Side effects: Updates the `property` member.
Errors: None

#### `Vertex::SetID(int)`
Sets the graph-defined ID of the vertex.
Params: `int` — The vertex ID.
Returns: None
Side effects: Updates the `id` member.
Errors: None

#### `Vertex::AddBond()`
Adds a new bond to the vertex.
Params: None
Returns: None
Side effects: Adds a bond (implementation detail).
Errors: None

#### `Vertex::CopyNBonds(PVertex)`
Copies the bonds from another vertex.
Params: `PVertex` — The source vertex.
Returns: None
Side effects: Copies bonds from the source vertex to this vertex.
Errors: None

#### `Vertex::SetUserID(int)`
Sets the user-defined ID of the vertex.
Params: `int` — The vertex ID.
Returns: None
Side effects: Updates the `user_id` member.
Errors: None

#### `Vertex::GetProperty() -> int`
Returns the property flag for the vertex.
Params: None
Returns: The property flag.
Side effects: None
Errors: None

#### `Vertex::GetID() -> int`
Returns the graph-defined ID of the vertex.
Params: None
Returns: The vertex ID.
Side effects: None
Errors: None

#### `Vertex::GetUserID() -> int`
Returns the user-defined ID of the vertex.
Params: None
Returns: The user-defined vertex ID.
Side effects: None
Errors: None

#### `Vertex::GetName() -> cpstr`
Returns the name of the vertex.
Params: None
Returns: The vertex name as a C-string.
Side effects: None
Errors: None

#### `Vertex::GetType() -> int`
Returns the type of the vertex.
Params: None
Returns: The vertex type.
Side effects: None
Errors: None

#### `Vertex::GetTypeExt() -> int`
Returns the vertex type extension.
Params: None
Returns: The vertex type extension.
Side effects: None
Errors: None

#### `Vertex::GetNBonds() -> int`
Returns the number of bonds attached to the vertex.
Params: None
Returns: The number of bonds.
Side effects: None
Errors: None

#### `Vertex::SaveType() -> void`
Saves the vertex type to the user ID.
Params: None
Returns: None
Side effects: Saves the vertex type to the user ID.
Errors: None

#### `Vertex::RestoreType() -> void`
Restores the vertex type from the user ID.
Params: None
Returns: None
Side effects: Restores the vertex type from the user ID.
Errors: None

#### `Vertex::CopyType(PVertex) -> void`
Copies the vertex type from another vertex.
Params: `PVertex` — The source vertex.
Returns: None
Side effects: Copies the vertex type from the source vertex to this vertex.
Errors: None

#### `Vertex::Print(int)`
Prints vertex information to the standard output.
Params: `int` — The PKey.
Returns: None
Side effects: Prints vertex information.
Errors: None

#### `Vertex::Copy(PVertex) -> void`
Copies vertex data from another vertex.
Params: `PVertex` — The source vertex.
Returns: None
Side effects: Copies vertex data from the source vertex.
Errors: None

#### `Vertex::read(io::RFile)`
Reads vertex data from a file.
Params: `io::RFile` — The file stream.
Returns: None
Side effects: Reads data from the file.
Errors: None

#### `Vertex::write(io::RFile)`
Writes vertex data to a file.
Params: `io::RFile` — The file stream.
Returns: None
Side effects: Writes data to the file.
Errors: None

#### `Vertex::mem_read(cpstr, int &)`
Reads data from a memory string.
Params: `cpstr` — The memory string, `int &` — The length of the string.
Returns: None
Side effects: Reads data from the memory string.
Errors: None

#### `Vertex::mem_write(pstr, int &)`
Writes data to a memory string.
Params: `pstr` — The memory string, `int &` — The length of the string.
Returns: None
Side effects: Writes data to the memory string.
Errors: None


---
