# MMDB2 API Reference (AST-driven)

---

## API Index

This is a comprehensive list of methods, properties, or attributes associated with various molecular modeling components, generally representing structure or chemical feature.

Since you did not ask a specific question, I will provide a summary of the context and suggest possible actions based on the provided data dump.

### Summary of Context
The provided set of components (`Atom`, `Bond`, `Residue`, `Chain`, `Protein`, `Ligand`, `System`) strongly suggests an object-oriented framework for handling chemical structures, likely in the context of biomolecular simulation or analysis.

## Per-Class Method Documentation

## `Alignment`
Manages and calculates the scoring matrix and alignment path for two input sequences.

#### `Alignment::Alignment() -> void`
Constructs a default, empty alignment object.
Receiver: instance — call as `obj->Alignment()`
Side effects: Initializes internal state variables.

#### `Alignment::Alignment(io::RPStream Object) -> void`
Constructs an alignment object from an input stream containing alignment data.
Receiver: instance — call as `obj->Alignment(io::RPStream Object)`
Side effects: Reads and parses alignment data from the provided stream.

#### `Alignment::~Alignment() -> void`
Destroys the alignment object, freeing allocated resources.
Receiver: instance — call as `obj->~Alignment()`
Side effects: Deallocates memory associated with alignment structures.

#### `Alignment::SetAffineModel(realtype WGap, realtype WSpace) -> void`
Sets the gap penalty weights for use in subsequent alignment calculations.
Receiver: instance — call as `obj->SetAffineModel(realtype WGap, realtype WSpace)`
Params: `WGap` — The weight applied to gap penalties.
Params: `WSpace` — The weight applied to spacing penalties.
Side effects: Updates internal gap penalty weights ($\text{Vg}$ and $\text{Ws}$).

#### `Alignment::Align(cpstr S, cpstr T, ALIGN_METHOD Method=ALIGN_GLOBAL) -> void`
Performs the sequence alignment between sequences S and T using a specified global or local alignment method.
Receiver: instance — call as `obj->Align(cpstr S, cpstr T, ALIGN_METHOD Method)`
Params: `S` — The identifier/string for the first sequence.
Params: `T` — The identifier/string for the second sequence.
Params: `Method` — The type of alignment strategy to use (default is global).
Side effects: Calculates and stores the full alignment score matrix, path, and final score.

#### `Alignment::GetAlignedS() -> pstr`
Retrieves the sequence string representing the aligned residues for the first sequence (S).
Receiver: instance — call as `obj->GetAlignedS()`
Returns: A pointer to the string containing the residues of the first sequence after alignment.

#### `Alignment::GetAlignedT() -> pstr`
Retrieves the sequence string representing the aligned residues for the second sequence (T).
Receiver: instance — call as `obj->GetAlignedT()`
Returns: A pointer to the string containing the residues of the second sequence after alignment.

#### `Alignment::GetScore() -> realtype`
Retrieves the maximum score achieved for the computed optimal alignment path.
Receiver: instance — call as `obj->GetScore()`
Returns: The optimal alignment score as a floating-point number.

#### `Alignment::GetNAlign() -> int`
Retrieves the total number of residues that were aligned across both sequences.
Receiver: instance — call as `obj->GetNAlign()`
Returns: The count of aligned positions in the resulting alignment.

#### `Alignment::InitAlignment() -> void`
Initializes all internal data structures and state variables required for alignment calculations.
Receiver: instance — call as `obj->InitAlignment()`
Side effects: Resets alignment parameters and data members.

#### `Alignment::SetScores(realtype SEqual, realtype SNEqual) -> void`
Sets the scoring parameters used for calculating residue matches and mismatches.
Receiver: instance — call as `obj->SetScores(realtype SEqual, realtype SNEqual)`
Params: `SEqual` — The score for identical residues.
Params: `SNEqual` — The score for non-identical, but non-gap, residues.
Side effects: Updates internal scoring parameters ($\text{SEq}$ and $\text{SNEq}$).

---

## `Atom`
Represents an individual atomic structure within a molecular model, holding coordinates, chemical identity, and connectivity information.

#### `Atom::SetResidue(PResidue res) -> void`
Sets the residue object to which this atom belongs.
Receiver: instance — call as `obj->SetResidue(...)`
Params: `res` — pointer to the residue to set.
Side effects: Updates the atom's associated residue.

#### `Atom::AddBond(PAtom bond_atom, int bond_order, int nAdd_bonds=1) -> int`
Adds a covalent bond connecting this atom to another specified atom.
Receiver: instance — call as `obj->AddBond(...)`
Params: `bond_atom` — pointer to the atom to bond to; `bond_order` — the bond order; `nAdd_bonds` — initial estimate of new bonds.
Returns: int — the current number of bonds (if > 0), or an error code (if $\le 0$).
Side effects: Modifies the atom's internal list of bonds.

#### `Atom::GetNBonds() -> int`
Retrieves the total number of covalent bonds attached to this atom.
Receiver: instance — call as `obj->GetNBonds()`
Returns: int — the number of bonds.

#### `Atom::GetBonds(RPAtomBond atomBond, int & nAtomBonds) -> void`
Returns a reference to the atom's internal bond structure array.
Receiver: instance — call as `obj->GetBonds(...)`
Params: `atomBond` — reference to the bond structure array; `nAtomBonds` — reference to store the count of bonds found.
Side effects: Fills the provided bond structure array and sets the bond count.

#### `Atom::GetBonds(RPAtomBondI atomBondI, int & nAtomBonds) -> void`
Returns a reference to the atom's internal bond index structure array.
Receiver: instance — call as `obj->GetBonds(...)`
Params: `atomBondI` — reference to the bond index structure array; `nAtomBonds` — reference to store the count of bonds found.
Side effects: Fills the provided bond index structure array and sets the bond count.

#### `Atom::GetBonds(PAtomBondI atomBondI, int & nAtomBonds, int maxlength) -> void`
Fills a user-provided bond index structure array without allocating or disposing memory.
Receiver: instance — call as `obj->GetBonds(...)`
Params: `atomBondI` — pointer to the user-provided bond index array; `nAtomBonds` — reference to store the count of bonds found; `maxlength` — maximum size of the provided array.
Side effects: Writes bond index data to the provided pointer and sets the bond count.

#### `Atom::SetAtomName(int ix, int sN, const AtomName aName, const AltLoc aLoc, const SegID sID, const Element eName) -> void`
Sets the comprehensive name, segment, and identification details for the atom.
Receiver: instance — call as `obj->SetAtomName(...)`
Params: `ix` — atom index; `sN` — serial number; `aName` — atom name; `aLoc` — alternative location; `sID` — segment ID; `eName` — element name.
Side effects: Overwrites the atom's name, segment, and element information.

#### `Atom::SetAtomName(const AtomName atomName) -> void`
Renames the atom using only the new atom name string.
Receiver: instance — call as `obj->SetAtomName(...)`
Params: `atomName` — the new atom name.
Side effects: Updates the atom's name.

#### `Atom::SetElementName(const Element elName) -> void`
Sets the chemical element symbol for the atom.
Receiver: instance — call as `obj->SetElementName(...)`
Params: `elName` — the element to set.
Side effects: Updates the atom's element symbol.

#### `Atom::MakeTer() -> void`
Marks this atom as a chain terminator.
Receiver: instance — call as `obj->MakeTer()`
Side effects: Sets the `Ter` flag.

#### `Atom::SetCoordinates(realtype xx, realtype yy, realtype zz, realtype occ, realtype tFac) -> void`
Sets the core structural parameters: coordinates, occupancy, and temperature factor.
Receiver: instance — call as `obj->SetCoordinates(...)`
Params: `xx, yy, zz` — coordinates in angstroms; `occ` — occupancy; `tFac` — temperature factor.
Side effects: Overwrites the atom's coordinates and associated metadata.

#### `Atom::GetModelNum() -> int`
Retrieves the model number associated with this atom.
Receiver: instance — call as `obj->GetModelNum()`
Returns: int — the model number.

#### `Atom::GetChainID() -> pstr`
Retrieves the chain identifier associated with this atom.
Receiver: instance — call as `obj->GetChainID()`
Returns: pstr — a string representing the chain ID.

#### `Atom::GetResName() -> pstr`
Retrieves the three-letter residue name associated with this atom.
Receiver: instance — call as `obj->GetResName()`
Returns: pstr — the residue name string.

#### `Atom::GetSeqNum() -> int`
Retrieves the sequence number of the residue containing this atom.
Receiver: instance — call as `obj->GetSeqNum()`
Returns: int — the residue sequence number.

#### `Atom::GetInsCode() -> pstr`
Retrieves the residue insertion code associated with this atom.
Receiver: instance — call as `obj->GetInsCode()`
Returns: pstr — the insertion code string.

#### `Atom::GetAtomName() -> pstr`
Retrieves the canonical atom name assigned to this atom.
Receiver: instance — call as `obj->GetAtomName()`
Returns: pstr — the atom name string.

#### `Atom::GetElementName() -> pstr`
Retrieves the element symbol associated with this atom.
Receiver: instance — call as `obj->GetElementName()`
Returns: pstr — the element symbol string.

#### `Atom::GetChainCalphas(PPAtom & Calphas, int & nCalphas, cpstr altLoc = "") -> void`
Populates an array with pointers to C-alpha atoms belonging to the same chain.
Receiver: instance — call as `obj->GetChainCalphas(...)`
Params: `Calphas` — pointer to an array to receive C-alpha pointers; `nCalphas` — reference to store the count of C-alphas; `altLoc` — optional alternative location indicator.
Side effects: Allocates and populates the `Calphas` array, which the caller must delete.

#### `Atom::isTer() -> bool`
Checks if this atom is designated as a chain terminator.
Receiver: instance — call as `obj->isTer()`
Returns: bool — true if the atom is a terminator.

#### `Atom::isMetal() -> bool`
Checks if this atom is designated as a metal atom.
Receiver: instance — call as `obj->isMetal()`
Returns: bool — true if the atom is a metal.

#### `Atom::isInSelection(int selHnd) -> bool`
Checks if this atom is currently included in a specified selection.
Receiver: instance — call as `obj->isInSelection(...)`
Params: `selHnd` — the selection handle ID.
Returns: bool — true if the atom belongs to the selection.

#### `Atom::GetResidue() -> PResidue`
Returns a pointer to the residue object containing this atom.
Receiver: instance — call as `obj->GetResidue()`
Returns: PResidue — pointer to the containing residue.

#### `Atom::GetChain() -> PChain`
Returns a pointer to the chain object containing this atom.
Receiver: instance — call as `obj->GetChain()`
Returns: PChain — pointer to the containing chain.

#### `Atom::GetModel() -> PModel`
Returns a pointer to the model object containing this atom.
Receiver: instance — call as `obj->GetModel()`
Returns: PModel — pointer to the containing model.

#### `Atom::GetResidueNo() -> int`
Retrieves the sequence number identifier of the residue.
Receiver: instance — call as `obj->GetResidueNo()`
Returns: int — the residue sequence number.

#### `Atom::GetAtomID(pstr AtomID) -> pstr`
Generates a unique, formatted identifier string for the atom based on its structural context.
Receiver: instance — call as `obj->GetAtomID(...)`
Params: `AtomID` — an optional initial string ID for formatting context.
Returns: pstr — the generated formatted atom ID string.

#### `Atom::GetAtomIDfmt(pstr AtomID) -> pstr`
Generates a formatted version of the atom ID string.
Receiver: instance — call as `obj->GetAtomIDfmt(...)`
Params: `AtomID` — the input string to format.
Returns: pstr — the formatted atom ID string.

#### `Atom::Transform(const mat33 & tm, vect3 & v) -> void`
Transforms the coordinates of the atom using a 3x3 rotation/scaling matrix.
Receiver: instance — call as `obj->Transform(...)`
Params: `tm` — the 3x3 transformation matrix; `v` — reference vector to store the transformed coordinates.
Side effects: Updates the atom's internal coordinates.

#### `Atom::Transform(const mat44 & tm) -> void`
Transforms the coordinates of the atom using a 4x4 transformation matrix.
Receiver: instance — call as `obj->Transform(...)`
Params: `tm` — the 4x4 transformation matrix.
Side effects: Updates the atom's internal coordinates.

#### `Atom::PutUDData(int UDDhandle, int iudd) -> int`
Stores an integer value into the atom's user-defined data block.
Receiver: instance — call as `obj->PutUDData(...)`
Params: `UDDhandle` — handle for the data block; `iudd` — integer value to store.
Returns: int — status code.
Side effects: Mutates the atom's internal UDData storage.

#### `Atom::PutUDData(int UDDhandle, realtype rudd) -> int`
Stores a real-valued floating-point number into the atom's user-defined data block.
Receiver: instance — call as `obj->PutUDData(...)`
Params: `UDDhandle` — handle for the data block; `rudd` — real value to store.
Returns: int — status code.
Side effects: Mutates the atom's internal UDData storage.

#### `Atom::PutUDData(int UDDhandle, cpstr sudd) -> int`
Stores a character string into the atom's user-defined data block.
Receiver: instance — call as `obj->PutUDData(...)`
Params: `UDDhandle` — handle for the data block; `sudd` — string value to store.
Returns: int — status code.
Side effects: Mutates the atom's internal UDData storage.

#### `Atom::GetUDData(int UDDhandle, int & iudd) -> int`
Retrieves an integer value from the atom's user-defined data block.
Receiver: instance — call as `obj->GetUDData(...)`
Params: `UDDhandle` — handle for the data block; `iudd` — reference to store the retrieved integer.
Returns: int — status code.
Side effects: None.

#### `Atom::GetUDData(int UDDhandle, realtype & rudd) -> int`
Retrieves a real-valued floating-point number from the atom's user-defined data block.
Receiver: instance — call as `obj->GetUDData(...)`
Params: `UDDhandle` — handle for the data block; `rudd` — reference to store the retrieved real value.
Returns: int — status code.
Side effects: None.

#### `Atom::GetUDData(int UDDhandle, pstr sudd, int maxLen) -> int`
Retrieves a string value from the atom's user-defined data block into a provided buffer.
Receiver: instance — call as `obj->GetUDData(...)`
Params: `UDDhandle` — handle for the data block; `sudd` — output string buffer; `maxLen` — maximum length of the output string.
Returns: int — status code.
Side effects: None.

#### `Atom::GetUDData(int UDDhandle, pstr & sudd) -> int`
Retrieves a string value from the atom's user-defined data block into an output string reference.
Receiver: instance — call as `obj->GetUDData(...)`
Params: `UDDhandle` — handle for the data block; `sudd` — reference to store the retrieved string.
Returns: int — status code.
Side effects: None.

#### `Atom::GetIndex() -> int`
Retrieves the internal file index assigned to the atom.
Receiver: instance — call as `obj->GetIndex()`
Returns: int — the internal index.

#### `Atom::Copy(PAtom atom) -> void`
Creates a deep copy of the atom's structural information from another atom.
Receiver: instance — call as `obj->Copy(...)`
Params: `atom` — pointer to the source atom to copy from.
Side effects: Copies coordinates, IDs, and attributes from the source atom.

---

## `Author`
Represents the author information associated with a structure file.

#### `Author::Author()`
Constructs a default, empty author object.
Receiver: instance — call as `obj->Author()`
Returns: None

#### `Author::Author(cpstr S)`
Constructs an author object from a string content.
Receiver: instance — call as `obj->Author(cpstr S)`
Params: `S` — the author string content.
Returns: None

#### `Author::Author(io::RPStream Object)`
Constructs an author object by reading from an input stream.
Receiver: instance — call as `obj->Author(io::RPStream Object)`
Params: `Object` — the readable stream containing the author data.
Returns: None

#### `Author::~Author()`
Destroys the author object.
Receiver: instance — call as `obj->~Author()`

#### `Author::PDBASCIIDump(pstr S, int N)`
Writes the author information into the specified ASCII string buffer.
Receiver: instance — call as `obj->PDBASCIIDump(pstr S, int N)`
Params: `S` — the buffer to write to.
`N` — the maximum size of the buffer.
Returns: None
Side effects: Populates the string buffer `S` with author details.

#### `Author::PDBASCIIDump1(io::RFile)`
Writes the author information to a specified file stream.
Receiver: instance — call as `obj->PDBASCIIDump1(io::RFile)`
Params: None
Returns: bool — True if the dump was successful, otherwise false.
Side effects: Writes author details to the provided file.

#### `Author::ConvertPDBASCII(cpstr S)`
Parses and converts author information found within a PDB ASCII string.
Receiver: instance — call as `obj->ConvertPDBASCII(cpstr S)`
Params: `S` — the PDB ASCII string to parse.
Returns: ERROR_CODE — Status code indicating success or failure of conversion.
Side effects: Updates the internal state of the object with parsed author details.

#### `Author::GetClassID()`
Retrieves the unique identifier for the Author class.
Receiver: instance — call as `obj->GetClassID()`
Returns: CLASS_ID — The constant ID representing the Author class.

#### `Author::write(io::RFile f)`
Writes the author object's data contents to a specified file stream.
Receiver: instance — call as `obj->write(io::RFile f)`
Params: `f` — the file stream to write to.
Returns: None
Side effects: Writes object data to the file stream `f`.

#### `Author::read(io::RFile f)`
Reads author object data from a specified file stream.
Receiver: instance — call as `obj->read(io::RFile f)`
Params: `f` — the file stream to read from.
Returns: None
Side effects: Updates the object's state from the data read from the file stream `f`.

---

## `BondManager`
Manages the creation and removal of covalent bonds within a molecular structure.

#### `BondManager::MakeBonds(bool calc_only)`
Creates or regenerates all covalent bonds present in the structure.
Receiver: instance — call as `obj->MakeBonds(...)` (or `obj.MakeBonds(...)`).
Params: `calc_only` — If true, only calculates bonds without modifying connectivity.
Returns: void.
Side effects: Updates the internal bond list and potentially the structure connectivity.

#### `BondManager::RemoveBonds()`
Clears all existing covalent bonds from the structure.
Receiver: instance — call as `obj->RemoveBonds()` (or `obj.RemoveBonds()`).
Returns: void.
Side effects: Removes all stored bond information, leaving the residue connectivity intact.

---

## `Category`
Represents a structural category within an mmCIF file, holding associated tags.

#### `Category::GetCategoryName() -> pstr`
Returns the name assigned to this category.
Receiver: instance — call as `obj->GetCategoryName()`
Returns: A pointer to the character string representing the category name, or NULL if no name was set.
Side effects: none.

#### `Category::GetTagNo(cpstr ttag) -> int`
Returns the zero-based serial number of a specified tag within the category.
Receiver: instance — call as `obj->GetTagNo(ttag)`
Params: `ttag` — The tag name (or field name) to check.
Returns: The tag's position index (>=0), or a negative number indicating insertion point (<0).
Side effects: none.

#### `Category::GetNofTags() -> int`
Returns the total count of tags currently stored in this category.
Receiver: instance — call as `obj->GetNofTags()`
Returns: The number of tags in the category.
Side effects: none.

#### `Category::GetTag(int tagNo) -> pstr`
Retrieves the tag string at a specified zero-based index.
Receiver: instance — call as `obj->GetTag(tagNo)`
Params: `tagNo` — The zero-based serial number of the tag to retrieve.
Returns: A pointer to the tag string at the given index, or NULL if the index is out of range.
Side effects: none.

#### `Category::CheckTags(cpstr * tagList) -> bool`
Checks if all provided tags are present within this category.
Receiver: instance — call as `obj->CheckTags(tagList)`
Params: `tagList` — A null-terminated array of tag names to check for presence.
Returns: True if every tag in the list is found in the category, false otherwise.
Side effects: none.

---

## `Chain`
Manages the collection of residues belonging to a single polypeptide chain within a molecular structure model.

#### `Chain::SetChain(const ChainID chID)`
Sets the chain ID for this chain object.
Receiver: instance — call as `obj->SetChain(...)`
Params: `chID` — The desired ChainID.
Returns: void
Side effects: Mutates `chainID` to reflect the new chain identifier.

#### `Chain::GetCoordHierarchy()`
Retrieves the coordinate hierarchy associated with this chain object.
Receiver: instance — call as `obj->GetCoordHierarchy()`
Params: none
Returns: PManager — A pointer to the coordinate hierarchy manager.
Side effects: None.

#### `Chain::GetNumberOfResidues()`
Returns the total count of residues currently stored in the chain.
Receiver: instance — call as `obj->GetNumberOfResidues()`
Params: none
Returns: int — The number of residues in the chain.
Side effects: None.

#### `Chain::GetResidue(int resNo)`
Retrieves the residue at the specified zero-based index within the chain's residue table.
Receiver: instance — call as `obj->GetResidue(...)`
Params: `resNo` — The zero-based residue number in the coordinate file.
Returns: PResidue — Pointer to the requested residue, or null if out of bounds.
Side effects: None.

#### `Chain::GetResidue(int seqNum, const InsCode insCode)`
Retrieves the residue whose sequence number and insertion code match the provided values.
Receiver: instance — call as `obj->GetResidue(...)`
Params: `seqNum` — The sequence number of the residue. `insCode` — The insertion code of the residue.
Returns: PResidue — Pointer to the requested residue, or null if absent.
Side effects: None.

#### `Chain::GetResidueNo(int seqNum, const InsCode insCode)`
Finds the index of a residue within the chain's internal residue table using its sequence identity.
Receiver: instance — call as `obj->GetResidueNo(...)`
Params: `seqNum` — The sequence number of the residue. `insCode` — The insertion code of the residue.
Returns: int — The zero-based index in the residue table (0..nres-1), or -1 if not found.
Side effects: None.

#### `Chain::GetResidueTable(PPResidue & resTable, int & NumberOfResidues)`
Populates the provided residue array pointer with all residues in the chain and sets the count.
Receiver: instance — call as `obj->GetResidueTable(...)`
Params: `resTable` — A pointer-to-pointer to the residue array that will receive the residues. `NumberOfResidues` — Reference to an integer that will receive the total count.
Returns: void
Side effects: Allocates and populates `*resTable` and sets `NumberOfResidues`.

#### `Chain::GetResidueCreate(const ResName resName, int seqNum, const InsCode insCode, bool Enforce)`
Retrieves an existing residue or creates a new one at the end of the chain matching the provided sequence identity and name.
Receiver: instance — call as `obj->GetResidueCreate(...)`
Params: `resName` — The desired name for the residue. `seqNum` — The sequence number for the residue. `insCode` — The insertion code for the residue. `Enforce` — If true, forces creation even if a mismatch exists.
Returns: PResidue — Pointer to the requested or newly created residue.
Side effects: Potentially adds a new residue to the chain structure.

#### `Chain::AddResidue(PResidue res)`
Adds a residue to the chain, shifting subsequent residues if the position requires it.
Receiver: instance — call as `obj->AddResidue(...)`
Params: `res` — The residue to be added to the chain.
Returns: int — The number of residues added (usually 1).
Side effects: Appends or inserts the residue into the chain's internal residue list.

#### `Chain::InsResidue(PResidue res, int pos)`
Inserts a residue at a specific zero-based position within the chain's residue table.
Receiver: instance — call as `obj->InsResidue(...)`
Params: `res` — The residue to be inserted. `pos` — The zero-based index where the residue should be placed.
Returns: int — The number of residues added (usually 1).
Side effects: Inserts the residue at the specified position, shifting residues after it.

#### `Chain::InsResidue(PResidue res, int seqNum, const InsCode insCode)`
Inserts a residue matching the given sequence identity at the position immediately preceding an existing residue.
Receiver: instance — call as `obj->InsResidue(...)`
Params: `res` — The residue to be inserted. `seqNum` — The sequence number of the residue before insertion. `insCode` — The insertion code of the residue before insertion.
Returns: int — The number of residues added (usually 1).
Side effects: Inserts the residue and updates sequence/coordinate links.

#### `Chain::DeleteResidue(int resNo)`
Removes the residue at the specified zero-based index from the chain.
Receiver: instance — call as `obj->DeleteResidue(...)`
Params: `resNo` — The zero-based residue number to delete.
Returns: int — The number of residues deleted (usually 1).
Side effects: Removes the residue from the chain; subsequent residues shift indices.

#### `Chain::DeleteResidue(int seqNum, const InsCode insCode)`
Removes a residue identified by its sequence number and insertion code.
Receiver: instance — call as `obj->DeleteResidue(...)`
Params: `seqNum` — The sequence number of the residue to delete. `insCode` — The insertion code of the residue to delete.
Returns: int — The number of residues deleted (usually 1).
Side effects: Removes the residue; subsequently shifts indices and updates tables.

#### `Chain::DeleteAllResidues()`
Removes all residues currently associated with this chain.
Receiver: instance — call as `obj->DeleteAllResidues()`
Params: none
Returns: int — The total number of residues deleted.
Side effects: Clears all residue data associated with the chain.

#### `Chain::TrimResidueTable()`
Cleans up the internal residue table after deletions, ensuring contiguous storage.
Receiver: instance — call as `obj->TrimResidueTable()`
Params: none
Returns: void
Side effects: Reallocates or marks the residue array to reflect only occupied space.

#### `Chain::Copy(PChain chain)`
Deep copies the contents of another chain object onto this chain.
Receiver: instance — call as `obj->Copy(...)`
Params: `chain` — The source chain object to copy from.
Returns: void
Side effects: Copies all residues and associated data from the source chain to this object.

#### `Chain::CopyAnnotations(PChain chain)`
Copies any annotations present in a source chain onto this chain object.
Receiver: instance — call as `obj->CopyAnnotations(...)`
Params: `chain` — The source chain object containing annotations to copy.
Returns: void
Side effects: Copies metadata/annotations from the source chain to this object.

#### `Chain::GetChainID()`
Retrieves the unique identifier assigned to this chain.
Receiver: instance — call as `obj->GetChainID()`
Params: none
Returns: cpstr — A string representation of the chain ID.
Side effects: None.

#### `Chain::SetChainID(const ChainID chID)`
Explicitly sets or updates the chain identifier associated with this object.
Receiver: instance — call as `obj->SetChainID(...)`
Params: `chID` — The new ChainID to assign.
Returns: void
Side effects: Mutates the internal `chainID`.

#### `Chain::GetChainID(pstr ChID)`
Retrieves the chain identifier represented by a given string.
Receiver: instance — call as `obj->GetChainID(...)`
Params: `ChID` — A string potentially containing the chain ID.
Returns: cpstr — The determined chain ID string.
Side effects: None.

#### `Chain::GetNumberOfDBRefs()`
Returns the count of unique database references associated with this chain.
Receiver: instance — call as `obj->GetNumberOfDBRefs()`
Params: none
Returns: int — The total number of stored database references.
Side effects: None.

#### `Chain::GetDBRef(int dbRefNo)`
Retrieves the database reference object at the specified zero-based index.
Receiver: instance — call as `obj->GetDBRef(...)`
Params: `dbRefNo` — The zero-based index of the database reference.
Returns: PDBReference — Pointer to the specified database reference.
Side effects: None.

#### `Chain::SortResidues()`
Reorders the residues within the chain according to established numbering conventions.
Receiver: instance — call as `obj->SortResidues()`
Params: none
Returns: void
Side effects: Reorders the internal residue array and updates associated indices.

#### `Chain::isSolventChain()`
Checks if this chain is classified as a solvent chain (e.g., water).
Receiver: instance — call as `obj->isSolventChain()`
Params: none
Returns: bool — True if the chain is solvent, false otherwise.
Side effects: None.

#### `Chain::isInSelection(int selHnd)`
Determines if the residues in this chain are currently marked as being part of a specific selection object.
Receiver: instance — call as `obj->isInSelection(...)`
Params: `selHnd` — The handle ID of the selection to check against.
Returns: bool — True if the chain is in the selection, false otherwise.
Side effects: None.

#### `Chain::isAminoacidChain()`
Checks if this chain consists exclusively of amino acid residues.
Receiver: instance — call as `obj->isAminoacidChain()`
Params: none
Returns: bool — True if the chain is composed of amino acids, false otherwise.
Side effects: None.

#### `Chain::isNucleotideChain()`
Checks if this chain consists exclusively of nucleotide residues.
Receiver: instance — call as `obj->isNucleotideChain()`
Params: none
Returns: bool — True if the chain is composed of nucleotides, false otherwise.
Side effects: None.

---

## `CisPep`
Represents a cis-peptide bond measurement structure element in a biological context.

#### `CisPep::PDBASCIIDump(pstr S, int N) -> void`
Writes the structural information of the cis-peptide bond to a PDB ASCII stream.
Receiver: instance — call as `obj->PDBASCIIDump(...)`
Params: `S` — The output string stream. `N` — The number of records to write.
Side effects: Writes structural data to the stream `S`.

#### `CisPep::ConvertPDBASCII(cpstr S) -> ERROR_CODE`
Parses and populates the cis-peptide structure object from a given PDB ASCII string.
Receiver: instance — call as `obj->ConvertPDBASCII(...)`
Params: `S` — The input PDB ASCII string.
Returns: An error code indicating success or failure.
Side effects: Populates the object's internal state from the string `S`.

#### `CisPep::GetClassID() -> CLASS_ID`
Retrieves the unique class identifier for a `CisPep` object.
Receiver: instance — call as `obj->GetClassID()`
Returns: The constant `ClassID_CisPep`.

#### `CisPep::Copy(PContainerClass cisPep) -> void`
Deep copies the content of another `PContainerClass` into this object.
Receiver: instance — call as `obj->Copy(...)`
Params: `cisPep` — The source container object.
Side effects: Copies the structural data from `cisPep` into this object.

#### `CisPep::write(io::RFile f) -> void`
Writes the object's structural data content to a specified file stream.
Receiver: instance — call as `obj->write(...)`
Params: `f` — The output file stream.
Side effects: Writes structural data to the file stream `f`.

#### `CisPep::read(io::RFile f) -> void`
Reads the object's structural data content from a specified file stream.
Receiver: instance — call as `obj->read(...)`
Params: `f` — The input file stream.
Side effects: Reads structural data from the file stream `f` and populates the object.

---

## `ClassContainer`
Manages and holds structural data for a collection of related structural units.

#### `ClassContainer::GetContainerClass(int ContClassNo) -> PContainerClass`
Retrieves the structural class object associated with a given container index.
Receiver: instance — call as `obj->GetContainerClass(...)`
Params: `ContClassNo` — The zero-based index of the container class to retrieve.
Returns: A pointer to the requested `PContainerClass` object.
Errors: Returns null if the index is out of bounds.

#### `ClassContainer::Length() -> int`
Returns the total number of structural units stored in the container.
Receiver: instance — call as `obj->Length()`
Returns: The integer length of the container.

---

## `Compound`
Represents a molecular compound, holding structural data compatible with PDB formats.

#### `Compound::Compound() -> void`
Constructs a default, empty Compound object.
Receiver: instance — call as `obj->Compound()` (or `obj.Compound()`).

#### `Compound::Compound(cpstr S) -> void`
Constructs a Compound object from a string containing PDB data.
Receiver: instance — call as `obj->Compound(S)` (or `obj.Compound(S)`).
Params: `S` — the string containing PDB data.

#### `Compound::Compound(io::RPStream Object) -> void`
Constructs a Compound object by reading data from a readable stream.
Receiver: instance — call as `obj->Compound(Object)` (or `obj.Compound(Object)`).
Params: `Object` — the input readable stream.

#### `Compound::~Compound() -> void`
Destroys the Compound object.
Receiver: instance — call as `obj->~Compound()` (or `obj.~Compound()`).

#### `Compound::PDBASCIIDump(pstr S, int N) -> void`
Writes the compound's structure to the specified ASCII string buffer.
Receiver: instance — call as `obj->PDBASCIIDump(S, N)` (or `obj.PDBASCIIDump(S, N)`).
Params: `S` — the destination string buffer.
Params: `N` — the maximum size of the buffer.
Side effects: Writes structural data to the provided string buffer.

#### `Compound::PDBASCIIDump1(io::RFile) -> bool`
Writes the compound's structure to a specified file handle using PDB ASCII format.
Receiver: instance — call as `obj->PDBASCIIDump1(f)` (or `obj.PDBASCIIDump1(f)`).
Params: `f` — the output file stream.
Returns: true if writing was successful, false otherwise.
Side effects: Writes structural data to the output file stream.

#### `Compound::ConvertPDBASCII(cpstr S) -> ERROR_CODE`
Parses and converts structural data from a PDB ASCII string into the internal compound representation.
Receiver: instance — call as `obj->ConvertPDBASCII(S)` (or `obj.ConvertPDBASCII(S)`).
Params: `S` — the string containing PDB data.
Returns: An error code indicating success or failure.
Side effects: Updates the internal compound structure based on the input string.

#### `Compound::GetClassID() -> CLASS_ID`
Returns the unique class ID associated with the Compound structure.
Receiver: instance — call as `obj->GetClassID()` (or `obj.GetClassID()`).
Returns: The constant `ClassID_Compound`.

#### `Compound::write(io::RFile f) -> void`
Writes the compound's entire structure contents to the specified file stream.
Receiver: instance — call as `obj->write(f)` (or `obj.write(f)`).
Params: `f` — the output file stream.
Side effects: Writes structural data to the file stream.

#### `Compound::read(io::RFile f) -> void`
Reads structural data from the specified file stream to initialize or update the compound object.
Receiver: instance — call as `obj->read(f)` (or `obj.read(f)`).
Params: `f` — the input file stream.
Side effects: Reads structural data from the file stream, modifying the object state.

---

## `CoorManager`
Manages coordinate data, including atoms, residues, chains, and models, within a structural biology context.

#### `CoorManager::GetNumberOfModels() -> int`
Returns the total number of models currently loaded in the manager.
Receiver: instance — call as `obj->GetNumberOfModels()`
Returns: The count of models (nModels).

#### `CoorManager::GetModel(int modelNo) -> PModel`
Retrieves a model by its serial number.
Receiver: instance — call as `obj->GetModel(...)`
Params: `modelNo` — The serial number of the model (1 <= modelNo <= nModels).
Returns: A pointer to the requested `Model` object.
Errors: Returns null if the model number is out of range.

#### `CoorManager::GetModel(cpstr CID) -> PModel`
Retrieves a model corresponding to a specific coordinate identifier string.
Receiver: instance — call as `obj->GetModel(...)`
Params: `CID` — The coordinate identifier string (e.g., "/mdl/1").
Returns: A pointer to the requested `Model` object.
Errors: Returns null if no model matches the identifier.

#### `CoorManager::GetModelTable(PPModel & modTable, int & NumberOfModels) -> void`
Fills the provided structure table with pointers to all loaded models and reports the total count.
Receiver: instance — call as `obj->GetModelTable(...)`
Params: `modTable` — Pointer to the model table to be filled.
`NumberOfModels` — Reference to an integer receiving the count of models.
Side effects: Populates `modTable` with pointers to all models.

#### `CoorManager::GetNumberOfChains(int modelNo) -> int`
Counts the number of chains present within a specified model.
Receiver: instance — call as `obj->GetNumberOfChains(...)`
Params: `modelNo` — The serial number of the model to check.
Returns: The total number of chains in the model.

#### `CoorManager::GetNumberOfChains(cpstr CID) -> int`
Counts the number of chains available based on a provided coordinate identifier string.
Receiver: instance — call as `obj->GetNumberOfChains(...)`
Params: `CID` — The coordinate identifier string.
Returns: The total number of chains matching the identifier.

#### `CoorManager::GetChain(int modelNo, const ChainID chainID) -> PChain`
Retrieves a specific chain object using its model number and chain identifier.
Receiver: instance — call as `obj->GetChain(...)`
Params: `modelNo` — The model serial number.
`chainID` — The identifier of the desired chain.
Returns: A pointer to the requested `Chain` object.
Errors: Returns null if the model or chain does not exist.

#### `CoorManager::GetChain(int modelNo, int chainNo) -> PChain`
Retrieves a specific chain object using its model number and sequential chain number.
Receiver: instance — call as `obj->GetChain(...)`
Params: `modelNo` — The model serial number.
`chainNo` — The sequential number of the chain (0-indexed).
Returns: A pointer to the requested `Chain` object.
Errors: Returns null if the model or chain does not exist.

#### `CoorManager::GetChain(cpstr CID) -> PChain`
Retrieves a chain object matching a specified coordinate identifier string.
Receiver: instance — call as `obj->GetChain(...)`
Params: `CID` — The coordinate identifier string.
Returns: A pointer to the requested `Chain` object.
Errors: Returns null if no chain matches the identifier.

#### `CoorManager::GetChainTable(int modelNo, PPChain & chainTable, int & NumberOfChains) -> void`
Fills the provided structure table with pointers to all chains in a model and reports the total count.
Receiver: instance — call as `obj->GetChainTable(...)`
Params: `modelNo` — The model serial number.
`chainTable` — Pointer to the chain table to be filled.
`NumberOfChains` — Reference to an integer receiving the count of chains.
Side effects: Populates `chainTable` with pointers to all chains.

#### `CoorManager::GetChainTable(cpstr CID, PPChain & chainTable, int & NumberOfChains) -> void`
Fills the provided structure table with pointers to all chains matching an identifier and reports the total count.
Receiver: instance — call as `obj->GetChainTable(...)`
Params: `CID` — The coordinate identifier string.
`chainTable` — Pointer to the chain table to be filled.
`NumberOfChains` — Reference to an integer receiving the count of chains.
Side effects: Populates `chainTable` with pointers to all chains.

#### `CoorManager::GetResidue(int modelNo, const ChainID chainID, int seqNo, const InsCode insCode) -> PResidue`
Retrieves a residue object by model number, chain ID, sequence number, and insertion code.
Receiver: instance — call as `obj->GetResidue(...)`
Params: `modelNo` — The model serial number.
`chainID` — The identifier of the chain.
`seqNo` — The residue sequence number.
`insCode` — The residue insertion code.
Returns: A pointer to the requested `Residue` object.

#### `CoorManager::GetResidue(int modelNo, int chainNo, int seqNo, const InsCode insCode) -> PResidue`
Retrieves a residue object by model number, sequential chain number, sequence number, and insertion code.
Receiver: instance — call as `obj->GetResidue(...)`
Params: `modelNo` — The model serial number.
`chainNo` — The sequential number of the chain.
`seqNo` — The residue sequence number.
`insCode` — The residue insertion code.
Returns: A pointer to the requested `Residue` object.

#### `CoorManager::GetResidue(int modelNo, const ChainID chainID, int resNo) -> PResidue`
Retrieves a residue object by model number, chain ID, and residue number.
Receiver: instance — call as `obj->GetResidue(...)`
Params: `modelNo` — The model serial number.
`chainID` — The identifier of the chain.
`resNo` — The residue number (0-indexed).
Returns: A pointer to the requested `Residue` object.

#### `CoorManager::GetResidue(int modelNo, int chainNo, int resNo) -> PResidue`
Retrieves a residue object by model number, sequential chain number, and residue number.
Receiver: instance — call as `obj->GetResidue(...)`
Params: `modelNo` — The model serial number.
`chainNo` — The sequential number of the chain.
`resNo` — The residue number (0-indexed).
Returns: A pointer to the requested `Residue` object.

#### `CoorManager::GetResidue(cpstr CID) -> PResidue`
Retrieves a residue object matching a specified coordinate identifier string.
Receiver: instance — call as `obj->GetResidue(...)`
Params: `CID` — The coordinate identifier string.
Returns: A pointer to the requested `Residue` object.

#### `CoorManager::GetResidueNo(int modelNo, const ChainID chainID, int seqNo, const InsCode insCode) -> int`
Retrieves the sequence number of a residue given model, chain ID, sequence number, and insertion code.
Receiver: instance — call as `obj->GetResidueNo(...)`
Params: `modelNo` — The model serial number.
`chainID` — The identifier of the chain.
`seqNo` — The residue sequence number.
`insCode` — The residue insertion code.
Returns: The residue sequence number.

#### `CoorManager::GetResidueNo(int modelNo, int chainNo, int seqNo, const InsCode insCode) -> int`
Retrieves the sequence number of a residue given model, sequential chain number, sequence number, and insertion code.
Receiver: instance — call as `obj->GetResidueNo(...)`
Params: `modelNo` — The model serial number.
`chainNo` — The sequential number of the chain.
`seqNo` — The residue sequence number.
`insCode` — The residue insertion code.
Returns: The residue sequence number.

#### `CoorManager::GetResidueTable(PPResidue & resTable, int & NumberOfResidues) -> void`
Fills the provided structure table with pointers to all residues in a model and reports the total count.
Receiver: instance — call as `obj->GetResidueTable(...)`
Params: `resTable` — Pointer to the residue table to be filled.
`NumberOfResidues` — Reference to an integer receiving the count of residues.
Side effects: Populates `resTable` with pointers to all residues.

#### `CoorManager::GetResidueTable(int modelNo, const ChainID chainID, PPResidue & resTable, int & NumberOfResidues) -> void`
Fills the provided structure table with residues belonging to a specific chain in a model and reports the total count.
Receiver: instance — call as `obj->GetResidueTable(...)`
Params: `modelNo` — The model serial number.
`chainID` — The identifier of the chain.
`resTable` — Pointer to the residue table to be filled.
`NumberOfResidues` — Reference to an integer receiving the count of residues.
Side effects: Populates `resTable` with pointers to residues on the chain.

#### `CoorManager::GetResidueTable(int modelNo, int chainNo, PPResidue & resTable, int & NumberOfResidues) -> void`
Fills the provided structure table with residues belonging to a specific chain in a model and reports the total count.
Receiver: instance — call as `obj->GetResidueTable(...)`
Params: `modelNo` — The model serial number.
`chainNo` — The sequential number of the chain.
`resTable` — Pointer to the residue table to be filled.
`NumberOfResidues` — Reference to an integer receiving the count of residues.
Side effects: Populates `resTable` with pointers to residues on the chain.

#### `CoorManager::GetResidueTable(cpstr CID, PPResidue & resTable, int & NumberOfResidues) -> void`
Fills the provided structure table with all residues matching a specified coordinate identifier string and reports the total count.
Receiver: instance — call as `obj->GetResidueTable(...)`
Params: `CID` — The coordinate identifier string.
`resTable` — Pointer to the residue table to be filled.
`NumberOfResidues` — Reference to an integer receiving the count of residues.
Side effects: Populates `resTable` with pointers to residues matching the identifier.

#### `CoorManager::GetAtomTable(PPAtom & atomTable, int & NumberOfAtoms) -> void`
Fills the provided structure table with pointers to all atoms in a model and reports the total count.
Receiver: instance — call as `obj->GetAtomTable(...)`
Params: `atomTable` — Pointer to the atom table to be filled.
`NumberOfAtoms` — Reference to an integer receiving the count of atoms.
Side effects: Populates `atomTable` with pointers to all atoms.

#### `CoorManager::GetAtomTable(int modelNo, const ChainID chainID, int seqNo, const InsCode insCode, PPAtom & atomTable, int & NumberOfAtoms) -> void`
Fills the provided structure table with atoms corresponding to a specific residue in a model and reports the total count.
Receiver: instance — call as `obj->GetAtomTable(...)`
Params: `modelNo` — The model serial number.
`chainID` — The identifier of the chain.
`seqNo` — The residue sequence number.
`insCode` — The residue insertion code.
`atomTable` — Pointer to the atom table to be filled.
`NumberOfAtoms` — Reference to an integer receiving the count of atoms.
Side effects: Populates `atomTable` with pointers to atoms belonging to the residue.

#### `CoorManager::GetAtomTable(int modelNo, int chainNo, int seqNo, const InsCode insCode, PPAtom & atomTable, int & NumberOfAtoms) -> void`
Fills the provided structure table with atoms corresponding to a specific residue in a model and reports the total count.
Receiver: instance — call as `obj->GetAtomTable(...)`
Params: `modelNo` — The model serial number.
`chainNo` — The sequential number of the chain.
`seqNo` — The residue sequence number.
`insCode` — The residue insertion code.
`atomTable` — Pointer to the atom table to be filled.
`NumberOfAtoms` — Reference to an integer receiving the count of atoms.
Side effects: Populates `atomTable` with pointers to atoms belonging to the residue.

#### `CoorManager::GetAtomTable(int modelNo, const ChainID chainID, int resNo, PPAtom & atomTable, int & NumberOfAtoms) -> void`
Fills the provided structure table with atoms belonging to a specific residue in a model and reports the total count.
Receiver: instance — call as `obj->GetAtomTable(...)`
Params: `modelNo` — The model serial number.
`chainID` — The identifier of the chain.
`resNo` — The residue number (0-indexed).
`atomTable` — Pointer to the atom table to be filled.
`NumberOfAtoms` — Reference to an integer receiving the count of atoms.
Side effects: Populates `atomTable` with pointers to atoms belonging to the residue.

#### `CoorManager::GetAtomTable(int modelNo, int chainNo, int resNo, PPAtom & atomTable, int & NumberOfAtoms) -> void`
Fills the provided structure table with atoms belonging to a specific residue in a model and reports the total count.
Receiver: instance — call as `obj->GetAtomTable(...)`
Params: `modelNo` — The model serial number.
`chainNo` — The sequential number of the chain.
`resNo` — The residue number (0-indexed).
`atomTable` — Pointer to the atom table to be filled.
`NumberOfAtoms` — Reference to an integer receiving the count of atoms.
Side effects: Populates `atomTable` with pointers to atoms belonging to the residue.

#### `CoorManager::GetAtomTable(cpstr CID, PPAtom & atomTable, int & NumberOfAtoms) -> void`
Fills the provided structure table with all atoms matching a specified coordinate identifier string and reports the total count.
Receiver: instance — call as `obj->GetAtomTable(...)`
Params: `CID` — The coordinate identifier string.
`atomTable` — Pointer to the atom table to be filled.
`NumberOfAtoms` — Reference to an integer receiving the count of atoms.
Side effects: Populates `atomTable` with pointers to atoms matching the identifier.

#### `CoorManager::GetAtomTable1(PPAtom & atomTable, int & NumberOfAtoms) -> void`
Fills the provided atom table with all non-TER atoms and reports the total count, avoiding NULL pointers.
Receiver: instance — call as `obj->GetAtomTable1(...)`
Params: `atomTable` — Pointer to the atom table to be filled.
`NumberOfAtoms` — Reference to an integer receiving the count of atoms.
Side effects: Populates `atomTable` with non-TER atoms.

#### `CoorManager::GetAtomTable1(int modelNo, const ChainID chainID, int seqNo, const InsCode insCode, PPAtom & atomTable, int & NumberOfAtoms) -> void`
Fills the provided atom table with non-TER atoms corresponding to a specific residue in a model and reports the total count.
Receiver: instance — call as `obj->GetAtomTable1(...)`
Params: `modelNo` — The model serial number.
`chainID` — The identifier of the chain.
`seqNo` — The residue sequence number.
`insCode` — The residue insertion code.
`atomTable` — Pointer to the atom table to be filled.
`NumberOfAtoms` — Reference to an integer receiving the count of atoms.
Side effects: Populates `atomTable` with pointers to non-TER atoms belonging to the residue.

#### `CoorManager::GetAtomTable1(int modelNo, int chainNo, int seqNo, const InsCode insCode, PPAtom & atomTable, int & NumberOfAtoms) -> void`
Fills the provided atom table with non-TER atoms corresponding to a specific residue in a model and reports the total count.
Receiver: instance — call as `obj->GetAtomTable1(...)`
Params: `modelNo` — The model serial number.
`chainNo` — The sequential number of the chain.
`seqNo` — The residue sequence number.
`insCode` — The residue insertion code.
`atomTable` — Pointer to the atom table to be filled.
`NumberOfAtoms` — Reference to an integer receiving the count of atoms.
Side effects: Populates `atomTable` with pointers to non-TER atoms belonging to the residue.

#### `CoorManager::GetAtomTable1(int modelNo, const ChainID chainID, int resNo, PPAtom & atomTable, int & NumberOfAtoms) -> void`
Fills the provided atom table with non-TER atoms belonging to a specific residue in a model and reports the total count.
Receiver: instance — call as `obj->GetAtomTable1(...)`
Params: `modelNo` — The model serial number.
`chainID` — The identifier of the chain.
`resNo` — The residue number (0-indexed).
`atomTable` — Pointer to the atom table to be filled.
`NumberOfAtoms` — Reference to an integer receiving the count of atoms.
Side effects: Populates `atomTable` with pointers to non-TER atoms belonging to the residue.

#### `CoorManager::GetAtomTable1(int modelNo, int chainNo, int resNo, PPAtom & atomTable, int & NumberOfAtoms) -> void`
Fills the provided atom table with non-TER atoms belonging to a specific residue in a model and reports the total count.
Receiver: instance — call as `obj->GetAtomTable1(...)`
Params: `modelNo` — The model serial number.
`chainNo` — The sequential number of the chain.
`resNo` — The residue number (0-indexed).
`atomTable` — Pointer to the atom table to be filled.
`NumberOfAtoms` — Reference to an integer receiving the count of atoms.
Side effects: Populates `atomTable` with pointers to non-TER atoms belonging to the residue.

#### `CoorManager::GetAtomTable1(cpstr CID, PPAtom & atomTable, int & NumberOfAtoms) -> void`
Fills the provided atom table with all non-TER atoms matching a specified coordinate identifier string and reports the total count.
Receiver: instance — call as `obj->GetAtomTable1(...)`
Params: `CID` — The coordinate identifier string.
`atomTable` — Pointer to the atom table to be filled.
`NumberOfAtoms` — Reference to an integer receiving the count of atoms.
Side effects: Populates `atomTable` with pointers to non-TER atoms matching the identifier.

#### `CoorManager::DeleteAtom(int modelNo, const ChainID chID, int seqNo, const InsCode insCode, const AtomName aname, const Element elmnt, const AltLoc aloc) -> int`
Deletes a specific atom specified by its detailed location in the structure.
Receiver: instance — call as `obj->DeleteAtom(...)`
Params: `modelNo` — The model serial number.
`chID` — The identifier of the chain.
`seqNo` — The residue sequence number.
`insCode` — The residue insertion code.
`aname` — The name of the atom to delete.
`elmnt` — The chemical element code.
`aloc` — The alternate location indicator.
Returns: The number of atoms deleted (0 or 1).
Side effects: Removes the specified atom from the structure.

#### `CoorManager::DeleteAtom(int modelNo, const ChainID chID, int seqNo, const InsCode insCode, int atomNo) -> int`
Deletes a specific atom identified by its model, chain ID, residue sequence number, insertion code, and sequential atom number.
Receiver: instance — call as `obj->DeleteAtom(...)`
Params: `modelNo` — The model serial number.
`chID` — The identifier of the chain.
`seqNo` — The residue sequence number.
`insCode` — The residue insertion code.
`atomNo` — The sequential atom number.
Returns: The number of atoms deleted (0 or 1).
Side effects: Removes the specified atom from the structure.

#### `CoorManager::DeleteAtom(int modelNo, const ChainID chID, int resNo, const AtomName aname, const Element elmnt, const AltLoc aloc) -> int`
Deletes a specific atom specified by its model, chain ID, residue number, atom name, element, and alternate location.
Receiver: instance — call as `obj->DeleteAtom(...)`
Params: `modelNo` — The model serial number.
`chainID` — The identifier of the chain.
`resNo` — The residue number (0-indexed).
`aname` — The name of the atom to delete.
`elmnt` — The chemical element code.
`aloc` — The alternate location indicator.
Returns: The number of atoms deleted (0 or 1).
Side effects: Removes the specified atom from the structure.

#### `CoorManager::DeleteAtom(int modelNo, int chainNo, int seqNo, const InsCode insCode, const AtomName aname, const Element elmnt, const AltLoc aloc) -> int`
Deletes a specific atom specified by model, sequential chain number, residue sequence number, insertion code, atom name, element, and alternate location.
Receiver: instance — call as `obj->DeleteAtom(...)`
Params: `modelNo` — The model serial number.
`chainNo` — The sequential number of the chain.
`seqNo` — The residue sequence number.
`insCode` — The residue insertion code.
`aname` — The name of the atom to delete.
`elmnt` — The chemical element code.
`aloc` — The alternate location indicator.
Returns: The number of atoms deleted (0 or 1).
Side effects: Removes the specified atom from the structure.

#### `CoorManager::DeleteAtom(int modelNo, int chainNo, int resNo, int atomNo) -> int`
Deletes a specific atom identified by model, sequential chain number, residue number, and sequential atom number.
Receiver: instance — call as `obj->DeleteAtom(...)`
Params: `modelNo` — The model serial number.
`chainNo` — The sequential number of the chain.
`resNo` — The residue number (0-indexed).
`atomNo` — The sequential atom number.
Returns: The number of atoms deleted (0 or 1).
Side effects: Removes the specified atom from the structure.

#### `CoorManager::DeleteAllAtoms(int modelNo, const ChainID chID, int seqNo, const InsCode insCode) -> int`
Deletes all atoms belonging to a specific residue defined by model, chain ID, sequence number, and insertion code.
Receiver: instance — call as `obj->DeleteAllAtoms(...)`
Params: `modelNo` — The model serial number.
`chID` — The identifier of the chain.
`seqNo` — The residue sequence number.
`insCode` — The residue insertion code.
Returns: The total number of atoms deleted.
Side effects: Removes all specified atoms from the structure.

#### `CoorManager::DeleteAllAtoms(int modelNo, const ChainID chID, int resNo) -> int`
Deletes all atoms belonging to a specific residue in a model, chain ID, and residue number.
Receiver: instance — call as `obj->DeleteAllAtoms(...)`
Params: `modelNo` — The model serial number.
`chID` — The identifier of the chain.
`resNo` — The residue number (0-indexed).
Returns: The total number of atoms deleted.
Side effects: Removes all specified atoms from the structure.

#### `CoorManager::DeleteAllAtoms(int modelNo, const ChainID chID) -> int`
Deletes all atoms belonging to a specific chain in a model.
Receiver: instance — call as `obj->DeleteAllAtoms(...)`
Params: `modelNo` — The model serial number.
`chID` — The identifier of the chain.
Returns: The total number of atoms deleted.
Side effects: Removes all specified atoms from the structure.

#### `CoorManager::DeleteAllAtoms(int modelNo, int chNo, int seqNo, const InsCode insCode) -> int`
Deletes all atoms belonging to a specific residue defined by model, sequential chain number, sequence number, and insertion code.
Receiver: instance — call as `obj->DeleteAllAtoms(...)`
Params: `modelNo` — The model serial number.
`chNo` — The sequential number of the chain.
`seqNo` — The residue sequence number.
`insCode` — The residue insertion code.
Returns: The total number of atoms deleted.
Side effects: Removes all specified atoms from the structure.

#### `CoorManager::DeleteAllAtoms(int modelNo, int chNo, int resNo) -> int`
Deletes all atoms belonging to a specific residue defined by model, sequential chain number, and residue number.
Receiver: instance — call as `obj->DeleteAllAtoms(...)`
Params: `modelNo` — The model serial number.
`chNo` — The sequential number of the chain.
`resNo` — The residue number (0-indexed).
Returns: The total number of atoms deleted.
Side effects: Removes all specified atoms from the structure.

#### `CoorManager::DeleteAllAtoms(int modelNo, int chNo) -> int`
Deletes all atoms belonging to a specific chain in a model.
Receiver: instance — call as `obj->DeleteAllAtoms(...)`
Params: `modelNo` — The model serial number.
`chNo` — The sequential number of the chain.
Returns: The total number of atoms deleted.
Side effects: Removes all specified atoms from the structure.

#### `CoorManager::DeleteAllAtoms(int modelNo) -> int`
Deletes all atoms contained within a specific model.
Receiver: instance — call as `obj->DeleteAllAtoms(...)`
Params: `modelNo` — The model serial number.
Returns: The total number of atoms deleted.
Side effects: Removes all atoms from the specified model.

#### `CoorManager::DeleteAllAtoms() -> int`
Deletes all atoms from the entire structure managed by the coor manager.
Receiver: instance — call as `obj->DeleteAllAtoms()`
Returns: The total number of atoms deleted.
Side effects: Clears all atomic coordinates from the structure.

#### `CoorManager::DeleteAltLocs() -> int`
Removes all alternative location atoms, retaining only those with maximal occupancy or those with the least alternative location indicator if occupancy is equal or unspecified.
Receiver: instance — call as `obj->DeleteAltLocs()`
Returns: The number of alternative location atoms deleted.
Side effects: Modifies the structure by deleting atoms.

#### `CoorManager::AddAtom(int modelNo, const ChainID chID, int seqNo, const InsCode insCode, PAtom atom) -> int`
Adds a specified atom to the structure at the location defined by model, chain ID, sequence number, and insertion code.
Receiver: instance — call as `obj->AddAtom(...)`
Params: `modelNo` — The model serial number.
`chID` — The identifier of the chain.
`seqNo` — The residue sequence number.
`insCode` — The residue insertion code.
`atom` — Pointer to the atom to be added.
Returns: The success status of the addition.
Side effects: Adds an atom to the structure.

#### `CoorManager::AddAtom(int modelNo, const ChainID chID, int resNo, PAtom atom) -> int`
Adds a specified atom to the structure at the location defined by model, chain ID, residue number, and the atom pointer.
Receiver: instance — call as `obj->AddAtom(...)`
Params: `modelNo` — The model serial number.
`chainID` — The identifier of the chain.
`resNo` — The residue number (0-indexed).
`atom` — Pointer to the atom to be added.
Returns: The success status of the addition.
Side effects: Adds an atom to the structure.

#### `CoorManager::AddAtom(int modelNo, int chNo, int seqNo, const InsCode insCode, PAtom atom) -> int`
Adds a specified atom to the structure using model number, sequential chain number, sequence number, insertion code, and atom pointer.
Receiver: instance — call as `obj->AddAtom(...)`
Params: `modelNo` — The model serial number.
`chNo` — The sequential number of the chain.
`seqNo` — The residue sequence number.
`insCode` — The residue insertion code.
`atom` — Pointer to the atom to be added.
Returns: The success status of the addition.
Side effects: Adds an atom to the structure.

#### `CoorManager::AddAtom(int modelNo, int chNo, int resNo, PAtom atom) -> int`
Adds a specified atom to the structure using model number, sequential chain number, residue number, and atom pointer.
Receiver: instance — call as `obj->AddAtom(...)`
Params: `modelNo` — The model serial number.
`chNo` — The sequential number of the chain.
`resNo` — The residue number (0-indexed).
`atom` — Pointer to the atom to be added.
Returns: The success status of the addition.
Side effects: Adds an atom to the structure.

#### `CoorManager::Frac2Orth(realtype xfrac, realtype yfrac, realtype zfrac, realtype & xorth, realtype & yorth, realtype & zorth) -> bool`
Transforms fractional coordinates (xfrac, yfrac, zfrac) into orthogonal coordinates (xorth, yorth, zorth) if transformation matrices are set.
Receiver: instance — call as `obj->Frac2Orth(...)`
Params: `xfrac, yfrac, zfrac` — The coordinates in fractional space.
`xorth, yorth, zorth` — References to receive the transformed coordinates in orthogonal space.
Returns: True if the transformation was performed using calculated matrices; False otherwise.
Side effects: Updates the output coordinate references.

#### `CoorManager::Orth2Frac(realtype xorth, realtype yorth, realtype zorth, realtype & xfrac, realtype & yfrac, realtype & zfrac) -> bool`
Transforms orthogonal coordinates (xorth, yorth, zorth) into fractional coordinates (xfrac, yfrac, zfrac) if transformation matrices are set.
Receiver: instance — call as `obj->Orth2Frac(...)`
Params: `xorth, yorth, zorth` — The coordinates in orthogonal space.
`xfrac, yfrac, zfrac` — References to receive the transformed coordinates in fractional space.
Returns: True if the transformation was performed using calculated matrices; False otherwise.
Side effects: Updates the output coordinate references.

#### `CoorManager::Frac2Orth(mat44 & F, mat44 & T) -> bool`
Performs a full transformation of coordinates from fractional space using matrices F and T.
Receiver: instance — call as `obj->Frac2Orth(...)`
Params: `F` — Reference to the fractional coordinate transformation matrix.
`T` — Reference to the orthogonal coordinate transformation matrix.
Returns: True if the transformation was performed; False otherwise.
Side effects: Updates the transformation matrices F and T.

#### `CoorManager::Orth2Frac(mat44 & T, mat44 & F) -> bool`
Performs a full transformation of coordinates from orthogonal space using matrices T and F.
Receiver: instance — call as `obj->Orth2Frac(...)`
Params: `T` — Reference to the orthogonal coordinate transformation matrix.
`F` — Reference to the fractional coordinate transformation matrix.
Returns: True if the transformation was performed; False otherwise.
Side effects: Updates the transformation matrices T and F.

#### `CoorManager::SeekContacts(PPAtom AIndex, int ilen, int atomNum, realtype dist1, realtype dist2, int seqDist, RPContact contact, int & ncontacts, int maxlen, long group) -> void`
Finds contacts between a reference atom index set and all atoms in the structure, starting the check from a specific atom index.
Receiver: instance — call as `obj->SeekContacts(...)`
Params: `AIndex` — Array of atom indices to check against.
`ilen` — Length of the atom index array.
`atomNum` — Index of the first contact atom within `AIndex`.
`dist1` — Minimal contact distance.
`dist2` — Maximal contact distance.
`seqDist` — Sequence distance restriction for contact checking.
`contact` — Array to store indices of contacting atoms.
`ncontacts` — Reference to store the number of contacts found.
`maxlen` — Maximum size of the contact array.
`group` — Contact group ID counter.
Side effects: Populates the `contact` array with contact pairs.

#### `CoorManager::SeekContacts(PAtom A, PPAtom AIndex, int ilen, realtype dist1, realtype dist2, int seqDist, RPContact contact, int & ncontacts, int maxlen, long group) -> void`
Finds contacts between a reference atom and a set of candidate atoms in the structure.
Receiver: instance — call as `obj->SeekContacts(...)`
Params: `A` — The first atom to check for contacts.
`AIndex` — Array of atoms to check against.
`ilen` — Length of the candidate atom array.
`dist1` — Minimal contact distance.
`dist2` — Maximal contact distance.
`seqDist` — Sequence distance restriction for contact checking.
`contact` — Array to store indices of contacting atoms.
`ncontacts` — Reference to store the number of contacts found.
`maxlen` — Maximum size of the contact array.
`group` — Contact group ID counter.
Side effects: Populates the `contact` array with contact pairs.

#### `CoorManager::SeekContacts(PPAtom AIndex1, int ilen1, PPAtom AIndex2, int ilen2, realtype dist1, realtype dist2, int seqDist, RPContact contact, int & ncontacts, int maxlen, mat44 * TMatrix, long group,

---

## `Cryst`
Represents the crystallographic structure information, including cell parameters and symmetry operations.

#### `Cryst::GetNumberOfSymOps()`
Returns the total number of symmetry operations stored for the crystal.
Receiver: instance — call as `obj->GetNumberOfSymOps()`
Returns: The count of symmetry operations.

#### `Cryst::GetTMatrix(mat44& TMatrix, int Nop, int cellshift_a, int cellshift_b, int cellshift_c, PSymOps symOpers)`
Calculates and returns the coordinate transformation matrix converting orthogonal coordinates based on the symmetry operation Nop, shifted by given cell parameters.
Receiver: instance — call as `obj->GetTMatrix(...)`
Params: `TMatrix` — The matrix to receive the result. `Nop` — The symmetry operation number. `cellshift_a`, `cellshift_b`, `cellshift_c` — The cell shifts applied to the result. `symOpers` — Optional list of symmetry operations.
Returns: 0 on success, 1 if Nop is undefined, 2 if transformations matrices are missing, or 3 if cell parameters are not set.
Side effects: Sets the value of `TMatrix`.

#### `Cryst::GetSymOp(int Nop)`
Retrieves the symmetry operation associated with the given operation number.
Receiver: instance — call as `obj->GetSymOp(...)`
Params: `Nop` — The symmetry operation number.
Returns: The PDB string representation of the symmetry operation.

#### `Cryst::GetNumberofNCSMatrices()`
Returns the count of Non-Crystallographic Symmetry matrices stored.
Receiver: instance — call as `obj->GetNumberOfNCSMatrices()`
Returns: The number of NCS matrices.

#### `Cryst::GetNCSMatrix(int NCSMatrixNo, mat33& ncs_m, vect3& ncs_v)`
Retrieves the Non-Crystallographic Symmetry matrix and translation vector for a given index.
Receiver: instance — call as `obj->GetNCSMatrix(...)`
Params: `NCSMatrixNo` — The index of the matrix to retrieve. `ncs_m` — The matrix to receive the result. `ncs_v` — The vector to receive the result.
Returns: True if the matrix and vector were successfully retrieved, false otherwise.
Side effects: Sets the values of `ncs_m` and `ncs_v`.

#### `Cryst::GetNCSMatrix(int NCSMatrixNo, mat44& ncs_m, int& iGiven)`
Retrieves the Non-Crystallographic Symmetry matrix and an indicator of if it was read from the file for a given index.
Receiver: instance — call as `obj->GetNCSMatrix(...)`
Params: `NCSMatrixNo` — The index of the matrix to retrieve. `ncs_m` — The matrix to receive the result. `iGiven` — Indicator of if the matrix was given in the file.
Returns: True if the matrix was successfully retrieved, false otherwise.
Side effects: Sets the value of `ncs_m` and `iGiven`.

#### `Cryst::GetUCTMatrix(mat44& TMatrix, int Nop, realtype x, realtype y, realtype z, int cellshift_a, int cellshift_b, int cellshift_c, PSymOps symOpers)`
Calculates and returns the coordinate transformation matrix for orthogonal coordinates based on the symmetry operation Nop, shifted by given cell parameters, positioning the point closest to the primary unit cell center.
Receiver: instance — call as `obj->GetUCTMatrix(...)`
Params: `TMatrix` — The matrix to receive the result. `Nop` — The symmetry operation number. `x`, `y`, `z` — Coordinates used for centering. `cellshift_a`, `cellshift_b`, `cellshift_c` — The cell shifts applied to the result. `symOpers` — Optional list of symmetry operations.
Returns: 0 on success, 1 if Nop is undefined, 2 if transformations matrices are missing, or 3 if cell parameters are not set.
Side effects: Sets the value of `TMatrix`.

#### `Cryst::GetFractMatrix(mat44& TMatrix, int Nop, int cellshift_a, int cellshift_b, int cellshift_c, PSymOps symOpers)`
Calculates and returns the coordinate transformation matrix for fractional coordinates based on the symmetry operation Nop, shifted by given cell parameters.
Receiver: instance — call as `obj->GetFractMatrix(...)`
Params: `TMatrix` — The matrix to receive the result. `Nop` — The symmetry operation number. `cellshift_a`, `cellshift_b`, `cellshift_c` — The cell shifts applied to the result. `symOpers` — Optional list of symmetry operations.
Returns: 0 on success, 1 if Nop is undefined, 2 if transformations matrices are missing, or 3 if cell parameters are not set.
Side effects: Sets the value of `TMatrix`.

#### `Cryst::GetSymOpMatrix(mat44& TMatrix, int Nop)`
Returns the transformation matrix for the Nop-th symmetry operator in the space group.
Receiver: instance — call as `obj->GetSymOpMatrix(...)`
Params: `TMatrix` — The matrix to receive the result. `Nop` — The symmetry operation number.
Returns: 0 on success, 1 if Nop is undefined, 2 if transformations matrices are missing, or 3 if cell parameters are not set.
Side effects: Sets the value of `TMatrix`.

---

## `DBReference`
A class representing a reference to a sequence segment within a structural database, handling I/O and structural information mapping.

#### `DBReference::DBReference()`
Constructs a default initialized `DBReference` object.
Receiver: instance — call as `obj->DBReference()`
Returns: None
Side effects: None

#### `DBReference::DBReference(PChain Chain_Owner)`
Constructs a `DBReference` associated with a given chain owner.
Receiver: instance — call as `obj->DBReference(Chain_Owner)`
Params: `Chain_Owner` — Pointer to the chain that owns this reference.
Returns: None
Side effects: None

#### `DBReference::DBReference(PChain Chain_Owner, cpstr S)`
Constructs a `DBReference` associated with a given chain owner and sequence string.
Receiver: instance — call as `obj->DBReference(Chain_Owner, S)`
Params: `Chain_Owner` — Pointer to the chain that owns this reference.
`S` — The sequence string to be referenced.
Returns: None
Side effects: None

#### `DBReference::DBReference(io::RPStream Object)`
Constructs a `DBReference` by reading from an input stream object.
Receiver: instance — call as `obj->DBReference(Object)`
Params: `Object` — An input stream object containing the reference data.
Returns: None
Side effects: None

#### `DBReference::~DBReference()`
Destructor for the `DBReference` object.
Receiver: instance — call as `obj->~DBReference()`
Returns: None
Side effects: None

#### `DBReference::PDBASCIIDump(pstr S, int N)`
Writes a textual representation of the PDB-formatted reference data to the specified string.
Receiver: instance — call as `obj->PDBASCIIDump(S, N)`
Params: `S` — The string buffer to write to.
`N` — The maximum number of characters to write.
Returns: None
Side effects: Writes data to string `S`.

#### `DBReference::MakeCIF(mmcif::PData CIF, int N)`
Writes the structural information of the reference into a Cif structure.
Receiver: instance — call as `obj->MakeCIF(CIF, N)`
Params: `CIF` — Pointer to the Cif data structure to populate.
`N` — The maximum number of atoms/records to process.
Returns: None
Side effects: Modifies the `CIF` data structure.

#### `DBReference::ConvertPDBASCII(cpstr S)`
Attempts to parse and convert PDB ASCII formatted data contained in `S` into the internal representation.
Receiver: instance — call as `obj->ConvertPDBASCII(S)`
Params: `S` — The PDB ASCII string to parse.
Returns: An error code indicating success or failure.
Side effects: Parses and populates the object's internal state.

#### `DBReference::GetCIF(mmcif::PData CIF, int & n)`
Extracts or generates Cif formatted data for the reference into the provided Cif structure.
Receiver: instance — call as `obj->GetCIF(CIF, n)`
Params: `CIF` — Pointer to the Cif data structure to receive data.
`n` — Reference parameter to receive the count of records written.
Returns: An error code indicating success or failure.
Side effects: Populates the `CIF` data structure.

#### `DBReference::GetClassID()`
Returns the specific class identifier for a `DBReference` object.
Receiver: instance — call as `obj->GetClassID()`
Returns: The constant class identifier `ClassID_DBReference`.
Side effects: None

#### `DBReference::Copy(PContainerClass DBRef)`
Copies the structural and reference information from another `DBReference` object.
Receiver: instance — call as `obj->Copy(DBRef)`
Params: `DBRef` — Pointer to the source `DBReference` object.
Returns: None
Side effects: Updates the object's internal state based on `DBRef`.

#### `DBReference::write(io::RFile f)`
Writes the structural data represented by the `DBReference` to a specified file stream.
Receiver: instance — call as `obj->write(f)`
Params: `f` — The file stream to write to.
Returns: None
Side effects: Performs file I/O to stream `f`.

#### `DBReference::read(io::RFile f)`
Reads structural data into the `DBReference` object from a specified file stream.
Receiver: instance — call as `obj->read(f)`
Params: `f` — The file stream to read from.
Returns: None
Side effects: Performs file I/O from stream `f` and updates object state.

---

## `Data`
Manages reading, storing, and retrieving structural data blocks, including structures and loops, from mmCIF files.

#### `Data::SetFlag(CIF_FLAG F)` -> void
Sets optional flags for reading mmCIF files, controlling warning output and suggesting categories/tags.
Receiver: instance — call as `obj->SetFlag(...)`
Params: `F` — flag or logical "or" of several flags to be set.
Returns: void
Side effects: Mutates internal reading flags.
Errors: None.

#### `Data::RemoveFlag(CIF_FLAG F)` -> void
Removes optional flags previously set for reading mmCIF files.
Receiver: instance — call as `obj->RemoveFlag(...)`
Params: `F` — flag or logical "or" of several flags to be removed.
Returns: void
Side effects: Mutates internal reading flags.
Errors: None.

#### `Data::ReadMMCIFData(cpstr FName, io::GZ_MODE gzipMode=io::GZM_CHECK)` -> int
Reads the first data block from the specified file, handling potential compression.
Receiver: instance — call as `obj->ReadMMCIFData(...)`
Params: `FName` — character null-terminated string with file name; `gzipMode` — flag for reading compressed files.
Returns: Returns status code (CIFRC_Ok, <0 error, >0 warning).
Side effects: Reads content from `FName` into the internal data structure.
Errors: Returns negative value if file reading fails.

#### `Data::ReadMMCIFData(io::RFile f, pstr S, int & lcount)` -> int
Reads sequential mmCIF data blocks from the current position of the file stream.
Receiver: instance — call as `obj->ReadMMCIFData(...)`
Params: `f` — reference to an opened file stream; `S` — buffer string for reading; `lcount` — global line counter passed by reference.
Returns: Returns status code (CIFRC_Ok, <0 error, >0 warning).
Side effects: Reads content from the file stream `f` into the internal data structure.
Errors: Returns negative value if reading fails.

#### `Data::PutDataName(cpstr dname)` -> void
Stores a name for the data block being written.
Receiver: instance — call as `obj->PutDataName(...)`
Params: `dname` — the name to associate with the data record.
Returns: void
Side effects: Sets the name for the next data block to be written.
Errors: None.

#### `Data::GetNumberOfCategories()` -> int
Returns the total count of categories (structures and loops) currently stored in the data block.
Receiver: instance — call as `obj->GetNumberOfCategories()`
Params: None.
Returns: The number of categories present.
Side effects: None.
Errors: None.

#### `Data::GetCategory(int categoryNo)` -> PCategory
Retrieves a pointer to the category (structure or loop) by its zero-based index number.
Receiver: instance — call as `obj->GetCategory(...)`
Params: `categoryNo` — the zero-based index of the category.
Returns: A pointer to the category, or `NULL` if the index is out of range.
Side effects: None.
Errors: Returns `NULL` if `categoryNo` is invalid.

#### `Data::GetStructure(cpstr CName)` -> PStruct
Retrieves a pointer to the mmCIF structure identified by its given name.
Receiver: instance — call as `obj->GetStructure(...)`
Params: `CName` — character string with the structure's name.
Returns: A pointer to the structure, or `NULL` if not found.
Side effects: None.
Errors: Returns `NULL` if the structure name is not found.

#### `Data::GetLoop(cpstr CName)` -> PLoop
Retrieves a pointer to the mmCIF loop identified by its given name.
Receiver: instance — call as `obj->GetLoop(...)`
Params: `CName` — character string with the loop's name.
Returns: A pointer to the loop, or `NULL` if not found.
Side effects: None.
Errors: Returns `NULL` if the loop name is not found.

#### `Data::FindLoop(cpstr *tagList)` -> PLoop
Searches for and returns a pointer to the first loop containing all specified tags from the provided list.
Receiver: instance — call as `obj->FindLoop(...)`
Params: `tagList` — list of tags to look for, terminated by an empty string `""`.
Returns: A pointer to the found loop, or `NULL` otherwise.
Side effects: None.
Errors: Returns `NULL` if no loop contains all the specified tags.

#### `Data::GetDataName(pstr & dname, bool Remove=false)` -> void
Retrieves the name of the data block into a provided string reference, optionally removing it from the block.
Receiver: instance — call as `obj->GetDataName(...)`
Params: `dname` — pointer reference to a string accepting the name; `Remove` — if true, the name is removed after retrieval.
Returns: void
Side effects: Populates the string referenced by `dname` and potentially modifies the data block.
Errors: None.

#### `Data::GetDataName()` -> pstr
Returns the name of the data block as a dynamically allocated string.
Receiver: instance — call as `obj->GetDataName()`
Params: None.
Returns: The name of the data block.
Side effects: None.
Errors: None.

#### `Data::GetString(pstr & Dest, cpstr CName, cpstr TName, bool Remove=false)` -> int
Retrieves a string field from the data block into a destination string, optionally removing the field.
Receiver: instance — call as `obj->GetString(...)`
Params: `Dest` — reference to the destination string; `CName` — category name; `TName` — tag name; `Remove` — if true, removes the field after reading.
Returns: An integer status code indicating success or failure reason.
Side effects: Populates `Dest` and potentially removes the field from the data block.
Errors: Returns codes like `CIFRC_WrongFormat` or `CIFRC_NoCategory` on failure.

#### `Data::GetString(cpstr CName, cpstr TName, int & RC)` -> pstr
Retrieves a string field from the data block, returning a dynamically allocated string.
Receiver: instance — call as `obj->GetString(...)`
Params: `CName` — category name; `TName` — tag name; `RC` — reference integer to receive return code.
Returns: A dynamically allocated string containing the value.
Side effects: None.
Errors: Returns `NULL` or handles conversion errors via `RC`.

#### `Data::GetReal(realtype & R, cpstr CName, cpstr TName, bool Remove=false)` -> int
Retrieves a real-valued field from the data block into a destination variable, optionally removing the field.
Receiver: instance — call as `obj->GetReal(...)`
Params: `R` — reference to the destination real value; `CName` — category name; `TName` — tag name; `Remove` — if true, removes the field after reading.
Returns: An integer status code indicating success or failure reason.
Side effects: Populates the reference `R` and potentially removes the field from the data block.
Errors: Returns codes like `CIFRC_WrongFormat` or `CIFRC_NoCategory` on failure.

#### `Data::GetInteger(int & I, cpstr CName, cpstr TName, bool Remove=false)` -> int
Retrieves an integer field from the data block into a destination variable, optionally removing the field.
Receiver: instance — call as `obj->GetInteger(...)`
Params: `I` — reference to the destination integer; `CName` — category name; `TName` — tag name; `Remove` — if true, removes the field after reading.
Returns: An integer status code indicating success or failure reason.
Side effects: Populates the reference `I` and potentially removes the field from the data block.
Errors: Returns codes like `CIFRC_WrongFormat` or `CIFRC_NoCategory` on failure.

#### `Data::GetLoopLength(cpstr CName)` -> int
Returns the total number of data entries (length) for a specified loop category.
Receiver: instance — call as `obj->GetLoopLength(...)`
Params: `CName` — character string with the loop's name.
Returns: The length of the loop data (0 if empty).
Side effects: None.
Errors: Returns `CIFRC_NotALoop` or `CIFRC_NoCategory` if the category is not a loop or does not exist.

#### `Data::GetLoopString(pstr & Dest, cpstr CName, cpstr TName, int nrow, bool Remove=false)` -> int
Retrieves the string data for a specific row number within a loop category.
Receiver: instance — call as `obj->GetLoopString(...)`
Params: `Dest` — reference to the destination string; `CName` — category name; `TName` — tag name; `nrow` — the zero-based index of the row; `Remove` — if true, removes the field after reading.
Returns: An integer status code indicating success or failure reason.
Side effects: Populates `Dest` and potentially removes the field from the data block.
Errors: Returns `CIFRC_WrongIndex` if `nrow` is out of range.

#### `Data::GetLoopString(cpstr CName, cpstr TName, int nrow, int & RC)` -> pstr
Retrieves the string data for a specific row number within a loop category, returning a dynamically allocated string.
Receiver: instance — call as `obj->GetLoopString(...)`
Params: `CName` — category name; `TName` — tag name; `nrow` — the zero-based index of the row; `RC` — reference integer to receive return code.
Returns: A dynamically allocated string containing the value.
Side effects: None.
Errors: Returns `NULL` or handles conversion errors via `RC`.

#### `Data::GetLoopReal(realtype & R, cpstr CName, cpstr TName, int nrow, bool Remove=false)` -> int
Retrieves a real-valued data entry for a specific row number within a loop category.
Receiver: instance — call as `obj->GetLoopReal(...)`
Params: `R` — reference to the destination real value; `CName` — category name; `TName` — tag name; `nrow` — the zero-based index of the row; `Remove` — if true, removes the field after reading.
Returns: An integer status code indicating success or failure reason.
Side effects: Populates the reference `R` and potentially removes the field from the data block.
Errors: Returns codes like `CIFRC_WrongIndex` or `CIFRC_NoCategory` on failure.

#### `Data::GetLoopInteger(int & I, cpstr CName, cpstr TName, int nrow, bool Remove=false)` -> int
Retrieves an integer data entry for a specific row number within a loop category.
Receiver: instance — call as `obj->GetLoopInteger(...)`
Params: `I` — reference to the destination integer; `CName` — category name; `TName` — tag name; `nrow` — the zero-based index of the row; `Remove` — if true, removes the field after reading.
Returns: An integer status code indicating success or failure reason.
Side effects: Populates the reference `I` and potentially removes the field from the data block.
Errors: Returns codes like `CIFRC_WrongIndex` or `CIFRC_NoCategory` on failure.

#### `Data::GetLoopSVector(psvector & S, cpstr CName, cpstr TName, int i1=0, int i2=MaxInt4, bool Remove=false)` -> int
Reads a range of string data entries from a loop category into a provided vector, optionally removing the data.
Receiver: instance — call as `obj->GetLoopSVector(...)`
Params: `S` — reference to the destination vector; `CName` — category name; `TName` — tag name; `i1` — starting index (inclusive); `i2` — ending index (inclusive); `Remove` — if true, removes the fields after reading.
Returns: An integer status code matching other vector reading functions.
Side effects: Populates the vector `S` and potentially removes the fields from the data block.
Errors: Returns status codes indicating data access failures.

#### `Data::GetLoopRVector(rvector & R, cpstr CName, cpstr TName, int i1=0, int i2=MaxInt4, bool Remove=false)` -> int
Reads a range of real-valued data entries from a loop category into a provided vector, optionally removing the data.
Receiver: instance — call as `obj->GetLoopRVector(...)`
Params: `R` — reference to the destination vector; `CName` — category name; `TName` — tag name; `i1` — starting index (inclusive); `i2` — ending index (inclusive); `Remove` — if true, removes the fields after reading.
Returns: An integer status code matching other vector reading functions.
Side effects: Populates the vector `R` and potentially removes the fields from the data block.
Errors: Returns status codes indicating data access failures.

#### `Data::GetLoopIVector(ivector & I, cpstr CName, cpstr TName, int i1=0, int i2=MaxInt4, bool Remove=false)` -> int
Reads a range of integer data entries from a loop category into a provided vector, optionally removing the data.
Receiver: instance — call as `obj->GetLoopIVector(...)`
Params: `I` — reference to the destination vector; `CName` — category name; `TName` — tag name; `i1` — starting index (inclusive); `i2` — ending index (inclusive); `Remove` — if true, removes the fields after reading.
Returns: An integer status code matching other vector reading functions.
Side effects: Populates the vector `I` and potentially removes the fields from the data block.
Errors: Returns status codes indicating data access failures.

#### `Data::PutDataName(cpstr dname)` -> int
Stores a data block name to be written to the file.
Receiver: instance — call as `obj->PutDataName(...)`
Params: `dname` — the name for the data block.
Returns: Status code (CIFRC_Ok or CIFRC_NotAStruct).
Side effects: Sets the name metadata for the upcoming data block write.
Errors: None.

#### `Data::GetCategoryNo(cpstr cname)` -> int
Searches for the index of a category by name within the stored data block.
Receiver: instance — call as `obj->GetCategoryNo(...)`
Params: `cname` — name of the category to search for.
Returns: The zero-based position of the category, or a negative number indicating the correct insertion point if not found.
Side effects: None.
Errors: None.

#### `Data::GetLoop(io::RFile f, pstr S, pstr & L, pstr & p, int & lcount, int & llen)` -> void
Reads a data block section intended to be a loop from the file stream.
Receiver: instance — call as `obj->GetLoop(...)`
Params: `f` — file stream reference; `S` — buffer string; `L` — output list of values; `p` — output pointer; `lcount` — line counter reference; `llen` — line length reference.
Returns: void
Side effects: Advances file stream position and populates the loop data structure.
Errors: None.

#### `Data::GetField(io::RFile f, pstr S, pstr & L, pstr & p, int & lcount, int & llen)` -> int
Reads a generic data field from the file stream, used for general data extraction.
Receiver: instance — call as `obj->GetField(...)`
Params: `f` — file stream reference; `S` — buffer string; `L` — output list of values; `p` — output pointer; `lcount` — line counter reference; `llen` — line length reference.
Returns: Status code indicating read success or failure.
Side effects: Advances file stream position and populates the field data structure.
Errors: Returns negative status code upon failure.

## `Data`
Manages reading, storing, and retrieving structural data blocks, including structures and loops, from mmCIF files.

---

## `Edge`
Represents a connection between two vertices in a graph structure.

#### `Edge::GetType() -> int`
Gets the type identifier associated with the edge.
Receiver: instance — call as `obj->GetType()`
Returns: The edge's type integer.
Side effects: None.

#### `Edge::GetVertex1() -> int`
Gets the index of the first vertex connected by this edge.
Receiver: instance — call as `obj->GetVertex1()`
Returns: The first vertex index.
Side effects: None.

#### `Edge::GetVertex2() -> int`
Gets the index of the second vertex connected by this edge.
Receiver: instance — call as `obj->GetVertex2()`
Returns: The second vertex index.
Side effects: None.

---

## `File`
Manages file I/O operations for macromolecular structure data, specifically supporting MMCIF and CIF formats.

#### `File::AddCIFData(cpstr DName) -> int`
Adds a new CIF dataset to the file structure.
Receiver: instance — call as `obj->AddCIFData(...)`
Params: `DName` — The name identifier for the dataset.
Returns: The index number assigned to the newly added dataset.
Side effects: The file structure dataset count is incremented.

#### `File::GetCIFData(int dataNo) -> PData`
Retrieves a pointer to a CIF dataset using its integer index.
Receiver: instance — call as `obj->GetCIFData(...)`
Params: `dataNo` — The zero-based index of the dataset to retrieve.
Returns: A pointer to the requested CIF dataset data structure.
Errors: Returns null if the index is out of bounds.

#### `File::GetCIFData(cpstr DName) -> PData`
Retrieves a pointer to a CIF dataset using its specified name.
Receiver: instance — call as `obj->GetCIFData(...)`
Params: `DName` — The name identifier of the dataset to retrieve.
Returns: A pointer to the requested CIF dataset data structure.
Errors: Returns null if no dataset with the given name exists.

#### `File::ReadMMCIFFile(cpstr FName, io::GZ_MODE gzipMode) -> int`
Reads structural data from a specified MMCIF file, handling optional gzip compression.
Receiver: instance — call as `obj->ReadMMCIFFile(...)`
Params: `FName` — The path to the MMCIF file.
`gzipMode` — The compression mode to use (e.g., compressed or uncompressed).
Returns: An integer status code indicating success or failure.
Side effects: The internal data structure is populated with the contents of the file.

#### `File::WriteMMCIFFile(cpstr FName, io::GZ_MODE gzipMode) -> int`
Writes the current structural data contained in the file object to a specified MMCIF file, handling optional gzip compression.
Receiver: instance — call as `obj->WriteMMCIFFile(...)`
Params: `FName` — The path where the MMCIF file will be written.
`gzipMode` — The compression mode to use (e.g., compressed or uncompressed).
Returns: An integer status code indicating success or failure.
Side effects: Creates or overwrites the file at `FName`.

#### `File::GetNofData() -> int`
Returns the total number of individual data datasets currently stored in the file object.
Receiver: instance — call as `obj->GetNofData()`
Returns: The count of datasets.

#### `File::GetCIFDataNo(cpstr DName) -> int`
Retrieves the zero-based integer index of a dataset given its name.
Receiver: instance — call as `obj->GetCIFDataNo(...)`
Params: `DName` — The name identifier of the dataset.
Returns: The index number, or -1 if the dataset name is not found.

#### `File::deleteCIFData(cpstr DName) -> int`
Deletes a CIF dataset identified by its name from the file object.
Receiver: instance — call as `obj->DeleteCIFData(...)`
Params: `DName` — The name identifier of the dataset to delete.
Returns: An integer status code indicating success or failure.

#### `File::deleteCIFData(int dataNo) -> int`
Deletes a CIF dataset identified by its zero-based index from the file object.
Receiver: instance — call as `obj->DeleteCIFData(...)`
Params: `dataNo` — The zero-based index of the dataset to delete.
Returns: An integer status code indicating success or failure.

---

## `Graph`
Represents a graph structure used to model connectivity relationships between atoms and residues in structural biology.

#### `Graph::SetName(cpstr gname) -> void`
Sets the user-defined name for this graph object.
Receiver: instance — call as `obj->SetName(...)`
Params: `gname` — The name to assign to the graph.
Side effects: Changes the internal name of the graph.

#### `Graph::AddVertex(PVertex V) -> void`
Adds a vertex to the graph, taking ownership of the provided vertex pointer.
Receiver: instance — call as `obj->AddVertex(...)`
Params: `V` — Pointer to the vertex to be added.
Side effects: The graph now manages the memory and structure of V.
Errors: The calling application must allocate V before calling this method.

#### `Graph::AddEdge(PEdge G) -> void`
Adds an edge to the graph, taking ownership of the provided edge pointer.
Receiver: instance — call as `obj->AddEdge(...)`
Params: `G` — Pointer to the edge to be added.
Side effects: The graph now manages the memory and structure of G.
Errors: The calling application must allocate G before calling this method.

#### `Graph::MakeGraph(PResidue R, cpstr altLoc=NULL) -> int`
Constructs a graph representing the connectivity of the residue R based on chemical bonds and distance cutoffs.
Receiver: instance — call as `obj->MakeGraph(...)`
Params: `R` — The residue whose connectivity defines the graph. `altLoc` — The specific conformation alternate location to use.
Returns: An integer code indicating success or the reason for graph construction.
Side effects: Populates the graph's vertices and edges based on R.
Errors: May return codes like `MKGRAPH_NoAtoms` if the residue lacks atoms.

#### `Graph::MakeGraph(PPAtom atom, int nAtoms) -> int`
Constructs a graph directly from a set of specified atoms.
Receiver: instance — call as `obj->MakeGraph(...)`
Params: `atom` — Pointer to the first atom. `nAtoms` — The count of atoms provided.
Returns: An integer indicating success or failure in constructing the graph.
Side effects: Populates the graph's vertices and edges based on the provided atoms.

#### `Graph::MakeSymmetryRelief(bool noCO2) -> void`
Calculates and sets up the graph structure to reflect inherent symmetry relationships within the molecule.
Receiver: instance — call as `obj->MakeSymmetryRelief(...)`
Params: `noCO2` — Flag indicating whether CO2 groups should be ignored during symmetry analysis.
Side effects: Modifies the graph structure to account for symmetry.

#### `Graph::Build(bool bondOrder) -> int`
Finalizes the graph construction by establishing chemical bonds between atoms based on geometric criteria.
Receiver: instance — call as `obj->Build(...)`
Params: `bondOrder` — Flag specifying whether bond order constraints should be enforced.
Returns: 0 if the graph building process was successful.
Side effects: Establishes chemical bonds represented by edges.

#### `Graph::MakeVertexIDs() -> void`
Assigns sequential, unique integer identifiers to all vertices currently in the graph.
Receiver: instance — call as `obj->MakeVertexIDs()`
Side effects: Updates vertex metadata with new IDs.

#### `Graph::GetVertexID(int vertexNo) -> int`
Retrieves the unique internal identifier assigned to a vertex given its sequential index.
Receiver: instance — call as `obj->GetVertexID(...)`
Params: `vertexNo` — The sequential number (1-based) of the vertex.
Returns: The internal vertex ID corresponding to the given sequential number.
Errors: Returns an error value if the vertex number is invalid.

#### `Graph::GetVertexNo(cpstr vname) -> int`
Finds the sequential number corresponding to a vertex based on its user-defined name.
Receiver: instance — call as `obj->GetVertexNo(...)`
Params: `vname` — The name of the vertex to look up.
Returns: The sequential number of the vertex.
Errors: Returns an error value if the vertex name is not found.

#### `Graph::GetNBondedVertices(int vertexNo) -> int`
Counts how many vertices are connected by a bond to the specified vertex.
Receiver: instance — call as `obj->GetNBondedVertices(...)`
Params: `vertexNo` — The sequential number of the vertex to check.
Returns: The count of bonded neighbors.
Errors: Returns 0 if the vertex number is invalid.

#### `Graph::GetBondedVertexID(int vertexNo, int bond_vx_type, int bondNo) -> int`
Retrieves the internal ID of a specific neighbor vertex connected to the input vertex via a defined bond.
Receiver: instance — call as `obj->GetBondedVertexID(...)`
Params: `vertexNo` — The sequential number of the starting vertex. `bond_vx_type` — The type of the bond connection. `bondNo` — The specific bond index.
Returns: The internal ID of the neighboring vertex.
Errors: Returns an error value if the specified bond does not exist or vertex numbers are invalid.

#### `Graph::GetVertex(int vertexNo) -> PVertex`
Retrieves a pointer to the vertex located at the specified sequential index.
Receiver: instance — call as `obj->GetVertex(...)`
Params: `vertexNo` — The sequential number (1-based) of the vertex to retrieve.
Returns: A pointer to the requested vertex, or nullptr if it does not exist.
Errors: Returns nullptr if `vertexNo` is out of bounds (1 to nVertices).

#### `Graph::GetVertices(PPVertex & V, int & nV) -> void`
Populates the output pointer `V` with all vertices belonging to the graph and sets `nV` to the total count.
Receiver: instance — call as `obj->GetVertices(...)`
Params: `V` — Reference to a pointer that will receive the array of vertices. `nV` — Reference to an integer that will receive the count.
Side effects: Populates the provided output arguments.

#### `Graph::GetEdge(int edgeNo) -> PEdge`
Retrieves a pointer to the edge located at the specified sequential index.
Receiver: instance — call as `obj->GetEdge(...)`
Params: `edgeNo` — The sequential number (1-based) of the edge.
Returns: A pointer to the requested edge, or nullptr if it does not exist.
Errors: Returns nullptr if `edgeNo` is out of bounds (1 to nEdges).

#### `Graph::GetEdges(PPEdge & E, int & nE) -> void`
Populates the output pointer `E` with all edges belonging to the graph and sets `nE` to the total count.
Receiver: instance — call as `obj->GetEdges(...)`
Params: `E` — Reference to a pointer that will receive the array of edges. `nE` — Reference to an integer that will receive the count.
Side effects: Populates the provided output arguments.

#### `Graph::Print() -> void`
Prints a detailed representation of the graph structure to standard output.
Receiver: instance — call as `obj->Print()`
Side effects: Writes descriptive information about the graph to the output stream.

#### `Graph::Print1() -> void`
Prints a simplified representation of the graph structure to standard output.
Receiver: instance — call as `obj->Print1()`
Side effects: Writes descriptive information about the graph to the output stream.

#### `Graph::Copy(PGraph G) -> void`
Performs a deep copy of the structure and content from another graph object `G` into this graph.
Receiver: instance — call as `obj->Copy(...)`
Params: `G` — Pointer to the graph object to copy from.
Side effects: Modifies the internal state of the current graph to mirror G.

#### `Graph::SetName(cpstr gname) -> void`
Sets the user-defined name for this graph object.
Receiver: instance — call as `obj->SetName(...)`
Params: `gname` — The name to assign to the graph.
Side effects: Changes the internal name of the graph.


---

## `GraphMatch`
Manages the process and results of finding maximal common subgraphs between two input structures.

#### `GraphMatch::SetTimeLimit(int maxTimeToRun=0) -> void`
Sets the maximum time allowed for the graph matching algorithm to run.
Receiver: instance — call as `obj->SetTimeLimit(...)`
Params: `maxTimeToRun` — The maximum time in seconds to allow the search.
Side effects: Sets the internal time limit for matching.

#### `GraphMatch::MatchGraphs(PGraph Gh1, PGraph Gh2, int minMatch, bool vertexType=true, VERTEX_EXT_TYPE vertexExt=EXTTYPE_Ignore) -> void`
Looks for maximal common subgraphs of size not less than `minMatch` between `Gh1` and `Gh2`.
Receiver: instance — call as `obj->MatchGraphs(...)`
Params: `Gh1` — The first input graph.
`Gh2` — The second input graph.
`minMatch` — The minimum size of the common subgraph to report.
`vertexType` — If true, vertex type is considered during matching.
`vertexExt` — Key used to consider extended vertex types.
Side effects: Computes and stores all found matches; populates internal match counters.

#### `GraphMatch::GetNofMatches() -> int`
Returns the total count of common subgraphs found during the last matching attempt.
Receiver: instance — call as `obj->GetNofMatches()`
Returns: The number of found matches.
Side effects: None.

#### `GraphMatch::GetMatch(int MatchNo, ivector & FV1, ivector & FV2, int & nv, realtype & p1, realtype & p2) -> void`
Retrieves the vertex indices and corresponding geometric coordinates for a specific match number.
Receiver: instance — call as `obj->GetMatch(...)`
Params: `MatchNo` — The 1-based index of the match to retrieve.
`FV1` — Output vector to store vertex indices from the first graph.
`FV2` — Output vector to store vertex indices from the second graph.
`nv` — Output reference for the number of vertices in the match.
`p1` — Output reference for the centroid coordinates in the first graph.
`p2` — Output reference for the centroid coordinates in the second graph.
Side effects: Populates the provided output containers with match data.

#### `GraphMatch::SetTimeLimit(int maxTimeToRun) -> void`
Sets the time limit for the graph matching algorithm.
Receiver: instance — call as `obj->SetTimeLimit(...)`
Params: `maxTimeToRun` — The maximum time allowed for the search.
Side effects: Sets the internal time limit for matching.

---

## `Helix`
Represents a structural helix motif within a protein structure.

#### `Helix::PDBASCIIDump(pstr S, int N)`
Writes the helix coordinates and information into a PDB ASCII formatted string.
Receiver: instance — call as `obj->PDBASCIIDump(...)`
Params: `S` — Output string buffer. `N` — Number of atoms to dump.
Returns: void
Side effects: Writes to `S`.

#### `Helix::MakeCIF(mmcif::PData CIF, int N)`
Populates the provided CIF structure with the coordinates and details of the helix.
Receiver: instance — call as `obj->MakeCIF(...)`
Params: `CIF` — Pointer to the output CIF data structure. `N` — Number of atoms included in the CIF.
Returns: void
Side effects: Mutates the CIF data structure.

#### `Helix::ConvertPDBASCII(cpstr S)`
Attempts to parse helix coordinates and structure from a PDB ASCII formatted string.
Receiver: instance — call as `obj->ConvertPDBASCII(...)`
Params: `S` — Input string containing PDB format data.
Returns: ERROR_CODE
Side effects: Updates the helix's internal structural state.

#### `Helix::GetCIF(mmcif::PData CIF, int & n)`
Writes the helix coordinates and information into the provided CIF structure.
Receiver: instance — call as `obj->GetCIF(...)`
Params: `CIF` — Pointer to the output CIF data structure. `n` — Reference to store the number of atoms written.
Returns: ERROR_CODE
Side effects: Writes data to `CIF` and updates `n`.

#### `Helix::GetClassID()`
Retrieves the unique identifier for the Helix class.
Receiver: instance — call as `obj->GetClassID()`
Returns: CLASS_ID
Side effects: none

#### `Helix::Copy(PContainerClass Helix)`
Copies the contents of another container object representing a helix structure.
Receiver: instance — call as `obj->Copy(...)`
Params: `Helix` — The source helix object to copy from.
Returns: void
Side effects: Copies data from the source helix to this object.

#### `Helix::write(io::RFile f)`
Writes the helix structure data to the specified output file stream.
Receiver: instance — call as `obj->write(...)`
Params: `f` — The file stream to write to.
Returns: void
Side effects: Writes to file stream `f`.

#### `Helix::read(io::RFile f)`
Reads helix structure data from the specified input file stream.
Receiver: instance — call as `obj->read(...)`
Params: `f` — The file stream to read from.
Returns: void
Side effects: Reads data from file stream `f` and updates object state.

---

## `Journal`
Manages the persistent storage and IO operations for a structural biology journal file.

#### `Journal::PDBASCIIDump(pstr S, int N) -> void`
Writes the journal contents to a specified file buffer in PDB ASCII format.
Receiver: instance — call as `obj->PDBASCIIDump(...)`
Params: `S` — The string buffer to write to. `N` — The number of characters to write.
Side effects: Writes data to the provided string buffer `S`.

#### `Journal::PDBASCIIDump1(io::RFile) -> bool`
Writes the journal contents to a specified file handle in PDB ASCII format.
Receiver: instance — call as `obj->PDBASCIIDump1(...)`
Params: None
Returns: A boolean indicating success or failure.
Side effects: Writes data to the specified file handle.

#### `Journal::ConvertPDBASCII(cpstr S) -> ERROR_CODE`
Parses and converts data from a PDB ASCII formatted string into the journal object's internal representation.
Receiver: instance — call as `obj->ConvertPDBASCII(...)`
Params: `S` — The source string containing PDB ASCII data.
Returns: An enumeration code indicating the success or type of conversion.
Side effects: Updates the object's internal state based on the parsed data.

#### `Journal::GetClassID() -> CLASS_ID`
Retrieves the unique identifier associated with the Journal class type.
Receiver: instance — call as `obj->GetClassID()`
Returns: The `CLASS_ID` constant for `Journal`.

#### `Journal::write(io::RFile f) -> void`
Writes the current journal object's state to the specified output file handle.
Receiver: instance — call as `obj->write(...)`
Params: None
Side effects: Writes the object's data to the file handle `f`.

#### `Journal::read(io::RFile f) -> void`
Reads journal data from the specified input file handle, populating the object's state.
Receiver: instance — call as `obj->read(...)`
Params: None
Side effects: Updates the object's internal state by reading from file handle `f`.

---

## `Link`
Represents a structural link connecting two specific atoms in a molecular structure, defined by their coordinates and connectivity information.

#### `Link::Copy(PContainerClass link) -> void`
Copies the contents of another container link into this link object.
Receiver: instance — call as `obj->Copy(...)` (or `obj.Copy(...)`).
Params: `link` — Pointer to the link object to copy from.
Side effects: Mutates the state of the current `Link` object.
Errors: None.

---

## `LinkContainer`
A container class used to manage and assemble structural components (links) within a larger model context.

#### `LinkContainer::MakeContainerClass(int ClassID) -> PContainerClass`
Creates and returns a new `PContainerClass` instance associated with the specified structural component ID.
Receiver: instance — call as `obj->MakeContainerClass(...)`
Params: `ClassID` — The unique identifier for the structural component.
Returns: A pointer to the newly created container class.
Side effects: None.
Errors: Returns null if the specified `ClassID` is invalid.

---

## `LinkR`
Represents a physical link connection between two atoms within a molecular structure, storing associated coordinates and metadata.

#### `LinkR::PDBASCIIDump(pstr S, int N)`
Dumps the link information into a PDB ASCII format string.
Receiver: instance — call as `obj->PDBASCIIDump(...)`
Params: `S` — The string buffer to write to. `N` — The maximum number of links to write.
Side effects: Writes link data to string `S`.

#### `LinkR::MakeCIF(mmcif::PData CIF, int N)`
Writes the link information into the provided MMCIF data structure.
Receiver: instance — call as `obj->MakeCIF(...)`
Params: `CIF` — Pointer to the target MMCIF data object. `N` — The number of links to write.
Side effects: Modifies the `CIF` object.

#### `LinkR::ConvertPDBASCII(cpstr S)`
Attempts to convert structural link information from a PDB ASCII string.
Receiver: instance — call as `obj->ConvertPDBASCII(...)`
Params: `S` — The PDB ASCII string containing link data.
Returns: An error code indicating conversion success or failure.
Side effects: Updates the link object's state from the string data.

#### `LinkR::GetCIF(mmcif::PData CIF, int & n)`
Writes the link information into the provided MMCIF data structure.
Receiver: instance — call as `obj->GetCIF(...)`
Params: `CIF` — Pointer to the target MMCIF data object. `n` — Reference to an integer to store the count of written links.
Returns: An error code indicating writing success or failure.
Side effects: Modifies the `CIF` object and updates `n`.

#### `LinkR::GetClassID()`
Retrieves the unique identifier for the LinkR class type.
Receiver: instance — call as `obj->GetClassID()`
Returns: The unique class identifier for LinkR.

#### `LinkR::Copy(PContainerClass LinkR)`
Copies the state of another LinkR object into this object.
Receiver: instance — call as `obj->Copy(...)`
Params: `LinkR` — Pointer to the source LinkR object to copy from.
Side effects: Updates the internal state of `obj` to match the source object.

#### `LinkR::write(io::RFile f)`
Writes the link information content to a specified file stream.
Receiver: instance — call as `obj->write(...)`
Params: `f` — The file stream to write to.
Side effects: Writes link data to the file stream `f`.

#### `LinkR::read(io::RFile f)`
Reads link information from a specified file stream, populating this object.
Receiver: instance — call as `obj->read(...)`
Params: `f` — The file stream to read from.
Side effects: Reads link data from the file stream `f` and populates the object's state.

---

## `Loop`
Manages structural data organized in a grid format, allowing tags (columns) and rows of values to be set, retrieved, and manipulated.

#### `Loop::GetInteger(int& I, cpstr TName, int nrow, bool Remove) -> int`
Fetches the integer value corresponding to the given tag in the specified row.
Receiver: instance — call as `obj->GetInteger(...)`
Params: `I` — Reference to accept the fetched integer value.
`TName` — Character string with the tag name.
`nrow` — Row number (0 to GetLoopLength()-1).
`Remove` — If true, removes the field after reading.
Returns: Returns a CIFRC code indicating success or failure (e.g., CIFRC_NoTag, CIFRC_WrongIndex).
Side effects: If `Remove` is true, the field is deleted from the loop.
Errors: Returns CIFRC_NoTag if the tag is missing, CIFRC_WrongIndex if the row is out of bounds, etc.

#### `Loop::GetLoopLength() -> int`
Returns the current number of data rows available in the loop.
Receiver: instance — call as `obj->GetLoopLength()`
Returns: The total number of data rows.

#### `Loop::GetString(pstr & S, cpstr TName, int nrow, bool Remove) -> int`
Fetches a string value associated with the given tag and row into the provided buffer.
Receiver: instance — call as `obj->GetString(...)`
Params: `S` — Pointer to a string buffer to receive the value.
`TName` — Character string with the tag name.
`nrow` — Row number (0 to GetLoopLength()-1).
`Remove` — If true, removes the field after reading.
Returns: A CIFRC code indicating success or failure.
Side effects: The string pointer `S` might be re-used and set to NULL if data is not found or available.
Errors: Returns codes indicating tag not found, wrong row index, or value not found.

#### `Loop::GetString(cpstr TName, int nrow, int & RC) -> pstr`
Fetches the string value associated with the given tag and row.
Receiver: instance — call as `obj->GetString(...)`
Params: `TName` — Character string with the tag name.
`nrow` — Row number (0 to GetLoopLength()-1).
`RC` — Reference to receive the return code.
Returns: A pointer to the string value, or NULL if the value is absent or data is missing.
Side effects: None.
Errors: Read return code `RC` to interpret NULL return (e.g., CIFRC_NoTag).

#### `Loop::CopyString(pstr buf, int maxlength, cpstr TName, int nrow, int & RC) -> void`
Copies the string value associated with the given tag and row into a specified buffer.
Receiver: instance — call as `obj->CopyString(...)`
Params: `buf` — Pointer to the character string buffer to accept the value.
`maxlength` — Maximum number of bytes to copy into the buffer.
`TName` — Character string with the tag name.
`nrow` — Row number (0 to GetLoopLength()-1).
`RC` — Reference to receive the return code.
Side effects: The buffer `buf` is modified with the value; the return code `RC` is set.
Errors: The return code `RC` indicates if the tag, row, or field was not found.

#### `Loop::DeleteField(cpstr TName, int nrow) -> int`
Deletes the field associated with the given tag in the specified row.
Receiver: instance — call as `obj->DeleteField(...)`
Params: `TName` — Character string with the tag name to delete.
`nrow` — Row number (0 to GetLoopLength()-1).
Returns: Returns >=0 if the field was deleted, or <0 if the field or tag was not found.
Side effects: Deletes the field data at the specified location.
Errors: Returns a negative value if the specified field or tag does not exist.

#### `Loop::DeleteRow(int nrow) -> int`
Deletes all fields contained within the specified row.
Receiver: instance — call as `obj->DeleteRow(...)`
Params: `nrow` — Row number (0 to GetLoopLength()-1).
Returns: Returns CIFRC_Ok on success, or CIFRC_WrongIndex if the row is not found.
Side effects: Clears all field data for the specified row.
Errors: Returns CIFRC_WrongIndex if the row index is invalid.

#### `Loop::GetReal(realtype & R, cpstr TName, int nrow, bool Remove) -> int`
Fetches the real number value corresponding to the given tag in the specified row.
Receiver: instance — call as `obj->GetReal(...)`
Params: `R` — Reference to receive the fetched real value.
`TName` — Character string with the tag name.
`nrow` — Row number (0 to GetLoopLength()-1).
`Remove` — If true, removes the field after reading.
Returns: A CIFRC code indicating success or failure.
Side effects: If `Remove` is true, the field is deleted from the loop.
Errors: Returns various CIFRC codes for missing data, wrong format, or indices.

#### `Loop::CopyReal(realtype & R, cpstr TName, int nrow, int & RC) -> void`
Copies the real number value associated with the given tag and row into a specified reference.
Receiver: instance — call as `obj->CopyReal(...)`
Params: `R` — Reference to receive the copied real value.
`TName` — Character string with the tag name.
`nrow` — Row number (0 to GetLoopLength()-1).
`RC` — Reference to receive the return code.
Side effects: Sets `R` to 0 if the operation fails.
Errors: The return code `RC` indicates if the tag, row, or field was not found.

#### `Loop::CopyInteger(int & I, cpstr TName, int nrow, int & RC) -> void`
Copies the integer value associated with the given tag and row into a specified reference.
Receiver: instance — call as `obj->CopyInteger(...)`
Params: `I` — Reference to receive the copied integer value.
`TName` — Character string with the tag name.
`nrow` — Row number (0 to GetLoopLength()-1).
`RC` — Reference to receive the return code.
Side effects: Sets `I` to 0 if the operation fails.
Errors: The return code `RC` indicates if the tag, row, or field was not found.

#### `Loop::GetInteger(int & I, cpstr TName, int nrow, bool Remove) -> int`
Fetches the integer value corresponding to the given tag in the specified row.
Receiver: instance — call as `obj->GetInteger(...)`
Params: `I` — Reference to receive the fetched integer value.
`TName` — Character string with the tag name.
`nrow` — Row number (0 to GetLoopLength()-1).
`Remove` — If true, removes the field after reading.
Returns: A CIFRC code indicating success or failure.
Side effects: If `Remove` is true, the field is deleted from the loop.
Errors: Returns various CIFRC codes for missing data, wrong format, or indices.

#### `Loop::GetSVector(psvector & S, cpstr TName, int i1, int i2, bool Remove) -> int`
Fetches a set of string values corresponding to the given tag in a specified row range.
Receiver: instance — call as `obj->GetSVector(...)`
Params: `S` — Reference to the string vector to accept the values.
`TName` — Character string with the tag name.
`i1` — Minimum row number to fetch.
`i2` — Maximum row number to fetch.
`Remove` — If true, removes fetched fields after reading.
Returns: A CIFRC code indicating success or failure.
Side effects: The contents of `S` are populated with fetched data.
Errors: Returns codes indicating tag not found or invalid row range.

#### `Loop::GetRVector(rvector & R, cpstr TName, int i1, int i2, bool Remove) -> int`
Fetches a set of real-point values corresponding to the given tag in a specified row range.
Receiver: instance — call as `obj->GetRVector(...)`
Params: `R` — Reference to the vector to accept the real values.
`TName` — Character string with the tag name.
`i1` — Minimum row number to fetch.
`i2` — Maximum row number to fetch.
`Remove` — If true, removes fetched fields after reading.
Returns: A CIFRC code indicating success or failure.
Side effects: The contents of `R` are populated with fetched data.
Errors: Returns codes indicating tag not found or invalid row range.

#### `Loop::GetIVector(ivector & I, cpstr TName, int i1, int i2, bool Remove) -> int`
Fetches a set of integer values corresponding to the given tag in a specified row range.
Receiver: instance — call as `obj->GetIVector(...)`
Params: `I` — Reference to the integer vector to accept the values.
`TName` — Character string with the tag name.
`i1` — Minimum row number to fetch.
`i2` — Maximum row number to fetch.
`Remove` — If true, removes fetched fields after reading.
Returns: A CIFRC code indicating success or failure.
Side effects: The contents of `I` are populated with fetched data.
Errors: Returns codes indicating tag not found or invalid row range.

#### `Loop::PutString(cpstr S, cpstr T, int nrow) -> void`
Sets a string value for a given tag and row, creating the tag or row if they do not exist.
Receiver: instance — call as `obj->PutString(...)`
Params: `S` — Character string with the value to set.
`T` — Character string with the tag name.
`nrow` — Row number to set the data.
Side effects: Mutates the loop structure by adding tags or rows if necessary, potentially reindexing the loop.
Errors: None specified.

#### `Loop::PutNoData(int NoDataType, cpstr T, int nrow) -> void`
Sets a specific data absence value for a given tag and row, creating the tag or row if they do not exist.
Receiver: instance — call as `obj->PutNoData(...)`
Params: `NoDataType` — Integer key specifying the absence type (CIF_NODATA_DOT or CIF_NODATA_QUESTION).
`T` — Character string with the tag name.
`nrow` — Row number to set the data.
Side effects: Mutates the loop structure by adding tags or rows if necessary, potentially reindexing the loop.
Errors: None specified.

#### `Loop::PutReal(realtype R, cpstr T, int nrow, int prec) -> void`
Sets a real-point value for a given tag and row, creating the tag or row if they do not exist.
Receiver: instance — call as `obj->PutReal(...)`
Params: `R` — Real number with value to be set.
`T` — Character string with the tag name.
`nrow` — Row number to set the data.
`prec` — Float-point precision to use in formatting.
Side effects: Mutates the loop structure by adding tags or rows if necessary, potentially reindexing the loop.
Errors: None specified.

#### `Loop::PutReal(realtype R, cpstr T, int nrow, cpstr format) -> void`
Sets a real-point value for a given tag and row using a specified format string, creating the tag or row if they do not exist.
Receiver: instance — call as `obj->PutReal(...)`
Params: `R` — Real number with value to be set.
`T` — Character string with the tag name.
`nrow` — Row number to set the data.
`format` — Format string to convert the real number.
Side effects: Mutates the loop structure by adding tags or rows if necessary, potentially reindexing the loop.
Errors: None specified.

#### `Loop::PutInteger(int I, cpstr T, int nrow) -> void`
Sets an integer value for a given tag and row, creating the tag or row if they do not exist.
Receiver: instance — call as `obj->PutInteger(...)`
Params: `I` — Integer number with value to be set.
`T` — Character string with the tag name.
`nrow` — Row number to set the data.
Side effects: Mutates the loop structure by adding tags or rows if necessary, potentially reindexing the loop.
Errors: None specified.

#### `Loop::PutSVector(psvector S, cpstr T, int i1, int i2) -> void`
Sets a set of string values for the given tag across a specified range of rows.
Receiver: instance — call as `obj->PutSVector(...)`
Params: `S` — String vector containing values to store.
`T` — Character string with the tag name.
`i1` — Minimum row index (inclusive).
`i2` — Maximum row index (inclusive).
Side effects: Mutates the loop structure by adding the tag if missing, and populates specified rows. Missing rows are created and initialized.
Errors: None specified.

#### `Loop::PutRVector(rvector R, cpstr T, int i1, int i2, int prec) -> void`
Sets a set of real-point values for the given tag across a specified range of rows.
Receiver: instance — call as `obj->PutRVector(...)`
Params: `R` — Vector of real numbers to store.
`T` — Character string with the tag name.
`i1` — Minimum row index (inclusive).
`i2` — Maximum row index (inclusive).
`prec` — Float-point precision to use in formatting.
Side effects: Mutates the loop structure by adding the tag if missing, and populates specified rows. Missing rows are created and initialized.
Errors: None specified.

#### `Loop::PutIVector(ivector I, cpstr T, int i1, int i2) -> void`
Sets a set of integer values for the given tag across a specified range of rows.
Receiver: instance — call as `obj->PutIVector(...)`
Params: `I` — Vector of integers to store.
`T` — Character string with the tag name.
`i1` — Minimum row index (inclusive).
`i2` — Maximum row index (inclusive).
Side effects: Mutates the loop structure by adding the tag if missing, and populates specified rows. Missing rows are created and initialized.
Errors: None specified.

#### `Loop::Optimize() -> void`
Optimizes the loop's internal data structures for better memory usage and data access speed.
Receiver: instance — call as `obj->Optimize()`
Side effects: Reorders and compacts the underlying data structures; the order is preserved based on appearance order.
Errors: None specified.

#### `Loop::WriteMMCIFLoop(cpstr FName, io::GZ_MODE gzipMode) -> bool`
Writes the loop's data content into a file in the mmCIF format.
Receiver: instance — call as `obj->WriteMMCIFLoop(...)`
Params: `FName` — Character string with the output file name.
`gzipMode` — Flag controlling compression of the output file.
Returns: True if writing was successful, false otherwise.
Side effects: Writes data to the specified file on disk.
Errors: Returns false if the file cannot be opened for writing.

#### `Loop::WriteMMCIF(io::RFile f) -> void`
Writes the loop's data content into the provided file stream.
Receiver: instance — call as `obj->WriteMMCIF(...)`
Params: `f` — Reference to the MMDB file stream object.
Side effects: Writes data to the provided file stream.
Errors: None specified.

#### `Loop::Copy(PCategory Loop) -> void`
Creates a deep copy of the loop's contents from another loop object.
Receiver: instance — call as `obj->Copy(...)`
Params: `Loop` — A pointer to the source loop object to copy data from.
Side effects: Duplicates all data and memory allocated within the loop structure.
Errors: None specified.

#### `Loop::write(io::RFile f) -> void`
Writes the loop's contents to the provided file stream using the MMDB stream mechanism.
Receiver: instance — call as `obj->write(...)`
Params: `f` — Reference to the MMDB file stream object.
Side effects: Writes data to the provided file stream.
Errors: None specified.

#### `Loop::read(io::RFile f) -> void`
Reads loop data from the provided file stream using the MMDB stream mechanism.
Receiver: instance — call as `obj->read(...)`
Params: `f` — Reference to the MMDB file stream object.
Side effects: Reads and populates the loop data from the stream.
Errors: None specified.

#### `Loop::InitLoop() -> void`
Initializes the loop's internal data structures.
Receiver: instance — call as `obj->InitLoop()`
Side effects: Sets up necessary internal data members (e.g., counts, memory).
Errors: None specified.

#### `Loop::FreeMemory() -> void`
Releases dynamically allocated memory associated with the loop's data fields.
Receiver: instance — call as `obj->FreeMemory()`
Side effects: Deallocates memory used by stored data strings and values.
Errors: None specified.

#### `Loop::DeleteFields() -> void`
Cleans up the stored field data within the loop structure.
Receiver: instance — call as `obj->DeleteFields()`
Side effects: Removes the actual data records from the loop without changing the row or tag structure.
Errors: None specified.

#### `Loop::ExpandRows(int nRowsNew) -> void`
Increases the number of rows supported by the loop structure.
Receiver: instance — call as `obj->ExpandRows(...)`
Params: `nRowsNew` — The new total number of rows required.
Side effects: Allocates and initializes memory for the specified number of rows.
Errors: None specified.

---

## `Manager`
Manages the structural data, coordination, and metadata for a molecular structure file.

#### `Manager::Copy(PManager MMDB, COPY_MASK CopyMask)` -> void
Transfers selected information between the current MMDB and the specified target MMDB.
Receiver: instance — call as `obj->Copy(...)` (or `obj.Copy(...)`).
Params: `MMDB` — The target MMDB to receive copied data. `CopyMask` — A bitmask defining which types of data are copied.
Returns: nothing.
Side effects: Replaces existing content in the target MMDB according to the copy mask.
Errors: None.

#### `Manager::Delete(word DelMask)` -> void
Removes structural information from the current MMDB based on the provided mask.
Receiver: instance — call as `obj->Delete(...)` (or `obj.Delete(...)`).
Params: `DelMask` — A bitmask defining which types of data should be deleted.
Returns: nothing.
Side effects: Mutates the internal state of the MMDB by deleting data.
Errors: None.

#### `Manager::GetRemarks()` -> PTitleContainer
Retrieves the title remarks metadata associated with the MMDB.
Receiver: instance — call as `obj->GetRemarks()` (or `obj.GetRemarks()`).
Params: None.
Returns: A pointer to a container holding structural remarks.
Side effects: None.
Errors: None.

---

## Model
Manages the structure and coordinates of a macromolecular model loaded from a PDB-like file format.

#### `Model::GetCoordHierarchy() -> PManager`
Returns the coordinate management structure associated with the model.
Receiver: instance — call as `obj->GetCoordHierarchy()`
Returns: A pointer to the coordination hierarchy manager.
Side effects: None.

#### `Model::GetChainCreate(const ChainID chID, bool enforceUniqueChainID) -> PChain`
Creates or retrieves a chain with the specified ID, handling potential uniqueness conflicts by appending serial numbers if necessary.
Receiver: instance — call as `obj->GetChainCreate(...)`
Params: `chID` — The desired chain identifier. `enforceUniqueChainID` — If true, guarantees uniqueness by modifying the ID if a conflict exists.
Returns: A pointer to the chain object.
Side effects: May modify the model by adding new chain components if required.

#### `Model::GetSerNum() -> int`
Retrieves the unique serial number assigned to this model within its containing file/context.
Receiver: instance — call as `obj->GetSerNum()`
Returns: The model's serial number.
Side effects: None.

#### `Model::GetModelID(pstr modelID) -> cpstr`
Returns the model's unique identifier string, typically derived from PDB header records.
Receiver: instance — call as `obj->GetModelID(...)`
Returns: The model's identifier as a string.
Side effects: None.

#### `Model::GetNumberOfChains() -> int`
Returns the total count of distinct chains present in the model.
Receiver: instance — call as `obj->GetNumberOfChains()`
Returns: The number of chains.
Side effects: None.

#### `Model::GetChain(const ChainID chID) -> PChain`
Retrieves a pointer to the chain identified by `chID`, returning NULL if no such chain exists in the model.
Receiver: instance — call as `obj->GetChain(...)`
Params: `chID` — The desired chain identifier.
Returns: A pointer to the chain object, or NULL if not found.
Side effects: None.

#### `Model::GetChain(int chainNo) -> PChain`
Retrieves a pointer to the chain located at the given zero-based index in the model's chain list.
Receiver: instance — call as `obj->GetChain(...)`
Params: `chainNo` — The zero-based index of the chain (0 <= chainNo < nChains).
Returns: A pointer to the chain object.
Side effects: None.

#### `Model::GetChainTable(PPChain & chainTable, int & NumberOfChains) -> void`
Fills the provided chain array and updates the count to list all chains currently in the model.
Receiver: instance — call as `obj->GetChainTable(...)`
Params: `chainTable` — A pointer to an array that will receive chain pointers. `NumberOfChains` — An integer reference updated with the chain count.
Returns: void
Side effects: Populates the `chainTable` array and sets `NumberOfChains`.

#### `Model::DeleteChain(const ChainID chID) -> int`
Removes the chain identified by `chID` from the model, returning the number of atoms deleted.
Receiver: instance — call as `obj->DeleteChain(...)`
Params: `chID` — The identifier of the chain to delete.
Returns: The number of atoms removed from the model.
Side effects: Removes all atoms belonging to the specified chain.

#### `Model::DeleteChain(int chainNo) -> int`
Removes the chain located at the given zero-based index from the model, returning the number of atoms deleted.
Receiver: instance — call as `obj->DeleteChain(...)`
Params: `chainNo` — The zero-based index of the chain to delete.
Returns: The number of atoms removed from the model.
Side effects: Removes all atoms belonging to the specified chain.

#### `Model::AddChain(PChain chn) -> int`
Adds a specified chain to the model, returning the number of atoms added.
Receiver: instance — call as `obj->AddChain(...)`
Params: `chn` — Pointer to the chain object to be added.
Returns: The number of atoms added to the model.
Side effects: Incorporates the chain's coordinates into the model.

#### `Model::SortChains(int sortKey) -> void`
Sorts the model's constituent chains based on the criteria specified by `sortKey`.
Receiver: instance — call as `obj->SortChains(...)`
Params: `sortKey` — An integer specifying the sorting key (e.g., SORT_CHAIN_NAME).
Returns: void
Side effects: Reorders the internal list of chains.

#### `Model::GetNumberOfHelices() -> int`
Returns the total count of defined helical secondary structure elements in the model.
Receiver: instance — call as `obj->GetNumberOfHelices()`
Returns: The number of helices.
Side effects: None.

#### `Model::GetNumberOfSheets() -> int`
Returns the total count of defined beta-sheet secondary structure elements in the model.
Receiver: instance — call as `obj->GetNumberOfSheets()`
Returns: The number of sheets.
Side effects: None.

#### `Model::GetHelix(int serialNum) -> PHelix`
Retrieves a pointer to the helix element at the specified one-based serial number.
Receiver: instance — call as `obj->GetHelix(...)`
Params: `serialNum` — The one-based index of the helix (1 <= serialNum <= NofHelices).
Returns: A pointer to the helix object.
Side effects: None.

#### `Model::GetSheetID(int serialNum, SheetID sheetID) -> void`
Retrieves the unique identifier associated with the sheet at the given serial number.
Receiver: instance — call as `obj->GetSheetID(...)`
Params: `serialNum` — The one-based index of the sheet. `sheetID` — An identifier object to store the result.
Returns: void
Side effects: Sets the provided `sheetID` parameter with the sheet's ID.

#### `Model::GetSheet(int serialNum) -> PSheet`
Retrieves a pointer to the sheet element at the specified one-based serial number.
Receiver: instance — call as `obj->GetSheet(...)`
Params: `serialNum` — The one-based index of the sheet (1 <= serialNum <= NofSheets).
Returns: A pointer to the sheet object.
Side effects: None.

#### `Model::GetSheet(const SheetID sheetID) -> PSheet`
Retrieves a pointer to the sheet element matching the provided unique ID.
Receiver: instance — call as `obj->GetSheet(...)`
Params: `sheetID` — The identifier of the sheet.
Returns: A pointer to the sheet object, or NULL if not found.
Side effects: None.

#### `Model::GetSheets() -> PSheetsContainer`
Returns a container object allowing access to all defined sheet secondary structure elements.
Receiver: instance — call as `obj->GetSheets()`
Returns: A pointer to the sheet container.
Side effects: None.

#### `Model::RemoveSecStructure() -> void`
Clears all pre-defined secondary structure assignments (helices, sheets) from the model.
Receiver: instance — call as `obj->RemoveSecStructure()`
Returns: void
Side effects: Clears secondary structure data members.

#### `Model::CalcSecStructure(bool flagBulge=true, int aminoSelHnd=-1) -> int`
Calculates and assigns secondary structure assignments (helix/sheet) to the model atoms based on structural analysis.
Receiver: instance — call as `obj->CalcSecStructure(...)`
Params: `flagBulge` — Whether to consider loop/bulge regions during calculation. `aminoSelHnd` — An amino acid selection handle to limit calculation scope.
Returns: An integer indicating success or failure.
Side effects: Populates secondary structure data structures.

#### `Model::GetNumberOfLinks() -> int`
Returns the total count of defined structural links connecting residues.
Receiver: instance — call as `obj->GetNumberOfLinks()`
Returns: The number of links.
Side effects: None.

#### `Model::GetLink(int serialNum) -> PLink`
Retrieves a pointer to the link element at the specified one-based serial number.
Receiver: instance — call as `obj->GetLink(...)`
Params: `serialNum` — The one-based index of the link (1 <= serialNum <= NofLinks).
Returns: A pointer to the link object.
Side effects: None.

#### `Model::GetLinks() -> PLinkContainer`
Returns a container object allowing access to all defined link secondary structure elements.
Receiver: instance — call as `obj->GetLinks()`
Returns: A pointer to the link container.
Side effects: None.

#### `Model::RemoveLinks() -> void`
Removes all structural links from the model.
Receiver: instance — call as `obj->RemoveLinks()`
Returns: void
Side effects: Clears link data members.

#### `Model::AddLink(PLink link) -> void`
Adds a specified structural link to the model.
Receiver: instance — call as `obj->AddLink(...)`
Params: `link` — Pointer to the link object to be added.
Returns: void
Side effects: Incorporates the link's coordinates into the model.

#### `Model::GetNumberOfLinkRs() -> int`
Returns the total count of Refmac link elements associated with the model.
Receiver: instance — call as `obj->GetNumberOfLinkRs()`
Returns: The number of Refmac links.
Side effects: None.

#### `Model::GetLinkR(int serialNum) -> PLinkR`
Retrieves a pointer to the Refmac link element at the specified one-based serial number.
Receiver: instance — call as `obj->GetLinkR(...)`
Params: `serialNum` — The one-based index of the linkR (1 <= serialNum <= NofLinks).
Returns: A pointer to the Refmac link object.
Side effects: None.

#### `Model::GetLinkRs() -> PLinkRContainer`
Returns a container object allowing access to all defined Refmac link secondary structure elements.
Receiver: instance — call as `obj->GetLinkRs()`
Returns: A pointer to the linkR container.
Side effects: None.

#### `Model::RemoveLinkRs() -> void`
Removes all Refmac links from the model.
Receiver: instance — call as `obj->RemoveLinkRs()`
Returns: void
Side effects: Clears Refmac link data members.

#### `Model::AddLinkR(PLinkR linkR) -> void`
Adds a specified Refmac link to the model.
Receiver: instance — call as `obj->AddLinkR(...)`
Params: `linkR` — Pointer to the Refmac link object to be added.
Returns: void
Side effects: Incorporates the linkR's coordinates into the model.

#### `Model::GetNumberOfCisPeps() -> int`
Returns the total count of defined cis-peptide bond features in the model.
Receiver: instance — call as `obj->GetNumberOfCisPeps()`
Returns: The number of cis-peptides.
Side effects: None.

#### `Model::GetCisPep(int CisPepNum) -> PCisPep`
Retrieves a pointer to the cis-peptide bond feature at the specified one-based serial number.
Receiver: instance — call as `obj->GetCisPep(...)`
Params: `CisPepNum` — The one-based index of the cis-peptide (1 <= CisPepNum <= NofCisPeps).
Returns: A pointer to the cis-peptide object.
Side effects: None.

#### `Model::GetCisPeps() -> PCisPepContainer`
Returns a container object allowing access to all defined cis-peptide bond features.
Receiver: instance — call as `obj->GetCisPeps()`
Returns: A pointer to the cis-peptide container.
Side effects: None.

#### `Model::RemoveCisPeps() -> void`
Removes all cis-peptide bond definitions from the model.
Receiver: instance — call as `obj->RemoveCisPeps()`
Returns: void
Side effects: Clears cis-peptide data members.

#### `Model::AddCisPep(PCisPep cisPep) -> void`
Adds a specified cis-peptide bond feature to the model.
Receiver: instance — call as `obj->AddCisPep(...)`
Params: `cisPep` — Pointer to the cis-peptide object to be added.
Returns: void
Side effects: Incorporates the cis-peptide's definition into the model.

#### `Model::Copy(PModel model) -> void`
Creates a deep copy of the structural data from another model object.
Receiver: instance — call as `obj->Copy(...)`
Params: `model` — Pointer to the model to be copied.
Returns: void
Side effects: Copies all coordinates and structural features from the source model.

#### `Model::CopyHets(PModel model) -> void`
Copies the heterocompound information from a source model to this model.
Receiver: instance — call as `obj->CopyHets(...)`
Params: `model` — Pointer to the model containing heterocompounds to copy.
Returns: void
Side effects: Updates the model's heterocompound data structure.

#### `Model::CopySecStructure(PModel model) -> void`
Copies the secondary structure assignments (helices, sheets) from a source model to this model.
Receiver: instance — call as `obj->CopySecStructure(...)`
Params: `model` — Pointer to the model containing secondary structure to copy.
Returns: void
Side effects: Updates the model's secondary structure data structures.

#### `Model::CopyLinks(PModel model) -> void`
Copies all structural links from a source model to this model.
Receiver: instance — call as `obj->CopyLinks(...)`
Params: `model` — Pointer to the model containing links to copy.
Returns: void
Side effects: Updates the model's link data structure.

#### `Model::CopyLinkRs(PModel model) -> void`
Copies all Refmac links from a source model to this model.
Receiver: instance — call as `obj->CopyLinkRs(...)`
Params: `model` — Pointer to the model containing Refmac links to copy.
Returns: void
Side effects: Updates the model's Refmac link data structure.

#### `Model::CopyCisPeps(PModel model) -> void`
Copies all cis-peptide bond definitions from a source model to this model.
Receiver: instance — call as `obj->CopyCisPeps(...)`
Params: `model` — Pointer to the model containing cis-peptides to copy.
Returns: void
Side effects: Updates the model's cis-peptide data structure.

#### `Model::PDBASCIIDumpPS(io::RFile f) -> void`
Writes the primary structure records (excluding cispeps) from the model in PDB ASCII format.
Receiver: instance — call as `obj->PDBASCIIDumpPS(...)`
Params: `f` — The output file stream for writing data.
Returns: void
Side effects: Writes content to the provided file stream `f`.

---

## `Remark`
Represents a single remark record found within a structural file format like PDB.

#### `Remark::PDBASCIIDump(pstr S, int N) -> void`
Writes the remark contents to a given ASCII string stream buffer.
Receiver: instance — call as `obj->PDBASCIIDump(...)`
Params: `S` — The string stream to write to. `N` — The number of records to process.
Side effects: Writes data to `S`.

#### `Remark::MakeCIF(mmcif::PData CIF, int N) -> void`
Populates the structure's remark section within a CIF data object.
Receiver: instance — call as `obj->MakeCIF(...)`
Params: `CIF` — The CIF data object to modify. `N` — The number of records to process.
Side effects: Mutates the passed `CIF` object.

#### `Remark::ConvertPDBASCII(cpstr S) -> ERROR_CODE`
Attempts to parse and populate the remark from a PDB ASCII string representation.
Receiver: instance — call as `obj->ConvertPDBASCII(...)`
Params: `S` — The PDB ASCII string containing the remark.
Returns: An error code indicating success or failure.
Side effects: Reads data from `S` and updates the object's state.

#### `Remark::GetCIF(mmcif::PData CIF, int & n) -> ERROR_CODE`
Retrieves the remark information from the structure and populates a CIF data object.
Receiver: instance — call as `obj->GetCIF(...)`
Params: `CIF` — The CIF data object to populate. `n` — A reference integer to store the count of processed records.
Returns: An error code indicating success or failure.
Side effects: Mutates the passed `CIF` object.

#### `Remark::GetClassID() -> CLASS_ID`
Retrieves the predefined structural class identifier for a Remark record.
Receiver: instance — call as `obj->GetClassID()`
Returns: The `CLASS_ID` associated with Remark records.

#### `Remark::Copy(PContainerClass RemarkClass) -> void`
Copies the content of another container class object into this one.
Receiver: instance — call as `obj->Copy(...)`
Params: `RemarkClass` — The container class object to copy from.
Side effects: Updates the internal state of the current object.

#### `Remark::write(io::RFile f) -> void`
Writes the object's representation to a specified output file stream.
Receiver: instance — call as `obj->write(...)`
Params: `f` — The file stream to write to.
Side effects: Performs file I/O to the stream `f`.

#### `Remark::read(io::RFile f) -> void`
Reads remark data from a specified input file stream and populates the object.
Receiver: instance — call as `obj->read(...)`
Params: `f` — The file stream to read from.
Side effects: Performs file I/O from the stream `f` and updates the object's state.

---

## `Residue`
Represents a single amino acid or nucleotide residue within a protein or nucleic acid structure.

#### `Residue::SetResID(const ResName resName, int sqNum, const InsCode ins) -> void`
Sets the residue's full identifier using name, sequence number, and insertion code.
Receiver: instance — call as `obj->SetResID(...)`
Params: `resName` — The residue name. `sqNum` — The sequence number. `ins` — The insertion code.
Side effects: Updates the residue's identity fields.

#### `Residue::GetChain() -> PChain`
Returns the chain object associated with this residue.
Receiver: instance — call as `obj->GetChain()`
Params: None
Returns: A pointer to the residue's associated `Chain` object.
Side effects: None.

#### `Residue::GetModel() -> PModel`
Returns the model object containing this residue.
Receiver: instance — call as `obj->GetModel()`
Params: None
Returns: A pointer to the residue's containing `Model` object.
Side effects: None.

#### `Residue::GetModelNum() -> int`
Retrieves the model number associated with this residue.
Receiver: instance — call as `obj->GetModelNum()`
Params: None
Returns: The integer model number.
Side effects: None.

#### `Residue::GetChainID() -> pstr`
Retrieves a string representation of the chain identifier.
Receiver: instance — call as `obj->GetChainID()`
Params: None
Returns: A string containing the chain ID.
Side effects: None.

#### `Residue::GetResName() -> pstr`
Retrieves a string representation of the residue name.
Receiver: instance — call as `obj->GetResName()`
Params: None
Returns: A string containing the residue name.
Side effects: None.

#### `Residue::SetResName(const ResName resName) -> void`
Sets the residue's name.
Receiver: instance — call as `obj->SetResName(...)`
Params: `resName` — The new residue name.
Side effects: Updates the residue name.

#### `Residue::GetSeqNum() -> int`
Retrieves the sequence number of the residue.
Receiver: instance — call as `obj->GetSeqNum()`
Params: None
Returns: The residue sequence number.
Side effects: None.

#### `Residue::GetInsCode() -> pstr`
Retrieves a string representation of the residue insertion code.
Receiver: instance — call as `obj->GetInsCode()`
Params: None
Returns: A string containing the insertion code.
Side effects: None.

#### `Residue::GetAtomStatistics(RAtomStat AS) -> void`
Populates atomic statistics structures based on the residue's atoms.
Receiver: instance — call as `obj->GetAtomStatistics(...)`
Params: `AS` — The structure to store atomic statistics.
Side effects: Populates the provided statistic structure.

#### `Residue::GetResidueID(pstr ResidueID) -> pstr`
Checks if the residue matches the provided residue ID string and returns the matching ID if successful.
Receiver: instance — call as `obj->GetResidueID(...)`
Params: `ResidueID` — The residue ID string to check.
Returns: A string representing the residue ID if matched, or potentially empty otherwise.
Side effects: None.

#### `Residue::isAminoacid() -> bool`
Checks if the residue is an amino acid type.
Receiver: instance — call as `obj->isAminoacid()`
Params: None
Returns: True if the residue is an amino acid.
Side effects: None.

#### `Residue::isNucleotide() -> bool`
Checks if the residue is a nucleotide type.
Receiver: instance — call as `obj->isNucleotide()`
Params: None
Returns: True if the residue is a nucleotide.
Side effects: None.

#### `Residue::isSolvent() -> bool`
Checks if the residue is modeled as solvent.
Receiver: instance — call as `obj->isSolvent()`
Params: None
Returns: True if the residue is solvent.
Side effects: None.

#### `Residue::isInSelection(int selHnd) -> bool`
Checks if the residue is included in a specified selection handle.
Receiver: instance — call as `obj->isInSelection(...)`
Params: `selHnd` — The selection handle ID.
Returns: True if the residue is in the selection.
Side effects: None.

#### `Residue::isNTerminus() -> bool`
Checks if the residue is the N-terminus of a chain.
Receiver: instance — call as `obj->isNTerminus()`
Params: None
Returns: True if the residue is the N-terminus.
Side effects: None.

#### `Residue::isCTerminus() -> bool`
Checks if the residue is the C-terminus of a chain.
Receiver: instance — call as `obj->isCTerminus()`
Params: None
Returns: True if the residue is the C-terminus.
Side effects: None.

#### `Residue::GetNumberOfAtoms() -> int`
Counts the total number of atoms associated with this residue.
Receiver: instance — call as `obj->GetNumberOfAtoms()`
Params: None
Returns: The total count of atoms.
Side effects: None.

#### `Residue::GetNumberOfAtoms(bool countTers) -> int`
Counts the number of atoms, optionally excluding terminal atoms.
Receiver: instance — call as `obj->GetNumberOfAtoms(...)`
Params: `countTers` — If true, terminal atoms are not counted.
Returns: The count of atoms matching criteria.
Side effects: None.

#### `Residue::GetAtom(const AtomName aname, const Element elname=NULL, const AltLoc aloc=cpstr("")) -> PAtom`
Retrieves a specific atom by name, element, and alternative location.
Receiver: instance — call as `obj->GetAtom(...)`
Params: `aname` — The desired atom name. `elname` — The desired element type. `aloc` — The alternative location to retrieve.
Returns: A pointer to the requested `Atom` object, or null if not found.
Side effects: None.

#### `Residue::GetAtom(int atomNo) -> PAtom`
Retrieves an atom by its index number within the residue.
Receiver: instance — call as `obj->GetAtom(...)`
Params: `atomNo` — The zero-based index of the atom.
Returns: A pointer to the requested `Atom` object, or null if the index is out of bounds.
Side effects: None.

#### `Residue::GetAtomTable(PPAtom & atomTable, int & NumberOfAtoms) -> void`
Fills the provided atom table pointer with all atoms in the residue and updates the count.
Receiver: instance — call as `obj->GetAtomTable(...)`
Params: `atomTable` — Pointer to allocate and receive the array of atom pointers. `NumberOfAtoms` — Reference to receive the count of atoms.
Returns: None.
Side effects: Populates `atomTable` and sets `NumberOfAtoms`.

#### `Residue::GetAtomTable1(PPAtom & atomTable, int & NumberOfAtoms) -> void`
Fills the provided atom table pointer with all atoms, excluding TER atoms and null pointers, and updates the count.
Receiver: instance — call as `obj->GetAtomTable1(...)`
Params: `atomTable` — Pointer to allocate and receive the array of atom pointers. `NumberOfAtoms` — Reference to receive the count of atoms.
Returns: None.
Side effects: Populates `atomTable` and sets `NumberOfAtoms`.

#### `Residue::DeleteAtom(const AtomName aname, const Element elname=NULL, const AltLoc aloc=cpstr("")) -> int`
Deletes an atom specified by name, element, and alternative location.
Receiver: instance — call as `obj->DeleteAtom(...)`
Params: `aname` — The atom name to delete. `elname` — The element type to match. `aloc` — The alternative location to target.
Returns: The number of atoms successfully deleted (1 or 0).
Side effects: Removes the atom from the residue structure.

#### `Residue::DeleteAtom(int atomNo) -> int`
Deletes an atom specified by its index number.
Receiver: instance — call as `obj->DeleteAtom(...)`
Params: `atomNo` — The zero-based index of the atom to delete.
Returns: The number of atoms successfully deleted (1 or 0).
Side effects: Removes the atom from the residue structure.

#### `Residue::DeleteAllAtoms() -> int`
Deletes every atom associated with this residue.
Receiver: instance — call as `obj->DeleteAllAtoms()`
Params: None
Returns: The total number of atoms deleted.
Side effects: Clears the atom table of the residue.

#### `Residue::DeleteAltLocs() -> int`
Removes all alternative location copies, keeping only the primary location (or based on occupancy rules).
Receiver: instance — call as `obj->DeleteAltLocs()`
Params: None
Returns: The total number of atoms that were deleted.
Side effects: Trims the atom table and removes alternative location data.

#### `Residue::TrimAtomTable() -> void`
Cleans up internal structures to accurately reflect the current number of atoms.
Receiver: instance — call as `obj->TrimAtomTable()`
Params: None
Returns: None.
Side effects: Resets internal atom count trackers.

#### `Residue::AddAtom(PAtom atm) -> int`
Adds an external atom pointer to the residue, handling necessary copying and checking.
Receiver: instance — call as `obj->AddAtom(...)`
Params: `atm` — Pointer to the atom to be added.
Returns: The count of atoms in the residue after addition.
Side effects: Increases atom count; potentially copies atom data.

#### `Residue::PutUDData(int UDDhandle, int iudd) -> int`
Stores an integer value associated with a unique user data handle.
Receiver: instance — call as `obj->PutUDData(...)`
Params: `UDDhandle` — Handle identifying the data slot. `iudd` — The integer value to store.
Returns: An error code or success status.
Side effects: Associates integer data with the residue.

#### `Residue::PutUDData(int UDDhandle, realtype rudd) -> int`
Stores a real-valued floating-point number associated with a unique user data handle.
Receiver: instance — call as `obj->PutUDData(...)`
Params: `UDDhandle` — Handle identifying the data slot. `rudd` — The floating-point value to store.
Returns: An error code or success status.
Side effects: Associates float data with the residue.

#### `Residue::PutUDData(int UDDhandle, cpstr sudd) -> int`
Stores a string value associated with a unique user data handle.
Receiver: instance — call as `obj->PutUDData(...)`
Params: `UDDhandle` — Handle identifying the data slot. `sudd` — The string value to store.
Returns: An error code or success status.
Side effects: Associates string data with the residue.

#### `Residue::GetUDData(int UDDhandle, int & iudd) -> int`
Retrieves an integer value previously stored using the specified handle.
Receiver: instance — call as `obj->GetUDData(...)`
Params: `UDDhandle` — Handle identifying the data slot. `iudd` — Reference where the retrieved integer will be stored.
Returns: An error code or success status.
Side effects: None.

#### `Residue::GetUDData(int UDDhandle, realtype & rudd) -> int`
Retrieves a real-valued floating-point number previously stored using the specified handle.
Receiver: instance — call as `obj->GetUDData(...)`
Params: `UDDhandle` — Handle identifying the data slot. `rudd` — Reference where the retrieved float will be stored.
Returns: An error code or success status.
Side effects: None.

#### `Residue::GetUDData(int UDDhandle, pstr sudd, int maxLen) -> int`
Retrieves a string up to a specified maximum length previously stored using the handle.
Receiver: instance — call as `obj->GetUDData(...)`
Params: `UDDhandle` — Handle identifying the data slot. `sudd` — Pointer to buffer for the result. `maxLen` — Maximum length of the string to read.
Returns: An error code or success status.
Side effects: Writes data into the provided string buffer.

#### `Residue::GetUDData(int UDDhandle, pstr & sudd) -> int`
Retrieves a string previously stored using the specified handle into a provided string reference.
Receiver: instance — call as `obj->GetUDData(...)`
Params: `UDDhandle` — Handle identifying the data slot. `sudd` — Reference where the retrieved string will be stored.
Returns: An error code or success status.
Side effects: Writes data into the provided string reference.

#### `Residue::_AddAtom(PAtom atm) -> int`
Internal method used to add an atom pointer to the residue's atom array.
Receiver: instance — call as `obj->_AddAtom(...)`
Params: `atm` — Pointer to the atom to add.
Returns: The updated atom count.
Side effects: Adds the atom to the internal list.

#### `Residue::_ExcludeAtom(int kndex) -> int`
Internal method used to mark an atom at a specific index as excluded.
Receiver: instance — call as `obj->_ExcludeAtom(...)`
Params: `kndex` — The index of the atom to exclude.
Returns: An integer status code.
Side effects: Marks the atom at `kndex` as excluded.

#### `Residue::Copy(PResidue res) -> void`
Creates a deep copy of another residue object.
Receiver: instance — call as `obj->Copy(...)`
Params: `res` — Pointer to the residue to copy from.
Returns: None.
Side effects: Copies data from `res` into `this` object.

#### `Residue::write(io::RFile f) -> void`
Writes the residue's data content to the specified file stream.
Receiver: instance — call as `obj->write(...)`
Params: `f` — Reference to the output file stream.
Returns: None.
Side effects: Writes data to `f`.

#### `Residue::read(io::RFile f) -> void`
Reads residue data content from the specified file stream.
Receiver: instance — call as `obj->read(...)`
Params: `f` — Reference to the input file stream.
Returns: None.
Side effects: Reads data from `f` and populates the object.

---

## `Root`
Manages the overall structure data, including coordinates, crystal, symmetry, and file I/O capabilities for structural biology assemblies.

#### `Root::SetFlag(word Flag) -> void`
Sets a specific internal flag controlling the root structure's behavior or state.
Receiver: instance — call as `obj->SetFlag(...)`
Params: `Flag` — the flag word to set.
Side effects: Mutates the internal `Flags` word.

#### `Root::ReadPDBASCII(cpstr PDBFileName, io::GZ_MODE gzipMode) -> ERROR_CODE`
Reads atomic coordinates and metadata from a standard PDB ASCII file.
Receiver: instance — call as `obj->ReadPDBASCII(...)`
Params: `PDBFileName` — filename of the PDB file; `gzipMode` — compression mode.
Returns: An error code indicating success or failure during file parsing.
Side effects: Populates the structure's atomic coordinates and metadata.
Errors: Failure reading or parsing the specified PDB file format.

#### `Root::ReadPDBASCII1(cpstr PDBLFName, io::GZ_MODE gzipMode) -> ERROR_CODE`
Reads atomic coordinates and metadata from a PDB file provided as a literal string.
Receiver: instance — call as `obj->ReadPDBASCII1(...)`
Params: `PDBLFName` — filename/string identifier of the PDB file; `gzipMode` — compression mode.
Returns: An error code indicating success or failure during file parsing.
Side effects: Populates the structure's atomic coordinates and metadata.
Errors: Failure reading or parsing the specified PDB file format.

#### `Root::ReadPDBASCII(io::RFile f) -> ERROR_CODE`
Reads atomic coordinates and metadata from a file stream.
Receiver: instance — call as `obj->ReadPDBASCII(...)`
Params: `f` — the input file stream.
Returns: An error code indicating success or failure during file parsing.
Side effects: Populates the structure's atomic coordinates and metadata.
Errors: Failure reading or parsing the file stream.

#### `Root::ReadCoorFile(cpstr LFName, io::GZ_MODE gzipMode) -> ERROR_CODE`
Reads coordinate information from a file specified by a name.
Receiver: instance — call as `obj->ReadCoorFile(...)`
Params: `LFName` — filename of the coordinate file; `gzipMode` — compression mode.
Returns: An error code indicating success or failure during file reading.
Side effects: Updates the structure's coordinates based on the input file.
Errors: Failure reading or parsing the coordinate file.

#### `Root::ReadCoorFile1(cpstr CFName, io::GZ_MODE gzipMode) -> ERROR_CODE`
Reads coordinate information from a file provided as a literal string.
Receiver: instance — call as `obj->ReadCoorFile1(...)`
Params: `CFName` — filename/string identifier of the coordinate file; `gzipMode` — compression mode.
Returns: An error code indicating success or failure during file reading.
Side effects: Updates the structure's coordinates based on the input file.
Errors: Failure reading or parsing the coordinate file.

#### `Root::ReadCoorFile(io::RFile f) -> ERROR_CODE`
Reads coordinate information from a file stream.
Receiver: instance — call as `obj->ReadCoorFile(...)`
Params: `f` — the input file stream.
Returns: An error code indicating success or failure during file reading.
Side effects: Updates the structure's coordinates based on the input file.
Errors: Failure reading or parsing the file stream.

#### `Root::WritePDBASCII(cpstr PDBFileName, io::GZ_MODE gzipMode) -> ERROR_CODE`
Writes the current structure's coordinates and metadata to a PDB ASCII file.
Receiver: instance — call as `obj->WritePDBASCII(...)`
Params: `PDBFileName` — filename for the output PDB file; `gzipMode` — compression mode.
Returns: An error code indicating success or failure during file writing.
Side effects: Writes coordinates and metadata to the specified file.
Errors: Failure writing to the specified file.

#### `Root::WritePDBASCII1(cpstr PDBLFName, io::GZ_MODE gzipMode) -> ERROR_CODE`
Writes the current structure's coordinates and metadata to a PDB file provided as a literal string.
Receiver: instance — call as `obj->WritePDBASCII1(...)`
Params: `PDBLFName` — filename/string identifier for the output PDB file; `gzipMode` — compression mode.
Returns: An error code indicating success or failure during file writing.
Side effects: Writes coordinates and metadata to the specified string output.
Errors: Failure writing to the specified string output.

#### `Root::WritePDBASCII(io::RFile f) -> void`
Writes the current structure's coordinates and metadata to a file stream.
Receiver: instance — call as `obj->WritePDBASCII(...)`
Params: `f` — the output file stream.
Side effects: Writes coordinates and metadata to the file stream.

#### `Root::WriteCIFASCII(cpstr CIFFileName, io::GZ_MODE gzipMode) -> ERROR_CODE`
Writes the current structure's coordinates and metadata to a CIF ASCII file.
Receiver: instance — call as `obj->WriteCIFASCII(...)`
Params: `CIFFileName` — filename for the output CIF file; `gzipMode` — compression mode.
Returns: An error code indicating success or failure during file writing.
Side effects: Writes coordinates and metadata to the specified file.
Errors: Failure writing to the specified file.

#### `Root::WriteCIFASCII1(cpstr CIFLFName, io::GZ_MODE gzipMode) -> ERROR_CODE`
Writes the current structure's coordinates and metadata to a CIF file provided as a literal string.
Receiver: instance — call as `obj->WriteCIFASCII1(...)`
Params: `CIFLFName` — filename/string identifier for the output CIF file; `gzipMode` — compression mode.
Returns: An error code indicating success or failure during file writing.
Side effects: Writes coordinates and metadata to the specified string output.
Errors: Failure writing to the specified string output.

#### `Root::GetInputBuffer(pstr Line, int & count) -> void`
Retrieves the last processed input line and counts how many lines were read for the current processing step.
Receiver: instance — call as `obj->GetInputBuffer(...)`
Params: `Line` — output string to receive the line content; `count` — reference to receive the line count.
Side effects: Populates `Line` and updates `count`.

#### `Root::PutPDBString(cpstr PDBString) -> ERROR_CODE`
Appends a plain keyworded string, typically a REMARK or JRNL card, to the structure's metadata, ensuring it is the last entry of its type.
Receiver: instance — call as `obj->PutPDBString(...)`
Params: `PDBString` — the string content to append.
Returns: An error code indicating success or failure.
Side effects: Appends the string to the structure's metadata record.
Errors: Failure appending the string.

#### `Root::PDBCleanup(word CleanKey) -> word`
Performs various cleanup operations on the structure's records (e.g., atom naming, residue renumbering) to conform to standard PDB formats.
Receiver: instance — call as `obj->PDBCleanup(...)`
Params: `CleanKey` — a bitmask specifying which cleanup operations to perform.
Returns: A word indicating the result, potentially detailing specific warnings or required subsequent cleanups.
Side effects: Mutates internal structures like atom names, chain IDs, and residue numbering.
Errors: Returns specific bit flags indicating if certain required cleanup steps could not be completed.

#### `Root::PutAtom(int index, int serNum, const AtomName atomName, const ResName resName, const ChainID chainID, int seqNum, const InsCode insCode, const AltLoc altLoc, const SegID segID, const Element element) -> int`
Adds a new atom to the structure at a specified position, potentially replacing an existing atom and updating associated residues.
Receiver: instance — call as `obj->PutAtom(...)`
Params: `index` — the position (1-based) where the atom is placed; `serNum` — the desired serial number; `atomName` — the atom name; `resName` — the residue name; `chainID` — the chain identifier; `seqNum` — the sequence number; `insCode` — insertion code; `altLoc` — alternative location; `segID` — segment identifier; `element` — the element type.
Returns: The resulting index/status of the operation.
Side effects: Increases atom count, populates atom records, and updates residue/chain metadata.
Errors: Fails if placement or required metadata population is impossible.

#### `Root::PutAtom(int index, PAtom A, int serNum) -> int`
Adds a new atom to the structure using a pre-existing `PAtom` object, allowing the user to specify index and serial number.
Receiver: instance — call as `obj->PutAtom(...)`
Params: `index` — the position (1-based) where the atom is placed; `A` — pointer to the atom object containing coordinates and names; `serNum` — the desired serial number (defaults to index if 0).
Returns: The resulting index/status of the operation.
Side effects: Increases atom count, populates atom records, and updates residue/chain metadata.
Errors: Fails if placement or required metadata population is impossible.

#### `Root::RemoveAtom(int index) -> void`
Deletes the atom located at the specified index in the primary atom array, and consequently cleans up any associated empty residues, chains, or models.
Receiver: instance — call as `obj->RemoveAtom(...)`
Params: `index` — the 1-based index of the atom to remove.
Side effects: Decreases atom count and may remove connected residue/chain/model records.
Errors: Attempting to remove an invalid index.

#### `Root::FinishStructEdit() -> int`
Finalizes all structural editing operations, ensuring that all derived numbering and links (e.g., residue sequence numbers) are consistent across the entire structure.
Receiver: instance — call as `obj->FinishStructEdit()`
Params: None.
Returns: An integer status code indicating completion success.
Side effects: Stabilizes the structural data model after modification.

#### `Root::DeleteAllModels() -> int`
Deletes all model records from the structure, resetting the model count and clearing model-specific data.
Receiver: instance — call as `obj->DeleteAllModels()`
Params: None.
Returns: An integer status code indicating the success of the deletion process.
Side effects: Resets model information in the structure.
Errors: None documented.

#### `Root::CrystReady() -> int`
Checks the completeness of crystallographic information required for structural refinement or analysis.
Receiver: instance — call as `obj->CrystReady()`
Params: None.
Returns: An integer flag combining status bits (e.g., `CRRDY_Complete`) indicating the completeness level of the crystal data.
Side effects: Reads and validates internal crystal structure parameters.
Errors: Returns specific flags indicating missing or ambiguous data (e.g., `CRRDY_NoCell`).

#### `Root::GetNumberOfNCSMatrices() -> int`
Returns the total count of Non-Crystallographic Symmetry (NCS) matrices available for the structure.
Receiver: instance — call as `obj->GetNumberOfNCSMatrices()`
Params: None.
Returns: The total integer number of available NCS matrices.
Side effects: None.

#### `Root::GetSymOp(int Nop) -> pstr`
Retrieves the symbolic name (e.g., XYZ card) for the symmetry operation defined by the given index.
Receiver: instance — call as `obj->GetSymOp(...)`
Params: `Nop` — the zero-based index of the symmetry operation.
Returns: A string containing the name of the symmetry operation.
Side effects: None.

#### `Root::RegisterUDInteger(UDR_TYPE udr_type, cpstr UDDataID) -> int`
Registers a unique identifier for a user-defined integer data type within the structure context.
Receiver: instance — call as `obj->RegisterUDInteger(...)`
Params: `udr_type` — the type of data being registered; `UDDataID` — the unique string identifier for this data.
Returns: An integer status code indicating success or failure.
Side effects: Updates the internal user-defined data registry.
Errors: Fails if the `UDDataID` is already registered.

#### `Root::RegisterUDReal(UDR_TYPE udr_type, cpstr UDDataID) -> int`
Registers a unique identifier for a user-defined real/floating-point data type within the structure context.
Receiver: instance — call as `obj->RegisterUDReal(...)`
Params: `udr_type` — the type of data being registered; `UDDataID` — the unique string identifier for this data.
Returns: An integer status code indicating success or failure.
Side effects: Updates the internal user-defined data registry.
Errors: Fails if the `UDDataID` is already registered.

#### `Root::RegisterUDString(UDR_TYPE udr_type, cpstr UDDataID) -> int`
Registers a unique identifier for a user-defined string data type within the structure context.
Receiver: instance — call as `obj->RegisterUDString(...)`
Params: `udr_type` — the type of data being registered; `UDDataID` — the unique string identifier for this data.
Returns: An integer status code indicating success or failure.
Side effects: Updates the internal user-defined data registry.
Errors: Fails if the `UDDataID` is already registered.

#### `Root::GetUDDHandle(UDR_TYPE udr_type, cpstr UDDataID) -> int`
Retrieves the internal handle ID corresponding to a registered user-defined data type.
Receiver: instance — call as `obj->GetUDDHandle(...)`
Params: `udr_type` — the type of data; `UDDataID` — the unique string identifier.
Returns: The integer handle ID associated with the data type.
Side effects: None.
Errors: Returns an error code if the data type is not registered.

#### `Root::SetSpaceGroup(cpstr spGroup) -> int`
Sets the space group identifier for the structure, defining the crystal symmetry group.
Receiver: instance — call as `obj->SetSpaceGroup(...)`
Params: `spGroup` — the string identifier of the space group (e.g., "P1").
Returns: An integer status code indicating success or failure.
Side effects: Updates the space group definition for the structure.
Errors: Fails if the space group string is invalid or unrecognized.

#### `Root::GetSpaceGroup() -> pstr`
Retrieves the string identifier of the currently set space group.
Receiver: instance — call as `obj->GetSpaceGroup()`
Params: None.
Returns: A string containing the space group identifier.
Side effects: None.

#### `Root::GetSpaceGroupFix() -> pstr`
Retrieves a fixed or standardized string representation of the space group, potentially correcting minor entries.
Receiver: instance — call as `obj->GetSpaceGroupFix()`
Params: None.
Returns: A string containing the fixed space group identifier.
Side effects: None.

#### `Root::SetCell(realtype cell_a, realtype cell_b, realtype cell_c, realtype cell_alpha, realtype cell_beta, realtype cell_gamma, int OrthCode) -> void`
Manually sets or updates the cell parameters (lattice constants and angles) for the structure.
Receiver: instance — call as `obj->SetCell(...)`
Params: `cell_a` — length of unit cell axis a; `cell_b` — length of unit cell axis b; `cell_c` — length of unit cell axis c; `cell_alpha` — angle $\alpha$; `cell_beta` — angle $\beta$; `cell_gamma` — angle $\gamma$; `OrthCode` — orthogonalization code.
Side effects: Updates the stored cell parameters.

#### `Root::GetCell(realtype & cell_a, realtype & cell_b, realtype & cell_c, realtype & cell_alpha, realtype & cell_beta, realtype & cell_gamma, realtype & vol, int & OrthCode) -> int`
Retrieves the current cell parameters, volume, and orthogonalization code from the structure.
Receiver: instance — call as `obj->GetCell(...)`
Params: `cell_a`, `cell_b`, `cell_c`, `cell_alpha`, `cell_beta`, `cell_gamma`, `vol`, `OrthCode` — reference parameters to populate with cell data.
Returns: An integer status code indicating success or failure.
Side effects: None.
Errors: Returns an error code if cell parameters are not set.

#### `Root::GetTMatrix(mat44 & TMatrix, int Nop, int cellshift_a, int cellshift_b, int cellshift_c) -> int`
Calculates and returns the coordinate transformation matrix required to map coordinates according to a symmetry operation.
Receiver: instance — call as `obj->GetTMatrix(...)`
Params: `TMatrix` — reference matrix to receive the resulting transformation matrix; `Nop` — the symmetry operation index; `cellshift_a`, `cellshift_b`, `cellshift_c` — cell shift parameters for the transformation.
Returns: An integer status code (0 for success, non-zero for failure/warning).
Side effects: Writes transformation data into `TMatrix`.
Errors: Returns non-zero codes if the symmetry operation is undefined, or if crystal/cell parameters are missing.

#### `Root::GetSymOpMatrix(mat44 & TMatrix, int Nop) -> int`
Returns the full transformation matrix that describes the Nop-th symmetry operator in the space group.
Receiver: instance — call as `obj->GetSymOpMatrix(...)`
Params: `TMatrix` — reference matrix to receive the resulting transformation matrix; `Nop` — the symmetry operation index.
Returns: An integer status code (0 for success, non-zero for failure/warning).
Side effects: Writes transformation data into `TMatrix`.
Errors: Returns non-zero codes if the symmetry operation is undefined, or if crystal/cell parameters are missing.

#### `Root::GetFileType() -> int`
Returns an enumeration value indicating the type of the file currently being read or processed.
Receiver: instance — call as `obj->GetFileType()`
Params: None.
Returns: An integer representing the file type (e.g., PDB, CIF, BIN).
Side effects: None.

#### `Root::PutUDData(int UDDhandle, int iudd) -> int`
Writes an integer value into the structure metadata, associated with a specific user-defined data handle.
Receiver: instance — call as `obj->PutUDData(...)`
Params: `UDDhandle` — the handle ID of the data; `iudd` — the integer value to store.
Returns: An integer status code indicating success or failure.
Side effects: Stores the integer data associated with the handle.
Errors: Fails if the handle ID is invalid.

#### `Root::PutUDData(int UDDhandle, realtype rudd) -> int`
Writes a real (floating-point) value into the structure metadata, associated with a specific user-defined data handle.
Receiver: instance — call as `obj->PutUDData(...)`
Params: `UDDhandle` — the handle ID of the data; `rudd` — the real value to store.
Returns: An integer status code indicating success or failure.
Side effects: Stores the real data associated with the handle.
Errors: Fails if the handle ID is invalid.

#### `Root::PutUDData(int UDDhandle, cpstr sudd) -> int`
Writes a string value into the structure metadata, associated with a specific user-defined data handle.
Receiver: instance — call as `obj->PutUDData(...)`
Params: `UDDhandle` — the handle ID of the data; `sudd` — the string value to store.
Returns: An integer status code indicating success or failure.
Side effects: Stores the string data associated with the handle.
Errors: Fails if the handle ID is invalid.

#### `Root::GetUDData(int UDDhandle, int & iudd) -> int`
Retrieves an integer value from the structure metadata associated with a given user-defined data handle.
Receiver: instance — call as `obj->GetUDData(...)`
Params: `UDDhandle` — the handle ID of the data; `iudd` — reference parameter to receive the integer value.
Returns: An integer status code indicating success or failure.
Side effects: None.
Errors: Returns an error code if the handle ID is invalid or the data is missing.

#### `Root::GetUDData(int UDDhandle, realtype & rudd) -> int`
Retrieves a real (floating-point) value from the structure metadata associated with a given user-defined data handle.
Receiver: instance — call as `obj->GetUDData(...)`
Params: `UDDhandle` — the handle ID of the data; `rudd` — reference parameter to receive the real value.
Returns: An integer status code indicating success or failure.
Side effects: None.
Errors: Returns an error code if the handle ID is invalid or the data is missing.

#### `Root::GetUDData(int UDDhandle, pstr sudd, int maxLen) -> int`
Retrieves a string value from the structure metadata, allocating space in a provided buffer.
Receiver: instance — call as `obj->GetUDData(...)`
Params: `UDDhandle` — the handle ID of the data; `sudd` — output string buffer; `maxLen` — maximum length of the string to retrieve.
Returns: An integer status code indicating success or failure.
Side effects: Writes the string data into `sudd`.
Errors: Fails if the handle ID is invalid or the data is missing.

#### `Root::GetUDData(int UDDhandle, pstr & sudd) -> int`
Retrieves a string value from the structure metadata by assigning it to a reference parameter.
Receiver: instance — call as `obj->GetUDData(...)`
Params: `UDDhandle` — the handle ID of the data; `sudd` — reference parameter to receive the string.
Returns: An integer status code indicating success or failure.
Side effects: Populates the string via the reference parameter.
Errors: Fails if the handle ID is invalid or the data is missing.

#### `Root::GetStructureTitle(pstr & L) -> pstr`
Retrieves the primary title string, attempting to resolve it from multiple record types in hierarchical order.
Receiver: instance — call as `obj->GetStructureTitle(...)`
Params: `L` — output string parameter receiving the final resolved title.
Returns: A string containing the structure title, or a default value if multiple records are missing.
Side effects: Reads metadata and resolves potential title sources.
Errors: If no title information is found, it returns a predefined fallback string.

#### `Root::GetCell(realtype & cell_a, realtype & cell_b, realtype & cell_c, realtype & cell_alpha, realtype & cell_beta, realtype & cell_gamma, realtype & vol, int & OrthCode) -> int`
Retrieves the current cell parameters, volume, and orthogonalization code from the structure.
Receiver: instance — call as `obj->GetCell(...)`
Params: `cell_a`, `cell_b`, `cell_c`, `cell_alpha`, `cell_beta`, `cell_gamma`, `vol`, `OrthCode` — reference parameters to receive the cell data.
Returns: An integer status code indicating success or failure.
Side effects: None.
Errors: Returns an error code if cell parameters are not set.

---

## `SelManager`
Manages the creation, manipulation, and retrieval of atom/object selections masks within the structural structure.

#### `SelManager::NewSelection() -> int`
Creates a new selection mask handle, which is required for all subsequent selection operations.
Receiver: instance — call as `obj->NewSelection()` (or `obj.NewSelection()`).
Returns: The new selection handle, which is a positive integer.
Side effects: A new selection mask is created; subsequent selection operations must use this returned handle.

#### `SelManager::DeleteSelection(int selHnd) -> void`
Deletes the specified selection mask handle and removes the corresponding selection attributes from all atoms.
Receiver: instance — call as `obj->DeleteSelection(selHnd)` (or `obj.DeleteSelection(selHnd)`).
Side effects: The specified mask handle (`selHnd`) becomes invalid.

#### `SelManager::SelectAtoms(int selHnd, int iSer1, int iSer2, SELECTION_KEY selKey) -> void`
Selects atoms within a specified serial number range by adding them to the atoms marked by the given mask.
Receiver: instance — call as `obj->SelectAtoms(selHnd, iSer1, iSer2, selKey)` (or `obj.SelectAtoms(selHnd, iSer1, iSer2, selKey)`).
Params: `selHnd` — The handle of the mask to modify. `iSer1` — The starting serial number of the range. `iSer2` — The ending serial number of the range. `selKey` — Determines how selection results are merged (default is OR).
Side effects: Atoms in the specified range are marked by the selection mask.

#### `SelManager::SelectAtoms(int selHnd, ivector asn, int nsn, SELECTION_KEY selKey) -> void`
Selects atoms whose serial numbers are listed in the provided vector.
Receiver: instance — call as `obj->SelectAtoms(selHnd, asn, nsn, selKey)` (or `obj.SelectAtoms(selHnd, asn, nsn, selKey)`).
Params: `selHnd` — The handle of the mask to modify. `asn` — Vector containing the serial numbers of atoms to select. `nsn` — The number of atoms in `asn`. `selKey` — Determines how selection results are merged (default is OR).
Side effects: Atoms listed in `asn` are marked by the selection mask.

#### `SelManager::UnselectAtoms(int selHnd, int iSer1, int iSer2) -> void`
Clears the specified selection mask for atoms within a given serial number range.
Receiver: instance — call as `obj->UnselectAtoms(selHnd, iSer1, iSer2)` (or `obj.UnselectAtoms(selHnd, iSer1, iSer2)`).
Side effects: Selection attributes corresponding to the mask are removed from atoms in the range.

#### `SelManager::SelectAtom(int selHnd, PAtom A, SELECTION_KEY selKey, bool makeIndex) -> void`
Selects a single atom identified by pointer A based on a selection key, optionally updating the internal selection index.
Receiver: instance — call as `obj->SelectAtom(selHnd, A, selKey, makeIndex)` (or `obj.SelectAtom(selHnd, A, selKey, makeIndex)`).
Params: `selHnd` — The handle of the mask to modify. `A` — Pointer to the specific atom to select. `selKey` — Determines how selection results are merged (default is OR). `makeIndex` — If false, selection index is not updated, improving performance.
Side effects: Atom A is marked by the selection mask.

#### `SelManager::SelectResidue(int selHnd, PResidue Res, SELECTION_TYPE sType, SELECTION_KEY sKey, bool makeIndex) -> void`
Selects a single residue (or its descendants) based on its type and properties.
Receiver: instance — call as `obj->SelectResidue(selHnd, Res, sType, sKey, makeIndex)` (or `obj.SelectResidue(selHnd, Res, sType, sKey, makeIndex)`).
Params: `selHnd` — The handle of the mask to modify. `Res` — Pointer to the specific residue to select. `sType` — The selection type defining what is being selected (e.g., residue, chain). `sKey` — Determines how selection results are merged (default is OR). `makeIndex` — If false, selection index is not updated.
Side effects: Residue Res and its associated atoms are marked by the selection mask.

#### `SelManager::SelectChain(int selHnd, PChain chain, SELECTION_TYPE sType, SELECTION_KEY sKey, bool makeIndex) -> void`
Selects a single chain (or its descendants) based on its type and properties.
Receiver: instance — call as `obj->SelectChain(selHnd, chain, sType, sKey, makeIndex)` (or `obj.SelectChain(selHnd, chain, sType, sKey, makeIndex)`).
Params: `selHnd` — The handle of the mask to modify. `chain` — Pointer to the specific chain to select. `sType` — The selection type defining what is being selected. `sKey` — Determines how selection results are merged (default is OR). `makeIndex` — If false, selection index is not updated.
Side effects: Chain objects associated with `chain` are marked by the selection mask.

#### `SelManager::SelectModel(int selHnd, PModel mdl, SELECTION_TYPE sType, SELECTION_KEY sKey, bool makeIndex) -> void`
Selects a single model (or its descendants) based on its type and properties.
Receiver: instance — call as `obj->SelectModel(selHnd, mdl, sType, sKey, makeIndex)` (or `obj.SelectModel(selHnd, mdl, sType, sKey, makeIndex)`).
Params: `selHnd` — The handle of the mask to modify. `mdl` — Pointer to the specific model to select. `sType` — The selection type defining what is being selected. `sKey` — Determines how selection results are merged (default is OR). `makeIndex` — If false, selection index is not updated.
Side effects: Model objects associated with `mdl` are marked by the selection mask.

#### `SelManager::MakeSelIndex(int selHnd) -> int`
Manually calculates the internal selection index count for the selection mask specified by `selHnd`.
Receiver: instance — call as `obj->MakeSelIndex(selHnd)` (or `obj.MakeSelIndex(selHnd)`).
Params: `selHnd` — The handle of the selection mask whose index needs calculation.
Returns: The final number of selected objects recorded by the mask.
Side effects: The internal selection index count is updated.

#### `SelManager::SelectAtoms(int selHnd, int iModel, cpstr Chains, int ResNo1, cpstr Ins1, int ResNo2, cpstr Ins2, cpstr RNames, cpstr ANames, cpstr Elements, cpstr altLocs, cpstr Segments, cpstr Charges, realtype occ1, realtype occ2, realtype x0, realtype y0, realtype z0, realtype d0, SELECTION_KEY selKey) -> void`
Selects atoms based on comprehensive criteria including model, chains, residue ranges, atom names, elemental composition, and geometric criteria (sphere/plane).
Receiver: instance — call as `obj->SelectAtoms(selHnd, iModel, Chains, ResNo1, Ins1, ResNo2, Ins2, RNames, ANames, Elements, altLocs, Segments, Charges, occ1, occ2, x0, y0, z0, d0, selKey)` (or `obj.SelectAtoms(selHnd, iModel, Chains, ResNo1, Ins1, ResNo2, Ins2, RNames, ANames, Elements, altLocs, Segments, Charges, occ1, occ2, x0, y0, z0, d0, selKey)`).
Params: `selHnd` — The handle of the mask to modify. `iModel` — Model number (0 for any). `Chains` — Comma-separated list of chains, or "*" for any. `ResNo1` — Starting residue sequence number. `Ins1` — Starting residue insertion code. `ResNo2` — Ending residue sequence number. `Ins2` — Ending residue insertion code. `RNames` — Comma-separated list of residue names, or "*" for any. `ANames` — Comma-separated list of atom names, or "*" for any. `Elements` — Comma-separated list of elements, or "*" for any. `altLocs` — Comma-separated list of alternative locations, or "*" for any. `Segments` — Comma-separated list of segment IDs, or "*" for any. `Charges` — Comma-separated list of charges, or "*" for any. `occ1` — Lowest occupancy for selection filtering. `occ2` — Highest occupancy for selection filtering. `x0` — Reference X-coordinate for geometric filtering. `y0` — Reference Y-coordinate for geometric filtering. `z0` — Reference Z-coordinate for geometric filtering. `d0` — Selection distance from the reference point (0.0 for any distance). `selKey` — Determines how selection results are merged (default is OR).
Side effects: Atoms matching all criteria are marked by the selection mask.

#### `SelManager::SelectAtoms(int selHnd, int iModel, cpstr Chains, int ResNo1, cpstr Ins1, int ResNo2, cpstr Ins2, cpstr RNames, cpstr ANames, cpstr Elements, cpstr altLocs, SELECTION_KEY selKey) -> void`
Selects atoms based on standard structural criteria excluding occupancy and geometric filtering.
Receiver: instance — call as `obj->SelectAtoms(selHnd, iModel, Chains, ResNo1, Ins1, ResNo2, Ins2, RNames, ANames, Elements, altLocs, selKey)` (or `obj.SelectAtoms(selHnd, iModel, Chains, ResNo1, Ins1, ResNo2, Ins2, RNames, ANames, Elements, altLocs, selKey)`).
Params: `selHnd` — The handle of the mask to modify. `iModel` — Model number (0 for any). `Chains` — Comma-separated list of chains, or "*" for any. `ResNo1` — Starting residue sequence number. `Ins1` — Starting residue insertion code. `ResNo2` — Ending residue sequence number. `Ins2` — Ending residue insertion code. `RNames` — Comma-separated list of residue names, or "*" for any. `ANames` — Comma-separated list of atom names, or "*" for any. `Elements` — Comma-separated list of elements, or "*" for any. `altLocs` — Comma-separated list of alternative locations, or "*" for any. `selKey` — Determines how selection results are merged (default is OR).
Side effects: Atoms matching all criteria are marked by the selection mask.

#### `SelManager::SelectUDD(int selHnd, SELECTION_TYPE sType, int UDDhandle, int selMin, int selMax, SELECTION_KEY sKey) -> void`
Selects atoms whose User-Defined Data (UDD) falls within a specified integer range.
Receiver: instance — call as `obj->SelectUDD(selHnd, sType, UDDhandle, selMin, selMax, sKey)` (or `obj.SelectUDD(selHnd, sType, UDDhandle, selMin, selMax, sKey)`).
Params: `selHnd` — The handle of the mask to modify. `sType` — The type of selection data being queried. `UDDhandle` — The handle of the UDD to check. `selMin` — The minimum inclusive boundary for the UDD. `selMax` — The maximum inclusive boundary for the UDD. `sKey` — Determines how selection results are merged (default is OR).
Side effects: Atoms meeting the UDD criteria are marked by the selection mask.

#### `SelManager::SelectUDD(int selHnd, SELECTION_TYPE sType, int UDDhandle, realtype selMin, realtype selMax, SELECTION_KEY sKey) -> void`
Selects atoms whose User-Defined Data (UDD) falls within a specified floating-point range.
Receiver: instance — call as `obj->SelectUDD(selHnd, sType, UDDhandle, selMin, selMax, sKey)` (or `obj.SelectUDD(selHnd, sType, UDDhandle, selMin, selMax, sKey)`).
Params: `selHnd` — The handle of the mask to modify. `sType` — The type of selection data being queried. `UDDhandle` — The handle of the UDD to check. `selMin` — The minimum inclusive boundary for the UDD. `selMax` — The maximum inclusive boundary for the UDD. `sKey` — Determines how selection results are merged (default is OR).
Side effects: Atoms meeting the UDD criteria are marked by the selection mask.

#### `SelManager::SelectUDD(int selHnd, SELECTION_TYPE sType, int UDDhandle, cpstr selStr, int cmpRule, SELECTION_KEY sKey) -> void`
Selects atoms based on matching a specified string pattern against a User-Defined Data handle using a comparison rule.
Receiver: instance — call as `obj->SelectUDD(selHnd, sType, UDDhandle, selStr, cmpRule, sKey)` (or `obj.SelectUDD(selHnd, sType, UDDhandle, selStr, cmpRule, sKey)`).
Params: `selHnd` — The handle of the mask to modify. `sType` — The type of selection data being queried. `UDDhandle` — The handle of the UDD to check. `selStr` — The string pattern to match. `cmpRule` — The rule defining how the comparison should occur. `sKey` — Determines how selection results are merged (default is OR).
Side effects: Atoms matching the UDD string pattern are marked by the selection mask.

#### `SelManager::SelectSphere(int selHnd, SELECTION_TYPE sType, realtype x, realtype y, realtype z, realtype r, SELECTION_KEY sKey) -> void`
Selects atoms located within a specified spherical region defined by center coordinates and a radius.
Receiver: instance — call as `obj->SelectSphere(selHnd, sType, x, y, z, r, selKey)` (or `obj.SelectSphere(selHnd, sType, x, y, z, r, selKey)`).
Params: `selHnd` — The handle of the mask to modify. `sType` — The type of selection data being queried. `x, y, z` — Coordinates of the sphere's center. `r` — The radius of the sphere. `selKey` — Determines how selection results are merged (default is OR).
Side effects: Atoms within the sphere are marked by the selection mask.

#### `SelManager::SelectCylinder(int selHnd, SELECTION_TYPE sType, realtype x1, realtype y1, realtype z1, realtype x2, realtype y2, realtype z2, realtype r, SELECTION_KEY sKey) -> void`
Selects atoms located within a specified cylindrical volume defined by two end points and a radius.
Receiver: instance — call as `obj->SelectCylinder(selHnd, sType, x1, y1, z1, x2, y2, z2, r, selKey)` (or `obj.SelectCylinder(selHnd, sType, x1, y1, z1, x2, y2, z2, r, selKey)`).
Params: `selHnd` — The handle of the mask to modify. `sType` — The type of selection data being queried. `x1, y1, z1` — Coordinates of the cylinder's first end. `x2, y2, z2` — Coordinates of the cylinder's second end. `r` — The radius of the cylinder. `selKey` — Determines how selection results are merged (default is OR).
Side effects: Atoms within the cylinder are marked by the selection mask.

#### `SelManager::SelectSlab(int selHnd, SELECTION_TYPE sType, realtype a, realtype b, realtype c, realtype d, realtype r, SELECTION_KEY sKey) -> void`
Selects atoms located within a specified slab (distance) from a plane defined by the equation ax+by+cz=d.
Receiver: instance — call as `obj->SelectSlab(selHnd, sType, a, b, c, d, r, selKey)` (or `obj.SelectSlab(selHnd, sType, a, b, c, d, r, selKey)`).
Params: `selHnd` — The handle of the mask to modify. `sType` — The type of selection data being queried. `a, b, c` — Coefficients defining the plane's normal vector. `d` — The distance constant for the plane equation. `r` — The maximum distance from the plane for selection. `selKey` — Determines how selection results are merged (default is OR).
Side effects: Atoms within the slab are marked by the selection mask.

#### `SelManager::SelectNeighbours(int selHnd, SELECTION_TYPE sType, PPAtom sA, int alen, realtype d1, realtype d2, SELECTION_KEY sKey) -> void`
Selects atoms located between minimum and maximum distances from a set of previously selected atoms (`sA`).
Receiver: instance — call as `obj->SelectNeighbours(selHnd, sType, sA, alen, d1, d2, sKey)` (or `obj.SelectNeighbours(selHnd, sType, sA, alen, d1, d2, sKey)`).
Params: `selHnd` — The handle of the mask to modify. `sType` — The type of selection data being queried. `sA` — Pointer to the array of atoms already selected (neighbors to check against). `alen` — The number of atoms in the `sA` array. `d1` — The minimum required distance to an atom in `sA`. `d2` — The maximum allowed distance to an atom in `sA`. `sKey` — Determines how selection results are merged (default is OR).
Side effects: Atoms satisfying the distance criteria are marked by the selection mask.

#### `SelManager::GetSelIndex(int selHnd, PPAtom & SelAtom, int & nSelAtoms) -> void`
Retrieves an indexed array of all atoms currently selected by the given mask handle.
Receiver: instance — call as `obj->GetSelIndex(selHnd, SelAtom, nSelAtoms)` (or `obj.GetSelIndex(selHnd, SelAtom, nSelAtoms)`).
Params: `selHnd` — The handle of the mask to check. `SelAtom` — Reference to the array where selected atom pointers will be stored. `nSelAtoms` — Reference to store the count of selected atoms.
Side effects: Fills the output index arrays/pointers; does not dispose of the data.

#### `SelManager::GetSelIndex(int selHnd, PPResidue & SelResidues, int & nSelResidues) -> void`
Retrieves an indexed array of all residues currently selected by the given mask handle.
Receiver: instance — call as `obj->GetSelIndex(selHnd, SelResidues, nSelResidues)` (or `obj.GetSelIndex(selHnd, SelResidues, nSelResidues)`).
Params: `selHnd` — The handle of the mask to check. `SelResidues` — Reference to the array where selected residue pointers will be stored. `nSelResidues` — Reference to store the count of selected residues.
Side effects: Fills the output index arrays/pointers; does not dispose of the data.

#### `SelManager::GetSelIndex(int selHnd, PPChain & SelChains, int & nSelChains) -> void`
Retrieves an indexed array of all chains currently selected by the given mask handle.
Receiver: instance — call as `obj->GetSelIndex(selHnd, SelChains, nSelChains)` (or `obj.GetSelIndex(selHnd, SelChains, nSelChains)`).
Params: `selHnd` — The handle of the mask to check. `SelChains` — Reference to the array where selected chain pointers will be stored. `nSelChains` — Reference to store the count of selected chains.
Side effects: Fills the output index arrays/pointers; does not dispose of the data.

#### `SelManager::GetSelIndex(int selHnd, PPModel & SelModels, int & nSelModels) -> void`
Retrieves an indexed array of all models currently selected by the given mask handle.
Receiver: instance — call as `obj->GetSelIndex(selHnd, SelModels, nSelModels)` (or `obj.GetSelIndex(selHnd, SelModels, nSelModels)`).
Params: `selHnd` — The handle of the mask to check. `SelModels` — Reference to the array where selected model pointers will be stored. `nSelModels` — Reference to store the count of selected models.
Side effects: Fills the output index arrays/pointers; does not dispose of the data.

#### `SelManager::Select(int selHnd, SELECTION_TYPE sType, cpstr CID, SELECTION_KEY sKey) -> void`
Selects atoms using a coordinate ID string (e.g., "/mdl/chn/s1.i1-s2.i2/at[el]:aloc").
Receiver: instance — call as `obj->Select(selHnd, sType, CID, sKey)` (or `obj.Select(selHnd, sType, CID, sKey)`).
Params: `selHnd` — The handle of the mask to modify. `sType` — The type of selection data being queried. `CID` — The coordinate ID string pattern. `sKey` — Determines how selection results are merged (default is OR).
Side effects: Atoms matching the specified coordinate ID pattern are marked by the selection mask.

#### `SelManager::Select(int selHnd1, SELECTION_TYPE sType, int selHnd2, SELECTION_KEY sKey) -> void`
Propagates selections between two selection masks.
Receiver: instance — call as `obj->Select(selHnd1, sType, selHnd2, sKey)` (or `obj.Select(selHnd1, sType, selHnd2, sKey)`).
Params: `selHnd1` — The destination mask handle (where selections will be added). `sType` — The type of selection data being queried. `selHnd2` — The source mask handle containing the selections to propagate. `sKey` — Determines how selection results are merged (default is OR).
Side effects: Objects selected by mask `selHnd2` are added to the selection set defined by `selHnd1`.

#### `SelManager::SelectProperty(int selHnd, SELECTION_PROPERTY propKey, SELECTION_TYPE sType, SELECTION_KEY sKey) -> void`
Selects atoms based on matching a specific molecular property key within the structure.
Receiver: instance — call as `obj->SelectProperty(selHnd, propKey, sType, sKey)` (or `obj.SelectProperty(selHnd, propKey, sType, sKey)`).
Params: `selHnd` — The handle of the mask to modify. `propKey` — The property key to filter by (e.g., charge, occupancy). `sType` — The type of selection data being queried. `sKey` — Determines how selection results are merged (default is OR).
Side effects: Atoms matching the property criteria are marked by the selection mask.

#### `SelManager::SelectDomain(int selHnd, cpstr domainRange, SELECTION_TYPE sType, SELECTION_KEY sKey, int modelNo) -> int`
Selects atoms based on complex domain range definitions, allowing for highly specific subset selection.
Receiver: instance — call as `obj->SelectDomain(selHnd, domainRange, sType, sKey, modelNo)` (or `obj.SelectDomain(selHnd, domainRange, sType, sKey, modelNo)`).
Params: `selHnd` — The handle of the mask to modify. `domainRange` — String defining the required structural domain subset. `sType` — The type of selection data being queried. `sKey` — Determines how selection results are merged (default is OR). `modelNo` — The model number to restrict the selection to (1 by default).
Returns: An integer status code indicating success or failure of the selection process.
Side effects: Atoms matching the domain criteria are marked by the selection mask.

#### `SelManager::SelectAtom(PAtom atm, int maskNo, SELECTION_KEY sKey, int & nsel) -> void`
(Protected implementation detail) Marks a specific atom within the selection list specified by `maskNo`.
Receiver: instance — call as `obj->SelectAtom(atm, maskNo, sKey, nsel)` (or `obj.SelectAtom(atm, maskNo, sKey, nsel)`).
Params: `atm` — Pointer to the atom to select. `maskNo` — The handle of the mask to modify. `sKey` — Determines how selection results are merged (default is OR). `nsel` — Reference to store the number of selected atoms.
Side effects: Atom `atm` is marked by the selection mask.

#### `SelManager::SelectObject(SELECTION_TYPE sType, PAtom atm, int maskNo, SELECTION_KEY sKey, int & nsel) -> void`
(Protected implementation detail) Selects an atom based on type and coordinates, updating the count.
Receiver: instance — call as `obj->SelectObject(sType, atm, maskNo, sKey, nsel)` (or `obj.SelectObject(sType, atm, maskNo, sKey, nsel)`).
Params: `sType` — The type of object to select. `atm` — Pointer to the atom defining the object. `maskNo` — The handle of the mask to modify. `sKey` — Determines how selection results are merged (default is OR). `nsel` — Reference to store the number of selected atoms.
Side effects: Atoms defining the object are marked by the selection mask.

#### `SelManager::SelectObject(PMask object, int maskNo, SELECTION_KEY sKey, int & nsel) -> void`
(Protected implementation detail) Selects a set of objects defined by an existing mask.
Receiver: instance — call as `obj->SelectObject(object, maskNo, sKey, nsel)` (or `obj.SelectObject(object, maskNo, sKey, nsel)`).
Params: `object` — Pointer to the mask containing objects to select. `maskNo` — The handle of the mask to modify. `sKey` — Determines how selection results are merged (default is OR). `nsel` — Reference to store the number of selected atoms.
Side effects: Atoms defined by the object mask are marked by the selection mask.

#### `SelManager::MakeSelIndex(int selHnd, SELECTION_TYPE sType, int nsel) -> void`
(Protected implementation detail) Manually calculates and updates the selection index count for a specific selection type and count.
Receiver: instance — call as `obj->MakeSelIndex(selHnd, sType, nsel)` (or `obj.MakeSelIndex(selHnd, sType, nsel)`).
Params: `selHnd` — The handle of the mask to process. `sType` — The type of selection data associated with the index. `nsel` — The number of selected items.
Side effects: The internal selection index count is updated for the given type.

## `SelManager`
Manages the creation, manipulation, and retrieval of atom/object selections masks within the structural structure.

---

## `Sheet`
Represents a biological sheet structure composed of multiple strands.

#### `Sheet::PDBASCIIDump(io::RFile f) -> void`
Writes the sheet structure data in PDB ASCII format to the specified file stream.
Receiver: instance — call as `obj->PDBASCIIDump(...)`
Params: `f` — The file stream to write the data to.
Side effects: Writes data to `f`.

#### `Sheet::MakeCIF(mmcif::PData CIF) -> void`
Generates the CIF representation of the sheet structure and populates the provided data structure.
Receiver: instance — call as `obj->MakeCIF(...)`
Params: `CIF` — Pointer to the CIF data structure to populate.
Side effects: Modifies the content of `CIF`.

#### `Sheet::ConvertPDBASCII(cpstr S) -> ERROR_CODE`
Converts the coordinates and structure information from a PDB ASCII string format into the internal sheet representation.
Receiver: instance — call as `obj->ConvertPDBASCII(...)`
Params: `S` — The PDB ASCII string containing the sheet data.
Returns: An error code indicating success or failure of the conversion.
Side effects: Populates internal structure based on the input string.

#### `Sheet::GetCIF(mmcif::PData CIF) -> int`
Extracts the sheet structure data into the CIF format and populates the provided data structure.
Receiver: instance — call as `obj->GetCIF(...)`
Params: `CIF` — Pointer to the CIF data structure to populate.
Returns: An integer indicating the number of successfully processed components.
Side effects: Modifies the content of `CIF`.

#### `Sheet::Copy(PSheet sheet) -> void`
Copies the contents and structure from another sheet object to this object.
Receiver: instance — call as `obj->Copy(...)`
Params: `sheet` — Pointer to the source sheet object to copy from.
Side effects: Sets the current object's internal state to match the source sheet.

#### `Sheet::write(io::RFile f) -> void`
Writes the entire sheet structure representation to the specified file stream using the stream's default format.
Receiver: instance — call as `obj->write(...)`
Params: `f` — The file stream to write to.
Side effects: Writes data to `f`.

#### `Sheet::read(io::RFile f) -> void`
Reads a sheet structure representation from the specified file stream and populates the object's state.
Receiver: instance — call as `obj->read(...)`
Params: `f` — The file stream to read from.
Side effects: Updates the object's internal state from the file stream.

---

## `Sheets`
Manages and provides I/O operations for data represented by one or more protein secondary structure sheets.

#### `Sheets::PDBASCIIDump(io::RFile f) -> void`
Writes the current sheet data structure content in PDB ASCII format to the specified file stream.
Receiver: instance — call as `obj->PDBASCIIDump(...)`
Params: `f` — The file stream where output should be written.
Side effects: Writes to the file stream `f`.

#### `Sheets::MakeCIF(mmcif::PData CIF) -> void`
Populates the given CIF data structure with the sheet information stored in the object.
Receiver: instance — call as `obj->MakeCIF(...)`
Params: `CIF` — Pointer to the CIF data structure to populate.
Side effects: Modifies the contents of `CIF`.

#### `Sheets::ConvertPDBASCII(cpstr S) -> ERROR_CODE`
Converts PDB ASCII format data from the string `S` into the internal sheet representation.
Receiver: instance — call as `obj->ConvertPDBASCII(...)`
Params: `S` — The string containing the PDB ASCII data.
Returns: An error code indicating success or failure of the conversion.
Side effects: Updates the internal state of the object upon success.

#### `Sheets::GetCIF(mmcif::PData CIF) -> int`
Extracts the sheet information and populates the provided CIF data structure, returning the number of sheets written.
Receiver: instance — call as `obj->GetCIF(...)`
Params: `CIF` — Pointer to the CIF data structure to populate.
Returns: The number of sheets successfully written to the CIF data structure.
Side effects: Modifies the contents of `CIF`.

#### `Sheets::Copy(PSheets Sheets) -> void`
Copies the sheet data and structure from another `Sheets` object pointed to by `Sheets`.
Receiver: instance — call as `obj->Copy(...)`
Params: `Sheets` — Pointer to the source `Sheets` object to copy from.
Side effects: Replaces the object's internal sheet data with the source data.

#### `Sheets::write(io::RFile f) -> void`
Writes the entire content of the object, including all associated sheets, to the specified file stream.
Receiver: instance — call as `obj->write(...)`
Params: `f` — The file stream to write the data to.
Side effects: Writes the structure data to the file stream `f`.

#### `Sheets::read(io::RFile f) -> void`
Reads sheet data from the specified file stream and populates the object's internal state.
Receiver: instance — call as `obj->read(...)`
Params: `f` — The file stream to read the data from.
Side effects: Reads and updates the object's internal sheet data.

#### `Sheets::InitSheets() -> void`
Initializes the internal data structures required to hold and manage sheet geometry information.
Receiver: instance — call as `obj->InitSheets()`
Side effects: Sets up the internal state of the object.

#### `Sheets::CIFFindSheets(mmcif::PData CIF, cpstr Category) -> void`
Performs a search within the provided CIF structure to find and populate sheet data corresponding to the specified category.
Receiver: instance — call as `obj->CIFFindSheets(...)`
Params: `CIF` — Pointer to the CIF data structure to search.
`Category` — The category string used to filter sheets during the search.
Side effects: Reads and populates the object's sheet data based on the search results in `CIF`.

---

## `Strand`
Represents a continuous segment of residues (a chain segment) within a macromolecular structure.

#### `Strand::PDBASCIIDump(pstr S) -> void`
Writes the structural information contained in the strand object to the provided PDB ASCII string buffer.
Receiver: instance — call as `obj->PDBASCIIDump(...)`
Params: `S` — The output string buffer to write to.
Side effects: Writes to the provided string `S`.

#### `Strand::MakeCIF(mmcif::PData CIF) -> void`
Writes the structural information contained in the strand object into a mutable CIF data structure.
Receiver: instance — call as `obj->MakeCIF(...)`
Params: `CIF` — The CIF data structure to populate.
Side effects: Populates the provided CIF data structure.

#### `Strand::ConvertPDBASCII(cpstr S) -> ERROR_CODE`
Parses and converts structural information from a PDB ASCII formatted string into the strand object's internal representation.
Receiver: instance — call as `obj->ConvertPDBASCII(...)`
Params: `S` — The input PDB ASCII string.
Returns: An error code indicating success or failure of the conversion.
Side effects: Updates the internal state of the strand object.

#### `Strand::GetCIF(mmcif::PData CIF, cpstr sheet_id) -> int`
Extracts the strand's structural coordinates and data into a specified CIF data structure, optionally filtering by a sheet ID.
Receiver: instance — call as `obj->GetCIF(...)`
Params: `CIF` — The CIF data structure to populate.
`sheet_id` — The optional sheet ID to restrict the data extraction to.
Returns: The number of atoms/residues successfully written to the CIF structure.
Side effects: Populates the provided CIF data structure.

#### `Strand::Copy(PStrand Strand) -> void`
Creates a deep copy of the data from another strand object into the current strand object.
Receiver: instance — call as `obj->Copy(...)`
Params: `Strand` — The source strand object to copy from.
Side effects: Updates the internal state of the strand object to match the source strand.

#### `Strand::write(io::RFile f) -> void`
Writes the structural information of the strand to a specified output file stream.
Receiver: instance — call as `obj->write(...)`
Params: `f` — The output file stream.
Side effects: Performs I/O write operations to the file stream `f`.

#### `Strand::read(io::RFile f) -> void`
Reads structural information for a strand object from a specified input file stream.
Receiver: instance — call as `obj->read(...)`
Params: `f` — The input file stream.
Side effects: Populates the internal state of the strand object from the file stream `f`.

---

## `Struct`
Represents a structural entity within the MMDB2 library, managing fields associated with a given molecular structure definition.

#### `Struct::GetField(int tagNo) -> pstr`
Returns the string value of the field located at the specified tag position number.
Receiver: instance — call as `obj->GetField(...)`
Params: `tagNo` — the zero-based index representing the tag's position in the structure.
Returns: A pointer to the string value of the field, or `NULL` if the tag does not exist, `CIF_NODATA_DOT_FIELD` if the field is "data not given", or `CIF_NODATA_QUESTION_FIELD` if the field is "data not available".
Side effects: None.
Errors: None.

#### `Struct::GetString(pstr & S, cpstr TName, bool Remove=false) -> int`
Fetches the string value associated with the given tag name into the provided buffer.
Receiver: instance — call as `obj->GetString(...)`
Params: `S` — a pointer to a string where the retrieved value will be stored. `TName` — the tag name to look up. `Remove` — if true, removes the tag from the structure after reading.
Returns: An integer return code (`CIFRC_NoTag`, `CIFRC_NoField`, or `CIFRC_Ok`).
Side effects: The content of `S` is overwritten, and if `Remove` is true, the tag is deleted.
Errors: Returns an error code if the tag or value cannot be found.

#### `Struct::GetString(cpstr TName, int & RC) -> pstr`
Retrieves a pointer to the string value associated with the given tag name.
Receiver: instance — call as `obj->GetString(...)`
Params: `TName` — the tag name to look up. `RC` — a reference integer where the return code will be stored.
Returns: A pointer to the value associated with `TName`, or `NULL` if the tag or value is not found or is a data-missing placeholder.
Side effects: None.
Errors: Returns `NULL` if the tag or value is missing or flagged as missing data.

---

## `SymOps`
Manages and applies symmetry operations to molecular coordinates.

#### `SymOps::SymOps() -> void`
Constructs a default, empty symmetry operations object.
Receiver: instance — call as `obj->SymOps()` (or `obj.SymOps()`).
Side effects: Initializes internal state.

#### `SymOps::SymOps(io::RPStream Object) -> void`
Constructs a symmetry operations object by reading from a specified stream.
Receiver: instance — call as `obj->SymOps(Object)` (or `obj.SymOps(Object)`).
Side effects: Populates symmetry operations from the input stream.

#### `SymOps::~SymOps() -> void`
Destroys the symmetry operations object and releases memory.
Receiver: instance — call as `obj->~SymOps()` (or `obj.~SymOps()`).

#### `SymOps::FreeMemory() -> void`
Deallocates all dynamically allocated memory associated with the object's symmetry operations.
Receiver: instance — call as `obj->FreeMemory()` (or `obj.FreeMemory()`).

#### `SymOps::SetGroupSymopLib(cpstr SpaceGroup, cpstr symop_lib=NULL) -> int`
Sets the symmetry operations from a specified library file associated with a space group.
Receiver: instance — call as `obj->SetGroupSymopLib(...)` (or `obj.SetGroupSymopLib(...)`).
Params: `SpaceGroup` — The name of the space group.
Params: `symop_lib` — Optional path to the symmetry operations library file.
Returns: An integer code indicating success (`SYMOP_Ok`) or failure.
Side effects: Updates the set of available symmetry operations.

#### `SymOps::SetGroup(cpstr SpaceGroup, cpstr syminfo_lib=NULL) -> int`
Sets the symmetry operations using information derived from a specified structure file.
Receiver: instance — call as `obj->SetGroup(...)` (or `obj.SetGroup(...)`).
Params: `SpaceGroup` — The name of the space group.
Params: `syminfo_lib` — Optional path to the structure information library file.
Returns: An integer code indicating success (`SYMOP_Ok`) or failure.
Side effects: Updates the set of available symmetry operations.

#### `SymOps::Reset() -> void`
Removes all currently loaded symmetry operations, resetting the object state.
Receiver: instance — call as `obj->Reset()` (or `obj.Reset()`).
Side effects: Clears internal symmetry operation data.

#### `SymOps::AddSymOp(cpstr XYZOperation) -> int`
Adds a new symmetry operation, potentially specified by its XYZ coordinates, to the object.
Receiver: instance — call as `obj->AddSymOp(...)` (or `obj.AddSymOp(...)`).
Params: `XYZOperation` — The definition of the symmetry operation.
Returns: A return code indicating success or failure.
Side effects: Increases the number of available symmetry operations.

#### `SymOps::PutGroupName(cpstr SpGroupName) -> void`
Sets or updates the name associated with the symmetry operations set.
Receiver: instance — call as `obj->PutGroupName(...)` (or `obj.PutGroupName(...)`).
Params: `SpGroupName` — The name to assign to the symmetry operation group.
Side effects: Updates the object's group name.

#### `SymOps::GetNofSymOps() -> int`
Retrieves the total number of symmetry operations currently stored within the object.
Receiver: instance — call as `obj->GetNofSymOps()` (or `obj.GetNofSymOps()`).
Returns: The count of symmetry operations.

#### `SymOps::GetSymOp(int Nop) -> pstr`
Retrieves a pointer to the symmetry operation located at the given index.
Receiver: instance — call as `obj->GetSymOp(...)` (or `obj.GetSymOp(...)`).
Params: `Nop` — The zero-based index of the desired symmetry operation.
Returns: A pointer to the requested symmetry operation structure.
Errors: Returns null if the provided index `Nop` is out of bounds (not from 0 to Nops-1).

#### `SymOps::Transform(realtype & x, realtype & y, realtype & z, int Nop) -> int`
Applies the symmetry operation specified by `Nop` to transform the input coordinates $(x, y, z)$.
Receiver: instance — call as `obj->Transform(...)` (or `obj.Transform(...)`).
Params: `x`, `y`, `z` — References to the coordinates to be transformed.
Params: `Nop` — The zero-based index of the symmetry operation to use.
Returns: Non-zero if `Nop` is invalid; zero on success.
Side effects: Mutates the values pointed to by $x, y,$ and $z$.
Errors: Returns non-zero if $Nop$ is outside the valid range $[0, Nops-1]$.

#### `SymOps::GetTMatrix(mat44 & TMatrix, int Nop) -> int`
Retrieves the $4\times4$ transformation matrix corresponding to the symmetry operation specified by `Nop`.
Receiver: instance — call as `obj->GetTMatrix(...)` (or `obj.GetTMatrix(...)`).
Params: `TMatrix` — Reference to the matrix variable that will receive the transformation matrix.
Params: `Nop` — The zero-based index of the symmetry operation.
Returns: Non-zero if `Nop` is invalid; zero on success.
Errors: Returns non-zero if $Nop$ is outside the valid range $[0, Nops-1]$.

#### `SymOps::Print() -> void`
Prints a textual representation of the symmetry operations object content to standard output.
Receiver: instance — call as `obj->Print()` (or `obj.Print()`).
Side effects: Writes data to the output stream.

#### `SymOps::Copy(PSymOps symOps) -> void`
Copies the state and data from another `PSymOps` object into the current object.
Receiver: instance — call as `obj->Copy(symOps)` (or `obj.Copy(symOps)`).
Params: `symOps` — The source symmetry operations object pointer.
Side effects: Replaces the object's current symmetry operations with those from `symOps`.

#### `SymOps::write(io::RFile f) -> void`
Writes the object's current symmetry operation data to the specified file stream.
Receiver: instance — call as `obj->write(...)` (or `obj.write(...)`).
Params: `f` — The output file stream.
Side effects: Writes object state to the file stream $f$.

#### `SymOps::read(io::RFile f) -> void`
Reads symmetry operation data from the specified file stream into the object.
Receiver: instance — call as `obj->read(...)` (or `obj.read(...)`).
Params: `f` — The input file stream.
Side effects: Overwrites object state from the file stream $f$.

---

## `Title`
Manages and provides access to metadata related to the deposition of the structural coordinates (e.g., journal, source, compound).

#### `Title::GetJournal() -> PTitleContainer`
Returns the container holding journal record details from the PDB file header.
Receiver: instance — call as `obj->GetJournal()`
Returns: A pointer to the `TitleContainer` holding journal data.
Side effects: None.
Errors: None.

---

## `TitleContainer`
Represents a container holding structural elements related to a title or collection of structures.

#### `TitleContainer::TitleContainer()`
Constructs a default TitleContainer instance.
Receiver: instance — call as `obj->TitleContainer()`
Returns: void

#### `TitleContainer::TitleContainer(io::RPStream Object)`
Constructs a TitleContainer from an input stream representing the structure data.
Receiver: instance — call as `obj->TitleContainer(Object)`
Params: `Object` — The input stream containing the container data.
Returns: void
Side effects: Initializes internal container state from the provided stream.

#### `TitleContainer::MakeContainerClass(int ClassID)`
Creates and returns a new container class instance matching the specified ID.
Receiver: instance — call as `obj->MakeContainerClass(ClassID)`
Params: `ClassID` — The unique identifier for the desired container type.
Returns: A pointer to the newly created container class object.
Side effects: None.

---

## `Vertex`
Represents a structural node within a graph structure, holding chemical and topological information.

#### `Vertex::GetUserID(void) -> int`
Retrieves the user-assigned unique identifier for this vertex.
Receiver: instance — call as `obj->GetUserID()` (or `obj.GetUserID()`).
Returns: The user ID integer.
Side effects: None.

#### `Vertex::GetName(void) -> cpstr`
Retrieves the chemical name associated with this vertex.
Receiver: instance — call as `obj->GetName()` (or `obj.GetName()`).
Returns: The vertex name string.
Side effects: None.

#### `Vertex::GetType(void) -> int`
Retrieves the fundamental structural type identifier for this vertex.
Receiver: instance — call as `obj->GetType()` (or `obj.GetType()`).
Returns: The vertex type integer.
Side effects: None.

#### `Vertex::GetTypeExt(void) -> int`
Retrieves the extended structural type identifier for this vertex.
Receiver: instance — call as `obj->GetTypeExt()` (or `obj.GetTypeExt()`).
Returns: The extended vertex type integer.
Side effects: None.

---
