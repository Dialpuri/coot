# API Reference: `/lmb/home/jdialpuri/Development/mmdb2/mmdb2`

---

## API Index

---

## `mmdb_matrix`

> This class represents a 4x4 matrix, providing methods for matrix operations and data manipulation.

### Methods

#### `Matrix()`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Constructs a new `Matrix` instance with all elements initialized to zero. |
| **Parameters** | `none`                                                              |
| **Returns**  | `void`                                                              |
| **Behaviour** | Initializes all the elements of the `Matrix` array to zero. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `Matrix(realtype data[16])`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Creates a new `Matrix` instance and initializes it with the provided data. |
| **Parameters** | `data` — A pointer to a 16-element array of `realtype` values.      |
| **Returns**  | `void`                                                              |
| **Behaviour** | Copies the data from the input array into the `Matrix` array. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `~Matrix()`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Destructor for the `Matrix` class, releasing allocated memory.     |
| **Parameters** | `none`                                                              |
| **Returns**  | `void`                                                              |
| **Behaviour** |  The destructor cleans up any resources associated with the matrix object. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `Matrix Multiply(Matrix other)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Multiplies this `Matrix` with another `Matrix` and returns the resulting matrix. |
| **Parameters** | `other` — A pointer to another `Matrix` object.                      |
| **Returns**  | `Matrix` — A new `Matrix` object representing the result of the multiplication. |
| **Behaviour** | Performs matrix multiplication on the two input matrices, using the standard matrix multiplication algorithm. |
| **Side effects** | Modifies the elements of the result matrix. |
| **Errors / exceptions** | `none`                                                              |

#### `Matrix Inverse()`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Calculates the inverse of this `Matrix` and returns the resulting matrix. |
| **Parameters** | `none`                                                              |
| **Returns**  | `Matrix` — A new `Matrix` object representing the inverse of this matrix. |
| **Behaviour** |  Performs the calculation of the matrix inverse. |
| **Side effects** | Modifies the elements of the inverse matrix.  |
| **Errors / exceptions** | `none`                                                              |

#### `void SetElement(int row, int col, realtype value)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Sets the value of the element at the specified row and column.       |
| **Parameters** | `row` — The row index (0-based). `col` — The column index (0-based). `value` — The value to set. |
| **Returns**  | `void`                                                              |
| **Behaviour** | Sets the value of the element at the specified `row` and `col` in the matrix. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `realtype GetElement(int row, int col)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Gets the value of the element at the specified row and column.       |
| **Parameters** | `row` — The row index (0-based). `col` — The column index (0-based). |
| **Returns**  | `realtype` — The value of the element at the specified `row` and `col`. |
| **Behaviour** | Retrieves the value of the element at the specified `row` and `col` in the matrix. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `void Print()`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Prints the contents of the `Matrix` to the standard output.           |
| **Parameters** | `none`                                                              |
| **Returns**  | `void`                                                              |
| **Behaviour** |  Prints the elements of the `Matrix` array to the console, row by row.  |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `void Copy(Matrix other)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Copies the state of another `Matrix` into this `Matrix`.              |
| **Parameters** | `other` — A pointer to another `Matrix` object whose state is to be copied. |
| **Returns**  | `void`                                                              |
| **Behaviour** | Copies the values from the other matrix to this matrix. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `void write(io::RFile f)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Writes the content of the `Matrix` to the output stream `f`.           |
| **Parameters** | `f` — An `io::RFile` object to write to. |
| **Returns**  | `void`                                                              |
| **Behaviour** | Writes the value of each element in the matrix to the file stream `f`. |
| **Side effects** | Writes the data to the output stream `f`. |
| **Errors / exceptions** | `none`                                                              |

#### `void read(io::RFile f)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Reads the content of the `Matrix` from the input stream `f`.           |
| **Parameters** | `f` — An `io::RFile` object to read from. |
| **Returns**  | `void`                                                              |
| **Behaviour** | Reads the value of each element in the matrix from the input file stream `f`. |
| **Side effects** | Reads the data from the input stream `f`. |
| **Errors / exceptions** | `none`                                                              |

---

## Free Functions



---

## Per-File Method Documentation

### `hybrid_36`

## `hy36encode`

> Converts an integer value to a string representation using the hybrid-36 encoding scheme, suitable for residue sequence numbers (4 characters) or atom serial numbers (5 characters).

### Methods

#### `const char* hy36encode(unsigned width, int value, char* result)`

| Field        | Detail                                                                                                                                                                                                                                         |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Converts an integer `value` to a string representation using the hybrid-36 encoding scheme, storing the result in the provided `result` buffer.  The result is guaranteed to be null-terminated.                                                        |
| **Parameters** | `width` — `unsigned` — The desired width of the output string (4 for residue sequence numbers, 5 for atom serial numbers).  If `width` is not 4 or 5, it returns an error message.  `value` — `int` — The integer value to be converted to a hybrid-36 string.  `result` — `char*` — A pointer to a character array where the resulting hybrid-36 string will be stored. This array must be large enough to hold the result plus the null terminator. |
| **Returns**  | `const char*` — A pointer to a string containing the hybrid-36 encoded representation of the input integer `value`. Returns `0` on success, and a pointer to an error message string (e.g., "value out of range.", "unsupported width.") on failure. |
| **Behaviour** | 1. **Handles different widths:** The function first checks the `width` parameter to determine the appropriate encoding scheme.  It handles both 4-character (residue sequence) and 5-character (atom serial) encodings. 2. **Handles ranges:**  Within each width, the function handles negative and positive values using a series of checks. 3. **`encode_pure` call:** For width 4, it performs the encoding by dividing the input value by 36 to determine the number of digits.  It then calculates the digits using the `digits_upper()` string. For width 5, similar logic is applied.  The result is built using the `encode_pure` helper function. 4. **Padding with asterisks:** If the resulting string is shorter than `width`, the remaining characters are filled with asterisks (`*`) using the `fill_with_stars` helper function. 5. **Error Handling:** Returns appropriate error messages via the function's return value when an invalid value is encountered or when the supplied width is invalid. |
| **Side effects** | None.                                                                                                                                                                                                                                           |
| **Errors / exceptions** | Returns a pointer to an error message string if:  - `width` is not 4 or 5.  - The input `value` is out of the representable range for the given `width`. -  An invalid character is encountered in the input string.   |

#### `const char* hy36decode(unsigned width, const char* s, unsigned s_size, int* result)`

| Field        | Detail                                                                                                                                                                                                                                                           |
|--------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Converts a string representation (in hybrid-36 format) to an integer value, storing the result in the provided `result` integer. The input string `s` is expected to be null-terminated, though it does not need to be. The result is guaranteed to be null-terminated. |
| **Parameters** | `width` — `unsigned` — The expected width of the hybrid-36 string (4 for residue sequence numbers, 5 for atom serial numbers).  This parameter is used to determine the length of the string to be parsed. `s` — `const char*` — A pointer to the string to be decoded.  This string should contain a hybrid-36 encoded representation of the integer value. `s_size` — `unsigned` — The size of the string pointed to by `s` (i.e., the number of characters in the string).  This *must* match the `width` parameter. `result` — `int*` — A pointer to an integer variable where the decoded integer value will be stored. |
| **Returns**  | `const char*` — A pointer to a string containing an error message, if any, or `0` on successful decoding. |
| **Behaviour** | 1. **Internal Error Handling:** The function contains internal error handling routines to check for internal inconsistencies. The 'first_call' flag and associated initialization are used to avoid a one-time expense. 2. **Initialization:** The digits_values_upper and digits_values_lower arrays are initialized during the first call to the function to store the corresponding values for the upper and lower case digits. 3. **Checks `s_size`:** The function first checks if `s_size` matches `width`. If they don't match, it returns an error message. 4. **Iteration & Parsing:** If `s_size` is correct, the function iterates through the string `s`.  For each character, it determines whether it's a digit or a minus sign.  If it's a digit, the corresponding value is added to the integer `value`. If it's a minus sign, the value is negated, and the 'have_minus' flag is set to true. 5. **`decode_pure` call:** The result of each iteration is passed to the `decode_pure` helper function. This function performs the actual conversion. 6. **Error Handling:** Returns appropriate error messages via the function's return value when an invalid character is encountered or when the supplied width is invalid.  |
| **Side effects** | None.                                                                                                                                                                                                                                           |
| **Errors / exceptions** | Returns a pointer to an error message string if:  - `width` is not 4 or 5.  - `s_size` does not equal `width`. - An invalid character is encountered in the input string (i.e., a character that's not a digit or a minus sign). - The value decoded is out of range. |

---

## `digits_upper`

> A constant string containing the upper-case hexadecimal digits used in the hybrid-36 encoding scheme.

### Methods

#### `const char* digits_upper()`

| Field        | Detail                               |
|--------------|--------------------------------------|
| **Summary**  | Returns a pointer to a string containing the uppercase hexadecimal digits ('0' through '9' and 'A' through 'Z').|
| **Parameters** | None                                 |
| **Returns**  | `const char*` — A pointer to the string "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ".|
| **Behaviour** | Simply returns the string literal.|
| **Side effects** | None.                              |
| **Errors / exceptions** | None.                              |
---

## `digits_lower`

> A constant string containing the lower-case hexadecimal digits used in the hybrid-36 encoding scheme.

### Methods

#### `const char* digits_lower()`

| Field        | Detail                               |
|--------------|--------------------------------------|
| **Summary**  | Returns a pointer to a string containing the lowercase hexadecimal digits ('0' through '9' and 'a' through 'z').|
| **Parameters** | None                                 |
| **Returns**  | `const char*` — A pointer to the string "0123456789abcdefghijklmnopqrstuvwxyz".|
| **Behaviour** | Simply returns the string literal.|
| **Side effects** | None.                              |
| **Errors / exceptions** | None.                              |
---

## `value_out_of_range`

> A constant string error message to be returned when the input integer value is out of the range representable by the hybrid-36 encoding scheme.

### Methods

#### `const char* value_out_of_range()`

| Field        | Detail                               |
|--------------|--------------------------------------|
| **Summary**  | Returns a pointer to a string indicating that the input integer value is out of range for the specified width. |
| **Parameters** | None                                 |
| **Returns**  | `const char*` — A pointer to the string "value out of range.".|
| **Behaviour** | Simply returns the string literal.|
| **Side effects** | None.                              |
| **Errors / exceptions** | None.                              |
---

## `invalid_number_literal`

> A constant string error message to be returned when an invalid character is encountered during the conversion process.

### Methods

#### `const char* invalid_number_literal()`

| Field        | Detail                               |
|--------------|--------------------------------------|
| **Summary**  | Returns a pointer to a string indicating that an invalid character was encountered during the conversion process. |
| **Parameters** | None                                 |
| **Returns**  | `const char*` — A pointer to the string "invalid number literal.".|
| **Behaviour** | Simply returns the string literal.|
| **Side effects** | None.                              |
| **Errors / exceptions** | None.                              |
---

## `fill_with_stars`

> Fills a character buffer with asterisks to the specified width.

### Methods

#### `void fill_with_stars(unsigned width, char* result)`

| Field        | Detail                                                                                             |
|--------------|-----------------------------------------------------------------------------------------------------|
| **Summary**  | Fills a character buffer `result` with asterisks until the buffer's width is reached.            |
| **Parameters** | `width` — `unsigned` — The width of the buffer to fill.  `result` — `char*` — A pointer to the buffer to be filled. |
| **Returns**  | `void`                                                                                             |
| **Behaviour** | Uses a `while` loop to iterate `width` times, copying the asterisk character (`*`) into the `result` buffer.  A null terminator is added at the end. |
| **Side effects** | Modifies the buffer pointed to by `result`.                                                     |
| **Errors / exceptions** | None.                                                                                          |
---

## `encode_pure`

> Encodes a integer into a hybrid-36 string.

### Methods

#### `void encode_pure(const char* digits, unsigned digits_size, unsigned width, int value, char* result)`

| Field        | Detail                                                                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Encodes an integer `value` into a hybrid-36 string, storing the result in the provided `result` buffer.                                                    |
| **Parameters** | `digits` — `const char*` — A string containing the hexadecimal digits to use for encoding. `digits_size` — `unsigned` — The size of the `digits` string. `width` — `unsigned` — The desired width of the output string (4 for residue sequence numbers, 5 for atom serial numbers). `value` — `int` — The integer value to be encoded. `result` — `char*` — A pointer to a character array where the resulting hybrid-36 string will be stored. |
| **Returns**  | `void`                                                                                                                                                    |
| **Behaviour** | 1. Divides the integer value `value` by `digits_size` to determine the number of digits. 2.  Converts the integer value to a string representation using the `digits` string. 3. Pads with asterisks if the number of digits is less than `width`.  4. Stores the result in the `result` buffer.|
| **Side effects** | Modifies the buffer pointed to by `result`.                                                                                                               |
| **Errors / exceptions** | None.                                                                                                                                                |
---

## `decode_pure`

> Decodes a string into an integer.

### Methods

#### `const char* decode_pure(const int* digits_values, unsigned digits_size, const char* s, unsigned s_size, int* result)`

| Field        | Detail                                                                                                                                                                                                                                  |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Decodes a hybrid-36 string (e.g., "A1T5") into an integer value, storing the result in the provided `result` integer.  The input string `s` is expected to be null-terminated.                                                   |
| **Parameters** | `digits_values` — `const int*` — An array of integers representing the hexadecimal digits (0-9 and a-z). `digits_size` — `unsigned` — The size of the `digits_values` array (36). `s` — `const char*` — A pointer to the string to be decoded. `s_size` — `unsigned` — The size of the string pointed to by `s`. `result` — `int*` — A pointer to an integer variable where the decoded integer value will be stored. |
| **Returns**  | `const char*` — A pointer to a string containing an error message, if any, or `0` on successful decoding.                                                                                                                            |
| **Behaviour** | 1. Handles different width parameters (4 or 5) to decode strings. 2. Handles the case where `s` has a length less than the `width` parameter. 3. Iterates through the characters in the input string `s`. 4. Checks if the character is a digit or a minus sign. 5. If it's a digit, the corresponding value is added to the integer `value`. 6. If it's a minus sign, the `value` is negated. 7. Returns the integer value.    |
| **Side effects** | Modifies the integer `value`.                                                                                                                                                                                                                 |
| **Errors / exceptions** | Returns a pointer to an error message string if:  - The input string `s` is not a valid hybrid-36 string. - The `digits_values` array is invalid. - A character in `s` is not a digit or a minus sign. - The width is invalid. |
---

## Free Functions

There are no free functions.


---

### `mmdb_atom`

Okay, here's the complete documentation output, incorporating all methods, internal behavior, and side effects, formatted as requested.  I've added detail where necessary, particularly regarding the `Atom` class, based on the code.

---
## `mmdb::Atom`

> Represents an atom in a protein or molecule, storing its coordinates, occupancy, temperature factor, and other relevant properties.

### Methods

#### `void Atom::InitAtom()`

| Field        | Detail                                                                                                                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Initializes the atom's internal data structures, setting all values to default or placeholder values.                                                                                                    |
| **Parameters** | none                                                                                                                                                                                                  |
| **Returns**  | `void`                                                                                                                                                                                                |
| **Behaviour** | Sets all member variables to their default values: Serial number (`serNum`) to -1, index to -1, name and label_atom_id to empty strings, altLoc to empty string, segID to empty string, element to "  ", occupancy to 0.0, and so on. Initializes the `WhatIsSet` flag to `ASET_All`. It also sets `finished` to `false`|
| **Side effects** | Resets all internal data fields of the `Atom` object, ensuring a clean state before data is read or set. |
| **Errors / exceptions** | None                                                                                                                                                                                                |

#### `void Atom::FreeMemory()`

| Field        | Detail                                                                                                                                                                                                                                                                                     |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Releases memory allocated for the `Atom` object's internal data structures, including the `Bond` array.                                                                                                                                                                            |
| **Parameters** | none                                                                                                                                                                                                                                                                                     |
| **Returns**  | `void`                                                                                                                                                                                                                                                                                     |
| **Behaviour** | Deletes the `Bond` array, setting all pointers to `NULL` to prevent dangling pointers.                                                                                                                                                                                                 |
| **Side effects** | Releases the memory occupied by the `Bond` array.  Also resets the  `nBonds` to 0.                                                                                                                                                                                                  |
| **Errors / exceptions** | None                                                                                                                                                                                                                                                                                     |

#### `int Atom::AddBond(PAtom bond_atom, int bond_order, int nAdd_bonds)`

| Field        | Detail                                                                                                                                                                                                                                                                    |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Adds a bond to the atom, connecting it to the specified `bond_atom` with the given `bond_order`.                                                                                                                                                                              |
| **Parameters** | `bond_atom` — A pointer to the `PAtom` that this atom is bonded to.  `bond_order` — The order of the bond (e.g., 1 for single, 2 for double, 3 for triple). `nAdd_bonds` — The number of bonds to add (defaulting to 1). |
| **Returns**  | `int` — The number of bonds that were added. Returns `-k` if the atom is already bonded to `bond_atom`, where `k` is the index of the already-existing bond (or -1 if no bonds exist). Returns 0 if the addition was successful.  |
| **Behaviour** | Checks if an atom with the same `bond_atom` already exists.  If so, returns -index. Otherwise, adds the new bond to the `Bond` array.  The `nBonds` count is incremented. If `nAdd_bonds` is greater than the current number of bonds, the `Bond` array is expanded to accommodate the new bonds. |
| **Side effects** | Modifies the `Bond` array by adding a new bond element. Increments the `nBonds` count.                                                                                                                                                                                |
| **Errors / exceptions** | None                                                                                                                                                                                                                                                                                     |

#### `int Atom::GetNBonds()`

| Field        | Detail                                                                                                                                     |
|--------------|--------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the number of bonds associated with the atom.                                                                                        |
| **Parameters** | none                                                                                                                                      |
| **Returns**  | `int` — The number of bonds.                                                                                                               |
| **Behaviour** | Returns the value of `nBonds & 0x000000FF`.                                                                                             |
| **Side effects** | None                                                                                                                                      |
| **Errors / exceptions** | None                                                                                                                                      |

#### `void Atom::GetBonds(RPAtomBond atomBond, int & nAtomBonds)`

| Field        | Detail                                                                                                                                                                                          |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns a pointer to the atom's internal bond structure.                                                                                                                                          |
| **Parameters** | `atomBond` — A pointer to the `RPAtomBond` object where the bond information will be stored.  `nAtomBonds` — An integer reference to store the number of bonds.                       |
| **Returns**  | `void`                                                                                                                                                                                        |
| **Behaviour** | Copies the contents of the `Bond` array into the `atomBond` object, and populates the `nAtomBonds` reference with the number of bonds.                                                           |
| **Side effects** | Modifies the `atomBond` object by populating it with bond data.                                                                                                                              |
| **Errors / exceptions** | None                                                                                                                                                                                        |

#### `void Atom::GetBonds(RPAtomBondI atomBondI, int & nAtomBonds)`

| Field        | Detail                                                                                                                                                                                          |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Disposes of the `AtomBondI` object and returns a pointer to an array of `AtomBondI` structures representing the atom's bonds.                                                              |
| **Parameters** | `atomBondI` — A pointer to the `RPAtomBondI` object where the bond information will be stored.  `nAtomBonds` — An integer reference to store the number of bonds.                                |
| **Returns**  | `void`                                                                                                                                                                                        |
| **Behaviour** | Deletes the `AtomBondI` object if it exists.  Allocates a new `AtomBondI` array of the appropriate size, and copies the bond information from the `Bond` array into it. Populates the `nAtomBonds` reference. |
| **Side effects** | Deallocates the memory for the `AtomBondI` object. Allocates memory for the `AtomBondI` array. Copies bond data.                                                                         |
| **Errors / exceptions** | None                                                                                                                                                                                        |

#### `void Atom::GetBonds(PAtomBondI AtomBondI, int & nAtomBonds, int maxlength)`

| Field        | Detail                                                                                                                                                                                          |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Does not dispose or allocate  `AtomBondI` structure. It assumes that the length of `AtomBondI` is sufficient to accommodate all bonds.                                                              |
| **Parameters** | `AtomBondI` — A pointer to the `PAtomBondI` object where the bond information will be stored.  `nAtomBonds` — An integer reference to store the number of bonds. `maxlength`— The maximum length of `AtomBondI`. |
| **Returns**  | `void`                                                                                                                                                                                        |
| **Behaviour** | Copies the contents of the `Bond` array into the `AtomBondI` object, and populates the `nAtomBonds` reference with the number of bonds.                                                           |
| **Side effects** | Modifies the `AtomBondI` object by populating it with bond data.                                                                                                                              |
| **Errors / exceptions** | None                                                                                                                                                                                        |

#### `int Atom::GetResidueNo()`

| Field        | Detail                                                                                                                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the residue number of the atom.                                                                                                                                                            |
| **Parameters** | none                                                                                                                                                                                              |
| **Returns**  | `int` — The residue number.                                                                                                                                                                      |
| **Behaviour** | Returns the value of `residue->seqNum`.                                                                                                                                                            |
| **Side effects** | None                                                                                                                                                                                              |
| **Errors / exceptions** | None                                                                                                                                                                                              |

#### `void Atom::SetResidue(PResidue res)`

| Field        | Detail                                                                                                                                                                                          |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the `residue` pointer to the specified `PResidue` object.                                                                                                                                       |
| **Parameters** | `res` — A pointer to the `PResidue` object to be associated with this atom.                                                                                                                       |
| **Returns**  | `void`                                                                                                                                                                                        |
| **Behaviour** | Assigns the `residue` pointer to the `res` parameter.                                                                                                                                               |
| **Side effects** | Modifies the `residue` pointer.                                                                                                                                                             |
| **Errors / exceptions** | None                                                                                                                                                                                        |

#### `void Atom::PDBASCIIDump(io::RFile f)`

| Field        | Detail                                                                                                                               |
|--------------|--------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Generates a PDB-ASCII style line for the atom, including its serial number, name, coordinates, occupancy, temperature factor, etc. |
| **Parameters** | `f` — An `io::RFile` object to write the line to.                                                                                    |
| **Returns**  | `void`                                                                                                                               |
| **Behaviour** |  Constructs a PDB-ASCII formatted string representing the atom's data and writes it to the output file.  It handles the serial number formatting and other required parameters. |
| **Side effects** | Writes to the output file.                                                                                                       |
| **Errors / exceptions** | None                                                                                                                               |

#### `void Atom::MakeAtomCIF(mmcif::PData CIF)`

| Field        | Detail                                                                                                                        |
|--------------|-------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Generates the CIF (Chemical Interchange Format) representation of the atom.                                                     |
| **Parameters** | `CIF` — A pointer to the `mmcif::PData` object to store the CIF data.                                                      |
| **Returns**  | `void`                                                                                                                         |
| **Behaviour** |  The implementation of this method is not provided.                                                                            |
| **Side effects** | None.                                                                                                                      |
| **Errors / exceptions** | None                                                                                                                      |

#### `void Atom::FreeBonds()`

| Field        | Detail                                                                                                                                 |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Releases the memory allocated for the `Bond` array.                                                                                  |
| **Parameters** | none                                                                                                                                 |
| **Returns**  | `void`                                                                                                                                |
| **Behaviour** | Sets the `Bond` pointer to `NULL` and resets `nBonds` to 0.                                                                            |
| **Side effects** | Deallocates the memory for the `Bond` array.                                                                                         |
| **Errors / exceptions** | None                                                                                                                                |

#### `void Atom::SetAtomName(int ix, int sN, const AtomName aName, const AltLoc aLoc, const SegID sID, const Element eName)`

| Field        | Detail                                                                                                                                                                                                    |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the atom's name, serial number, and alternative location.                                                                                                                                            |
| **Parameters** | `ix` — The index of the atom in the `atom` array. `sN` — The serial number of the atom. `aName` — The atom name. `aLoc` — The alternative location indicator. `sID` — The segment identifier. `eName` — The chemical element name.|
| **Returns**  | `void`                                                                                                                                                                                                |
| **Behaviour** | Updates the `index`, `serNum`, `name`, `label_atom_id`, `altLoc`, `segID`, and `element` member variables with the given values.                                                                                |
| **Side effects** | Modifies the `index`, `serNum`, `name`, `label_atom_id`, `altLoc`, `segID`, and `element` member variables.                                                                                                |
| **Errors / exceptions** | None                                                                                                                                                                                                |

#### `void Atom::ConvertPDBATOM(int ix, cpstr S)`

| Field        | Detail                                                                                                                                                                       |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads atom data from a PDB ASCII "ATOM" record.                                                                                                                              |
| **Parameters** | `ix` — The index of the atom in the `atom` array. `S` — A pointer to the PDB ASCII string representing the "ATOM" record.                                                       |
| **Returns**  | `void`                                                                                                                                                                    |
| **Behaviour** | Reads the coordinates (x, y, z), occupancy, temperature factor, and other related data from the PDB ASCII "ATOM" record, storing them in the `Atom` object's member variables. |
| **Side effects** | Modifies the `x`, `y`, `z`, `occupancy`, `tempFactor`, and `sigX`, `sigY`, `sigZ` member variables.  Also sets the `WhatIsSet` flag.                                           |
| **Errors / exceptions** | None                                                                                                                                                                    |

#### `void Atom::ConvertPDBSIGATM(int ix, cpstr S)`

| Field        | Detail                                                                                                                                                                                   |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads atom data from a PDB ASCII "SIGATM" record.                                                                                                                                       |
| **Parameters** | `ix` — The index of the atom in the `atom` array. `S` — A pointer to the PDB ASCII string representing the "SIGATM" record.                                                              |
| **Returns**  | `void`                                                                                                                                                                             |
| **Behaviour** | Reads the standard deviation values (sigX, sigY, sigZ), occupancy standard deviation, and temperature factor standard deviation from the PDB ASCII "SIGATM" record, storing them in the `Atom` object's member variables. |
| **Side effects** | Modifies the `sigX`, `sigY`, `sigZ`, `sigOcc`, `sigTemp` member variables. Also sets the `WhatIsSet` flag.                                                                              |
| **Errors / exceptions** | None                                                                                                                                                                             |

#### `void Atom::ConvertPDBANISOU(int ix, cpstr S)`

| Field        | Detail                                                                                                                                                                                         |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads atom data from a PDB ASCII "ANISOU" record.                                                                                                                                                |
| **Parameters** | `ix` — The index of the atom in the `atom` array. `S` — A pointer to the PDB ASCII string representing the "ANISOU" record.                                                                |
| **Returns**  | `void`                                                                                                                                                                                      |
| **Behaviour** | Reads the anisotropic temperature factors (u11, u22, u33, u12, u13, u23) from the PDB ASCII "ANISOU" record, storing them in the `Atom` object's member variables.                                  |
| **Side effects** | Modifies the `u11`, `u22`, `u33`, `u12`, `u13`, and `u23` member variables. Also sets the `WhatIsSet` flag.                                                                                   |
| **Errors / exceptions** | None                                                                                                                                                                                      |

#### `void Atom::ConvertPDBSIGUIJ(int ix, cpstr S)`

| Field        | Detail                                                                                                                                                                   |
|--------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads atom data from a PDB ASCII "SIGUIJ" record.                                                                                                                          |
| **Parameters** | `ix` — The index of the atom in the `atom` array. `S` — A pointer to the PDB ASCII string representing the "SIGUIJ" record.                                                |
| **Returns**  | `void`                                                                                                                                                                  |
| **Behaviour** | Reads the sigma values for the anisotropic temperature factors (su11, su22, su33, su12, su13, su23) from the PDB ASCII "SIGUIJ" record, storing them in the `Atom` object's member variables.|
| **Side effects** | Modifies the `su11`, `su22`, `su33`, `su12`, `su13`, and `su23` member variables. Also sets the `WhatIsSet` flag.                                                        |
| **Errors / exceptions** | None                                                                                                                                                                  |

#### `ERROR_CODE Atom::ConvertPDBTER(int ix, cpstr S)`

| Field        | Detail                                                                                                                                                                 |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads atom data from a PDB ASCII "TER" record.                                                                                                                          |
| **Parameters** | `ix` — The index of the atom in the `atom` array. `S` — A pointer to the PDB ASCII string representing the "TER" record.                                              |
| **Returns**  | `ERROR_CODE` — An error code indicating the success or failure of the operation.                                                                                      |
| **Behaviour** | Reads the serial number, occupancy, temperature factor, and other related data from the PDB ASCII "TER" record, storing them in the `Atom` object's member variables. |
| **Side effects** | Modifies the `serNum`, `occupancy`, `tempFactor`, and `WhatIsSet` member variables.                                                                                   |
| **Errors / exceptions** |  Return `Error_ATOM_AlreadySet` if atom is already set. Return `Error_ATOM_Unrecognized` if reading value fails.                                                              |

#### `ERROR_CODE Atom::GetAtomSiteCIF(int ix, mmcif::PLoop Loop)`

| Field        | Detail                                                                                                                                                                                                                                                                               |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads the atom site information from the CIF file.                                                                                                                                                                                                                                |
| **Parameters** | `ix` — The index of the atom in the CIF file. `Loop` — A pointer to the CIF loop object.                                                                                                                                                                                           |
| **Returns**  | `ERROR_CODE` — An error code indicating the success or failure of the operation.                                                                                                                                                                                                   |
| **Behaviour** | Retrieves the atom ID from the CIF file using `Loop->GetString`. The retrieved atom ID is assigned to the `atom_site_id` member variable of the atom object.                                                                                                                            |
| **Side effects** | Modifies the `atom_site_id` member variable.                                                                                                                                                                                                                                    |
| **Errors / exceptions** | Returns `Error_EmptyCIFLoop` if `Loop` is empty.                                                                                                                                                                                                                               |

#### `ERROR_CODE Atom::GetremapCIF (int ix, mmcif::PLoop Loop)`

| Field        | Detail                                                                                                                                                                                                                                                                             |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Read CIF data that remaps related atom ID, chain and residue.                                                                                                                                                                                                                         |
| **Parameters** | `ix` — The index of the atom in the CIF file. `Loop` — A pointer to the CIF loop object.                                                                                                                                                                                           |
| **Returns**  | `ERROR_CODE` — An error code indicating the success or failure of the operation.                                                                                                                                                                                                  |
| **Behaviour** | Retrieves orig_label_asym_id,  label_comp_id, label_asym_id from the CIF loop, and assigns them to respective variables.                                                                                                                                                            |
| **Side effects** | Modifies `orig_label_asym_id`,`label_comp_id`, `label_asym_id`                                                                                                                                                                                                              |
| **Errors / exceptions** | Return `Error_EmptyCIFLoop` if `Loop` is empty.                                                                                                                                                                                                                             |

#### `ERROR_CODE Atom::GetAtomSiteCIF (int ix, mmcif::PLoop Loop)`

| Field        | Detail                                                                                                                                                                                                                                                                                  |
|--------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads atom site information from CIF.                                                                                                                                                                                                                                               |
| **Parameters** | `ix` — The index of the atom in the CIF file. `Loop` — A pointer to the CIF loop object.                                                                                                                                                                                          |
| **Returns**  | `ERROR_CODE` — An error code indicating the success or failure of the operation.                                                                                                                                                                                                  |
| **Behaviour** | Retrieves atom site ID from the CIF file using `Loop->GetString`. The retrieved atom site ID is assigned to the `atom_site_id` member variable of the atom object.                                                                                                                     |
| **Side effects** | Modifies the `atom_site_id` member variable.                                                                                                                                                                                                                                   |
| **Errors / exceptions** | Returns `Error_EmptyCIFLoop` if `Loop` is empty.                                                                                                                                                                                                                              |

#### `void Atom::MakeremapCIF (mmcif::PData CIF)`

| Field        | Detail                                                                                                                                                                                                        |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | This method is used to apply the changes of the atom from the CIF file.                                                                                                                                     |
| **Parameters** | `CIF` — A pointer to the mmcif::PData object, which represents the CIF file.                                                                                                                            |
| **Returns**  | `void`                                                                                                                                                                                                      |
| **Behaviour** | The method calls `mmcif::CIF->AddLoop` to create the loop, and sets labels from different entities.                                                                                                            |
| **Side effects** | Modifies the data stored in CIF.                                                                                                                                                                           |
| **Errors / exceptions** | None                                                                                                                                                                                                      |

#### `void Atom::SetAtomName(int ix, int sN, const AtomName aName, const AltLoc aLoc, const SegID sID, const Element eName)`

| Field        | Detail                                                                                                                                            |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the atom's name, serial number, and alternative location.                                                                                        |
| **Parameters** | `ix` — The index of the atom in the `atom` array. `sN` — The serial number of the atom. `aName` — The atom name. `aLoc` — The alternative location indicator. `sID` — The segment identifier. `eName` — The chemical element name. |
| **Returns**  | `void`                                                                                                                                               |
| **Behaviour** | Updates the `index`, `serNum`, `name`, `label_atom_id`, `altLoc`, `segID`, and `element` member variables with the given values.                     |
| **Side effects** | Modifies the `index`, `serNum`, `name`, `label_atom_id`, `altLoc`, `segID`, and `element` member variables.                                       |
| **Errors / exceptions** | None                                                                                                                                               |

#### `void Atom::SetElementName(cpstr chrg)`

| Field        | Detail                                                                                                                                     |
|--------------|--------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the element name.                                                                                                                        |
| **Parameters** | `chrg` — A pointer to the character string representing the element name.                                                                   |
| **Returns**  | `void`                                                                                                                                      |
| **Behaviour** | Sets the `element` member variable to the given element name.                                                                               |
| **Side effects** | Modifies the `element` member variable.                                                                                                    |
| **Errors / exceptions** | None                                                                                                                                      |

#### `void Atom::SetCharge(cpstr chrg)`

| Field        | Detail                                                                                                                                                                              |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the atom's charge.                                                                                                                                                           |
| **Parameters** | `chrg` — A pointer to the character string representing the charge.                                                                                                                |
| **Returns**  | `void`                                                                                                                                                                             |
| **Behaviour** | Formats the charge value as a positive string (e.g., "+1") and assigns it to the `charge` member variable.                                                                        |
| **Side effects** | Modifies the `charge` member variable.                                                                                                                                         |
| **Errors / exceptions** | None                                                                                                                                                                             |

#### `void Atom::SetAtomIndex(int ix)`

| Field        | Detail                                                                                                                                                           |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the index of the atom.                                                                                                                                        |
| **Parameters** | `ix` — The index of the atom.                                                                                                                                    |
| **Returns**  | `void`                                                                                                                                                            |
| **Behaviour** | Sets the `index` member variable to the given index.                                                                                                               |
| **Side effects** | Modifies the `index` member variable.                                                                                                                             |
| **Errors / exceptions** | None                                                                                                                                                            |

#### `void Atom::MakeTer()`

| Field        | Detail                                                                                                                                                                         |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the `Ter` flag to `true`, indicating that the atom is a chain terminator.                                                                                                  |
| **Parameters** | none                                                                                                                                                                         |
| **Returns**  | `void`                                                                                                                                                                       |
| **Behaviour** | Sets the `Ter` member variable to `true`.                                                                                                                                      |
| **Side effects** | Modifies the `Ter` member variable.                                                                                                                                       |
| **Errors / exceptions** | None                                                                                                                                                                       |

#### `void Atom::SetCoordinates(realtype xx, realtype yy, realtype zz, realtype occ, realtype tFac)`

| Field        | Detail                                                                                                                                                                                                   |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the atom's orthogonal coordinates, occupancy, and temperature factor.                                                                                                                            |
| **Parameters** | `xx` — The x-coordinate. `yy` — The y-coordinate. `zz` — The z-coordinate. `occ` — The occupancy. `tFac` — The temperature factor.                                                               |
| **Returns**  | `void`                                                                                                                                                                                                 |
| **Behaviour** | Sets the `x`, `y`, `z`, `occupancy`, and `tempFactor` member variables to the given values.                                                                                                               |
| **Side effects** | Modifies the `x`, `y`, `z`, `occupancy`, and `tempFactor` member variables.                                                                                                                             |
| **Errors / exceptions** | None                                                                                                                                                                                                 |

#### `int Atom::GetModelNum()`

| Field        | Detail                                                                                                                                 |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the model number of the atom.                                                                                                   |
| **Parameters** | none                                                                                                                                   |
| **Returns**  | `int` — The model number.                                                                                                              |
| **Behaviour** | Returns the value of `residue->model->GetSerNum()`.                                                                                 |
| **Side effects** | None                                                                                                                                   |
| **Errors / exceptions** | None                                                                                                                                   |

#### `pstr Atom::GetChainID()`

| Field        | Detail                                                                                                                       |
|--------------|-----------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the chain ID of the atom.                                                                                             |
| **Parameters** | none                                                                                                                        |
| **Returns**  | `pstr` — The chain ID.                                                                                                        |
| **Behaviour** | Returns the value of `residue->chain->chainID`.                                                                              |
| **Side effects** | None                                                                                                                        |
| **Errors / exceptions** | None                                                                                                                        |

#### `pstr Atom::GetLabelAsymID()`

| Field        | Detail                                                                                                                      |
|--------------|----------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the chain ID of the atom, as it is presented in the asymmetric unit.                                                   |
| **Parameters** | none                                                                                                                        |
| **Returns**  | `pstr` — The chain ID.                                                                                                        |
| **Behaviour** | Returns the value of `residue->label_asym_id`.                                                                               |
| **Side effects** | None                                                                                                                        |
| **Errors / exceptions** | None                                                                                                                        |

#### `pstr Atom::GetOrigLabelAsymID()`

| Field        | Detail                                                                                                                                                                                          |
|--------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the original chain ID of the atom.                                                                                                                                                |
| **Parameters** | none                                                                                                                                                                                           |
| **Returns**  | `pstr` — The chain ID.                                                                                                                                                                         |
| **Behaviour** | Returns the value of `residue->orig_label_asym_id`.                                                                                                                                            |
| **Side effects** | None                                                                                                                                                                                           |
| **Errors / exceptions** | None                                                                                                                                                                                           |

#### `pstr Atom::GetUniprotnum()`

| Field        | Detail                                                                                                                                                                                          |
|--------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the uniprot number of the atom.                                                                                                                                                  |
| **Parameters** | none                                                                                                                                                                                           |
| **Returns**  | `pstr` — The uniprot number.                                                                                                                                                                   |
| **Behaviour** | Returns the value of `residue->pdbx_sifts_xref_db_num`.                                                                                                                                        |
| **Side effects** | None                                                                                                                                                                                           |
| **Errors / exceptions** | None                                                                                                                                                                                           |

#### `pstr Atom::GetUniprotName()`

| Field        | Detail                                                                                                                                            |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the uniprot name of the atom.                                                                                                         |
| **Parameters** | none                                                                                                                                         |
| **Returns**  | `pstr` — The uniprot name.                                                                                                                   |
| **Behaviour** | Returns the value of `residue->pdbx_sifts_xref_db_acc`.                                                                                     |
| **Side effects** | None                                                                                                                                         |
| **Errors / exceptions** | None                                                                                                                                         |

#### `pstr Atom::GetUniprotAcc()`

| Field        | Detail                                                                                                                                                                      |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the uniprot accession number of the atom.                                                                                                                      |
| **Parameters** | none                                                                                                                                                                      |
| **Returns**  | `pstr` — The uniprot accession number.                                                                                                                                  |
| **Behaviour** | Returns the value of `residue->pdbx_sifts_xref_db_acc`.                                                                                                                  |
| **Side effects** | None                                                                                                                                                                      |
| **Errors / exceptions** | None                                                                                                                                                                      |

#### `pstr Atom::GetUniprotnum()`

| Field        | Detail                                                                                                                                                                      |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the uniprot number of the atom.                                                                                                                      |
| **Parameters** | none                                                                                                                                                                      |
| **Returns**  | `pstr` — The uniprot number.                                                                                                                                  |
| **Behaviour** | Returns the value of `residue->pdbx_sifts_xref_db_num`.                                                                                                                  |
| **Side effects** | None                                                                                                                                                                      |
| **Errors / exceptions** | None                                                                                                                                                                      |

#### `pstr Atom::GetAtomSiteID()`

| Field        | Detail                                                                                                                            |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the atom site ID.                                                                                                           |
| **Parameters** | none                                                                                                                             |
| **Returns**  | `pstr` — The atom site ID.                                                                                                      |
| **Behaviour** | Returns the value of `atom_site_id`.                                                                                              |
| **Side effects** | None                                                                                                                             |
| **Errors / exceptions** | None                                                                                                                             |

#### `int Atom::GetModelNum()`

| Field        | Detail                                                                                                                                                                                                                                |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the model number of the atom.                                                                                                                                                                                                 |
| **Parameters** | none                                                                                                                                                                                                                                 |
| **Returns**  | `int` — The model number.                                                                                                                                                                                                          |
| **Behaviour** | Returns the value of `residue->model->GetSerNum()`.                                                                                                                                                                                      |
| **Side effects** | None                                                                                                                                                                                                                                 |
| **Errors / exceptions** | None                                                                                                                                                                                                                                 |

#### `void Atom::ApplyTransform(const mat44 & TMatrix)`

| Field        | Detail                                                                                                                                                                                  |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------

---

### `mmdb_bondmngr`

---
## `mmdb::BondManager`

> This class is responsible for creating and removing chemical bonds within a macromolecular model based on atomic coordinates.

### Methods

#### `BondManager::BondManager()`

| Field        | Detail                                                           |
|--------------|------------------------------------------------------------------|
| **Summary**  | Constructs a `BondManager` object.                               |
| **Parameters** | none                                                             |
| **Returns**  | `void`                                                          |
| **Behaviour** | Initializes the `BondManager` object by calling the base class constructor, `SelManager()`. |
| **Side effects** | none                                                          |
| **Errors / exceptions** | none                                                          |

#### `BondManager::BondManager(io::RPStream Object)`

| Field        | Detail                                                            |
|--------------|------------------------------------------------------------------|
| **Summary**  | Constructs a `BondManager` object, initialized from an `io::RPStream`. |
| **Parameters** | `Object` — An `io::RPStream` object used for stream-based object initialization. |
| **Returns**  | `void`                                                          |
| **Behaviour** | Initializes the `BondManager` object by calling the base class constructor, `SelManager(Object)`. |
| **Side effects** | none                                                          |
| **Errors / exceptions** | none                                                          |

#### `BondManager::~BondManager()`

| Field        | Detail                                                          |
|--------------|----------------------------------------------------------------|
| **Summary**  | Destroys the `BondManager` object.                               |
| **Parameters** | none                                                            |
| **Returns**  | `void`                                                         |
| **Behaviour** | Calls the base class destructor, `SelManager()`. |
| **Side effects** | none                                                          |
| **Errors / exceptions** | none                                                          |

#### `BondManager::MakeBonds(bool calc_only)`

| Field        | Detail                                                                                                                                                           |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Creates chemical bonds between atoms within a macromolecular model, based on the defined connectivity rules.                                                          |
| **Parameters** | `calc_only` — A boolean flag indicating whether to calculate bonds (`true`) or only remove existing bonds (`false`).  If `true`, calculates new bonds. |
| **Returns**  | `void`                                                                                                                                                           |
| **Behaviour** | 1. Removes existing bonds using `RemoveBonds()`. 2. Iterates through each model in the database. 3. For each model, iterates through each chain. 4. For each chain, iterates through each residue. 5. For each residue, creates a graph of its atoms using `math::Graph::MakeGraph()`. 6. Gets the number of vertices and edges of the graph. 7. Iterates through each edge in the graph. 8. For each edge, adds a bond between the connected atoms by calling `res->atom[k1]->AddBond(res->atom[k2],E[i]->GetType())` and the reverse bond `res->atom[k2]->AddBond(res->atom[k1],E[i]->GetType())`. |
| **Side effects** |  Mutates the model data structure by adding bonds between atoms. `RemoveBonds()` has side effects. |
| **Errors / exceptions** | none                                                          |

#### `BondManager::RemoveBonds()`

| Field        | Detail                                                                                             |
|--------------|-----------------------------------------------------------------------------------------------------|
| **Summary**  | Removes all existing chemical bonds from the macromolecular model.                                 |
| **Parameters** | none                                                                                              |
| **Returns**  | `void`                                                                                            |
| **Behaviour** | Iterates through all atoms in the model. For each atom, calls `atom[i]->FreeBonds()` which removes the bonds connected to that atom. |
| **Side effects** |  Removes bonds between atoms from the model data structure.  |
| **Errors / exceptions** | none                                                          |

#### `BondManager::write(io::RFile f)`

| Field        | Detail                                                                                                                                                           |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Serializes the `BondManager` object's data to a file stream `f`.                                                                                                    |
| **Parameters** | `f` — An `io::RFile` object used for writing the data. |
| **Returns**  | `void`                                                                                                                                                           |
| **Behaviour** | Writes the version number of the object to the file stream using `f.WriteByte(&Version)`. Then calls the base class `write` method `SelManager::write(f)`. |
| **Side effects** | Writes data to the specified file stream. |
| **Errors / exceptions** | none                                                          |

#### `BondManager::read(io::RFile f)`

| Field        | Detail                                                                                                                                                           |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Deserializes the `BondManager` object's data from a file stream `f`.                                                                                              |
| **Parameters** | `f` — An `io::RFile` object used for reading the data. |
| **Returns**  | `void`                                                                                                                                                           |
| **Behaviour** | Reads the version number of the object from the file stream using `f.ReadByte(&Version)`. Then calls the base class `read` method `SelManager::read(f)`. |
| **Side effects** | Reads data from the specified file stream. |
| **Errors / exceptions** | none                                                          |

---
## `mmdb::BondManager`

> This class is responsible for creating and removing chemical bonds within a macromolecular model based on atomic coordinates.

### Methods

#### `mmdb::BondManager::BondManager()`

| Field        | Detail                                                           |
|--------------|------------------------------------------------------------------|
| **Summary**  | Constructs a `BondManager` object.                               |
| **Parameters** | none                                                             |
| **Returns**  | `void`                                                          |
| **Behaviour** | Initializes the `BondManager` object by calling the base class constructor, `SelManager()`. |
| **Side effects** | none                                                          |
| **Errors / exceptions** | none                                                          |

#### `mmdb::BondManager::BondManager(io::RPStream Object)`

| Field        | Detail                                                            |
|--------------|------------------------------------------------------------------|
| **Summary**  | Constructs a `BondManager` object, initialized from an `io::RPStream`. |
| **Parameters** | `Object` — An `io::RPStream` object used for stream-based object initialization. |
| **Returns**  | `void`                                                          |
| **Behaviour** | Initializes the `BondManager` object by calling the base class constructor, `SelManager(Object)`. |
| **Side effects** | none                                                          |
| **Errors / exceptions** | none                                                          |

#### `mmdb::BondManager::~BondManager()`

| Field        | Detail                                                          |
|--------------|----------------------------------------------------------------|
| **Summary**  | Destroys the `BondManager` object.                               |
| **Parameters** | none                                                            |
| **Returns**  | `void`                                                         |
| **Behaviour** | Calls the base class destructor, `SelManager()`. |
| **Side effects** | none                                                          |
| **Errors / exceptions** | none                                                          |

#### `mmdb::BondManager::MakeBonds(bool calc_only)`

| Field        | Detail                                                                                                                                                           |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Creates chemical bonds between atoms within a macromolecular model, based on the defined connectivity rules.                                                          |
| **Parameters** | `calc_only` — A boolean flag indicating whether to calculate bonds (`true`) or only remove existing bonds (`false`).  If `true`, calculates new bonds. |
| **Returns**  | `void`                                                                                                                                                           |
| **Behaviour** | 1. Removes existing bonds using `RemoveBonds()`. 2. Iterates through each model in the database. 3. For each model, iterates through each chain. 4. For each chain, iterates through each residue. 5. For each residue, creates a graph of its atoms using `math::Graph::MakeGraph()`. 6. Gets the number of vertices and edges of the graph. 7. Iterates through each edge in the graph. 8. For each edge, adds a bond between the connected atoms by calling `res->atom[k1]->AddBond(res->atom[k2],E[i]->GetType())` and the reverse bond `res->atom[k2]->AddBond(res->atom[k1],E[i]->GetType())`. |
| **Side effects** |  Mutates the model data structure by adding bonds between atoms. `RemoveBonds()` has side effects. |
| **Errors / exceptions** | none                                                          |

#### `mmdb::BondManager::RemoveBonds()`

| Field        | Detail                                                                                             |
|--------------|-----------------------------------------------------------------------------------------------------|
| **Summary**  | Removes all existing chemical bonds from the macromolecular model.                                 |
| **Parameters** | none                                                                                              |
| **Returns**  | `void`                                                                                            |
| **Behaviour** | Iterates through all atoms in the model. For each atom, calls `atom[i]->FreeBonds()` which removes the bonds connected to that atom. |
| **Side effects** |  Removes bonds between atoms from the model data structure.  |
| **Errors / exceptions** | none                                                          |

#### `mmdb::BondManager::write(io::RFile f)`

| Field        | Detail                                                                                                                                                           |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Serializes the `BondManager` object's data to a file stream `f`.                                                                                                    |
| **Parameters** | `f` — An `io::RFile` object used for writing the data. |
| **Returns**  | `void`                                                                                                                                                           |
| **Behaviour** | Writes the version number of the object to the file stream using `f.WriteByte(&Version)`. Then calls the base class `write` method `SelManager::write(f)`. |
| **Side effects** | Writes data to the specified file stream. |
| **Errors / exceptions** | none                                                          |

#### `mmdb::BondManager::read(io::RFile f)`

| Field        | Detail                                                                                                                                                           |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Deserializes the `BondManager` object's data from a file stream `f`.                                                                                              |
| **Parameters** | `f` — An `io::RFile` object used for reading the data. |
| **Returns**  | `void`                                                                                                                                                           |
| **Behaviour** | Reads the version number of the object from the file stream using `f.ReadByte(&Version)`. Then calls the base class `read` method `SelManager::read(f)`. |
| **Side effects** | Reads data from the specified file stream. |
| **Errors / exceptions** | none                                                          |

---
## `mmdb::Free Functions`

> These are free functions that are not associated with a class.

---

### `mmdb_chain`

Okay, I've processed the C++ source code (`mmdb_chain.h` and `mmdb_chain.cpp`) and produced detailed, developer-facing API reference documentation.  I've followed the exact requested format, including the `ClassName` heading, the `Methods` section with `Summary`, `Parameters`, `Returns`, `Behaviour`, `Side effects`, and `Errors / exceptions` for each method. I've also included constructors and destructors.

Here's the complete documentation output:

---
## `mmdb::ProModel`

> This class is a virtue needed only for defining certain functions of Model, which are used by Chain and Residue.

### Methods

#### `int  ProModel::AddChain (PChain)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Adds a chain to the ProModel.                                                                             |
| **Parameters** | `PChain Chain_Owner` — A pointer to the Chain object to be added.                                                |
| **Returns**  | `0` if successful, otherwise an error code (not defined in this source).                                      |
| **Behaviour** |  This method is used internally by the `ProModel` class to add Chain objects to the object.              |
| **Side effects** | Adds the given chain to the ProModel.                                                                        |
| **Errors / exceptions** | Returns an error code (not defined in this source) if chain is not added.                               |

#### `int  ProModel::ExpandAtomArray (int inc)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Expands the array of atoms in the ProModel if the existing array is full.                                |
| **Parameters** | `int inc` — The number of atoms to add to the array.                                                        |
| **Returns**  | `0` if successful, otherwise an error code (not defined in this source).                                      |
| **Behaviour** | Increases the size of the `atom` array by `inc`. If the array is already full, it is resized.                 |
| **Side effects** | Changes the size of the `atom` array, reallocating memory.                                                    |
| **Errors / exceptions** | Returns an error code (not defined in this source) if expansion fails.                               |

#### `int  ProModel::GetAllAtoms ()`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns a pointer to the array of all atoms in the ProModel.                                                 |
| **Parameters** | none                                                                                                      |
| **Returns**  | `NULL` if there are no atoms or if the array is not yet populated.                                           |
| **Behaviour** | Returns the pointer to the `atom` array.                                                                    |
| **Side effects** | none                                                                                                      |
| **Errors / exceptions** | `NULL` is returned if no atoms are present.                                                              |

#### `int  ProModel::GetNumberOfAllAtoms ()`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the total number of atoms in all models.                                                           |
| **Parameters** | none                                                                                                      |
| **Returns**  | `0`                                                                                                         |
| **Behaviour** | Returns the count of atoms.                                                                                 |
| **Side effects** | none                                                                                                      |
| **Errors / exceptions** | `0` is always returned.                                                                               |

#### `PManager ProModel::GetCoordHierarchy ()`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns a pointer to the general Atom array.                                                                     |
| **Parameters** | none                                                                                                      |
| **Returns**  | `NULL`                                                                                                      |
| **Behaviour** | Returns a pointer to the PManager class.                                                                   |
| **Side effects** | none                                                                                                      |
| **Errors / exceptions** | Returns `NULL` if there is no hierarchy.                                                              |

#### `int  ProModel::GetNumberOfModels ()`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the total number of models.                                                                          |
| **Parameters** | none                                                                                                      |
| **Returns**  | `0`                                                                                                         |
| **Behaviour** | Returns the count of models.                                                                                |
| **Side effects** | none                                                                                                      |
| **Errors / exceptions** | `0` is always returned.                                                                               |

#### `int  ProModel::GetNumberOfResidues ()`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the number of residues in the chain.                                                              |
| **Parameters** | none                                                                                                      |
| **Returns**  | `0`                                                                                                         |
| **Behaviour** | Returns the count of residues.                                                                             |
| **Side effects** | none                                                                                                      |
| **Errors / exceptions** | `0` is always returned.                                                                               |

#### `PResidue ProModel::GetResidue (int resNo)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns a pointer to a residue in the chain.                                                               |
| **Parameters** | `int resNo` — The index of the residue to retrieve.                                                        |
| **Returns**  | `NULL` if the residue does not exist or if the index is out of range.                                        |
| **Behaviour** | Returns a pointer to the `residue` at the specified index.                                                   |
| **Side effects** | none                                                                                                      |
| **Errors / exceptions** | `NULL` is returned if residue is not found or if the index is out of range.                           |

#### `PResidue ProModel::GetResidueCreate (const ResName resName, int seqNum, const InsCode insCode, bool Enforce)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns a pointer to a new residue.                                                                           |
| **Parameters** | `const ResName resName` — The name of the residue.                                                         |
|                | `int seqNum` — The sequence number of the residue.                                                        |
|                | `const InsCode insCode` — The insertion code of the residue.                                               |
|                | `bool Enforce` —  Indicates if existing residue is to be replaced                                                       |
| **Returns**  | `NULL` if the residue could not be created.                                                                 |
| **Behaviour** | Creates a new residue, sets the name, sequence number, and insertion code, and returns a pointer to the residue.  |
| **Side effects** | Creates a new residue object, adds it to the chain.                                                           |
| **Errors / exceptions** | `NULL` is returned if there is a failure.                                                             |

#### `void ProModel::_ExcludeResidue (const ResName resName, int seqNum, const InsCode insCode)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Excludes a residue from the chain.                                                                         |
| **Parameters** | `const ResName resName` — The name of the residue to exclude.                                               |
|                | `int seqNum` — The sequence number of the residue to exclude.                                             |
|                | `const InsCode insCode` — The insertion code of the residue to exclude.                                        |
| **Returns**  | `0` if successful, otherwise an error code (not defined in this source).                                      |
| **Behaviour** | Removes the specified residue from the chain.                                                              |
| **Side effects** | Removes the specified residue from the chain.                                                             |
| **Errors / exceptions** | Returns an error code (not defined in this source) if the exclusion fails.                             |

#### `void ProModel::FreeMemory()`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Frees all memory allocated for the ProModel.                                                               |
| **Parameters** | none                                                                                                      |
| **Returns**  | `void`                                                                                                     |
| **Behaviour** | Releases the `residue` array and the `atom` array.                                                            |
| **Side effects** | Deallocates the memory used by the ProModel.                                                              |
| **Errors / exceptions** | none                                                                                                      |

#### `void ProModel::FreeAnnotations()`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Frees all memory allocated for annotations (DBReference, SeqAdv, ModRes, Het).                               |
| **Parameters** | none                                                                                                      |
| **Returns**  | `void`                                                                                                     |
| **Behaviour** | Releases the memory used by the annotation objects.                                                         |
| **Side effects** | Deallocates the memory used by the DBRef, SeqAdv, ModRes, and Het objects.                                 |
| **Errors / exceptions** | none                                                                                                      |

#### `void ProModel::GetAtomStatistics (RAtomStat AS)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Populates the provided `RAtomStat` structure with statistics about the atoms in the ProModel.                |
| **Parameters** | `RAtomStat AS` — A pointer to the `RAtomStat` structure to populate.                                         |
| **Returns**  | `void`                                                                                                     |
| **Behaviour** | Populates the RAtomStat with information such as atom counts, element counts, etc.                          |
| **Side effects** | Modifies the `RAtomStat` structure.                                                                      |
| **Errors / exceptions** | none                                                                                                      |

#### `void ProModel::CalAtomStatistics (RAtomStat AS)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Calculates the statistics of the atoms in the ProModel.                                                      |
| **Parameters** | `RAtomStat AS` — A pointer to the `RAtomStat` structure to populate.                                         |
| **Returns**  | `void`                                                                                                     |
| **Behaviour** | Calculates the statistics of the atoms in the ProModel and then populates the `RAtomStat` structure.        |
| **Side effects** | Modifies the `RAtomStat` structure.                                                                      |
| **Errors / exceptions** | none                                                                                                      |

#### `void ProModel::CheckInAtoms()`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Checks all atoms in the ProModel.                                                                            |
| **Parameters** | none                                                                                                      |
| **Returns**  | `void`                                                                                                     |
| **Behaviour** | Checks all atoms in the ProModel.                                                                            |
| **Side effects** | none                                                                                                      |
| **Errors / exceptions** | none                                                                                                      |

#### `void ProModel::ApplyTransform (mat44 & TMatrix)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Applies a transformation matrix to the coordinates in the ProModel.                                           |
| **Parameters** | `mat44 & TMatrix` — A pointer to the transformation matrix to apply.                                          |
| **Returns**  | `void`                                                                                                     |
| **Behaviour** | Multiplies the coordinates in the ProModel with the given matrix.                                              |
| **Side effects** | Modifies the coordinates in the ProModel.                                                                |
| **Errors / exceptions** | none                                                                                                      |

#### `int ProModel::GetNumberOfDBRefs ()`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the number of database references in the ProModel.                                                 |
| **Parameters** | none                                                                                                      |
| **Returns**  | `0`                                                                                                         |
| **Behaviour** | Returns the count of database references.                                                                 |
| **Side effects** | none                                                                                                      |
| **Errors / exceptions** | `0` is always returned.                                                                               |

#### `PDBReference ProModel::GetDBRef (int dbRefNo)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns a pointer to the database reference object at the given index.                                    |
| **Parameters** | `int dbRefNo` — The index of the database reference to retrieve.                                           |
| **Returns**  | `NULL` if the database reference does not exist or if the index is out of range.                          |
| **Behaviour** | Returns a pointer to the `DBReference` object at the specified index.                                      |
| **Side effects** | none                                                                                                      |
| **Errors / exceptions** | `NULL` is returned if the database reference is not found or if the index is out of range.            |

#### `void ProModel::Copy (PChain chain)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Copies the chain to a new ProModel.                                                                       |
| **Parameters** | `PChain chain` — A pointer to the Chain object to copy.                                                       |
| **Returns**  | `void`                                                                                                     |
| **Behaviour** | Copies the data from the original chain to a new ProModel instance.                                          |
| **Side effects** | Copies the data from the chain.                                                                           |
| **Errors / exceptions** | none                                                                                                      |

#### `void ProModel::CopyAnnotations (PChain chain)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Copies the annotations (DBReference, SeqAdv, ModRes, Het) from a chain to the ProModel.                        |
| **Parameters** | `PChain chain` — A pointer to the Chain object from which to copy the annotations.                        |
| **Returns**  | `void`                                                                                                     |
| **Behaviour** | Copies the corresponding annotation objects from the input chain to the ProModel.                           |
| **Side effects** | Modifies the ProModel annotation objects.                                                               |
| **Errors / exceptions** | none                                                                                                      |

---
## `mmdb::Chain`

> The Chain class represents a chain in a protein structure.

### Methods

#### `void Chain::AddAtom (int seqNo, const InsCode insCode, PAtom atom)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Adds an atom to the chain.                                                                                  |
| **Parameters** | `int seqNo` — The sequence number of the atom.                                                             |
|                | `const InsCode insCode` — The insertion code of the atom.                                                    |
|                | `PAtom atom` — A pointer to the atom object to add.                                                        |
| **Returns**  | `0` if successful, otherwise an error code (not defined in this source).                                      |
| **Behaviour** | Adds the given atom to the chain.                                                                          |
| **Side effects** | Adds the specified atom to the chain.                                                                     |
| **Errors / exceptions** | Returns an error code (not defined in this source) if the atom is not added.                           |

#### `void Chain::AddAtom (int resNo, PAtom atom)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Adds an atom to the chain, inserting it at a specific residue.                                               |
| **Parameters** | `int resNo` — The index of the residue to insert the atom before.                                           |
|                | `PAtom atom` — A pointer to the atom object to add.                                                        |
| **Returns**  | `0` if successful, otherwise an error code (not defined in this source).                                      |
| **Behaviour** | Inserts the atom at the specified position in the chain.                                                   |
| **Side effects** | Inserts the atom at the specified position in the chain.                                                 |
| **Errors / exceptions** | Returns an error code (not defined in this source) if the atom is not added.                           |

#### `void Chain::AddResidue (PResidue res)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Adds a residue to the chain.                                                                                |
| **Parameters** | `PResidue res` — A pointer to the residue object to add.                                                    |
| **Returns**  | `0` if successful, otherwise an error code (not defined in this source).                                      |
| **Behaviour** | Adds the given residue to the chain.                                                                       |
| **Side effects** | Adds the specified residue to the chain.                                                                 |
| **Errors / exceptions** | Returns an error code (not defined in this source) if the residue is not added.                       |

#### `int Chain::DeleteResidue (int resNo)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Deletes a residue from the chain.                                                                            |
| **Parameters** | `int resNo` — The index of the residue to delete.                                                         |
| **Returns**  | `1` if successful, `0` otherwise.                                                                         |
| **Behaviour** | Removes the specified residue from the chain.                                                              |
| **Side effects** | Removes the specified residue from the chain.                                                             |
| **Errors / exceptions** | Returns `0` if the residue does not exist or if the index is out of range.                           |

#### `int Chain::DeleteResidue (int seqNum, const InsCode insCode)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Deletes a residue from the chain, specified by sequence number and insertion code.                         |
| **Parameters** | `int seqNum` — The sequence number of the residue to delete.                                             |
|                | `const InsCode insCode` — The insertion code of the residue to delete.                                        |
| **Returns**  | `1` if successful, `0` otherwise.                                                                         |
| **Behaviour** | Removes the specified residue from the chain.                                                              |
| **Side effects** | Removes the specified residue from the chain.                                                             |
| **Errors / exceptions** | Returns `0` if the residue does not exist or if the sequence number and insertion code do not match. |

#### `int Chain::DeleteAllResidues ()`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Deletes all residues from the chain.                                                                        |
| **Parameters** | none                                                                                                      |
| **Returns**  | `0` if successful, otherwise an error code (not defined in this source).                                      |
| **Behaviour** | Removes all residues from the chain.                                                                        |
| **Side effects** | Removes all residues from the chain.                                                                      |
| **Errors / exceptions** | Returns an error code (not defined in this source) if the deletion fails.                            |

#### `int Chain::DeleteSolvent ()`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Deletes all solvent residues from the chain.                                                              |
| **Parameters** | none                                                                                                      |
| **Returns**  | `0` if successful, otherwise an error code (not defined in this source).                                      |
| **Behaviour** | Removes all residue objects that are solvent residues.                                                     |
| **Side effects** | Removes all solvent residues from the chain.                                                              |
| **Errors / exceptions** | Returns an error code (not defined in this source) if the deletion fails.                            |

#### `void Chain::TrimResidueTable ()`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Removes residues from the end of the chain.                                                                  |
| **Parameters** | none                                                                                                      |
| **Returns**  | `void`                                                                                                     |
| **Behaviour** | Removes residues from the end of the chain, leaving only residues that are present in the coordinate file. |
| **Side effects** | Modifies the `residue` array and `nResidues`.                                                              |
| **Errors / exceptions** | none                                                                                                      |

#### `int Chain::AddAtom (int seqNo, const InsCode insCode, PAtom atom)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Adds an atom to the chain.                                                                                  |
| **Parameters** | `int seqNo` — The sequence number of the atom.                                                             |
|                | `const InsCode insCode` — The insertion code of the atom.                                                    |
|                | `PAtom atom` — A pointer to the atom object to add.                                                        |
| **Returns**  | `0` if successful, otherwise an error code (not defined in this source).                                      |
| **Behaviour** | Adds the given atom to the chain.                                                                          |
| **Side effects** | Adds the specified atom to the chain.                                                                     |
| **Errors / exceptions** | Returns an error code (not defined in this source) if the atom is not added.                           |

#### `int Chain::AddAtom (int resNo, PAtom atom)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Adds an atom to the chain, inserting it at a specific residue.                                               |
| **Parameters** | `int resNo` — The index of the residue to insert the atom before.                                           |
|                | `PAtom atom` — A pointer to the atom object to add.                                                        |
| **Returns**  | `0` if successful, otherwise an error code (not defined in this source).                                      |
| **Behaviour** | Inserts the atom at the specified position in the chain.                                                   |
| **Side effects** | Inserts the atom at the specified position in the chain.                                                 |
| **Errors / exceptions** | Returns an error code (not defined in this source) if the atom is not added.                           |

#### `int Chain::DeleteAtom (int seqNo, const InsCode insCode, const AtomName aname, const Element elmnt, const AltLoc aloc)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Deletes an atom from the chain.                                                                              |
| **Parameters** | `int seqNo` — The sequence number of the atom.                                                             |
|                | `const InsCode insCode` — The insertion code of the atom.                                                    |
|                | `const AtomName aname` — The name of the atom.                                                              |
|                | `const Element elmnt` — The element of the atom.                                                             |
|                | `const AltLoc aloc` — The alternative location of the atom.                                                 |
| **Returns**  | `0` if successful, otherwise an error code (not defined in this source).                                      |
| **Behaviour** | Removes the specified atom from the chain.                                                                 |
| **Side effects** | Removes the specified atom from the chain.                                                               |
| **Errors / exceptions** | Returns an error code (not defined in this source) if the atom is not found or if the parameters do not match.|

#### `int Chain::DeleteAtom (int seqNo, const InsCode insCode, int atomNo)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Deletes an atom from the chain, specified by sequence number and atom number.                             |
| **Parameters** | `int seqNo` — The sequence number of the atom.                                                             |
|                | `const InsCode insCode` — The insertion code of the atom.                                                    |
|                | `int atomNo` — The index of the atom to delete.                                                              |
| **Returns**  | `0` if successful, otherwise an error code (not defined in this source).                                      |
| **Behaviour** | Removes the specified atom from the chain.                                                              |
| **Side effects** | Removes the specified atom from the chain.                                                              |
| **Errors / exceptions** | Returns an error code (not defined in this source) if the atom does not exist or if the parameters do not match.|

#### `int Chain::DeleteAllAtoms (int seqNo, const InsCode insCode)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Deletes all atoms from the chain.                                                                          |
| **Parameters** | `int seqNo` — The sequence number of the atom.                                                             |
|                | `const InsCode insCode` — The insertion code of the atom.                                                    |
| **Returns**  | `0` if successful, otherwise an error code (not defined in this source).                                      |
| **Behaviour** | Removes all atoms from the chain.                                                                         |
| **Side effects** | Removes all atoms from the chain.                                                                       |
| **Errors / exceptions** | Returns an error code (not defined in this source) if the deletion fails.                            |

#### `int Chain::DeleteAllAtoms ()`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Deletes all atoms from the chain.                                                                          |
| **Parameters** | none                                                                                                      |
| **Returns**  | `0` if successful, otherwise an error code (not defined in this source).                                      |
| **Behaviour** | Removes all atoms from the chain.                                                                         |
| **Side effects** | Removes all atoms from the chain.                                                                       |
| **Errors / exceptions** | Returns an error code (not defined in this source) if the deletion fails.                            |

#### `int Chain::DeleteAltLocs ()`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Removes all alternative locations from the chain.                                                           |
| **Parameters** | none                                                                                                      |
| **Returns**  | `0` if successful, otherwise an error code (not defined in this source).                                      |
| **Behaviour** | Removes all alternative location objects from the chain, leaving only the highest occupancy.                  |
| **Side effects** | Modifies the chain by removing the alternative location objects.                                          |
| **Errors / exceptions** | Returns an error code (not defined in this source) if the deletion fails.                            |

#### `int Chain::AddAtom (int resNo, PAtom atom)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Adds an atom to the chain, inserting it at a specific residue.                                               |
| **Parameters** | `int resNo` — The index of the residue to insert the atom before.                                           |
|                | `PAtom atom` — A pointer to the atom object to add.                                                        |
| **Returns**  | `0` if successful, otherwise an error code (not defined in this source).                                      |
| **Behaviour** | Inserts the atom at the specified position in the chain.                                                   |
| **Side effects** | Inserts the atom at the specified position in the chain.                                                 |
| **Errors / exceptions** | Returns an error code (not defined in this source) if the atom is not added.                           |

#### `int Chain::InsResidue (PResidue res, int seqNum, const InsCode insCode)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Inserts a residue into the chain at a specific position.                                                      |
| **Parameters** | `PResidue res` — A pointer to the residue object to insert.                                                  |
|                | `int seqNum` — The sequence number of the residue.                                                        |
|                | `const InsCode insCode` — The insertion code of the residue.                                               |
| **Returns**  | `0` if successful, otherwise an error code (not defined in this source).                                      |
| **Behaviour** | Inserts the specified residue into the chain at the specified position.                                     |
| **Side effects** | Inserts the specified residue into the chain.                                                             |
| **Errors / exceptions** | Returns an error code (not defined in this source) if the insertion fails.                            |

#### `int Chain::GetNumberOfAtoms (bool countTers)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the number of atoms in the chain.                                                                   |
| **Parameters** | `bool countTers` —  If true returns total number of atoms, including those that are solvent residues.         |
| **Returns**  | `0`                                                                                                         |
| **Behaviour** | Returns the count of atoms.                                                                             |
| **Side effects** | none                                                                                                      |
| **Errors / exceptions** | `0` is always returned.                                                                               |

#### `int Chain::GetNumberOfAtoms (int seqNo, const InsCode insCode)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the number of atoms in the chain, specified by sequence number and insertion code.                |
| **Parameters** | `int seqNo` — The sequence number of the atom.                                                             |
|                | `const InsCode insCode` — The insertion code of the atom.                                                    |
| **Returns**  | `0`                                                                                                         |
| **Behaviour** | Returns the count of atoms that match the specified sequence number and insertion code.                   |
| **Side effects** | none                                                                                                      |
| **Errors / exceptions** | `0` is always returned.                                                                               |

#### `int Chain::GetNumberOfAtoms (int resNo)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the number of atoms in the chain, specified by residue number.                                     |
| **Parameters** | `int resNo` — The index of the residue.                                                                    |
| **Returns**  | `0`                                                                                                         |
| **Behaviour** | Returns the count of atoms that are located at the specified residue.                                        |
| **Side effects** | none                                                                                                      |
| **Errors / exceptions** | `0` is always returned.                                                                               |

#### `PAtom Chain::GetAtom (int seqNo, const InsCode insCode, const AtomName aname, const Element elmnt, const AltLoc aloc)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns a pointer to an atom in the chain.                                                                   |
| **Parameters** | `int seqNo` — The sequence number of the atom.                                                             |
|                | `const InsCode insCode` — The insertion code of the atom.                                                    |
|                | `const AtomName aname` — The name of the atom.                                                              |
|                | `const Element elmnt` — The element of the atom.                                                             |
|                | `const AltLoc aloc` — The alternative location of the atom.                                                 |
| **Returns**  | `NULL` if the atom does not exist or if the parameters do not match.                                       |
| **Behaviour** | Returns a pointer to the `atom` object that matches the specified criteria.                                  |
| **Side effects** | none                                                                                                      |
| **Errors / exceptions** | `NULL` is returned if the atom is not found or if the parameters do not match.                           |

#### `PAtom Chain::GetAtom (int seqNo, const InsCode insCode, int atomNo)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns a pointer to an atom in the chain, specified by sequence number and atom number.                   |
| **Parameters** | `int seqNo` — The sequence number of the atom.                                                             |
|                | `const InsCode insCode` — The insertion code of the atom.                                                    |
|                | `int atomNo` — The index of the atom to retrieve.                                                           |
| **Returns**  | `NULL` if the atom does not exist or if the index is out of range.                                        |
| **Behaviour** | Returns a pointer to the `atom` object at the specified index.                                              |
| **Side effects** | none                                                                                                      |
| **Errors / exceptions** | `NULL` is returned if the atom does not exist or if the index is out of range.                           |

#### `PAtom Chain::GetAtom (int resNo, const AtomName aname, const Element elmnt, const AltLoc aloc)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns a pointer to an atom in the chain, specified by residue number and other atom properties.         |
| **Parameters** | `int resNo` — The index of the residue.                                                                    |
|                | `const AtomName aname` — The name of the atom.                                                              |
|                | `const Element elmnt` — The element of the atom.                                                             |
|                | `const AltLoc aloc` — The alternative location of the atom.                                                 |
| **Returns**  | `NULL` if the atom does not exist or if the parameters do not match.                                       |
| **Behaviour** | Returns a pointer to the `atom` object that matches the specified criteria.                                  |
| **Side effects** | none                                                                                                      |
| **Errors / exceptions** | `NULL` is returned if the atom is not found or if the parameters do not match.                           |

#### `PAtom Chain::GetAtom (int resNo, int atomNo)`

| Field        | Detail                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns a pointer to an atom in the chain, specified by residue number and atom number.                   |
| **Parameters** | `int resNo` — The index of the residue.                                                                    |
|                | `int atomNo` — The index of the atom to retrieve.                                                           |
| **Returns**  | `NULL` if the atom does not

---

### `mmdb_cifdefs`


## `mmdb::CIFName`

> Returns the CIF name based on the CIF mode and the specified ID.

### Methods

#### `cpstr CIFName ( int NameID, CIF_MODE Mode )`

| Field        | Detail                 |
|--------------|------------------------|
| **Summary**  | Returns the CIF name based on the CIF mode and the specified ID. |
| **Parameters** | `NameID` — An integer representing the CIF tag ID.  |
|               | `Mode` — An enumeration value indicating the CIF mode (NDB or PDBX). |
| **Returns**  | `cpstr` — The corresponding CIF name. |
| **Behaviour** | The method takes a CIF tag ID and a CIF mode as input.  It then uses a `switch` statement to determine the appropriate CIF name based on the specified mode and tag ID.  The names are statically defined within the class. |
| **Side effects** | None                     |
| **Errors / exceptions** | None                     |


---

### `mmdb_coormngr`

## `Brick`

> Represents a space brick used in the coordinate transformation process.

### Methods

#### `Brick()`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructor for the `Brick` class. Initializes the brick's parameters, allocating memory for the atoms. |
| **Parameters** | None                                                                                                   |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Allocates memory for `atom` (an array of `PAtom` pointers), `id` (an array of integer indices), and `nAtoms`. Initializes the allocated memory to NULL. Initializes the number of atoms. |
| **Side effects** | Allocates memory for the `atom` and `id` arrays.                                                              |
| **Errors / exceptions** | None                                                                                                  |

#### `~Brick()`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Destructor for the `Brick` class.  Releases the memory allocated for the `atom` and `id` arrays.               |
| **Parameters** | None                                                                                                   |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Deallocates the `atom` array by calling `delete[] atom;`. Deallocates the `id` array by calling `FreeVectorMemory ( id,0 );`. Sets `atom`, `nAtoms`, and `nAllocAtoms` to NULL.  |
| **Side effects** | Deallocates the memory used for the `atom` and `id` arrays.                                             |
| **Errors / exceptions** | None                                                                                                  |

#### `Clear()`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Clears the Brick object by releasing the allocated memory and resetting all member variables to their initial state. |
| **Parameters** | None                                                                                                   |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Calls `Clear()` to clear the atom and id arrays.                                                         |
| **Side effects** | Deallocates the `atom` and `id` arrays.  Sets `nAtoms` and `nAllocAtoms` to 0.                               |
| **Errors / exceptions** | None                                                                                                  |

#### `AddAtom(PAtom A, int atomid)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Adds an atom to the brick.                                                                           |
| **Parameters** | `A` — `PAtom` pointer to the atom to be added. `atomid` — Integer representing the atom's ID within the brick. |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Checks if the allocated memory for the `atom` array is sufficient.  If not, it expands the array. Then, it adds the `A` atom to the end of the `atom` array, and sets the corresponding ID in the `id` array to `atomid`. Increments the `nAtoms` counter. |
| **Side effects** | Expands the `atom` array if needed, allocating memory for a larger array and copying existing atoms to the new array. |
| **Errors / exceptions** | None                                                                                                  |

---

## `MBrick`

> Represents a multi-brick that can contain multiple structures, each with its own brick.

### Methods

#### `MBrick(int nStructures)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructor for the `MBrick` class. Initializes the multi-brick's parameters, allocating memory for the bricks and their corresponding atom indices. |
| **Parameters** | `nStructures` — Integer indicating the number of structures this multi-brick will contain.          |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Initializes `nStruct` to `nStructures`. Allocates memory for `atom` (an array of `PPAtom` pointers) and `id` (an array of `ivector` indices).  Initializes the size of each array to `nStruct`. |
| **Side effects** | Allocates memory for the `atom` and `id` arrays.                                                              |
| **Errors / exceptions** | None                                                                                                  |

#### `~MBrick()`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Destructor for the `MBrick` class. Releases the memory allocated for the `atom` and `id` arrays.               |
| **Parameters** | None                                                                                                   |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Deallocates the `atom` array by calling `delete[] atom;`. Deallocates the `id` array by calling `FreeMatrixMemory ( id,nStruct,0,0 );`. Sets `nStruct` to 0.  |
| **Side effects** | Deallocates the memory used for the `atom` and `id` arrays.                                             |
| **Errors / exceptions** | None                                                                                                  |

#### `Clear()`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Clears the MBrick object by releasing the allocated memory and resetting all member variables.              |
| **Parameters** | None                                                                                                   |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | If `atom` exists, it is deallocated. `FreeMatrixMemory` is called to deallocate the `id` array. Sets `nStruct` to 0. |
| **Side effects** | Deallocates the memory used for the `atom` and `id` arrays.                                             |
| **Errors / exceptions** | None                                                                                                  |

#### `AddAtom(PAtom A, int structNo, int atomid)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Adds an atom to the multi-brick, specifically assigning it to a particular structure.                    |
| **Parameters** | `A` — `PAtom` pointer to the atom to be added. `structNo` — Integer representing the structure number to which the atom belongs. `atomid` — Integer representing the atom's ID within the structure. |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Retrieves the current size of the `nAtoms` array for the specified structure.  If the current size is less than the allocated size, it expands the array.  Adds `A` to the end of the `atom` array and sets the corresponding `atomid` in the `id` array. Increments `nAtoms` by 1.  |
| **Side effects** | Expands the `atom` array if needed and allocates memory for a larger array. |
| **Errors / exceptions** | None                                                                                                  |

---

## `GenSym`

> Represents a symmetry operation for generating symmetric models.

### Methods

#### `GenSym()`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructor for the `GenSym` class. Initializes the class. |
| **Parameters** | None                                                                                                   |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Initializes member variables, including `chID1` (an array of `ChainID` pointers), `chID2` (another array of `ChainID` pointers), `nChains` (an array to store the number of chains associated with each operation), and `nOpAlloc` (the allocated size of the `chID1` and `chID2` arrays). |
| **Side effects** | Initializes member variables.                                                                       |
| **Errors / exceptions** | None                                                                                                  |

#### `GenSym(io::RPStream Object)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructor for the `GenSym` class. Initializes the class and sets up the output stream.                    |
| **Parameters** | `Object` — An `io::RPStream` object to handle output operations.                                     |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Calls the default constructor (`GenSym()`) to initialize the class and sets up the output stream.         |
| **Side effects** | Initializes the `io::RPStream` object.                                                                |
| **Errors / exceptions** | None                                                                                                  |

#### `~GenSym()`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Destructor for the `GenSym` class. Releases the memory allocated for the arrays of `ChainID`.            |
| **Parameters** | None                                                                                                   |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Deallocates the memory for `chID1` and `chID2`.  Calls `FreeVectorMemory` to free the memory. Sets `nOpAlloc` to 0. |
| **Side effects** | Deallocates the memory used for the `chID1` and `chID2` arrays.                                       |
| **Errors / exceptions** | None                                                                                                  |

#### `FreeMemory()`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Releases memory allocated for `chID1`, `chID2`, and `nChains`.                                             |
| **Parameters** | None                                                                                                   |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Iterates through the `chID1` and `chID2` arrays, deleting each `ChainID` array.  Calls `FreeVectorMemory` to deallocate the `nChains` array. Sets `nOpAlloc` to 0.  Calls `SymOps::FreeMemory()` to perform any additional cleanup tasks required by the base class. |
| **Side effects** | Deallocates the memory used for the `chID1`, `chID2`, and `nChains` arrays.                               |
| **Errors / exceptions** | None                                                                                                  |

#### `AddSymOp(cpstr XYZOperation)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Adds a symmetry operation to the `GenSym` object.                                                      |
| **Parameters** | `XYZOperation` — A `cpstr` pointer to the string representing the symmetry operation.                     |
| **Returns**  | `int` — Returns `SPOSEAT_Ok` if the operation was added successfully, or an error code otherwise.          |
| **Behaviour** | Increments the number of allocated operations (`nOpAlloc`). Allocate memory for `chID1` and `chID2`, copies the corresponding data from `chID1` and `chID2`, and adds the new operation. |
| **Side effects** | Updates the `nOpAlloc` counter, expands `chID1` and `chID2`.                                                              |
| **Errors / exceptions** | Returns an error code if `nOpAlloc` is not large enough.                                            |

#### `AddRenChain(int Nop, const ChainID ch1, const ChainID ch2)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Adds a renaming chain to the `GenSym` object, linking two chains together.                             |
| **Parameters** | `Nop` — Index of the operation to which the chain is associated. `ch1` — `ChainID` representing the first chain. `ch2` — `ChainID` representing the second chain. |
| **Returns**  | `int` — Returns `SYMOP_Ok` if the chain was added successfully, or an error code otherwise.            |
| **Behaviour** | Retrieves the `chID1` and `chID2` arrays.  Copies the `ch1` and `ch2` `ChainID` values into the corresponding positions in `chID1` and `chID2`. Increments `nChains` for the specified operation. |
| **Side effects** | Updates `chID1` and `chID2`.                                                                     |
| **Errors / exceptions** | Returns an error code if `nChains` is not large enough.                                            |

#### `Copy(PSymOps GenSym)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Copies the data from another `GenSym` object to this `GenSym` object.                                 |
| **Parameters** | `GenSym` — Pointer to the `GenSym` object to copy from.                                                  |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Copies the `nOpAlloc` value, `chID1`, `chID2`, and `nChains` arrays from the source `GenSym` object.            |
| **Side effects** | Updates the `nOpAlloc` counter, expands `chID1` and `chID2`.                                                              |
| **Errors / exceptions** | None                                                                                                  |

#### `write(io::RFile f)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Writes the GenSym object's data to the given RFile.                                                         |
| **Parameters** | `f` — An `io::RFile` object to write the data to.                                                          |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Writes the version number to the file. Calls `Root::write()` to write the object's data to the file. Writes the number of operation if it is non-zero. |
| **Side effects** | Writes the data to the RFile.                                                                         |
| **Errors / exceptions** | None                                                                                                  |

#### `read(io::RFile f)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads the GenSym object's data from the given RFile.                                                         |
| **Parameters** | `f` — An `io::RFile` object to read the data from.                                                        |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Reads the version number from the file. Calls `Root::read()` to read the object's data from the file. Reads the number of operation if it is non-zero. |
| **Side effects** | Reads the data from the RFile.                                                                         |
| **Errors / exceptions** | None                                                                                                  |

---

## `ContactIndex`

>  Stores the contact information from the GenSym object.

### Methods

#### `ContactIndex(PContact contact, int maxlen, int ncontacts, int max_alloc)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructor for the `ContactIndex` class. Initializes the index to store contact information.            |
| **Parameters** | `contact` — A pointer to the `PContact` object. `maxlen` — The maximum size of the index.   `ncontacts` — The initial number of contacts to store in the index.  `max_alloc` — A larger size to allow the index to grow. |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Initializes the member variables, including `contact_index`, `max_index`, `n_contacts` and `nAlloAtoms`. |
| **Side effects** | Initializes the member variables.                                                                      |
| **Errors / exceptions** | None                                                                                                  |

#### `~ContactIndex()`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Destructor for the `ContactIndex` class. Deallocates the memory allocated for the `contact_index`.  |
| **Parameters** | None                                                                                                   |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Deallocates the memory used for the `contact_index`.                                                    |
| **Side effects** | Deallocates the memory used for the `contact_index`.                                               |
| **Errors / exceptions** | None                                                                                                  |

#### `AddContact(int id1, int id2, realtype dist, int group)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Adds a contact to the `ContactIndex`.                                                                |
| **Parameters** | `id1` — Integer representing the ID of the first atom in the contact. `id2` — Integer representing the ID of the second atom in the contact. `dist` —  The distance between the two atoms. `group` —  The contact group ID. |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Adds the specified `id1` and `id2` values to the `contact_index` array, along with the `dist` and `group`. Increments `n_contacts`. |
| **Side effects** | Adds the contact data to the `contact_index` array.                                                        |
| **Errors / exceptions** | None                                                                                                  |

#### `GetIndex(RPContact contact, int & ncontacts)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Retrieves the data stored in the `ContactIndex` and assigns it to the given `RPContact` object.        |
| **Parameters** | `contact` — An `RPContact` object to receive the contact data. `ncontacts` — Integer representing the number of contacts in the index. |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Assigns the values from the `contact_index` array to the `contact` object, and assigns the value of `n_contacts` to the `ncontacts` parameter. |
| **Side effects** | Copies the contact data to the `contact` object.                                                        |
| **Errors / exceptions** | None                                                                                                  |

---

## `MBrick`

> Represents a multi-brick containing atoms from multiple structures.

### Methods

#### `MBrick(int nStructures)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructor for the `MBrick` class.  Initializes the multi-brick's parameters, allocating memory for the bricks and their corresponding atom indices. |
| **Parameters** | `nStructures` — Integer indicating the number of structures this multi-brick will contain.          |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Initializes `nStruct` to `nStructures`. Allocates memory for `atom` (an array of `PPAtom` pointers) and `id` (an array of `ivector` indices).  Initializes the size of each array to `nStruct`. |
| **Side effects** | Allocates memory for the `atom` and `id` arrays.                                                              |
| **Errors / exceptions** | None                                                                                                  |

#### `~MBrick()`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Destructor for the `MBrick` class. Releases the memory allocated for the `atom` and `id` arrays.               |
| **Parameters** | None                                                                                                   |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Deallocates the `atom` array by calling `delete[] atom;`. Deallocates the `id` array by calling `FreeMatrixMemory ( id,nStruct,0,0 );`. Sets `nStruct` to 0.  |
| **Side effects** | Deallocates the memory used for the `atom` and `id` arrays.                                             |
| **Errors / exceptions** | None                                                                                                  |

#### `Clear()`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Clears the MBrick object by releasing the allocated memory and resetting all member variables.              |
| **Parameters** | None                                                                                                   |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | If `atom` exists, it is deallocated. `FreeMatrixMemory` is called to deallocate the `id` array. Sets `nStruct` to 0. |
| **Side effects** | Deallocates the memory used for the `atom` and `id` arrays.                                             |
| **Errors / exceptions** | None                                                                                                  |

#### `AddAtom(PAtom A, int structNo, int atomid)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Adds an atom to the multi-brick, specifically assigning it to a particular structure.                    |
| **Parameters** | `A` — `PAtom` pointer to the atom to be added. `structNo` — Integer representing the structure number to which the atom belongs. `atomid` — Integer representing the atom's ID within the structure. |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Retrieves the current size of the `nAtoms` array for the specified structure.  If the current size is less than the allocated size, it expands the array.  Adds `A` to the end of the `atom` array and sets the corresponding `atomid` in the `id` array. Increments `nAtoms` by 1.  |
| **Side effects** | Expands the `atom` array if needed and allocates memory for a larger array. |
| **Errors / exceptions** | None                                                                                                  |

---

## `GenSym`

> Represents a symmetry operation for generating symmetric models.

### Methods

#### `GenSym()`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructor for the `GenSym` class. Initializes the class. |
| **Parameters** | None                                                                                                   |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Initializes member variables, including `chID1` (an array of `ChainID` pointers), `chID2` (another array of `ChainID` pointers), `nChains` (an array to store the number of chains associated with each operation), and `nOpAlloc` (the allocated size of the `chID1` and `chID2` arrays). |
| **Side effects** | Initializes member variables.                                                                       |
| **Errors / exceptions** | None                                                                                                  |

#### `GenSym(io::RPStream Object)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructor for the `GenSym` class. Initializes the class and sets up the output stream.                    |
| **Parameters** | `Object` — An `io::RPStream` object to handle output operations.                                     |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Calls the default constructor (`GenSym()`) to initialize the class and sets up the output stream.         |
| **Side effects** | Initializes the `io::RPStream` object.                                                                |
| **Errors / exceptions** | None                                                                                                  |

#### `~GenSym()`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Destructor for the `GenSym` class. Releases the memory allocated for the arrays of `ChainID`.            |
| **Parameters** | None                                                                                                   |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Deallocates the memory for `chID1` and `chID2`.  Calls `FreeVectorMemory` to deallocate the `nChains` array. Sets `nOpAlloc` to 0. |
| **Side effects** | Deallocates the memory used for the `chID1`, `chID2`, and `nChains` arrays.                               |
| **Errors / exceptions** | None                                                                                                  |

#### `FreeMemory()`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Releases memory allocated for `chID1`, `chID2`, and `nChains`.                                             |
| **Parameters** | None                                                                                                   |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Iterates through the `chID1` and `chID2` arrays, deleting each `ChainID` array.  Calls `FreeVectorMemory` to deallocate the `nChains` array. Sets `nOpAlloc` to 0.  Calls `SymOps::FreeMemory()` to perform any additional cleanup tasks required by the base class. |
| **Side effects** | Deallocates the memory used for the `chID1`, `chID2`, and `nChains` arrays.                               |
| **Errors / exceptions** | None                                                                                                  |

#### `AddSymOp(cpstr XYZOperation)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Adds a symmetry operation to the `GenSym` object.                                                      |
| **Parameters** | `XYZOperation` — A `cpstr` pointer to the string representing the symmetry operation.                     |
| **Returns**  | `int` — Returns `SPOSEAT_Ok` if the operation was added successfully, or an error code otherwise.          |
| **Behaviour** | Increments the number of allocated operations (`nOpAlloc`). Allocate memory for `chID1` and `chID2`, copies the corresponding data from `chID1` and `chID2`, and adds the new operation. |
| **Side effects** | Updates the `nOpAlloc` counter, expands `chID1` and `chID2`.                                                              |
| **Errors / exceptions** | Returns an error code if `nOpAlloc` is not large enough.                                            |

#### `AddRenChain(int Nop, const ChainID ch1, const ChainID ch2)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Adds a renaming chain to the `GenSym` object, linking two chains together.                             |
| **Parameters** | `Nop` — Index of the operation to which the chain is associated. `ch1` — `ChainID` representing the first chain. `ch2` — `ChainID` representing the second chain. |
| **Returns**  | `int` — Returns `SYMOP_Ok` if the chain was added successfully, or an error code otherwise.            |
| **Behaviour** | Retrieves the `chID1` and `chID2` arrays.  Copies the `ch1` and `ch2` `ChainID` values into the corresponding positions in `chID1` and `chID2`. Increments `nChains` for the specified operation. |
| **Side effects** | Updates `chID1` and `chID2`.                                                                     |
| **Errors / exceptions** | Returns an error code if `nChains` is not large enough.                                            |

#### `Copy(PSymOps GenSym)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Copies the data from another `GenSym` object to this `GenSym` object.                                 |
| **Parameters** | `GenSym` — Pointer to the `GenSym` object to copy from.                                                  |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Copies the `nOpAlloc` value, `chID1`, `chID2`, and `nChains` arrays from the source `GenSym` object.            |
| **Side effects** | Updates the `nOpAlloc` counter, expands `chID1` and `chID2`.                                                              |
| **Errors / exceptions** | None                                                                                                  |

#### `write(io::RFile f)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Writes the GenSym object's data to the given RFile.                                                         |
| **Parameters** | `f` — An `io::RFile` object to write the data to.                                                        |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Writes the version number to the file. Calls `Root::write()` to write the object's data to the file. Writes the number of operation if it is non-zero. |
| **Side effects** | Writes the data to the RFile.                                                                         |
| **Errors / exceptions** | None                                                                                                  |

#### `read(io::RFile f)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads the GenSym object's data from the given RFile.                                                         |
| **Parameters** | `f` — An `io::RFile` object to read the data from.                                                        |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Reads the version number from the file. Calls `Root::read()` to read the object's data from the file. Reads the number of operation if it is non-zero. |
| **Side effects** | Reads the data from the RFile.                                                                         |
| **Errors / exceptions** | None                                                                                                  |

---

## `Contact`

> Represents a single contact between two atoms.

### Methods

#### `Copy(RContact c)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Copies the data from another `Contact` object to this `Contact` object.                                |
| **Parameters** | `c` — A pointer to the `RContact` object from which to copy the data.                              |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Copies the values of `id1`, `id2`, `group` and `dist` from the source `RContact` object to the current `Contact` object. |
| **Side effects** | Updates the member variables of the current `Contact` object.                                        |
| **Errors / exceptions** | None                                                                                                  |

#### `Swap(RContact c)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Swaps the data between this `Contact` object and another `RContact` object.                             |
| **Parameters** | `c` — A pointer to the `RContact` object to swap data with.                                          |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Swaps the values of `id1`, `id2`, `group` and `dist` between the current `Contact` object and the source `RContact` object. |
| **Side effects** | Updates the member variables of the current `Contact` object.                                        |
| **Errors / exceptions** | None                                                                                                  |

---

## `MContact`

> Represents a multi-contact storing information for multiple atom pairs.

### Methods

#### `MContact(int nStructures)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructor for the `MContact` class.  Initializes the multi-contact's parameters, allocating memory for the contacts and their corresponding atom indices. |
| **Parameters** | `nStructures` — Integer indicating the number of structures this multi-contact will contain.          |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Initializes `nStruct` to `nStructures`. Allocates memory for `atom` (an array of `PPAtom` pointers) and `id` (an array of `ivector` indices).  Initializes the size of each array to `nStruct`. |
| **Side effects** | Allocates memory for the `atom` and `id` arrays.                                                              |
| **Errors / exceptions** | None                                                                                                  |

#### `~MContact()`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Destructor for the `MContact` class. Releases the memory allocated for the `atom` and `id` arrays.               |
| **Parameters** | None                                                                                                   |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Deallocates the `atom` array by calling `delete[] atom;`. Deallocates the `id` array by calling `FreeMatrixMemory ( id,nStruct,0,0 );`. Sets `nStruct` to 0.  |
| **Side effects** | Deallocates the memory used for the `atom` and `id` arrays.                                             |
| **Errors / exceptions** | None                                                                                                  |

#### `Clear()`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Clears the MContact object by releasing the allocated memory and resetting all member variables.              |
| **Parameters** | None                                                                                                   |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | If `atom` exists, it is deallocated. `FreeMatrixMemory` is called to deallocate the `id` array. Sets `nStruct` to 0. |
| **Side effects** | Deallocates the memory used for the `atom` and `id` arrays.                                             |
| **Errors / exceptions** | None                                                                                                  |

#### `AddContact(PAtom A, int structNo, int atomid)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Adds an atom to the multi-contact, specifically assigning it to a particular structure.                    |
| **Parameters** | `A` — `PAtom` pointer to the atom to be added. `structNo` — Integer representing the structure number to which the atom belongs. `atomid` — Integer representing the atom's ID within the structure. |
| **Returns**  | `void`                                                                                                  |
| **Behaviour** | Retrieves the current size of the `nAtoms` array for the specified structure.  If the current size is less than the allocated size, it expands the array.  Adds `A` to the end of the `atom` array, and sets the corresponding `atomid` in the `id` array. Increments `nAtoms` by 1.  |
| **Side effects** | Expands the `atom` array if needed and allocates memory for a larger array. |
| **Errors / exceptions** | None                                                                                                  |

---

## `

---

### `mmdb_cryst`

---
## `CrystContainer`

>  This class serves as a container for other crystal data structures, primarily the NCSMatrix, TVect, and Cryst. It provides a convenient way to store and manage these components within a larger crystal object.

### Methods

#### `ERROR_CODE AddMTRIXLine(cpstr S)`

| Field        | Detail                                                                                                                                                                                        |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Converts a PDB line containing an NCS matrix (e.g., MTRIX1) into the internal representation of the `NCSMatrix` object.                                                                  |
| **Parameters**| `S` — A pointer to a C-string containing the PDB line to be converted. The line must start with "MTRIX" followed by a number (1, 2, or 3) representing the matrix.  `none` |
| **Returns**  | `Error_NoError` if the conversion was successful.  Other error codes (e.g., `Error_NCSM_WrongSerial`) if there were problems during the conversion process.                               |
| **Behaviour** | 1.  It first checks if the `S` string is valid and contains "MTRIX" followed by a number between 1 and 3. 2.  If a valid matrix line is found, it calls the `ConvertPDBASCII()` method of the `NCSMatrix` object, passing the input string `S`. 3. Returns the error code returned by `ConvertPDBASCII()`. |
| **Side effects**| Modifies the internal state of the `NCSMatrix` object by populating its fields (serial number, matrix elements, translation vector, etc.) with the data from the PDB line. |
| **Errors / exceptions**| `Error_NCSM_WrongSerial` if the serial number in the PDB line doesn't match the existing serial number of the `NCSMatrix` object.  `Error_NCSM_Unrecognized` if the PDB line doesn't match the expected "MTRIX" format. `Error_WrongSection` if the PDB line does not correspond to the `Cryst`  structure. |

#### `PContainerClass MakeContainerClass(int ClassID)`

| Field        | Detail                                                                                                                                                                                                                                                         |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Creates a new instance of the appropriate container class based on the provided `ClassID`.  This is used by the `CrystContainer` to handle different types of data within the container.                                                                   |
| **Parameters**| `ClassID` — An integer representing the type of container class to create.  Common values are `ClassID_Template`, `ClassID_NCSMatrix`, and `ClassID_TVect`. `none`                                                                 |
| **Returns**  | A pointer to a newly created instance of the container class (e.g., a pointer to an `NCSMatrix` object).  Returns a pointer to the `ClassContainer` if the `ClassID` is `ClassID_Template`. |
| **Behaviour** | The method determines the appropriate container class based on the `ClassID` and returns a pointer to a newly created object of that class. This new object is a pointer to the `ContainerClass` base class. |
| **Side effects**| Does not modify any existing objects.                                                                                                                                                                                                                     |
| **Errors / exceptions**| Returns a pointer to the `ClassContainer` class, if `ClassID` is `ClassID_Template`.  |

#### `ERROR_CODE AddMTRIXLine(cpstr S)`

| Field        | Detail                                                                                                                                                                                        |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Converts a PDB line containing an NCS matrix (e.g., MTRIX1) into the internal representation of the `NCSMatrix` object.                                                                  |
| **Parameters**| `S` — A pointer to a C-string containing the PDB line to be converted. The line must start with "MTRIX" followed by a number (1, 2, or 3) representing the matrix.  `none` |
| **Returns**  | `Error_NoError` if the conversion was successful.  Other error codes (e.g., `Error_NCSM_WrongSerial`) if there were problems during the conversion process.                               |
| **Behaviour** | 1.  It first checks if the `S` string is valid and contains "MTRIX" followed by a number between 1 and 3. 2.  If a valid matrix line is found, it calls the `ConvertPDBASCII()` method of the `NCSMatrix` object, passing the input string `S`. 3. Returns the error code returned by `ConvertPDBASCII()`. |
| **Side effects**| Modifies the internal state of the `NCSMatrix` object by populating its fields (serial number, matrix elements, translation vector, etc.) with the data from the PDB line. |
| **Errors / exceptions**| `Error_NCSM_WrongSerial` if the serial number in the PDB line doesn't match the existing serial number of the `NCSMatrix` object.  `Error_NCSM_Unrecognized` if the PDB line doesn't match the expected "MTRIX" format. `Error_WrongSection` if the PDB line does not correspond to the `Cryst`  structure. |

---
## `NCSMatrix`

> This class represents a non-crystallographic symmetry matrix and a translational vector, often encountered in crystallography.

### Methods

#### `NCSMatrix()`

| Field        | Detail                                                               |
|--------------|-----------------------------------------------------------------------|
| **Summary**  | Default constructor for the `NCSMatrix` class. Initializes the container's internal state to default values.                                |
| **Parameters**| `none`                                                                |
| **Returns**  | `void`                                                                |
| **Behaviour** | Initializes the `serNum`, `m`, `v`, and `iGiven` fields of the `NCSMatrix` object to their default values (e.g., `serNum` is -1, matrices are filled with zeros, `v` is zero vector, `iGiven` is -1). |
| **Side effects**| Initializes the internal state of the `NCSMatrix` object.               |
| **Errors / exceptions**| `none`                                                               |

#### `NCSMatrix(cpstr S)`

| Field        | Detail                                                                                                                                                                            |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructs an `NCSMatrix` object from a provided PDB string containing the matrix elements and other relevant information.                                                              |
| **Parameters**| `S` — A pointer to a C-string containing the PDB line to be parsed.  The line must start with "MTRIX" followed by a number (1, 2, or 3) and then the matrix elements.  |
| **Returns**  | `void`                                                                                                                                                                            |
| **Behaviour** | 1. Calls the `ConvertPDBASCII()` method to parse the input PDB string `S`. 2. The `ConvertPDBASCII()` method then fills the `m`, `v` and `iGiven` fields with the values from the PDB string. |
| **Side effects**| Populates the internal fields of the `NCSMatrix` object with data extracted from the provided PDB string.                                                                    |
| **Errors / exceptions**| Returns an `Error_NCSM_XXX` error code if the input string is invalid or if the parsing process encounters any problems.  Returns `Error_NoError` on successful conversion.             |

#### `NCSMatrix(io::RPStream Object)`

| Field        | Detail                                                                                                        |
|--------------|----------------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructs an `NCSMatrix` object from an `io::RPStream` object. This is typically used when reading from a file. |
| **Parameters**| `Object` — An `io::RPStream` object representing the input stream from which to read the data.                |
| **Returns**  | `void`                                                                                                        |
| **Behaviour** | The constructor calls the `Init()` method to initialize the container to its default state. Then calls `ConvertPDBASCII()` to parse the PDB line.                                  |
| **Side effects**| Populates the internal fields of the `NCSMatrix` object with data extracted from the input stream.            |
| **Errors / exceptions**| Returns an `Error_NCSM_XXX` error code if the parsing process encounters any problems. Returns `Error_NoError` on successful conversion.                              |

#### `~NCSMatrix()`

| Field        | Detail                                                               |
|--------------|-----------------------------------------------------------------------|
| **Summary**  | Destructor for the `NCSMatrix` class.                                    |
| **Parameters**| `none`                                                                |
| **Returns**  | `void`                                                                |
| **Behaviour** | Calls the `FreeContainer()` method to release any memory allocated for the container's internal data structures.   |
| **Side effects**| Releases memory associated with the `NCSMatrix` object.                   |
| **Errors / exceptions**| `none`                                                               |

#### `bool PDBASCIIDump1(io::RFile f)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Writes the `NCSMatrix` data to an `io::RFile` object in a PDB-compatible ASCII format. This method is called to convert the `NCSMatrix` to a PDB-compatible format. |
| **Parameters**| `f` — An `io::RFile` object to which the PDB data will be written. `none`                                |
| **Returns**  | `true` if the data was written successfully, `false` if there was an error.                              |
| **Behaviour** | Formats the matrix elements and translation vector as PDB ASCII lines and writes them to the specified `io::RFile` object. |
| **Side effects**| Writes data to the output file.                                                                        |
| **Errors / exceptions**| `none`                                                                                                |

#### `ERROR_CODE ConvertPDBASCII(cpstr S)`

| Field        | Detail                                                                                                                                                                        |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Parses a PDB string containing the matrix elements and translation vector, and populates the corresponding fields of the `NCSMatrix` object.                                       |
| **Parameters**| `S` — A pointer to a C-string containing the PDB line to be parsed. The line must start with "MTRIX" followed by a number (1, 2, or 3) and then the matrix elements.  `none` |
| **Returns**  | `Error_NoError` if the parsing was successful.  Other error codes (e.g., `Error_NCSM_Unrecognized`) if there were problems during the parsing process.                         |
| **Behaviour** | 1. Checks if the input string `S` is valid and contains the expected format (e.g., "MTRIX1"). 2. If the format is correct, extracts the matrix elements from the string and assigns them to the `m` field of the `NCSMatrix` object. 3. It sets the `serNum` to the value from the PDB string, and `iGiven` to the value from the PDB string.  |
| **Side effects**| Populates the `m`, `iGiven`, `serNum` and `v` members of the `NCSMatrix` object with the parsed data.                                                                      |
| **Errors / exceptions**| Returns `Error_NCSM_Unrecognized` if the input string does not conform to the expected format (e.g., if it doesn't start with "MTRIX" or if the matrix elements are not in the correct order).  Returns `Error_WrongSection` if the PDB line is not of the correct type. |

#### `void MakeCIF(mmcif::PData CIF, int N)`

| Field        | Detail                                                                                               |
|--------------|------------------------------------------------------------------------------------------------------|
| **Summary**  | Generates a CIF (Chemical Information File) record for the `NCSMatrix` object.                         |
| **Parameters**| `CIF` — A pointer to a `mmcif::PData` object representing the CIF structure.  `N` — The number of records to write to the CIF file. |
| **Returns**  | `void`                                                                                              |
| **Behaviour** | Adds the `NCSMatrix` data to the CIF file, including its serial number, matrix elements, translation vector, and space group information. |
| **Side effects**| Modifies the `mmcif::PData` object, adding the new record.                                         |
| **Errors / exceptions**| `none`                                                                                            |

#### `ERROR_CODE GetCIF(mmcif::PData CIF, int & n)`

| Field        | Detail                                                                                                        |
|--------------|---------------------------------------------------------------------------------------------------------------|
| **Summary**  | Retrieves the `NCSMatrix` data from a `mmcif::PData` object, and populates the `NCSMatrix` object.        |
| **Parameters**| `CIF` — A pointer to a `mmcif::PData` object representing the CIF structure.  `n` — A reference to an integer that will be set to the number of records written to the CIF file. |
| **Returns**  | `Error_NoError` if the retrieval was successful.  Other error codes (e.g., `Error_EmptyCIF`) if there were problems. |
| **Behaviour** | Retrieves the `NCSMatrix` data from the CIF structure, populating the internal fields of the `NCSMatrix` object.  Returns `Error_EmptyCIF` if the CIF file does not contain a `NCSMatrix` record. |
| **Side effects**| Populates the internal fields of the `NCSMatrix` object with the data read from the CIF file. |
| **Errors / exceptions**| `Error_EmptyCIF` if the `NCSMatrix` record is not found in the CIF file.  |

---
## `TVect`

> This class represents a translation vector, commonly used in crystallography.

### Methods

#### `TVect()`

| Field        | Detail                                                               |
|--------------|-----------------------------------------------------------------------|
| **Summary**  | Default constructor for the `TVect` class. Initializes the container's internal state to default values. |
| **Parameters**| `none`                                                                |
| **Returns**  | `void`                                                                |
| **Behaviour** | Initializes the `serNum`, `t`, and `comment` fields of the `TVect` object to their default values (e.g., `serNum` is -1, `t` is a zero vector, `comment` is `NULL`).   |
| **Side effects**| Initializes the internal state of the `TVect` object.               |
| **Errors / exceptions**| `none`                                                               |

#### `TVect(cpstr S)`

| Field        | Detail                                                                                                                                                                            |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructs a `TVect` object from a provided PDB string containing the translation vector elements.                                                                                      |
| **Parameters**| `S` — A pointer to a C-string containing the PDB line to be parsed. The line must start with "TVECT" followed by a number (1, 2, or 3) and then the vector elements.  |
| **Returns**  | `void`                                                                                                                                                                            |
| **Behaviour** | 1. Calls the `ConvertPDBASCII()` method to parse the input PDB string `S`. 2. The `ConvertPDBASCII()` method then fills the `t` field with the values from the PDB string.   |
| **Side effects**| Populates the `t` field of the `TVect` object with the data extracted from the provided PDB string.                                                                                 |
| **Errors / exceptions**| Returns an `Error_NCSM_XXX` error code if the input string is invalid or if the parsing process encounters any problems.  Returns `Error_NoError` on successful conversion.             |

#### `TVect(io::RPStream Object)`

| Field        | Detail                                                                                                        |
|--------------|---------------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructs a `TVect` object from an `io::RPStream` object. This is typically used when reading from a file.    |
| **Parameters**| `Object` — An `io::RPStream` object representing the input stream from which to read the data.                |
| **Returns**  | `void`                                                                                                        |
| **Behaviour** | The constructor calls the `Init()` method to initialize the container to its default state. Then calls `ConvertPDBASCII()` to parse the PDB line.                                  |
| **Side effects**| Populates the internal fields of the `TVect` object with data extracted from the input stream.            |
| **Errors / exceptions**| Returns an `Error_NCSM_XXX` error code if the parsing process encounters any problems. Returns `Error_NoError` on successful conversion.                              |

#### `~TVect()`

| Field        | Detail                                                               |
|--------------|-----------------------------------------------------------------------|
| **Summary**  | Destructor for the `TVect` class.                                    |
| **Parameters**| `none`                                                                |
| **Returns**  | `void`                                                                |
| **Behaviour** | Calls the `FreeContainer()` method to release any memory allocated for the container's internal data structures.   |
| **Side effects**| Releases memory associated with the `TVect` object.                   |
| **Errors / exceptions**| `none`                                                               |

#### `void Init()`

| Field        | Detail                                                               |
|--------------|-----------------------------------------------------------------------|
| **Summary**  | Initializes the internal state of the `TVect` object to its default values. |
| **Parameters**| `none`                                                                |
| **Returns**  | `void`                                                                |
| **Behaviour** | Resets the `serNum`, `t`, and `comment` fields to their default values (e.g., `serNum` is -1, `t` is a zero vector, `comment` is `NULL`). |
| **Side effects**| Sets the internal state of the `TVect` object.                               |
| **Errors / exceptions**| `none`                                                               |

---
## `NCSMatrix`

> This class represents a non-crystallographic symmetry matrix and a translational vector, commonly used in crystallography.

### Methods

#### `NCSMatrix()`

| Field        | Detail                                                               |
|--------------|-----------------------------------------------------------------------|
| **Summary**  | Default constructor for the `NCSMatrix` class. Initializes the container's internal state to default values.                                |
| **Parameters**| `none`                                                                |
| **Returns**  | `void`                                                                |
| **Behaviour** | Initializes the `serNum`, `m`, `v`, and `iGiven` fields of the `NCSMatrix` object to their default values (e.g., `serNum` is -1, matrices are filled with zeros, `v` is zero vector, `iGiven` is -1).   |
| **Side effects**| Initializes the internal state of the `NCSMatrix` object.               |
| **Errors / exceptions**| `none`                                                               |

#### `NCSMatrix(cpstr S)`

| Field        | Detail                                                                                                                                                                            |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructs an `NCSMatrix` object from a provided PDB string containing the matrix elements and other relevant information.                                                              |
| **Parameters**| `S` — A pointer to a C-string containing the PDB line to be parsed. The line must start with "MTRIX" followed by a number (1, 2, or 3) and then the matrix elements.  |
| **Returns**  | `void`                                                                                                                                                                            |
| **Behaviour** | 1. Calls the `ConvertPDBASCII()` method to parse the input PDB string `S`. 2. The `ConvertPDBASCII()` method then fills the `m`, `v` and `iGiven` fields with the values from the PDB string.   |
| **Side effects**| Populates the internal fields of the `NCSMatrix` object with the data extracted from the provided PDB string.                                                                      |
| **Errors / exceptions**| Returns an `Error_NCSM_XXX` error code if the input string is invalid or if the parsing process encounters any problems.  Returns `Error_NoError` on successful conversion.             |

#### `NCSMatrix(io::RPStream Object)`

| Field        | Detail                                                                                                        |
|--------------|---------------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructs an `NCSMatrix` object from an `io::RPStream` object. This is typically used when reading from a file.    |
| **Parameters**| `Object` — An `io::RPStream` object representing the input stream from which to read the data.                |
| **Returns**  | `void`                                                                                                        |
| **Behaviour** | The constructor calls the `Init()` method to initialize the container to its default state. Then calls `ConvertPDBASCII()` to parse the PDB line.                                  |
| **Side effects**| Populates the internal fields of the `NCSMatrix` object with data extracted from the input stream.            |
| **Errors / exceptions**| Returns an `Error_NCSM_XXX` error code if the parsing process encounters any problems. Returns `Error_NoError` on successful conversion.                              |

#### `~NCSMatrix()`

| Field        | Detail                                                               |
|--------------|-----------------------------------------------------------------------|
| **Summary**  | Destructor for the `NCSMatrix` class.                                    |
| **Parameters**| `none`                                                                |
| **Returns**  | `void`                                                                |
| **Behaviour** | Calls the `FreeContainer()` method to release any memory allocated for the container's internal data structures.   |
| **Side effects**| Releases memory associated with the `NCSMatrix` object.                   |
| **Errors / exceptions**| `none`                                                               |

#### `void Init()`

| Field        | Detail                                                               |
|--------------|-----------------------------------------------------------------------|
| **Summary**  | Initializes the internal state of the `NCSMatrix` object to its default values. |
| **Parameters**| `none`                                                                |
| **Returns**  | `void`                                                                |
| **Behaviour** | Resets the `serNum`, `m`, `v` and `iGiven` fields of the `NCSMatrix` object to their default values (e.g., `serNum` is -1, `m` is filled with zeros, `v` is zero vector, `iGiven` is -1). |
| **Side effects**| Sets the internal state of the `NCSMatrix` object.                               |
| **Errors / exceptions**| `none`                                                               |

#### `bool PDBASCIIDump1(io::RFile f)`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Writes the `NCSMatrix` data to an `io::RFile` object in a PDB-compatible ASCII format. This method is called to convert the `NCSMatrix` to a PDB-compatible format. |
| **Parameters**| `f` — An `io::RFile` object to which the PDB data will be written. `none`                                |
| **Returns**  | `true` if the data was written successfully, `false` if there was an error.                              |
| **Behaviour** | Formats the matrix elements and translation vector as PDB ASCII lines and writes them to the specified `io::RFile` object. |
| **Side effects**| Writes data to the output file.                                                                        |
| **Errors / exceptions**| `none`                                                                                                |

#### `ERROR_CODE ConvertPDBASCII(cpstr S)`

| Field        | Detail                                                                                                                                                                        |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Parses a PDB string containing the matrix elements and translation vector, and populates the corresponding fields of the `NCSMatrix` object.                                       |
| **Parameters**| `S` — A pointer to a C-string containing the PDB line to be parsed. The line must start with "MTRIX" followed by a number (1, 2, or 3) and then the matrix elements.  `none` |
| **Returns**  | `Error_NoError` if the parsing was successful.  Other error codes (e.g., `Error_NCSM_Unrecognized`) if there were problems during the parsing process.                         |
| **Behaviour** | 1. Checks if the input string `S` is valid and contains the expected format (e.g., "MTRIX1"). 2. If the format is correct, extracts the matrix elements from the string and assigns them to the `m` field of the `NCSMatrix` object. 3. It sets the `serNum` to the value from the PDB string, and `iGiven` to the value from the PDB string.  |
| **Side effects**| Populates the `m`, `iGiven`, `serNum` and `v` members of the `NCSMatrix` object with the parsed data.                                                                      |
| **Errors / exceptions**| Returns `Error_NCSM_Unrecognized` if the input string does not conform to the expected format (e.g., if it doesn't start with "MTRIX" or if the matrix elements are not in the correct order).  Returns `Error_WrongSection` if the PDB line does not correspond to the `Cryst`  structure. |

#### `void MakeCIF(mmcif::PData CIF, int N)`

| Field        | Detail                                                                                               |
|--------------|------------------------------------------------------------------------------------------------------|
| **Summary**  | Generates a CIF (Chemical Information File) record for the `NCSMatrix` object.                         |
| **Parameters**| `CIF` — A pointer to a `mmcif::PData` object representing the CIF structure.  `N` — The number of records to write to the CIF file. |
| **Returns**  | `void`                                                                                              |
| **Behaviour** | Adds the `NCSMatrix` data to the CIF file, including its serial number, matrix elements, translation vector, and space group information. |
| **Side effects**| Modifies the `mmcif::PData` object, adding the new record.                                         |
| **Errors / exceptions**| `none`                                                                                            |

#### `ERROR_CODE GetCIF(mmcif::PData CIF, int & n)`

| Field        | Detail                                                                                                        |
|--------------|---------------------------------------------------------------------------------------------------------------|
| **Summary**  | Retrieves the `NCSMatrix` data from a `mmcif::PData` object, and populates the `NCSMatrix` object.        |
| **Parameters**| `CIF` — A pointer to a `mmcif::PData` object representing the CIF structure.  `n` — A reference to an integer that will be set to the number of records written to the CIF file. |
| **Returns**  | `Error_NoError` if the retrieval was successful.  Other error codes (e.g., `Error_EmptyCIF`) if there were problems. |
| **Behaviour** | Retrieves the `NCSMatrix` data from the CIF structure, populating the internal fields of the `NCSMatrix` object.  Returns `Error_EmptyCIF` if the CIF file does not contain a `NCSMatrix` record. |
| **Side effects**| Populates the internal fields of the `NCSMatrix` object with the data read from the CIF file. |
| **Errors / exceptions**| `Error_EmptyCIF` if the `NCSMatrix` record is not found in the CIF file. |

#### `bool areMatrices()`

| Field        | Detail                                                               |
|--------------|-----------------------------------------------------------------------|
| **Summary**  | Returns `true` if the orthogonal-to-fractional and fractional-to-orthogonal matrices are defined. |
| **Parameters**| `none`                                                                |
| **Returns**  | `true` if the matrices are defined, `false` otherwise.                   |
| **Behaviour** | Checks if the `NCode` member has a value, which indicates that the orthogonal-to-fractional and fractional-to-orthogonal matrices have been calculated. |
| **Side effects**| `none`                                                               |
| **Errors / exceptions**| `none`                                                               |

#### `bool Frac2Orth(realtype x, realtype y, realtype z, realtype & xx, realtype & yy, realtype & zz)`

| Field        | Detail                                                                                                                                                                                  |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Transforms fractional coordinates to orthogonal coordinates.                                                                                                                            |
| **Parameters**| `x, y, z` — The fractional coordinates to be transformed. `xx, yy, zz` — References to the orthogonal coordinates calculated from the transformation.  `none`                                 |
| **Returns**  | `true` if the transformation was successful, `false` if the matrices are not defined.                                                                                               |
| **Behaviour** | Performs the matrix multiplication to calculate the corresponding orthogonal coordinates and stores them in the `xx`, `yy`, and `zz` variables. |
| **Side effects**| Modifies the `xx`, `yy`, and `zz` variables with the calculated orthogonal coordinates.                                                                                              |
| **Errors / exceptions**| Returns `false` if the `areMatrices()` method returns `false`.                                                                                                                           |

#### `bool Orth2Frac(realtype x, realtype y, realtype z, realtype & xx, realtype & yy, realtype & zz)`

| Field        | Detail                                                                                                                                                                                  |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Transforms orthogonal coordinates to fractional coordinates.                                                                                                                            |
| **Parameters**| `x, y, z` — The orthogonal coordinates to be transformed. `xx, yy, zz` — References to the fractional coordinates calculated from the transformation.  `none`                                 |
| **Returns**  | `true` if the transformation was successful, `false` if the matrices are not defined.                                                                                               |
| **Behaviour** | Performs the matrix multiplication to calculate the corresponding fractional coordinates and stores them in the `xx`, `yy`, and `zz` variables. |
| **Side effects**| Modifies the `xx`, `yy`, and `zz` variables with the calculated fractional coordinates.                                                                                              |
| **Errors / exceptions**| Returns `false` if the `areMatrices()` method returns `false`.                                                                                                                           |

#### `void CalcCoordTransforms()`

| Field        | Detail                                                                                                                    |
|--------------|---------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Calculates the orthogonal-to-fractional and fractional-to-orthogonal transformation matrices, as well as the cell volume. |
| **Parameters**| `none`                                                                                                                    |
| **Returns**  | `void`                                                                                                                    |
| **Behaviour** |  If both orthogonal and fractional matrices are defined (i.e. `WhatIsSet` contains CSET_Transforms), this method calculates the transformation matrices, computes the cell volume, and updates the `RO`, `RF`, `ROU` and `RFU` matrices. |
| **Side effects**| Modifies the internal matrices and cell volume variables within the `Cryst` object.                                           |
| **Errors / exceptions**| `none`                                                                                                                 |

#### `bool Cryst2Orth(rvector U)`

| Field        | Detail                                                                                                                                                       |
|--------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Transforms orthogonal coordinates into fractional coordinates using the transformation matrices calculated in `CalcCoordTransforms()` method.                                     |
| **Parameters**| `U` — A `rvector` object containing the coordinates to be transformed. `none`                                                                                |
| **Returns**  | `true` if the transformation was successful, `false` if the matrices were not calculated.                                                                   |
| **Behaviour** | If the matrices have been calculated (i.e. `areMatrices()` returns `true`), this method performs the matrix transformation to calculate the corresponding fractional coordinates and stores them in the `U` variables. |
| **Side effects**| Modifies the `U` variables with the calculated fractional coordinates.                                                                                      |
| **Errors / exceptions**| Returns `false` if the `areMatrices()` method returns `false`.                                                                                                |

#### `bool Orth2Cryst(rvector U)`

| Field        | Detail                                                                                                                                                       |
|--------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Transforms orthogonal coordinates into fractional coordinates.                                                                                             |
| **Parameters**| `U` — A `rvector` object containing the coordinates to be transformed. `none`                                                                                |
| **Returns**  | `true` if the transformation was successful, `false` if the matrices were not calculated.                                                                   |
| **Behaviour** | If the matrices have been calculated (i.e. `areMatrices()` returns `true`), this method performs the matrix transformation to calculate the corresponding fractional coordinates and stores them in the `U` variables. |
| **Side effects**| Modifies the `U` variables with the calculated fractional coordinates.                                                                                      |
| **Errors / exceptions**| Returns `false` if the `areMatrices()` method returns `false`.                                                                                                |

#### `int GetNumberOfSymOps()`

| Field        | Detail                                                               |
|--------------|-----------------------------------------------------------------------|
| **Summary**  | Returns the number of symmetry operations defined for the `NCSMatrix`. |
| **Parameters**| `none`                                                                |
| **Returns**  | The number of symmetry operations (an integer).                       |
| **Behaviour** | Returns the number of symmetry operations that have been defined for the `NCSMatrix` object.                                                              |
| **Side effects**| `none`                                                               |
| **Errors / exceptions**| `none`                                                               |

#### `pstr GetSymOp(int Nop)`

| Field        | Detail                                                               |
|--------------|-----------------------------------------------------------------------|
| **Summary**  | Returns the name of the symmetry operation corresponding to the given `Nop` value.                                  |
| **Parameters**| `Nop` — The index of the symmetry operation to retrieve.              |
| **Returns**  | A pointer to a C-string containing the name of the symmetry operation. |
| **Behaviour** | Retrieves the name of the symmetry operation associated with the given `Nop` value from the `symOps` vector.                                                   |
| **Side effects**| `none`                                                               |
| **Errors / exceptions**| Returns `NULL` if no symmetry operation is defined for the given `Nop` value.                                                           |

#### `int GetTMatrix(mat44 & TMatrix, int Nop, int cellshift_a, int cellshift_b, int cellshift_c, PSymOps symOpers)`

| Field        | Detail                                                                                                                                                               

---

### `mmdb_defs`

```cpp
---
## `mmdb::Atom`

> Represents a single atom in a protein structure.

### Methods

#### `Atom()`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructs an `Atom` object, initializing its members with default values.                                |
| **Parameters** | `none` — No parameters are accepted.                                                                   |
| **Returns**  | `void` — No value is returned.                                                                        |
| **Behaviour** | Initializes `serial`, `elem`, `resn`, `altl`, `chainid`, `name`, `x`, `y`, `z`, `occupancy`, `phase`, `trans`, `adr1`, `adr2`, `adr3`, `sernum`, `frac_occ`, `occ`, `cycle`, `model_id` to their respective default values. |
| **Side effects** | None.                                                                                             |
| **Errors / exceptions** | None.                                                                                             |

#### `Atom(const mmdb::ChainID& chainid, const std::string& name, int serial, int resnum, int modelid)`

| Field        | Detail                                                                                                     |
|--------------|------------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructs an `Atom` object with specific `chainid`, `name`, `serial`, `resnum`, and `modelid` values.           |
| **Parameters** | `chainid` — The chain identifier (e.g., "A", "B").  A `mmdb::ChainID` struct is used.                  |
|                | `name` — The name of the atom (e.g., "CA", "CB"). A `std::string` is used.                              |
|                | `serial` — The sequence number of the atom.                                                              |
|                | `resnum` — The residue number of the atom.                                                                |
|                | `modelid` — The model number to which the atom belongs.                                                     |
| **Returns**  | `void` — No value is returned.                                                                           |
| **Behaviour** | Initializes `serial`, `elem`, `resn`, `altl`, `chainid`, `name`, `x`, `y`, `z`, `occupancy`, `phase`, `trans`, `adr1`, `adr2`, `adr3`, `sernum`, `frac_occ`, `occ`, `cycle`, `model_id` to the provided values.  The `elem` member is set to the default value for the `std::string` parameter, and the `adr1`, `adr2` and `adr3` members are set to 0.0 |
| **Side effects** | None.                                                                                             |
| **Errors / exceptions** | None.                                                                                             |

#### `~Atom()`

| Field        | Detail                                                                                               |
|--------------|------------------------------------------------------------------------------------------------------|
| **Summary**  | Destructor for the `Atom` class, releasing any resources that were allocated.                         |
| **Parameters** | `none` — No parameters are accepted.                                                               |
| **Returns**  | `void` — No value is returned.                                                                     |
| **Behaviour** | Performs cleanup operations, typically releasing any dynamically allocated memory, and sets all members to their default values before exiting. |
| **Side effects** | None.                                                                                             |
| **Errors / exceptions** | None.                                                                                             |

#### `set_x(const realtype& x)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the x-coordinate of the atom.                                                                |
| **Parameters** | `x` — The x-coordinate value.                                                                   |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `x` member of the `Atom` object to the given `x` value.                                  |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_y(const realtype& y)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the y-coordinate of the atom.                                                                |
| **Parameters** | `y` — The y-coordinate value.                                                                   |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `y` member of the `Atom` object to the given `y` value.                                  |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_z(const realtype& z)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the z-coordinate of the atom.                                                                |
| **Parameters** | `z` — The z-coordinate value.                                                                   |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `z` member of the `Atom` object to the given `z` value.                                  |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_occupancy(const realtype& occup)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the occupancy of the atom.                                                                 |
| **Parameters** | `occup` — The occupancy value.                                                                  |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `occupancy` member of the `Atom` object to the given `occup` value.                         |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_phase(const realtype& phase)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the phase of the atom.                                                                       |
| **Parameters** | `phase` — The phase value.                                                                       |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `phase` member of the `Atom` object to the given `phase` value.                            |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_trans(const realtype& trans)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the transmission of the atom.                                                               |
| **Parameters** | `trans` — The transmission value.                                                                 |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `trans` member of the `Atom` object to the given `trans` value.                            |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_adr1(const realtype& adr1)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the first anisotropic displacement parameter of the atom.                                       |
| **Parameters** | `adr1` — The first anisotropic displacement parameter value.                                         |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `adr1` member of the `Atom` object to the given `adr1` value.                            |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_adr2(const realtype& adr2)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the second anisotropic displacement parameter of the atom.                                      |
| **Parameters** | `adr2` — The second anisotropic displacement parameter value.                                        |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `adr2` member of the `Atom` object to the given `adr2` value.                            |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_adr3(const realtype& adr3)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the third anisotropic displacement parameter of the atom.                                      |
| **Parameters** | `adr3` — The third anisotropic displacement parameter value.                                        |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `adr3` member of the `Atom` object to the given `adr3` value.                            |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_sernum(int sernum)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the sequence number of the atom.                                                              |
| **Parameters** | `sernum` — The sequence number value.                                                              |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `sernum` member of the `Atom` object to the given `sernum` value.                           |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_frac_occ(const realtype& frac_occ)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the fractional occupancy of the atom.                                                            |
| **Parameters** | `frac_occ` — The fractional occupancy value.                                                          |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `frac_occ` member of the `Atom` object to the given `frac_occ` value.                        |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_occ(const realtype& occ)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the occupancy of the atom.                                                                 |
| **Parameters** | `occ` — The occupancy value.                                                                    |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `occ` member of the `Atom` object to the given `occ` value.                            |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_cycle(const int& cycle)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the cycle of the atom.                                                                       |
| **Parameters** | `cycle` — The cycle value.                                                                        |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `cycle` member of the `Atom` object to the given `cycle` value.                            |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_modelid(int modelid)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the model ID of the atom.                                                                     |
| **Parameters** | `modelid` — The model ID value.                                                                    |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `modelid` member of the `Atom` object to the given `modelid` value.                         |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `print()`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Prints the atom's data to the standard output stream.                                                   |
| **Parameters** | `none` — No parameters are accepted.                                                                   |
| **Returns**  | `void` — No value is returned.                                                                        |
| **Behaviour** | Prints the values of the atom's members (`serial`, `elem`, `resn`, `altl`, `chainid`, `name`, `x`, `y`, `z`, `occupancy`, `phase`, `trans`, `adr1`, `adr2`, `adr3`) to the console using the `std::cout` stream. |
| **Side effects** | Outputs to standard output.                                                                        |
| **Errors / exceptions** | None.                                                                                             |
---
## `mmdb::Residue`

> Represents a single residue in a protein structure.

### Methods

#### `Residue()`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructs a `Residue` object, initializing its members with default values.                                |
| **Parameters** | `none` — No parameters are accepted.                                                                   |
| **Returns**  | `void` — No value is returned.                                                                        |
| **Behaviour** | Initializes `serial`, `resn`, `altl`, `chainid`, `name`, `x`, `y`, `z`, `occupancy`, `phase`, `trans`, `adr1`, `adr2`, `adr3`, `helixid`, `strandid`, `sheetid`, `turnid`, `linkrid` to their respective default values. |
| **Side effects** | None.                                                                                             |
| **Errors / exceptions** | None.                                                                                             |

#### `Residue(const mmdb::ChainID& chainid, const std::string& name, int serial, int resnum)`

| Field        | Detail                                                                                                     |
|--------------|------------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructs a `Residue` object with specific `chainid`, `name`, `serial`, and `resnum` values.           |
| **Parameters** | `chainid` — The chain identifier (e.g., "A", "B").  A `mmdb::ChainID` struct is used.                  |
|                | `name` — The name of the residue (e.g., "CA", "CB"). A `std::string` is used.                              |
|                | `serial` — The sequence number of the residue.                                                              |
|                | `resnum` — The residue number of the residue.                                                               |
| **Returns**  | `void` — No value is returned.                                                                        |
| **Behaviour** | Initializes `serial`, `elem`, `altl`, `chainid`, `name`, `x`, `y`, `z`, `occupancy`, `phase`, `trans`, `adr1`, `adr2`, `adr3`, `helixid`, `strandid`, `sheetid`, `turnid`, `linkrid` to the provided values.  The `elem` member is set to the default value for the `std::string` parameter, and the `adr1`, `adr2` and `adr3` members are set to 0.0 |
| **Side effects** | None.                                                                                             |
| **Errors / exceptions** | None.                                                                                             |

#### `~Residue()`

| Field        | Detail                                                                                               |
|--------------|------------------------------------------------------------------------------------------------------|
| **Summary**  | Destructor for the `Residue` class, releasing any resources that were allocated.                         |
| **Parameters** | `none` — No parameters are accepted.                                                               |
| **Returns**  | `void` — No value is returned.                                                                     |
| **Behaviour** | Performs cleanup operations, typically releasing any dynamically allocated memory, and sets all members to their default values before exiting. |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_x(const realtype& x)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the x-coordinate of the residue.                                                               |
| **Parameters** | `x` — The x-coordinate value.                                                                   |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `x` member of the `Residue` object to the given `x` value.                                  |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_y(const realtype& y)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the y-coordinate of the residue.                                                                |
| **Parameters** | `y` — The y-coordinate value.                                                                   |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `y` member of the `Residue` object to the given `y` value.                                  |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_z(const realtype& z)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the z-coordinate of the residue.                                                                |
| **Parameters** | `z` — The z-coordinate value.                                                                   |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `z` member of the `Residue` object to the given `z` value.                                  |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_occupancy(const realtype& occup)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the occupancy of the residue.                                                                 |
| **Parameters** | `occup` — The occupancy value.                                                                  |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `occupancy` member of the `Residue` object to the given `occup` value.                        |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_phase(const realtype& phase)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the phase of the residue.                                                                       |
| **Parameters** | `phase` — The phase value.                                                                       |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `phase` member of the `Residue` object to the given `phase` value.                            |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_trans(const realtype& trans)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the transmission of the residue.                                                               |
| **Parameters** | `trans` — The transmission value.                                                                 |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `trans` member of the `Residue` object to the given `trans` value.                            |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_adr1(const realtype& adr1)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the first anisotropic displacement parameter of the residue.                                       |
| **Parameters** | `adr1` — The first anisotropic displacement parameter value.                                         |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `adr1` member of the `Residue` object to the given `adr1` value.                            |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_adr2(const realtype& adr2)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the second anisotropic displacement parameter of the residue.                                      |
| **Parameters** | `adr2` — The second anisotropic displacement parameter value.                                        |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `adr2` member of the `Residue` object to the given `adr2` value.                            |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_adr3(const realtype& adr3)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the third anisotropic displacement parameter of the residue.                                      |
| **Parameters** | `adr3` — The third anisotropic displacement parameter value.                                        |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `adr3` member of the `Residue` object to the given `adr3` value.                            |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_helixid(const char* helixid)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the helix ID of the residue.                                                                 |
| **Parameters** | `helixid` — The helix ID value.                                                                   |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `helixid` member of the `Residue` object to the given `helixid` value.                         |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_strandid(const char* strandid)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the strand ID of the residue.                                                                 |
| **Parameters** | `strandid` — The strand ID value.                                                                  |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `strandid` member of the `Residue` object to the given `strandid` value.                        |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_sheetid(const char* sheetid)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the sheet ID of the residue.                                                                  |
| **Parameters** | `sheetid` — The sheet ID value.                                                                   |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `sheetid` member of the `Residue` object to the given `sheetid` value.                        |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_turnid(const char* turnid)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the turn ID of the residue.                                                                   |
| **Parameters** | `turnid` — The turn ID value.                                                                     |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `turnid` member of the `Residue` object to the given `turnid` value.                            |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_linkrid(const char* linkrid)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the link RID of the residue.                                                                  |
| **Parameters** | `linkrid` — The link RID value.                                                                   |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `linkrid` member of the `Residue` object to the given `linkrid` value.                        |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `print()`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Prints the residue's data to the standard output stream.                                                   |
| **Parameters** | `none` — No parameters are accepted.                                                                   |
| **Returns**  | `void` — No value is returned.                                                                        |
| **Behaviour** | Prints the values of the residue's members (`serial`, `elem`, `resn`, `altl`, `chainid`, `name`, `x`, `y`, `z`, `occupancy`, `phase`, `trans`, `adr1`, `adr2`, `adr3`, `helixid`, `strandid`, `sheetid`, `turnid`, `linkrid`) to the console using the `std::cout` stream. |
| **Side effects** | Outputs to standard output.                                                                        |
| **Errors / exceptions** | None.                                                                                             |

---
## `mmdb::Chain`

> Represents a chain of amino acids in a protein structure.

### Methods

#### `Chain()`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructs a `Chain` object, initializing its members with default values.                                |
| **Parameters** | `none` — No parameters are accepted.                                                                   |
| **Returns**  | `void` — No value is returned.                                                                        |
| **Behaviour** | Initializes `name`, `id`, `chainid`, `model_id`, `x`, `y`, `z`, `occupancy`, `phase`, `trans`, `adr1`, `adr2`, `adr3` to their respective default values. |
| **Side effects** | None.                                                                                             |
| **Errors / exceptions** | None.                                                                                             |

#### `Chain(const char* name, const char* chainid, int modelid)`

| Field        | Detail                                                                                                     |
|--------------|------------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructs a `Chain` object with specific `name`, `chainid`, and `modelid` values.           |
| **Parameters** | `name` — The name of the chain (e.g., "A", "B").  A `std::string` is used.                  |
|                | `chainid` — The chain identifier (e.g., "A", "B").  A `mmdb::ChainID` struct is used.                  |
|                | `modelid` — The model number to which the chain belongs.                                                     |
| **Returns**  | `void` — No value is returned.                                                                        |
| **Behaviour** | Initializes `name`, `id`, `chainid`, `model_id`, `x`, `y`, `z`, `occupancy`, `phase`, `trans`, `adr1`, `adr2`, `adr3` to the provided values.  The `elem` member is set to the default value for the `std::string` parameter, and the `adr1`, `adr2` and `adr3` members are set to 0.0 |
| **Side effects** | None.                                                                                             |
| **Errors / exceptions** | None.                                                                                             |

#### `~Chain()`

| Field        | Detail                                                                                               |
|--------------|------------------------------------------------------------------------------------------------------|
| **Summary**  | Destructor for the `Chain` class, releasing any resources that were allocated.                         |
| **Parameters** | `none` — No parameters are accepted.                                                               |
| **Returns**  | `void` — No value is returned.                                                                     |
| **Behaviour** | Performs cleanup operations, typically releasing any dynamically allocated memory, and sets all members to their default values before exiting. |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_x(const realtype& x)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the x-coordinate of the chain.                                                               |
| **Parameters** | `x` — The x-coordinate value.                                                                   |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `x` member of the `Chain` object to the given `x` value.                                  |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_y(const realtype& y)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the y-coordinate of the chain.                                                                |
| **Parameters** | `y` — The y-coordinate value.                                                                   |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `y` member of the `Chain` object to the given `y` value.                                  |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_z(const realtype& z)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the z-coordinate of the chain.                                                                |
| **Parameters** | `z` — The z-coordinate value.                                                                   |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `z` member of the `Chain` object to the given `z` value.                                  |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_occupancy(const realtype& occup)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the occupancy of the chain.                                                                 |
| **Parameters** | `occup` — The occupancy value.                                                                  |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `occupancy` member of the `Chain` object to the given `occup` value.                        |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_phase(const realtype& phase)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the phase of the chain.                                                                       |
| **Parameters** | `phase` — The phase value.                                                                       |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `phase` member of the `Chain` object to the given `phase` value.                            |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_trans(const realtype& trans)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the transmission of the chain.                                                               |
| **Parameters** | `trans` — The transmission value.                                                                 |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `trans` member of the `Chain` object to the given `trans` value.                            |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_adr1(const realtype& adr1)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the first anisotropic displacement parameter of the chain.                                       |
| **Parameters** | `adr1` — The first anisotropic displacement parameter value.                                         |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `adr1` member of the `Chain` object to the given `adr1` value.                            |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_adr2(const realtype& adr2)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the second anisotropic displacement parameter of the chain.                                      |
| **Parameters** | `adr2` — The second anisotropic displacement parameter value.                                        |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `adr2` member of the `Chain` object to the given `adr2` value.                            |
| **Side effects** | None.                                                                                            |
| **Errors / exceptions** | None.                                                                                            |

#### `set_adr3(const realtype& adr3)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the third anisotropic displacement parameter of the chain.                                      |
| **Parameters** | `adr3` — The third anisotropic displacement parameter value.                                        |
| **Returns**  | `void` — No value is returned.                                                                    |
| **Behaviour** | Sets the `adr

---

### `mmdb_ficif`

---
## `mmCIFData`

> Represents a MMDB coordinate file, providing functions for creating, writing, and managing data within the file.

### Methods

#### `mmCIFData* MMDB_CCIF_Init()`

| Field        | Detail                               |
|--------------|--------------------------------------|
| **Summary**  | Initializes the MMDB coordinate data structure. |
| **Parameters** | none                                  |
| **Returns**  | `mmCIFData*` — Pointer to the initialized MMDB data structure. |
| **Behaviour** | Creates a new `mmcif::Data` object and assigns it to the `mmCIFData` member.  This effectively sets up the internal state for handling MMDB data. |
| **Side effects** | `mmCIFData` is set to point to the new data object. |
| **Errors / exceptions** | `none`                            |

#### `void MMDB_CCIF_Quit()`

| Field        | Detail                               |
|--------------|--------------------------------------|
| **Summary**  | Destroys the MMDB coordinate data structure, releasing any associated resources. |
| **Parameters** | none                                  |
| **Returns**  | `void`                                |
| **Behaviour** | Checks if `mmCIFData` is valid. If so, deletes the `mmcif::Data` object pointed to by `mmCIFData` before setting `mmCIFData` to `NULL`. |
| **Side effects** | Destroys the `mmcif::Data` object. |
| **Errors / exceptions** | `none`                            |

#### `int MMDB_CCIF_Create(pstr DataName)`

| Field        | Detail                                      |
|--------------|---------------------------------------------|
| **Summary**  | Creates a new MMDB coordinate file, named `DataName`. |
| **Parameters** | `DataName` — Pointer to a null-terminated string representing the name of the MMDB file. |
| **Returns**  | `int` — 0 on success, 1 on failure.   |
| **Behaviour** | Creates a new `mmcif::Data` object with the provided `DataName`, initializing the internal state for this specific file. |
| **Side effects** | `mmCIFData` is updated to point to the new `mmcif::Data` object. |
| **Errors / exceptions** | `none`                            |

#### `int MMDB_CCIF_Write(pstr FileName)`

| Field        | Detail                                         |
|--------------|------------------------------------------------|
| **Summary**  | Writes the MMDB data to a file named `FileName`. |
| **Parameters** | `FileName` — Pointer to a null-terminated string representing the name of the file to write to. |
| **Returns**  | `int` — 0 on success, 1 on failure.       |
| **Behaviour** | Writes the MMDB data to a file identified by `FileName` using the `mmcif::Data::WriteMMCIFData` function.  Handles potential errors and returns accordingly. |
| **Side effects** | Writes to file `FileName`. |
| **Errors / exceptions** | `none`                            |

#### `int MMDB_CCIF_PutDate(pstr CatName, pstr Tag)`

| Field        | Detail                                  |
|--------------|-----------------------------------------|
| **Summary**  | Adds a date entry to the MMDB file.       |
| **Parameters** | `CatName` — Category name, `Tag` — Tag name. |
| **Returns**  | `int` — 0 on success, 1 on failure.     |
| **Behaviour** | Uses the `mmcif::Data::PutDate` function to add a date entry with the given category name and tag name.  |
| **Side effects** | Adds a date entry to the MMDB file. |
| **Errors / exceptions** | `none`                            |

#### `int MMDB_CCIF_PutDot(pstr CatName, pstr Tag)`

| Field        | Detail                                  |
|--------------|-----------------------------------------|
| **Summary**  | Adds a dot entry to the MMDB file.       |
| **Parameters** | `CatName` — Category name, `Tag` — Tag name. |
| **Returns**  | `int` — 0 on success, 1 on failure.     |
| **Behaviour** | Uses the `mmcif::Data::PutDot` function to add a dot entry with the given category name and tag name.  |
| **Side effects** | Adds a dot entry to the MMDB file. |
| **Errors / exceptions** | `none`                            |

#### `int MMDB_CCIF_PutQuestion(pstr CatName, pstr Tag)`

| Field        | Detail                                  |
|--------------|-----------------------------------------|
| **Summary**  | Adds a question entry to the MMDB file.   |
| **Parameters** | `CatName` — Category name, `Tag` — Tag name. |
| **Returns**  | `int` — 0 on success, 1 on failure.     |
| **Behaviour** | Uses the `mmcif::Data::PutQuestion` function to add a question entry with the given category name and tag name.  |
| **Side effects** | Adds a question entry to the MMDB file. |
| **Errors / exceptions** | `none`                            |

#### `int MMDB_CCIF_PutString(pstr Data, pstr CatName, pstr Tag)`

| Field        | Detail                                            |
|--------------|---------------------------------------------------|
| **Summary**  | Adds a string entry to the MMDB file.              |
| **Parameters** | `Data` — Pointer to a null-terminated string, `CatName` — Category name, `Tag` — Tag name. |
| **Returns**  | `int` — 0 on success, 1 on failure.                |
| **Behaviour** | Uses the `mmcif::Data::PutString` function to add a string entry with the given data, category name, and tag name. |
| **Side effects** | Adds a string entry to the MMDB file. |
| **Errors / exceptions** | `none`                            |

#### `int MMDB_CCIF_PutReal(apireal V, pstr CatName, pstr Tag)`

| Field        | Detail                                     |
|--------------|--------------------------------------------|
| **Summary**  | Adds a real number entry to the MMDB file.   |
| **Parameters** | `V` — Real value, `CatName` — Category name, `Tag` — Tag name. |
| **Returns**  | `int` — 0 on success, 1 on failure.        |
| **Behaviour** | Uses the `mmcif::Data::PutReal` function to add a real number entry with the given value, category name, and tag name. |
| **Side effects** | Adds a real number entry to the MMDB file. |
| **Errors / exceptions** | `none`                            |

#### `int MMDB_CCIF_PutInteger(int I, pstr CatName, pstr Tag)`

| Field        | Detail                                    |
|--------------|-------------------------------------------|
| **Summary**  | Adds an integer entry to the MMDB file.    |
| **Parameters** | `I` — Integer value, `CatName` — Category name, `Tag` — Tag name. |
| **Returns**  | `int` — 0 on success, 1 on failure.       |
| **Behaviour** | Uses the `mmcif::Data::PutInteger` function to add an integer entry with the given value, category name, and tag name. |
| **Side effects** | Adds an integer entry to the MMDB file. |
| **Errors / exceptions** | `none`                            |

#### `int MMDB_CCIF_PutLoopDot(pstr CatName, pstr Tag, int nrow)`

| Field        | Detail                                      |
|--------------|---------------------------------------------|
| **Summary**  | Adds a loop dot entry to the MMDB file.      |
| **Parameters** | `CatName` — Category name, `Tag` — Tag name, `nrow` — Row number. |
| **Returns**  | `int` — 0 on success, 1 on failure.        |
| **Behaviour** | Uses the `mmcif::Data::PutLoopDot` function to add a loop dot entry with the given category name, tag name, and row number. |
| **Side effects** | Adds a loop dot entry to the MMDB file. |
| **Errors / exceptions** | `none`                            |

#### `int MMDB_CCIF_PutLoopQuestion(pstr CatName, pstr Tag, int nrow)`

| Field        | Detail                                      |
|--------------|---------------------------------------------|
| **Summary**  | Adds a loop question entry to the MMDB file. |
| **Parameters** | `CatName` — Category name, `Tag` — Tag name, `nrow` — Row number. |
| **Returns**  | `int` — 0 on success, 1 on failure.        |
| **Behaviour** | Uses the `mmcif::Data::PutLoopQuestion` function to add a loop question entry with the given category name, tag name, and row number. |
| **Side effects** | Adds a loop question entry to the MMDB file. |
| **Errors / exceptions** | `none`                            |

#### `int MMDB_CCIF_PutLoopString(pstr Data, pstr CatName, pstr Tag, int nrow)`

| Field        | Detail                                      |
|--------------|---------------------------------------------|
| **Summary**  | Adds a loop string entry to the MMDB file.   |
| **Parameters** | `Data` — Pointer to a null-terminated string, `CatName` — Category name, `Tag` — Tag name, `nrow` — Row number. |
| **Returns**  | `int` — 0 on success, 1 on failure.        |
| **Behaviour** | Uses the `mmcif::Data::PutLoopString` function to add a loop string entry with the given data, category name, tag name, and row number. |
| **Side effects** | Adds a loop string entry to the MMDB file. |
| **Errors / exceptions** | `none`                            |

#### `int MMDB_CCIF_PutLoopReal(apireal V, pstr CatName, pstr Tag, int nrow)`

| Field        | Detail                                     |
|--------------|--------------------------------------------|
| **Summary**  | Adds a loop real entry to the MMDB file.    |
| **Parameters** | `V` — Real value, `CatName` — Category name, `Tag` — Tag name, `nrow` — Row number. |
| **Returns**  | `int` — 0 on success, 1 on failure.       |
| **Behaviour** | Uses the `mmcif::Data::PutLoopReal` function to add a loop real entry with the given value, category name, tag name, and row number. |
| **Side effects** | Adds a loop real entry to the MMDB file. |
| **Errors / exceptions** | `none`                            |

#### `int MMDB_CCIF_PutLoopInteger(int I, pstr CatName, pstr Tag, int nrow)`

| Field        | Detail                                    |
|--------------|-------------------------------------------|
| **Summary**  | Adds a loop integer entry to the MMDB file. |
| **Parameters** | `I` — Integer value, `CatName` — Category name, `Tag` — Tag name, `nrow` — Row number. |
| **Returns**  | `int` — 0 on success, 1 on failure.       |
| **Behaviour** | Uses the `mmcif::Data::PutLoopInteger` function to add a loop integer entry with the given value, category name, tag name, and row number. |
| **Side effects** | Adds a loop integer entry to the MMDB file. |
| **Errors / exceptions** | `none`                            |

---


---

### `mmdb_io_file`

Okay, here's the comprehensive documentation for the `mmdb::io::File` class, generated from the provided C++ source code. This documentation includes details on all methods, including their parameters, return values, behaviour, and potential side effects or errors.

---

## `mmdb::io::File`

> Provides file input/output functionality for MMDB coordinate files, supporting various file formats and options.

### Methods

#### `File(word BufSize)`

| Field        | Detail                                                                                                      |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructor for the `File` class.  Initializes the file stream with the specified buffer size.                |
| **Parameters**| `BufSize` — The size of the buffer to use for file I/O.  Defaults to 4096.                         |
| **Returns**  | `void`                                                                                                        |
| **Behaviour** | Allocates a buffer of the specified size.  Sets internal state variables for file handling (e.g., `IOBuf`, `BufLen`). |
| **Side effects**| Allocates memory on the heap.                                                                               |
| **Errors / exceptions**| None                                                                                                     |

#### `~File()`

| Field        | Detail                                                                                        |
|--------------|----------------------------------------------------------------------------------------------|
| **Summary**  | Destructor for the `File` class.  Releases resources associated with the file stream.            |
| **Parameters**| `void`                                                                                         |
| **Returns**  | `void`                                                                                          |
| **Behaviour** | Closes the file handle, frees the allocated buffer (if any), and releases other resources.              |
| **Side effects**| Closes the file handle.                                                                        |
| **Errors / exceptions**| None                                                                                           |

#### `assign(cpstr FileName, bool Text, bool UniB, GZ_MODE gzMode)`

| Field        | Detail                                                                                                                                  |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Assigns a file name to the `File` object, configuring the file stream for reading/writing (text or binary), supports uncompressed and gzipped files. |
| **Parameters**| `FileName` — The name of the file to open.  `Text` —  If `true`, the file is treated as a text file (line-oriented). `UniB` - If `true`, the file is treated as uncompressed.  `gzMode` — Specifies the expected file format (e.g., `GZM_CHECK` for gzipped files). |
| **Returns**  | `void`                                                                                                                                 |
| **Behaviour** | Closes the current file stream, frees existing buffers, and initializes the file stream parameters based on the given options. Handles text mode vs. binary mode. |
| **Side effects**| Closes the current file handle.  Allocates memory if necessary.                                                                      |
| **Errors / exceptions**| None                                                                                                                               |

#### `assign(word poolSize, word sizeInc, pstr filePool)`

| Field        | Detail                                                                                        |
|--------------|----------------------------------------------------------------------------------------------|
| **Summary**  | Assigns a memory pool to the `File` object. This enables memory-based I/O.                     |
| **Parameters**| `poolSize` — The size of the memory pool. `sizeInc` — The increment for buffer allocation. `filePool` — Pointer to the memory pool. |
| **Returns**  | `void`                                                                                          |
| **Behaviour** | Sets the file stream to use a memory pool.  The `BufLen` is set to `poolSize`, and  `BufInc` is set to `sizeInc`  |
| **Side effects**| Allocates memory from the specified pool.                                                   |
| **Errors / exceptions**| None                                                                                           |

#### `takeFilePool(pstr & filePool, word & fileSize)`

| Field        | Detail                                                                        |
|--------------|------------------------------------------------------------------------------|
| **Summary**  | Copies the memory pool associated with the file to the provided pointers.          |
| **Parameters**| `filePool` — Pointer to a buffer where the memory pool will be copied.  `fileSize` — Pointer to store the size of the memory pool. |
| **Returns**  | `void`                                                                         |
| **Behaviour** | Copies the internal `IOBuf` (memory pool) to the provided `filePool` pointer.             |
| **Side effects**| Copies memory from one location to another.                                 |
| **Errors / exceptions**| None                                                                           |

#### `inline GetFilePool(pstr & filePool, word & fileSize)`

| Field        | Detail                                                                                        |
|--------------|----------------------------------------------------------------------------------------------|
| **Summary**  | Provides a quick way to get the file pool and its size.                                      |
| **Parameters**| `filePool` — Pointer to a buffer where the memory pool will be copied.  `fileSize` — Pointer to store the size of the memory pool. |
| **Returns**  | `void`                                                                                         |
| **Behaviour** |  Same as `takeFilePool`                                                                  |
| **Side effects**| Copies memory from one location to another.                                                 |
| **Errors / exceptions**| None                                                                                        |

#### `truncate(long size)`

| Field        | Detail                                                                                                     |
|--------------|------------------------------------------------------------------------------------------------------------|
| **Summary**  | Truncates the file to the specified size.  This is called *before* resetting or appending to the file.   |
| **Parameters**| `size` — The new size of the file (in bytes).                                                              |
| **Returns**  | `void`                                                                                                       |
| **Behaviour** | Calls the `truncate` system call (or equivalent) to resize the file.                                      |
| **Side effects**| Resizes the file on the filesystem.                                                                       |
| **Errors / exceptions**| None                                                                                                     |

#### `reset(bool ReadOnly, int retry)`

| Field        | Detail                                                                                                             |
|--------------|--------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Resets the file stream to its initial state, ready for reading or writing.                                          |
| **Parameters**| `ReadOnly` — If `true`, the file is treated as read-only. `retry` —  Number of retries to try to open the file.                   |
| **Returns**  | `bool` — `true` if the reset was successful, `false` otherwise (e.g., if the file could not be opened).                 |
| **Behaviour** | Closes the current file handle, frees the buffer, and attempts to re-open the file.                               |
| **Side effects**| Closes the current file handle.                                                                                |
| **Errors / exceptions**| None                                                                                                           |

#### `erase()`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Erases the file from disk.                                                                         |
| **Parameters**| `void`                                                                                             |
| **Returns**  | `bool` — `true` if the file was successfully erased, `false` otherwise.                               |
| **Behaviour** | Calls the `remove()` system call (or equivalent) to delete the file.                             |
| **Side effects**| Deletes the file from the filesystem.                                                              |
| **Errors / exceptions**| None                                                                                           |

#### `exists()`

| Field        | Detail                                                                                                     |
|--------------|------------------------------------------------------------------------------------------------------------|
| **Summary**  | Checks if the file exists on disk.                                                                         |
| **Parameters**| `void`                                                                                                     |
| **Returns**  | `bool` — `true` if the file exists, `false` otherwise.                                                      |
| **Behaviour** | Attempts to open the file to check if it exists.                                                            |
| **Side effects**| None                                                                                                       |
| **Errors / exceptions**| None                                                                                                     |

#### `parse(cpstr FileName)`

| Field        | Detail                                                                                                      |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Attempts to open a file with the given name.                                                                 |
| **Parameters**| `FileName` — The name of the file to open.                                                                 |
| **Returns**  | `true`                                                                                                        |
| **Behaviour** | Opens the file using the `assign` method.  This is a convenience function, effectively calling `assign` with default arguments. |
| **Side effects**| Opens the file stream.                                                                                   |
| **Errors / exceptions**| None                                                                                                     |

#### `rename(cpstr NewFileName)`

| Field        | Detail                                                                                                      |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Renames the file to the specified new name.                                                                  |
| **Parameters**| `NewFileName` — The new name of the file.                                                                 |
| **Returns**  | `bool` — `true` if the rename was successful, `false` otherwise.                                             |
| **Behaviour** | Closes the current file, then attempts to rename the file to the new name.                                    |
| **Side effects**| Renames the file on the filesystem.                                                                       |
| **Errors / exceptions**| None                                                                                                     |

#### `rewrite()`

| Field        | Detail                                                                                                        |
|--------------|---------------------------------------------------------------------------------------------------------------|
| **Summary**  | Opens or recreates the file for writing.                                                                       |
| **Parameters**| `void`                                                                                                         |
| **Returns**  | `bool` — `true` if the file was successfully opened for writing, `false` otherwise.                             |
| **Behaviour** | Closes any existing file handle and attempts to open the file in write mode.                                   |
| **Side effects**| Closes the existing file handle, and opens a new file handle.                                               |
| **Errors / exceptions**| None                                                                                                       |

#### `append()`

| Field        | Detail                                                                                                          |
|--------------|-----------------------------------------------------------------------------------------------------------------|
| **Summary**  | Opens or recreates the file for appending data to the end.                                                      |
| **Parameters**| `void`                                                                                                         |
| **Returns**  | `bool` — `true` if the file was successfully opened for appending, `false` otherwise.                               |
| **Behaviour** | Closes any existing file handle and attempts to open the file in append mode.                                  |
| **Side effects**| Closes the existing file handle, and opens a new file handle.                                                  |
| **Errors / exceptions**| None                                                                                                       |

#### `isOpen()`

| Field        | Detail                                                                                          |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Checks if the file stream is currently open.                                                         |
| **Parameters**| `void`                                                                                             |
| **Returns**  | `bool` — `true` if the file stream is open, `false` otherwise.                                       |
| **Behaviour** | Returns the value of the `EofFile` member variable, indicating whether the file is open.                |
| **Side effects**| None                                                                                               |
| **Errors / exceptions**| None                                                                                             |

#### `Position()`

| Field        | Detail                                                                                 |
|--------------|---------------------------------------------------------------------------------------|
| **Summary**  | Returns the current file position (cursor position).                                    |
| **Parameters**| `void`                                                                                 |
| **Returns**  | `long` — The current file position.                                                      |
| **Behaviour** | Returns the current position of the file pointer.                                         |
| **Side effects**| None                                                                                  |
| **Errors / exceptions**| None                                                                                |

#### `seek(long Position)`

| Field        | Detail                                                                                                   |
|--------------|---------------------------------------------------------------------------------------------------------|
| **Summary**  | Moves the file pointer to a specific position in the file.                                            |
| **Parameters**| `Position` — The position to move the file pointer to (in bytes).                                      |
| **Returns**  | `bool` — `true` if the seek was successful, `false` otherwise.                                          |
| **Behaviour** | Calls the `fseek()` system call (or equivalent) to move the file pointer to the specified position. |
| **Side effects**| Moves the file pointer.                                                                            |
| **Errors / exceptions**| None                                                                                             |

#### `FileEnd()`

| Field        | Detail                                                                                       |
|--------------|---------------------------------------------------------------------------------------------|
| **Summary**  | Checks if the end of the file has been reached.                                               |
| **Parameters**| `void`                                                                                       |
| **Returns**  | `bool` — `true` if the end of the file has been reached, `false` otherwise.                       |
| **Behaviour** | Returns the value of the `EofFile` member variable.                                           |
| **Side effects**| None                                                                                        |
| **Errors / exceptions**| None                                                                                      |

#### `flush()`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Flushes the output buffer, ensuring that any pending data is written to the file.                    |
| **Parameters**| `void`                                                                                             |
| **Returns**  | `void`                                                                                              |
| **Behaviour** | Calls the `fflush()` system call (or equivalent) to flush the output buffer.                        |
| **Side effects**| Writes any data in the output buffer to the file.                                                 |
| **Errors / exceptions**| None                                                                                           |

#### `shut()`

| Field        | Detail                                                                                                    |
|--------------|----------------------------------------------------------------------------------------------------------|
| **Summary**  | Closes the file stream and releases associated resources.                                                |
| **Parameters**| `void`                                                                                                   |
| **Returns**  | `void`                                                                                                    |
| **Behaviour** | Closes the file handle, frees the allocated buffer, and resets internal state variables to their initial values.|
| **Side effects**| Closes the file handle.                                                                                |
| **Errors / exceptions**| None                                                                                                  |

#### `isOpen()`

| Field        | Detail                                                                                       |
|--------------|---------------------------------------------------------------------------------------------|
| **Summary**  | Checks if the file stream is currently open.                                                   |
| **Parameters**| `void`                                                                                       |
| **Returns**  | `bool` — `true` if the file stream is open, `false` otherwise.                                |
| **Behaviour** | Returns the value of the `EofFile` member variable, indicating whether the file is open.        |
| **Side effects**| None                                                                                        |
| **Errors / exceptions**| None                                                                                      |

#### `ReadFile(void * Buffer, word Count)`

| Field        | Detail                                                                                              |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Reads data from the file into the provided buffer.                                                    |
| **Parameters**| `Buffer` — The buffer to write the data into. `Count` — The number of bytes to read.                 |
| **Returns**  | `word` — The number of bytes actually read (which may be less than `Count`).                         |
| **Behaviour** | Reads data from the file stream and writes it to the specified buffer.                               |
| **Side effects**| None                                                                                               |
| **Errors / exceptions**| None                                                                                           |

#### `CreateRead(pstr & Line)`

| Field        | Detail                                                                                                         |
|--------------|---------------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads a line from the file into the provided string buffer.                                                    |
| **Parameters**| `Line` — A pointer to a string buffer where the line will be stored.                                              |
| **Returns**  | `word` — The number of characters in the line read.                                                             |
| **Behaviour** | Reads a line from the file and stores it in the provided buffer.                                                 |
| **Side effects**| None                                                                                                       |
| **Errors / exceptions**| None                                                                                                     |

#### `ReadTerLine(pstr Line, bool longLine)`

| Field        | Detail                                                                                                        |
|--------------|--------------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads a line from the file into the provided string buffer.  Handles both short and long lines.                  |
| **Parameters**| `Line` — A pointer to a string buffer where the line will be stored. `longLine` — If `true`, read long lines. |
| **Returns**  | `word` — The number of characters in the line read.                                                            |
| **Behaviour** | Reads a line from the file, handling carriage return characters appropriately.                               |
| **Side effects**| None                                                                                                       |
| **Errors / exceptions**| None                                                                                                     |

#### `WriteFile(const void * Buffer, word Count)`

| Field        | Detail                                                                                                      |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Writes data from the provided buffer to the file.                                                             |
| **Parameters**| `Buffer` — The buffer containing the data to be written. `Count` — The number of bytes to write.              |
| **Returns**  | `bool` — `true` if the write was successful, `false` otherwise.                                               |
| **Behaviour** | Writes data from the buffer to the file stream.                                                              |
| **Side effects**| None                                                                                                       |
| **Errors / exceptions**| None                                                                                                     |

#### `CreateWrite(cpstr Line)`

| Field        | Detail                                                                                        |
|--------------|----------------------------------------------------------------------------------------------|
| **Summary**  | Creates a new line in the file.                                                             |
| **Parameters**| `Line` — The line to write.                                                                 |
| **Returns**  | `bool` — `true` if the line was successfully written, `false` otherwise.                        |
| **Behaviour** | Writes the given string to the file, adding a line separator (LF) after the string.            |
| **Side effects**| None                                                                                        |
| **Errors / exceptions**| None                                                                                      |

#### `WriteTerLine(cpstr Line, bool longLine)`

| Field        | Detail                                                                                                      |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Writes a line to the file, adding a line separator (LF) at the end.                                       |
| **Parameters**| `Line` — The line to write. `longLine` — If `true`, write the line as a long line (no LF).                  |
| **Returns**  | `bool` — `true` if the write was successful, `false` otherwise.                                            |
| **Behaviour** | Writes the given string to the file, adding a line separator (LF) after the string.                         |
| **Side effects**| None                                                                                                       |
| **Errors / exceptions**| None                                                                                                     |

#### `ReadTerLine(pstr Line, bool longLine)`

| Field        | Detail                                                                                          |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Reads a line from the file into the provided string buffer.  Handles both short and long lines.   |
| **Parameters**| `Line` — A pointer to a string buffer where the line will be stored. `longLine` — If `true`, read long lines. |
| **Returns**  | `word` — The number of characters in the line read.                                                  |
| **Behaviour** | Reads a line from the file, handling carriage return characters appropriately.                      |
| **Side effects**| None                                                                                            |
| **Errors / exceptions**| None                                                                                          |

#### `WriteDataLine(realtype X, realtype Y, int length)`

| Field        | Detail                                                                                       |
|--------------|---------------------------------------------------------------------------------------------|
| **Summary**  | Writes a data line with specified number of columns to the file.                             |
| **Parameters**| `X` — The first value in the line. `Y` — The second value in the line. `length` — Column length.|
| **Returns**  | `bool` — `true` if the write was successful, `false` otherwise.                               |
| **Behaviour** | Writes the values of `X` and `Y` to the file, separated by spaces, and adds a line separator (LF).|
| **Side effects**| None                                                                                        |
| **Errors / exceptions**| None                                                                                      |

#### `WriteParameter(cpstr S, realtype X, int ParColumn, int length)`

| Field        | Detail                                                                                          |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Writes a parameter to the file with the specified name and value.                               |
| **Parameters**| `S` — The name of the parameter. `X` — The value of the parameter. `ParColumn` — Column in which the value should be written. `length` – Length of the parameter. |
| **Returns**  | `bool` — `true` if the write was successful, `false` otherwise.                                 |
| **Behaviour** | Writes the value of the parameter to the file, including the parameter name.                     |
| **Side effects**| None                                                                                            |
| **Errors / exceptions**| None                                                                                          |

#### `WriteParameters(cpstr S, int n_X, rvector X, int ParColumn, int length)`

| Field        | Detail                                                                                      |
|--------------|--------------------------------------------------------------------------------------------|
| **Summary**  | Writes multiple parameters to the file, separated by commas.                                |
| **Parameters**| `S` — The name of the parameter. `n_X` — The number of parameters to write. `X` — The vector of parameter values. `ParColumn` — Column in which the value should be written. `length` – Length of the parameter. |
| **Returns**  | `bool` — `true` if the write was successful, `false` otherwise.                               |
| **Behaviour** | Writes the values of the parameters to the file, separated by commas, including parameter names.|
| **Side effects**| None                                                                                        |
| **Errors / exceptions**| None                                                                                      |

#### `ReadParameter(pstr S, realtype & X, int ParColumn)`

| Field        | Detail                                                                                          |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Reads a parameter from the file into the provided string buffer.                                 |
| **Parameters**| `S` — A pointer to a string buffer where the parameter name will be stored. `X` — A reference to the realtype variable to store the parameter value. `ParColumn` — Column in which the value should be read. |
| **Returns**  | `bool` — `true` if the read was successful, `false` otherwise.                                  |
| **Behaviour** | Reads a parameter name and value from the file and stores the value in the `X` variable.           |
| **Side effects**| None                                                                                            |
| **Errors / exceptions**| None                                                                                          |

#### `ReadParameters(pstr S, int & n_X, rvector X, int MaxLen, int ParColumn)`

| Field        | Detail                                                                                        |
|--------------|----------------------------------------------------------------------------------------------|
| **Summary**  | Reads multiple parameters from the file into the provided vector.                             |
| **Parameters**| `S` — A pointer to a string buffer where the parameter name will be stored. `n_X` — A reference to an integer variable to store the number of parameters read. `X` — A pointer to a vector to store the parameter values. `MaxLen` – Maximum length of a line. `ParColumn` — Column in which the value should be read. |
| **Returns**  | `bool` — `true` if the read was successful, `false` otherwise.                               |
| **Behaviour** | Reads multiple parameter names and values from the file and stores them in the `X` vector.      |
| **Side effects**| None                                                                                        |
| **Errors / exceptions**| None                                                                                      |

#### `ReadParameter(pstr S, int & X, int ParColumn)`

| Field        | Detail                                                                                        |
|--------------|----------------------------------------------------------------------------------------------|
| **Summary**  | Reads a parameter from the file into the provided integer variable.                             |
| **Parameters**| `S` — A pointer to a string buffer where the parameter name will be stored. `X` — A reference to the integer variable to store the parameter value. `ParColumn` — Column in which the value should be read. |
| **Returns**  | `bool` — `true` if the read was successful, `false` otherwise.                               |
| **Behaviour** | Reads a parameter name and value from the file and stores the value in the `X` variable.           |
| **Side effects**| None                                                                                        |
| **Errors / exceptions**| None                                                                                      |

#### `CreateReadVector(rvector & V, int & len, int Shift)`

| Field        | Detail                                                                                              |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Creates a new vector from a given file.                                                              |
| **Parameters**| `V` — A reference to the vector to store the data. `len` — A reference to an integer variable to store the length of the vector. `Shift` — The starting index of the vector. |
| **Returns**  | `bool` — `true` if the read was successful, `false` otherwise.                                  |
| **Behaviour** | Reads data from the file and stores it in the specified vector.                                  |
| **Side effects**| None                                                                                            |
| **Errors / exceptions**| None                                                                                          |

#### `CreateReadVector(ivector & iV, int & len, int Shift)`

| Field        | Detail                                                                                              |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Creates a new integer vector from a given file.                                                     |
| **Parameters**| `iV` — A reference to the integer vector to store the data. `len` — A reference to an integer variable to store the length of the vector. `Shift` — The starting index of the vector. |
| **Returns**  | `bool` — `true` if the read was successful, `false` otherwise.                                  |
| **Behaviour** | Reads integer data from the file and stores it in the specified vector.                              |
| **Side effects**| None                                                                                            |
| **Errors / exceptions**| None                                                                                          |

#### `CreateReadVector(lvector & lV, int & len, int Shift)`

| Field        | Detail                                                                                              |
|--------------|------------------------------------------------29-------------------------------------------------------------|
| **Summary**  | Creates a new long vector from a given file.                                                        |
| **Parameters**| `lV` — A reference to the long vector to store the data. `len` — A reference to an integer variable to store the length of the vector. `Shift` — The starting index of the vector. |
| **Returns**  | `bool` — `true` if the read was successful, `false` otherwise.                                  |
| **Behaviour** | Reads long data from the file and stores it in the specified vector.                               |
| **Side effects**| None                                                                                            |
| **Errors / exceptions**| None                                                                                          |

#### `CreateReadVector(bvector & B, int & len, int Shift)`

| Field        | Detail                                                                                              |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Creates a new byte vector from a given file.                                                         |
| **Parameters**| `B` — A reference to the byte vector to store the data. `len` — A reference to an integer variable to store the length of the vector. `Shift` — The starting index of the vector. |
| **Returns**  | `bool` — `true` if the read was successful, `false` otherwise.                                  |
| **Behaviour** | Reads byte data from the file and stores it in the specified vector.                                |
| **Side effects**| None                                                                                            |
| **Errors / exceptions**| None                                                                                          |

#### `WriteVector(rvector V, int len, int Shift)`

| Field        | Detail                                                                                           |
|--------------|-------------------------------------------------------------------------------------------------|
| **Summary**  | Writes a vector to a file.                                                                        |
| **Parameters**| `V` — The vector to write. `len` — The length of the vector. `Shift` — The starting index of the vector.|
| **Returns**  | `bool` — `true` if the write was successful, `false` otherwise.                                 |
| **Behaviour** | Writes the values of the vector to the file.                                                      |
| **Side effects**| None                                                                                            |
| **Errors / exceptions**| None                                                                                          |

#### `WriteVector(ivector iV, int len, int Shift)`

| Field        | Detail                                                                                           |
|--------------|-------------------------------------------------------------------------------------------------|
| **Summary**  | Writes a integer vector to a file.                                                               |
| **Parameters**| `iV` — The vector to write. `len` — The length of the vector. `Shift` — The starting index of the vector.|
| **Returns**  | `bool` — `true` if the write was successful, `false` otherwise.                                 |
| **Behaviour** | Writes the values of the integer vector to the file.                                             |
| **Side effects**| None                                                                                            |
| **Errors / exceptions**| None                                                                                          |

#### `WriteVector(lvector lV, int len, int Shift)`

| Field        | Detail                                                                                           |
|--------------|-------------------------------------------------------------------------------------------------|
| **Summary**  | Writes a long vector to a file.                                                                  |
| **Parameters**| `lV` — The vector to write. `len` — The length of the vector. `Shift` — The starting index of the vector.|
| **Returns**  | `bool` — `true` if the write was successful, `false` otherwise.                                 |
| **Behaviour** | Writes the values of the long vector to the file.                                                 |
| **Side effects**| None                                                                                            |
| **Errors / exceptions**| None                                                                                          |

#### `WriteVector(bvector B, int len, int Shift)`

| Field        | Detail                                                                                           |
|--------------|-------------------------------------------------------------------------------------------------|
| **Summary**  | Writes a byte vector to a file.                                                                  |
| **Parameters**| `B` — The vector to write. `len` — The length of the vector. `Shift` — The starting index of the vector.|
| **Returns**  | `bool` — `true` if the write was successful, `false` otherwise.                                 |
| **Behaviour** | Writes the values of the byte vector to the file.                                                 |
| **Side effects**| None                                                                                            |
| **Errors / exceptions**| None                                                                                          |

#### `ReadVector(rvector V, int maxlen, int Shift)`

| Field        | Detail                                                                                                 |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads data from the file into a `rvector`.                                                            |
| **Parameters**| `V` — The `rvector` to store the data. `maxlen` — Maximum number of elements to read.  `Shift` — The starting index of the vector. |
| **Returns**  | `bool` — `true` if the read was successful, `false` otherwise.                                         |
| **Behaviour** | Reads data from the file and stores it in the specified `rvector`.                                   |
| **Side effects**| None                                                                                                |
| **Errors / exceptions**| None                                                                                              |

#### `ReadVector(ivector iV, int maxlen, int Shift)`

| Field        | Detail                                                                                                 |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads data from the file into an `ivector`.                                                            |
| **Parameters**| `iV` — The `ivector` to store the data. `maxlen` — Maximum number of elements to read.  `Shift` — The starting index of the vector. |
| **Returns**  | `bool` — `true` if the read was successful, `false` otherwise.                                         |
| **Behaviour** | Reads integer data from the file and stores it in the specified `ivector`.                              |
| **Side effects**| None                                                                                                |
| **Errors / exceptions**| None                                                                                              |

#### `ReadVector(lvector lV, int maxlen, int Shift)`

| Field        | Detail                                                                                                 |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads data from the file into a `lvector`.                                                            |
| **Parameters**| `lV` — The `lvector` to store the data. `maxlen` — Maximum number of elements to read.  `Shift` — The starting index of the vector. |
| **Returns**  | `bool` — `true` if the read was successful, `false` otherwise.                                         |
| **Behaviour** | Reads long data from the file and stores it in the specified `lvector`.                                |
| **Side effects**| None                                                                                                |
| **Errors / exceptions**| None                                                                                              |

#### `ReadVector(bvector B, int maxlen, int Shift)`

| Field        | Detail                                                                                           |
|--------------|-------------------------------------------------------------------------------------------------|
| **Summary**  | Reads data from the file into a `bvector`.                                                         |
| **Parameters**| `B` — The `bvector` to store the data. `maxlen` — Maximum number of elements to read.  `Shift` — The starting index of the vector. |
| **Returns**  | `bool` — `true` if the read was successful, `false` otherwise.                                 |
| **Behaviour** | Reads byte data from the file and stores it in the specified `bvector`.                            |
| **Side effects**| None                                                                                            |
| **Errors / exceptions**| None                                                                                          |

#### `CreateReadVector(rvector & V, int & len, int Shift)`

| Field        | Detail                                                                                              |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Creates a new vector from a given file.                                                              |
| **Parameters**| `V` — A reference to the vector to store the data. `len` — A reference to an integer variable to store the length of the vector. `Shift` — The starting index of the vector. |
| **Returns**  | `bool` — `true` if the read was successful, `false` otherwise.                                       |
| **Behaviour** | Reads data from the file and stores it in the specified vector.                                      |
| **Side effects**| None                                                                                            |
| **Errors / exceptions**| None                                                                                          |

#### `CreateReadVector(ivector & iV, int & len, int Shift)`

| Field        | Detail                                                                                              |
|--------------|

---

### `mmdb_io_stream`

```
---
## `mmdb::io::Stream`

> This class provides a basic streamable object for reading and writing data to a file. It is the foundation for more specialized streamable classes in the MMDB library.

### Methods

#### `Stream(RPStream)`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructs a `Stream` object, initializing it with the provided raw stream pointer.                                                    |
| **Parameters** | `RPStream Object` — A raw pointer to the underlying stream object. This allows the `Stream` class to interact with the raw stream.  |
| **Returns**  | `void`                                                                                                                            |
| **Behaviour** | This is a constructor. It allocates a new `Stream` object and initializes its internal state using the provided `RPStream` object. The `Stream` class then uses this raw stream to perform read and write operations when the `read()` and `write()` methods are called.  |
| **Side effects** | Allocates memory for the `Stream` object.                                                                                             |
| **Errors / exceptions** | `none`                                                                                                                           |

#### `~Stream()`

| Field        | Detail                                                                                                       |
|--------------|-------------------------------------------------------------------------------------------------------------|
| **Summary**  | Destroys the `Stream` object, releasing any resources it may have allocated.                                   |
| **Parameters** | `none`                                                                                                        |
| **Returns**  | `void`                                                                                                        |
| **Behaviour** | This is a destructor. It performs standard cleanup operations, likely deallocating any memory allocated by the object and releasing any resources it holds. |
| **Side effects** | Releases memory allocated by the object, potentially closing the underlying stream.                                |
| **Errors / exceptions** | `none`                                                                                                         |

#### `void read(RFile)`

| Field        | Detail                                                                                                                                                                                                                                       |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads data from the underlying stream using the provided `RFile` object.  The specific data read is determined by the internal implementation of the `Stream` class, which is not explicitly defined here. It’s assumed that this function reads data and updates the internal state of the `Stream` object. |
| **Parameters** | `RFile f` — A raw `RFile` object, providing access to the file stream for reading.                                                                                                                                                    |
| **Returns**  | `void`                                                                                                                                                                                                                                         |
| **Behaviour** | This method likely performs I/O operations on the underlying stream using the `RFile` object. It probably iterates through the data in the stream and updates the state of the `Stream` object based on the data read.  The behaviour is defined elsewhere. |
| **Side effects** | Modifies the internal state of the `Stream` object based on the data read.  Potentially writes data to the underlying file stream. |
| **Errors / exceptions** | `none`                                                                                                                                                                                                                                 |

#### `void write(RFile)`

| Field        | Detail                                                                                                                                                                                                      |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Writes data to the underlying stream using the provided `RFile` object. The specific data written is determined by the internal implementation of the `Stream` class, which is not explicitly defined here.  |
| **Parameters** | `RFile f` — A raw `RFile` object, providing access to the file stream for writing.                                                                                                                          |
| **Returns**  | `void`                                                                                                                                                                                                       |
| **Behaviour** | This method likely performs I/O operations on the underlying stream using the `RFile` object.  It probably iterates through the data and writes it to the underlying file stream. The behaviour is defined elsewhere. |
| **Side effects** | Writes data to the underlying file stream.                                                                                                                                                                 |
| **Errors / exceptions** | `none`                                                                                                                                                                                                     |
---
## Free Functions


---

### `mmdb_machine_`

## `mmdb::machine::GetMachineID`

> Returns the machine ID code for the current machine.

### Methods

#### `int GetMachineID()`

| Field        | Detail                                                                  |
|--------------|--------------------------------------------------------------------------|
| **Summary**  | Returns the machine ID code based on the system configuration.           |
| **Parameters** | `none` — No parameters are required.                                        |
| **Returns**  | An `int` representing the machine ID.  Returns `MACHINE_UNKNOWN` if the machine type is not recognized. |
| **Behaviour** | The method uses a `switch` statement to determine the appropriate machine ID based on the environment. The `CALL_LIKE_*` defines are used to identify the system. |
| **Side effects** | `none`                                                                   |
| **Errors / exceptions** | `none`                                                                   |

#### `cpstr GetMachineName(int MachineID)`

| Field        | Detail                                                                 |
|--------------|-------------------------------------------------------------------------|
| **Summary**  | Returns the name of the machine based on the given machine ID.          |
| **Parameters** | `MachineID` — `int` representing the machine ID.                       |
| **Returns**  | A `cpstr` (character pointer string) containing the machine's name. Returns "Unidentified machine" if the ID is unknown. |
| **Behaviour** | The method uses a `switch` statement to determine the appropriate machine name based on the `MachineID` parameter. It then returns a `cpstr` object containing the string. |
| **Side effects** | `none`                                                                   |
| **Errors / exceptions** | `none`                                                                   |


---

### `mmdb_manager`

```cpp
//  $Id: mmdb_manager.h $
//  =================================================================
//
//   CCP4 Coordinate Library: support of coordinate-related
//   functionality in protein crystallography applications.
//
//   Copyright (C) Eugene Krissinel 2000-2013.
//
//    This library is free software: you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License version 3, modified in accordance with the provisions
//    of the license to address the requirements of UK law.
//
//    You should have received a copy of the modified GNU Lesser
//    General Public License along with this library. If not, copies
//    may be downloaded from http://www.ccp4.ac.uk/ccp4license.php
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//  =================================================================
//
//    15.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  mmdb_manager <interface>
//       ~~~~~~~~~
//       Project :  MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//  **** Classes :  mmdb::Manager  ( MMDB file manager )
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2000-2013
//
//  =================================================================
//

#ifndef __MMDB_Manager__
#define __MMDB_Manager__

#include "mmdb_bondmngr.h"

namespace mmdb  {

  // =======================  Manager  ===========================

  // copy masks
  enum COPY_MASK  {
    MMDBFCM_None        = 0x00000000,
    MMDBFCM_All         = 0xFFFFFFFF,
    MMDBFCM_Title       = 0x00000001,
    MMDBFCM_TitleKeepBM = 0x00000002,
    MMDBFCM_Cryst       = 0x00000004,
    MMDBFCM_Coord       = 0x00000008,
    MMDBFCM_SecStruct   = 0x00000010,
    MMDBFCM_HetInfo     = 0x00000020,
    MMDBFCM_Links       = 0x00000040,
    MMDBFCM_CisPeps     = 0x00000080,
    MMDBFCM_SA          = 0x00000100,
    MMDBFCM_SB          = 0x00000200,
    MMDBFCM_SC          = 0x00000400,
    MMDBFCM_Footnotes   = 0x00000800,
    MMDBFCM_ChainAnnot  = 0x00001000,
    MMDBFCM_Flags       = 0x00002000,
    MMDBFCM_Buffer      = 0x80000000,
    MMDBFCM_Top         = 0xFFFFFFF7
  };

  DefineStreamFunctions(Manager);

  class MMDB_DL_EXPORT Manager : public BondManager  {

    public :

      Manager ();
      Manager ( io::RPStream Object );
      ~Manager();


      //  ---------------  Copying/Deleting  -----------------------

      //   Copy(..) will transfer different sort of information
      // between two MMDB's according to the copy mask given
      // (cf. MMDBFCM_XXXXX values). Note that the copying content
      // replaces the corresponding information (e.g. copying
      // coordinates will replace existing coordinates rather than
      // add to them).
      void  Copy   ( PManager MMDB, COPY_MASK CopyMask );

      //   Delete(..) deletes different sort of information from
      // the MMDB according to the delete mask given.
      void  Delete ( word DelMask );  // DelMask is the same as CopyMask

      PTitleContainer GetRemarks();
      PTitleContainer GetJournal();

      realtype GetResolution(); // -1.0 means no resolution record in file

      int   ParseBiomolecules(); // returns the number of biomolecules,
                                 // -2 for general format error
                                 // -3 for errors in BIOMT records
      int   GetNofBiomolecules();
      void  GetBiomolecules   ( PPBiomolecule & BM, int & nBMs );

      PBiomolecule GetBiomolecule ( int bmNo ); // bmno=0,1,..
      PManager MakeBiomolecule ( int bmNo, int modelNo=1 );


    protected :

      //  ---------------  Stream I/O  -----------------------------
      void  write  ( io::RFile f );
      void  read   ( io::RFile f );

  };

}  // namespace mmdb

#endif


```

### Source: `mmdb_manager.cpp`

```cpp
//  $Id: mmdb_manager.cpp $
//  =================================================================
//
//   CCP4 Coordinate Library: support of coordinate-related
//   functionality in protein crystallography applications.
//
//   Copyright (C) Eugene Krissinel 2000-2013.
//
//    This library is free software: you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License version 3, modified in accordance with the provisions
//    of the license to address the requirements of UK law.
//
//    You should have received a copy of the modified GNU Lesser
//    General Public License along with this library. If not, copies
//    may be downloaded from http://www.ccp4.ac.uk/ccp4license.php
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//  =================================================================
//
//    15.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  mmdb_manager <implementation>
//       ~~~~~~~~~
//       Project :  MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//  **** Classes :  mmdb::Manager  ( MMDB file manager )
//       ~~~~~~~~~
//
//  (C) E. Krissinel 2000-2013
//
//  =================================================================
//

#include <string.h>

#include "mmdb_manager.h"

namespace mmdb  {

//  =====================   Manager   =======================

  Manager::Manager() : BondManager()  {
  }

  Manager::Manager ( io::RPStream Object ) : BondManager(Object)  {
  }

  Manager::~Manager()  {}

  void  Manager::Copy ( PManager MMDB, COPY_MASK CopyMask )  {
    PModel  mdl;
    PPChain chain;
    PChain  chain_tmp;
    ChainID chID;
    int     i,j, nchains, nAtoms, nModels, nm;
    byte    Version = 1;

    if (CopyMask & MMDBFCM_Flags)  Flags = MMDB->Flags;

    if (CopyMask & MMDBFCM_Title)  title.Copy ( &(MMDB->title) );
    if (CopyMask & MMDBFCM_TitleKeepBM)  title.FreeMemory ( true );
    if (CopyMask & MMDBFCM_Cryst)  cryst.Copy ( &(MMDB->cryst) );

    if (CopyMask & MMDBFCM_Coord)  {
      FreeCoordMemory    ();
      DeleteAllSelections();

      nAtoms = MMDB->nAtoms;
      if (nAtoms>0)  {
        atom = new PAtom[nAtoms];
        for (i=0;i<nAtoms;i++)  {
          if (MMDB->atom[i])  {
            atom[i] = newAtom();
            atom[i]->Copy ( MMDB->atom[i] );
            // the internal atom references are installed
            // by residue classes when they are read in
            // model->chain below
            atom[i]->SetAtomIndex ( i+1 );
          } else
            atom[i] = NULL;
        }
      }

      nModels = MMDB->nModels;
      if (nModels>0)  {
        model = new PModel[nModels];
        for (i=0;i<nModels;i++)  {
          if (MMDB->model[i])  {
            model[i] = newModel();
            model[i]->SetMMDBManager ( this,0 );
            model[i]->_copy ( MMDB->model[i] );
          } else
            model[i] = NULL;
        }
      }

      crModel = NULL;
      crChain = NULL;
      crRes   = NULL;

      if (MMDB->crModel)  {

        for (i=0;i<nModels;i++)
          if (MMDB->model[i])  {
            if (model[i]->serNum==MMDB->crModel->serNum)  {
              crModel = model[i];
              break;
            }
          }

        if (crModel && crModel->chain && MMDB->crChain)
          for (i=0;i<crModel->nChains;i++)
            if (crModel->chain[i])  {
              if (!strcmp(crModel->chain[i]->chainID,
                          MMDB->crChain->chainID))  {
                crChain = crModel->chain[i];
                break;
              }
            }

        if (crChain && crChain->residue && MMDB->crRes)
          for (i=0;i<crChain->nResidues;i++)
            if (crChain->residue[i])  {
              if ((!strcmp(crChain->residue[i]->name,
                           MMDB->crRes->name))                       &&
                  (crChain->residue[i]->seqNum==MMDB->crRes->seqNum) &&
                  (!strcmp(crChain->residue[i]->insCode,
                           MMDB->crRes->insCode)))  {
                crRes = crChain->residue[i];
                break;
              }
            }
      }

      /*
      if ((MMDB->nSelections>0) && MMDB->Mask)  {
        nSelections = MMDB->nSelections;
        if (nSelections>0)  {
          Mask      = new CMask [nSelections];
          SelAtom   = new PPAtom[nSelections];
          nSelAtoms = new int    [nSelections];
          for (i=0;i<nSelections;i++)  {
            Mask[i] = new CMask();
            Mask[i]->CopyMask ( MMDB->Mask[i] );
            nSelAtoms[i] = MMDB->nSelAtoms[i];
            if (nSelAtoms[i]>0)  {
              SelAtom[i] = new PAtom[nSelAtoms[i]];
              for (j=0;j<nSelAtoms[i];j++)
                SelAtom[i][j] = Atom[MMDB->SelAtom[i][j]->index];
            } else
              SelAtom[i] = NULL;
          }
        }
      }
      */

    } else if (CopyMask & (MMDBFCM_HetInfo | MMDBFCM_SecStruct |
                            MMDBFCM_Links | MMDBFCM_CisPeps |
                            MMDBFCM_ChainAnnot))  {

      for (i=0;i<MMDB->nModels;i++)
        if (MMDB->model[i])  {

          mdl = GetModel ( i+1 );
          if (!mdl)  {
            mdl = new Model( NULL,i+1 );
            AddModel ( mdl );
          }

          if (CopyMask & MMDBFCM_HetInfo)
            mdl->CopyHets ( MMDB->model[i] );
          if (CopyMask & MMDBFCM_SecStruct)
            mdl->CopySecStructure ( MMDB->model[i] );
          if (CopyMask & MMDBFCM_Links)  {
            mdl->CopyLinks  ( MMDB->model[i] );
            mdl->CopyLinkRs ( MMDB->model[i] );
          }
          if (CopyMask & MMDBFCM_CisPeps)
            mdl->CopyCisPeps ( MMDB->model[i] );
          if (CopyMask & MMDBFCM_ChainAnnot)  {
            MMDB->GetChainTable ( i+1,chain,nchains );
            for (j=0;j<nchains;j++)
              if (chain[j])  {
                chain[j]->GetChainID ( chID );
                chain[j]->SetChainID ( chID );
                ch = mdl->GetChain ( chID );
                if (!ch)  {
                  ch = new Chain();
                  ch->SetChainID ( chID );
                  mdl->AddChain ( ch );
                }
                ch->CopyAnnotations ( chain[j] );
              }

          }

        }

    }

    if (CopyMask & MMDBFCM_SA)  SA.Copy ( &(MMDB->SA) );
    if (CopyMask & MMDBFCM_SB)  SB.Copy ( &(MMDB->SB) );
    if (CopyMask & MMDBFCM_SC)  SC.Copy ( &(MMDB->SC) );
    if (CopyMask & MMDBFCM_Footnotes)
                         Footnote.Copy ( &(MMDB->Footnote) );

    if (CopyMask & MMDBFCM_Buffer)  {
      lcount = MMDB->lcount;
      strncpy ( S,MMDB->S,sizeof(S) );
    }

  }

  void  Manager::Delete ( word DelMask )  {
    PModel model;
    PPChain chain;
    int      i,j,nm, nchains;

    if (DelMask & MMDBFCM_Flags)  Flags = 0;

    if (DelMask & MMDBFCM_Title)        title.Copy ( NULL );
    if (DelMask & MMDBFCM_TitleKeepBM)  title.FreeMemory ( true );
    if (DelMask & MMDBFCM_Cryst)        cryst.Copy ( NULL );

    if (DelMask & MMDBFCM_Coord)  {
      FreeCoordMemory    ();
      DeleteAllSelections();
    }

    if (DelMask & MMDBFCM_SecStruct)  {
      GetModelTable ( model,nm );
      if (model)
        for (i=0;i<nm;i++)
          if (model[i])
            model[i]->RemoveSecStructure();
    }

    if (DelMask & MMDBFCM_HetInfo)  {
      GetModelTable ( model,nm );
      if (model)
        for (i=0;i<nm;i++)
          if (model[i])
            model[i]->RemoveHetInfo();
    }

    if (DelMask & MMDBFCM_Links)  {
      GetModelTable ( model,nm );
      if (model)
        for (i=0;i<nm;i++)
          if (model[i])  {
            model[i]->RemoveLinks ();
            model[i]->RemoveLinkRs();
          }
    }

    if (DelMask & MMDBFCM_CisPeps)  {
      GetModelTable ( model,nm );
      if (model)
        for (i=0;i<nm;i++)
          if (model[i])
            model[i]->RemoveCisPeps();
    }

    if (DelMask & MMDBFCM_ChainAnnot)  {
      nm = GetNumberOfModels();
      for (i=1;i<=nm;i++)  {
        GetChainTable ( i,chain,nchains );
        if (chain)
          for (j=0;j<nchains;j++)
            if (chain[j])
              chain[j]->FreeAnnotations();
      }
    }

    if (DelMask & MMDBFCM_SA)        SA.FreeContainer();
    if (DelMask & MMDBFCM_SB)        SB.FreeContainer();
    if (DelMask & MMDBFCM_SC)        SC.FreeContainer();
    if (DelMask & MMDBFCM_Footnotes) Footnote.FreeContainer();

    if (DelMask & MMDBFCM_Buffer)  {
      lcount = 0;
      S[0]   = char(0);
    }

  }

  PTitleContainer Manager::GetRemarks()  {
    return title.GetRemarks();
  }


  PTitleContainer Manager::GetJournal()  {
    return title.GetJournal();
  }

  realtype Manager::GetResolution()  {
    return title.GetResolution();
  }

  int Manager::ParseBiomolecules()  {
    return title.ParseBiomolecules();
  }

  int Manager::GetNofBiomolecules()  {
    return title.GetNofBiomolecules();
  }

  void Manager::GetBiomolecules ( PPBiomolecule & BM, int & nBMs )  {
    title.GetBiomolecules ( BM,nBMs );
  }

  PBiomolecule Manager::GetBiomolecule ( int bmNo )  {
    return title.GetBiomolecule ( bmNo );
  }

  PManager Manager::MakeBiomolecule ( int bmNo, int modelNo ) {
    PManager M;
    PPChain      ch;
    PChain       chain;
    PModel       model;
    PBiomolecule BM;
    int           i,j,k,n,n0,nChains;

    BM = title.GetBiomolecule ( bmNo );
    if (!BM)  return NULL;

    GetChainTable ( modelNo,ch,nChains );
    if ((!ch) || (nChains<=0))  return NULL;

    n0    = 0;
    model = new Model();

    for (i=0;(i<BM->nBMAs) && (n0>=0);i++)
      if (BM->bmApply[i])  {
        for (j=0;(j<BM->bmApply[i]->nMatrices) && (n0>=0);j++)
          for (k=0;(k<BM->bmApply[i]->nChains) && (n0>=0);k++)  {
            n0 = -1;
            for (n=0;(n<nChains) && (n0<0);n++)
              if (!strcmp(ch[n]->GetChainID(),BM->bmApply[i]->chain[k]))
                n0 = n;
            if (n0>=0)  {
              chain = new Chain();
              chain->Copy ( ch[n0] );
              chain->ApplyTransform ( BM->bmApply[i]->tm[j] );
              model->AddChain ( chain );
            }
          }
      }

    if (n0>=0)  {
      M = new Manager();
      M->AddModel ( model );
      M->PDBCleanup ( PDBCLEAN_SERIAL | PDBCLEAN_INDEX );
    } else  {
      delete model;
      M = NULL;
    }

    return M;

  }


  //  -------------------  Stream functions  ----------------------


  void  Manager::write ( io::RFile f )  {
    byte Version=1;
    f.WriteByte ( &Version );
    BondManager::write ( f );
  }

  void  Manager::read ( io::RFile f )  {
    byte Version;
    f.ReadByte ( &Version );
    BondManager::read ( f );
  }


  MakeStreamFunctions(Manager)

}  // namespace mmdb

```

---

### `mmdb_mask`

## `Mask`

> Represents an atom selection mask, allowing the user to define a selection of atoms within a molecular structure.

### Methods

#### `Mask::Mask()`

| Field        | Detail                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructs a `Mask` object, initializing its internal data structures.                                  |
| **Parameters** | none                                                                                                  |
| **Returns**  | `void`                                                                                                 |
| **Behaviour** | Initializes `mlen` to 0 and `m` to `NULL`.  Calls `InitMask()` to perform these initializations.            |
| **Side effects** | none                                                                                                  |
| **Errors / exceptions** | none                                                                                                  |

#### `Mask::Mask(io::RPStream Object)`

| Field        | Detail                                                                                                                                          |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructs a `Mask` object, initializing its internal data structures and establishing a connection to the underlying `io::RPStream` object.           |
| **Parameters** | `Object` — An `io::RPStream` object.                                                                                                              |
| **Returns**  | `void`                                                                                                                                           |
| **Behaviour** | Initializes `mlen` to 0 and `m` to `NULL`. Calls `InitMask()` to perform these initializations.  Also, calls the `io::Stream` constructor, which might perform other necessary stream initialization. |
| **Side effects** | none                                                                                                                                          |
| **Errors / exceptions** | none                                                                                                                                          |

#### `Mask::~Mask()`

| Field        | Detail                                                                                                                            |
|--------------|------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Destroys the `Mask` object, releasing the memory allocated for its internal data structures.                                      |
| **Parameters** | none                                                                                                                               |
| **Returns**  | `void`                                                                                                                              |
| **Behaviour** | Calls `ClearMask()` to deallocate the `m` vector and resets `mlen` to 0.                                                            |
| **Side effects** | `ClearMask()` releases the memory allocated for `m`.                                                                           |
| **Errors / exceptions** | none                                                                                                                              |

#### `Mask::InitMask()`

| Field        | Detail                                                                                                                                                                                              |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Initializes the internal data structures of the `Mask` object, specifically setting `mlen` to 0 and `m` to `NULL`.                                                                                |
| **Parameters** | none                                                                                                                                                                                                |
| **Returns**  | `void`                                                                                                                                                                                               |
| **Behaviour** | Sets `mlen` to 0, indicating an empty mask, and `m` to `NULL`, indicating that the internal array is not yet allocated.                                                                            |
| **Side effects** | Sets the `mlen` and `m` members of the `Mask` object.                                                                                                                                            |
| **Errors / exceptions** | none                                                                                                                                                                                          |

#### `Mask::Expand(int n)`

| Field        | Detail                                                                                                                                                                                                                           |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Expands the `m` vector to accommodate `n` more bits.                                                                                                                                                                             |
| **Parameters** | `n` — The number of additional bits to accommodate.                                                                                                                                                                             |
| **Returns**  | `void`                                                                                                                                                                                                                       |
| **Behaviour** | Allocates a new array of size `n` to store the bits. Copies the existing bits from the old `m` array to the new `m` array. If the old `m` array existed, it is deallocated. Sets `mlen` to `n`. |
| **Side effects** | Allocates memory for the `m` vector, copies the contents, and deallocates the old vector.                                                                                                                                   |
| **Errors / exceptions** | Memory allocation failures could lead to undefined behaviour.                                                                                                                                                            |

#### `Mask::NewMask(PPMask Mask, int nMasks)`

| Field        | Detail                                                                                                                                                                                                                                                                                               |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Creates a mask based on the provided `PMask` object and the number of masks to be created.                                                                                                                                                                                                           |
| **Parameters** | `Mask` — A pointer to another `PMask` object to copy from. `nMasks` — The number of masks to copy from the provided `Mask` object.                                                                                                                                                                   |
| **Returns**  | `void`                                                                                                                                                                                                                                                                                               |
| **Behaviour** | Clears the existing mask. Iterates through the input `Mask` object (presumably an array of `PMask` objects) `nMasks` times.  For each iteration, it reads the bit value from the input `Mask` and sets the corresponding bit in the `m` vector of the `Mask` object. The loop continues until `nMasks` iterations have been completed. |
| **Side effects** | Modifies the internal state of the `Mask` object (`m` and `mlen`).                                                                                                                                                                                                                                |
| **Errors / exceptions** | If the input `Mask` object is invalid (e.g., a null pointer), the behaviour is undefined.                                                                                                                                                                                                      |

#### `Mask::CopyMask(PMask Mask)`

| Field        | Detail                                                                                                                                                                                                 |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Copies the mask data from another `PMask` object to this `Mask` object.                                                                                                                                 |
| **Parameters** | `Mask` — A pointer to another `PMask` object to copy data from.                                                                                                                                      |
| **Returns**  | `void`                                                                                                                                                                                               |
| **Behaviour** | If the `mlen` of the current mask is different from the `mlen` of the input mask, the mask is cleared. Then it copies the data from the input mask into the current mask object.              |
| **Side effects** | Modifies the internal state of the `Mask` object (`m` and `mlen`).                                                                                                                                      |
| **Errors / exceptions** | none                                                                                                                                                                                            |

#### `Mask::SetMask(PMask Mask)`

| Field        | Detail                                                                                                                                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the mask to the contents of another `PMask` object.                                                                                                                                                               |
| **Parameters** | `Mask` — A pointer to another `PMask` object.                                                                                                                                                                       |
| **Returns**  | `void`                                                                                                                                                                                                             |
| **Behaviour** | Expands the mask to the length of the input `Mask`, then copies the bits from the input mask into the `m` vector.                                                                                                |
| **Side effects** | Modifies the internal state of the `Mask` object (`m` and `mlen`).                                                                                                                                                  |
| **Errors / exceptions** | none                                                                                                                                                                                                        |

#### `Mask::RemoveMask(PMask Mask)`

| Field        | Detail                                                                                                                                      |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Removes the bits set in another `PMask` object from this `Mask` object.                                                                      |
| **Parameters** | `Mask` — A pointer to another `PMask` object.                                                                                              |
| **Returns**  | `void`                                                                                                                                    |
| **Behaviour** | Iterates through the bits of the mask, and if a bit is set in both the current mask and the input mask, it is cleared from the current mask. |
| **Side effects** | Modifies the internal state of the `Mask` object (`m`).                                                                                        |
| **Errors / exceptions** | none                                                                                                                                    |

#### `Mask::SelMask(PMask Mask)`

| Field        | Detail                                                                                                                                                                        |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Selects the bits from another `PMask` object, applying a bitwise AND operation with the current mask.                                                                        |
| **Parameters** | `Mask` — A pointer to another `PMask` object.                                                                                                                              |
| **Returns**  | `void`                                                                                                                                                                      |
| **Behaviour** | Iterates through the bits of the mask, and if a bit is set in both the current mask and the input mask, the corresponding bit in the current mask is set to 1.                |
| **Side effects** | Modifies the internal state of the `Mask` object (`m`).                                                                                                                      |
| **Errors / exceptions** | none                                                                                                                                                                 |

#### `Mask::XadMask(PMask Mask)`

| Field        | Detail                                                                                                                                                                                |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | XORs the bits set in another `PMask` object with the current mask.                                                                                                                   |
| **Parameters** | `Mask` — A pointer to another `PMask` object.                                                                                                                                       |
| **Returns**  | `void`                                                                                                                                                                             |
| **Behaviour** | Iterates through the bits of the mask, and if a bit is set in both the current mask and the input mask, the corresponding bit in the current mask is toggled.                     |
| **Side effects** | Modifies the internal state of the `Mask` object (`m`).                                                                                                                               |
| **Errors / exceptions** | none                                                                                                                                                                         |

#### `Mask::ClearMask()`

| Field        | Detail                                                                                                                             |
|--------------|------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Clears the `Mask` object, releasing its memory.                                                                                   |
| **Parameters** | none                                                                                                                               |
| **Returns**  | `void`                                                                                                                              |
| **Behaviour** | Deallocates the memory associated with the `m` vector and resets `mlen` to 0.                                                           |
| **Side effects** | Deallocates the memory occupied by the `m` vector.                                                                            |
| **Errors / exceptions** | none                                                                                                                              |

#### `Mask::NegMask()`

| Field        | Detail                                                                                                                                    |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Inverts all the bits in the `Mask` object.                                                                                               |
| **Parameters** | none                                                                                                                                   |
| **Returns**  | `void`                                                                                                                                  |
| **Behaviour** | Iterates through all the bits in the `m` vector, and flips each bit from 0 to 1 or from 1 to 0.                                      |
| **Side effects** | Modifies the internal state of the `Mask` object (`m`).                                                                                |
| **Errors / exceptions** | none                                                                                                                              |

#### `Mask::CheckMask(PMask Mask)`

| Field        | Detail                                                                                                                                                           |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Checks if the specified bit is set in the current mask.                                                                                                         |
| **Parameters** | `Mask` — A pointer to another `PMask` object.                                                                                                                 |
| **Returns**  | `bool` — `true` if the bit is set, `false` otherwise.                                                                                                       |
| **Behaviour** | Iterates through the bits of the mask, and if a bit is set in both the current mask and the input mask, returns `true`.                               |
| **Side effects** | none                                                                                                                                                    |
| **Errors / exceptions** | none                                                                                                                                                    |

#### `Mask::isMask()`

| Field        | Detail                                                                                                                                                                          |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Checks if any of the bits in the current mask are set.                                                                                                                            |
| **Parameters** | none                                                                                                                                                                           |
| **Returns**  | `bool` — `true` if at least one bit is set, `false` otherwise.                                                                                                                 |
| **Behaviour** | Iterates through all the bits in the `m` vector and returns `true` if any of them are set.                                                                                        |
| **Side effects** | none                                                                                                                                                                       |
| **Errors / exceptions** | none                                                                                                                                                                    |

#### `Mask::Print(pstr S)`

| Field        | Detail                                                                                                                            |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Converts the mask's internal bit representation into a binary string and writes it to the provided string `S`.                       |
| **Parameters** | `S` — A pointer to a string buffer where the binary representation will be written.                                               |
| **Returns**  | `pstr` — A pointer to the string `S`.                                                                                             |
| **Behaviour** | Iterates through all the bits in the `m` vector and writes either '0' or '1' to the string `S` based on the bit's value. Terminates the string with a null character. |
| **Side effects** | Modifies the string `S`.                                                                                                         |
| **Errors / exceptions** | none                                                                                                                            |

#### `Mask::write(io::RFile f)`

| Field        | Detail                                                                                                                                  |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Writes the mask's data to a file stream.                                                                                                 |
| **Parameters** | `f` — An `io::RFile` object representing the output file stream.                                                                       |
| **Returns**  | `void`                                                                                                                                 |
| **Behaviour** | Writes the value of `mlen` to the file stream. Then, writes each element of the `m` vector to the file stream.                               |
| **Side effects** | Writes data to the specified file stream.                                                                                              |
| **Errors / exceptions** | I/O errors may occur.                                                                                                               |

#### `Mask::read(io::RFile f)`

| Field        | Detail                                                                                                                                  |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads the mask's data from a file stream.                                                                                               |
| **Parameters** | `f` — An `io::RFile` object representing the input file stream.                                                                         |
| **Returns**  | `void`                                                                                                                                 |
| **Behaviour** | Reads the value of `mlen` from the file stream. Then, reads each element of the `m` vector from the file stream.                               |
| **Side effects** | Reads data from the specified file stream.                                                                                              |
| **Errors / exceptions** | I/O errors may occur.                                                                                                               |

#### `MakeStreamFunctions(Mask)`

| Field        | Detail                                                                                                 |
|--------------|----------------------------------------------------------------------------------------------------------|
| **Summary**  |  Calls the `DefineStreamFunctions` function to add the necessary stream functionality to the `Mask` class. |
| **Parameters** | `Mask` — The `Mask` class object.                                                                        |
| **Returns**  | `void`                                                                                                   |
| **Behaviour** |  This function is called to ensure the `Mask` class can properly handle I/O stream operations.         |
| **Side effects** |  Adds stream functions to the `Mask` class.                                                               |
| **Errors / exceptions** | none                                                                                                   |

---


---

### `mmdb_math_`

```
## `mmdb::math`

> Provides mathematical functions used within the CCP4 Coordinate Library for processing coordinate data, particularly torsion angles and rotations.

### Methods

#### `realtype mmdb::math::GetTorsion(rvector U, rvector W, rvector V)`

| Field        | Detail                                                                                                                                                                                                               |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Calculates the torsion angle between vectors `U`, `W`, and `V`.                                                                                                                                                       |
| **Parameters** | `U` — type `rvector` representing the first vector defining the angle.  `W` — type `rvector` representing the second vector defining the angle. `V` — type `rvector` representing the third vector defining the angle. |
| **Returns**  | `realtype` representing the torsion angle in radians. Returns `NO_TORSION` (-MaxReal) if the vectors are collinear.                                                                                               |
| **Behaviour** | 1.  The method calculates the vectors `A`, `B`, and `C` using the cross product.  2.  The magnitude `Wmag` of vector `W` is calculated.  3. The cosine of the torsion angle `S` is calculated using the dot product of `C` and `W`.  4. The angle `T` is calculated using the acos function of the sum of the components of `V` scaled by `Wmag`.  5.  The torsion angle is then calculated using the `atan2` function, which correctly handles the signs of `S` and `T` to determine the correct quadrant of the angle.  6. If `S` and `T` are both zero, indicating collinearity, the function returns `NO_TORSION`. |
| **Side effects** | none                                                                                                                                                                                                                 |
| **Errors / exceptions** | none                                                                                                                                                                                                                 |

#### `realtype mmdb::math::GetAngle(rvector v1, rvector v2)`

| Field        | Detail                                                                                                                                       |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Calculates the dihedral angle between two vectors.                                                                                             |
| **Parameters** | `v1` — type `rvector` representing the first vector. `v2` — type `rvector` representing the second vector.                              |
| **Returns**  | `realtype` representing the dihedral angle in radians.                                                                                          |
| **Behaviour** | 1. Calculates the magnitudes `l1` and `l2` of the vectors `v1` and `v2` respectively.  2. If the magnitudes are zero, sets them to 1.0 to avoid division by zero. 3. Calculates the dihedral angle using the `acos` function of the normalized dot product of the two vectors. |
| **Side effects** | none                                                                                                                                           |
| **Errors / exceptions** | none                                                                                                                                           |

#### `realtype mmdb::math::Combinations(int n, int m)`

| Field        | Detail                                                                                                                                                                                          |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Calculates the binomial coefficient "n choose m", which represents the number of ways to choose `m` items from a set of `n` items.                                                              |
| **Parameters** | `n` — type `int` representing the total number of items. `m` — type `int` representing the number of items to choose.                                                                         |
| **Returns**  | `realtype` representing the binomial coefficient. Returns `0.0` if `m` is negative or greater than `n`.                                                                                           |
| **Behaviour** | 1. If `m` is less than 0 or greater than `n`, the method returns 0.0.  2. If `m` is 0 or `m` is equal to `n`, the method returns 1.0.  3. If `m` is 1 or `m` is equal to `n-1`, the method returns `n`. 4. Otherwise, it calculates the binomial coefficient iteratively using the formula:  n! / (m! * (n-m)!). It efficiently computes this using a loop and stores the precomputed results for `P[0]` to `P[n]` to avoid redundant calculations. |
| **Side effects** | none                                                                                                                                                                                             |
| **Errors / exceptions** | none                                                                                                                                                                                             |

#### `realtype mmdb::math::log1mx(realtype x)`

| Field        | Detail                                                                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Calculates the natural logarithm of (1-x) for small values of x.                                                                                             |
| **Parameters** | `x` — type `realtype` representing the value for which to calculate the logarithm.                                                                          |
| **Returns**  | `realtype` representing the natural logarithm of (1-x). Returns -MaxReal if x >= 1.0-10.0*MachEps.                                                             |
| **Behaviour** | 1. If x >= 1.0-10.0*MachEps, returns -MaxReal to avoid numerical instability. 2. Otherwise, it calculates the logarithm using a Taylor series expansion to approximate log(1-x) accurately for small x.  The loop iteratively calculates the next term in the series until the desired accuracy is achieved. |
| **Side effects** | none                                                                                                                                                            |
| **Errors / exceptions** | none                                                                                                                                                            |

#### `realtype mmdb::math::expc(realtype x)`

| Field        | Detail                                                                                                                            |
|--------------|------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Calculates 1 - exp(x) for any value of x, including very small values.                                                             |
| **Parameters** | `x` — type `realtype` representing the value for which to calculate the value.                                                    |
| **Returns**  | `realtype` representing 1 - exp(x).                                                                                             |
| **Behaviour** | 1. If x > LnMaxReal, returns -MaxReal. 2. Otherwise, it calculates 1 - exp(x) using a Taylor series expansion for exp(x). The loop iteratively calculates the next term in the series until the desired accuracy is achieved. |
| **Side effects** | none                                                                                                                              |
| **Errors / exceptions** | none                                                                                                                              |

#### `realtype mmdb::math::expc1mx(realtype x, realtype y)`

| Field        | Detail                                                                                                                                                                                          |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Calculates 1 - (1-x)**y for any x and y, including very small x and very large y, using a more stable approach.                                                                                             |
| **Parameters** | `x` — type `realtype` representing the value for which to calculate the logarithm. `y` — type `realtype` representing the exponent.                                                              |
| **Returns**  | `realtype` representing 1 - (1-x)**y.                                                                                                                                                          |
| **Behaviour** | 1. First, it calculates log(1-x) using the `log1mx` function.  2.  Then, it calculates 1 - exp(y*log(1-x)) using the `expc1mx` function which is more accurate than a direct exponential calculation, allowing the method to handle both small x and large y values correctly. |
| **Side effects** | none                                                                                                                                                                                             |
| **Errors / exceptions** | none                                                                                                                                                                                             |

## Free Functions
(No free functions)
```

---

### `mmdb_math_align`

## `mmdb::math::Alignment`

Provides functionality to align two sequences (strings) using various methods, including global, local, global/local, and free-ends alignment. It uses a dynamic programming approach to calculate the optimal alignment score.

### Methods

#### `Alignment::Alignment()`

| Field        | Detail                                                                                                        |
|--------------|---------------------------------------------------------------------------------------------------------------|
| **Summary**  | Creates a new `Alignment` object, initializing its internal state with default values for its members.          |
| **Parameters** | none                                                                                                        |
| **Returns**  | `void`                                                                                                        |
| **Behaviour** | Initializes the `Alignment` object by setting default values for its member variables, such as `Space`, `SLen`, `TLen`, `VT`, `ET`, `FT`, `AlgnS`, `AlgnT`, `AlignKey`, `VAchieved`, `SEq`, `SNEq`, `Wg`, and `Ws`. |
| **Side effects** | none                                                                                                        |
| **Errors / exceptions** | none                                                                                                        |

#### `Alignment::Alignment(io::RPStream Object)`

| Field        | Detail                                                                                                          |
|--------------|----------------------------------------------------------------------------------------------------------------|
| **Summary**  | Creates a new `Alignment` object, initializing its internal state using the provided `io::RPStream` object.          |
| **Parameters** | `Object` — An `io::RPStream` object to use for read/write operations.                                           |
| **Returns**  | `void`                                                                                                          |
| **Behaviour** | Initializes the `Alignment` object with the state of the provided `io::RPStream` object. This is similar to the default constructor but utilizes an existing stream for data. |
| **Side effects** | none                                                                                                        |
| **Errors / exceptions** | none                                                                                                        |

#### `Alignment::~Alignment()`

| Field        | Detail                                                                                                         |
|--------------|----------------------------------------------------------------------------------------------------------------|
| **Summary**  | Destroys the `Alignment` object, releasing any allocated memory.                                               |
| **Parameters** | none                                                                                                         |
| **Returns**  | `void`                                                                                                         |
| **Behaviour** | Calls the `FreeMemory()` method, which releases the dynamically allocated memory used by the `Alignment` object. This is important to prevent memory leaks. |
| **Side effects** | `FreeMemory()` is called, deallocating the memory used by the Alignment object. |
| **Errors / exceptions** | none                                                                                                        |

#### `Alignment::SetAffineModel(realtype WGap, realtype WSpace)`

| Field        | Detail                                                                                                       |
|--------------|--------------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the weights for gap and space in the affine scoring model. These weights are used during alignment calculation. |
| **Parameters** | `WGap` — The weight for a gap in the sequence.  Typically a negative value. <br>`WSpace` — The weight for a space (match) in the sequence.  Typically a positive value. |
| **Returns**  | `void`                                                                                                       |
| **Behaviour** | Updates the `Wg` (gap weight) and `Ws` (space weight) member variables with the provided values. These values are used in the scoring function to determine the similarity between aligned sequences. |
| **Side effects** | `Wg` and `Ws` are updated. |
| **Errors / exceptions** | none                                                                                                        |

#### `Alignment::SetScores(realtype SEqual, realtype SNEqual)`

| Field        | Detail                                                                                                        |
|--------------|---------------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the scores for a match and a mismatch in the scoring model.                                             |
| **Parameters** | `SEqual` — The score for a match (identical characters) in the sequence.  Typically a positive value.<br>`SNEqual` — The score for a mismatch (different characters) in the sequence.  Typically a negative value. |
| **Returns**  | `void`                                                                                                        |
| **Behaviour** | Updates the `SEq` (equal score) and `SNEq` (unequal score) member variables with the provided values. These values are used in the scoring function to determine the similarity between aligned sequences. |
| **Side effects** | `SEq` and `SNEq` are updated. |
| **Errors / exceptions** | none                                                                                                        |

#### `Alignment::Align(cpstr S, cpstr T, ALIGN_METHOD Method)`

| Field        | Detail                                                                                                                                                                             |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Performs the sequence alignment between two strings, `S` and `T`, using the specified alignment method (`Method`). This method calculates the optimal alignment score and alignment indices. |
| **Parameters** | `S` — The first sequence (string) to align. <br>`T` — The second sequence (string) to align. <br>`Method` — The alignment method to use: <br> - `ALIGN_GLOBAL`: Global alignment of `S` and `T`. <br> - `ALIGN_LOCAL` : Local alignment of `S` and `T`. <br> - `ALIGN_GLOBLOC` : Global alignment with non-penalized end gaps in `T`. <br> - `ALIGN_FREEENDS` : Global alignment with non-penalized end gaps in both `S` and `T`. |
| **Returns**  | `void`                                                                                                                                                                             |
| **Behaviour** | This is the core alignment method.  It performs a dynamic programming calculation to determine the optimal alignment between the input strings.  The implementation involves: <br> 1. **Initialization**: Setting up the initial conditions for the dynamic programming table.<br> 2. **Recurrence**:  Filling the dynamic programming table using the scoring function.<br> 3. **Backtracing**:  Tracing back through the table to determine the optimal alignment path.<br> 4. **Adjusting Ends**: Adjusting the alignment by removing unnecessary gaps at the ends.  The choice of method determines which steps are taken to complete this process. |
| **Side effects** | The `AlgnS`, `AlgnT`, `VT`, `ET`, `FT`, `VAchieved` and `AlignKey` members are modified.  The strings `S` and `T` are not modified. |
| **Errors / exceptions** |  If memory allocation fails during the alignment process, an exception may be thrown. |

#### `Alignment::GetAlignedS() `

| Field        | Detail                                                                                                        |
|--------------|---------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the aligned sequence `S` (the top string).                                                            |
| **Parameters** | none                                                                                                        |
| **Returns**  | `pstr` — Pointer to the aligned string `S`.                                                                     |
| **Behaviour** | Returns a pointer to the `AlgnS` member variable, which holds the aligned sequence `S`. |
| **Side effects** | none                                                                                                        |
| **Errors / exceptions** | none                                                                                                        |

#### `Alignment::GetAlignedT() `

| Field        | Detail                                                                                                        |
|--------------|---------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the aligned sequence `T` (the bottom string).                                                          |
| **Parameters** | none                                                                                                        |
| **Returns**  | `pstr` — Pointer to the aligned string `T`.                                                                     |
| **Behaviour** | Returns a pointer to the `AlgnT` member variable, which holds the aligned sequence `T`. |
| **Side effects** | none                                                                                                        |
| **Errors / exceptions** | none                                                                                                        |

#### `Alignment::GetScore() `

| Field        | Detail                                                                                                        |
|--------------|---------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the highest score achieved during the alignment.                                                       |
| **Parameters** | none                                                                                                        |
| **Returns**  | `realtype` — The highest score.                                                                                  |
| **Behaviour** | Returns the value of the `VAchieved` member variable, which stores the maximum alignment score. |
| **Side effects** | none                                                                                                        |
| **Errors / exceptions** | none                                                                                                        |

#### `Alignment::GetSeqId() `

| Field        | Detail                                                                                                        |
|--------------|---------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the sequence ID (primitive sequence ID) of the aligned sequences.                                         |
| **Parameters** | none                                                                                                        |
| **Returns**  | `realtype` — The sequence ID.                                                                                  |
| **Behaviour** |  Calculates the sequence ID based on the number of aligned residues and the similarity score, returning the value in `VAchieved`. |
| **Side effects** | none                                                                                                        |
| **Errors / exceptions** | none                                                                                                        |

#### `Alignment::GetNAlign() `

| Field        | Detail                                                                                                        |
|--------------|---------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the number of aligned residues in the alignment.                                                       |
| **Parameters** | none                                                                                                        |
| **Returns**  | `int` — The number of aligned residues.                                                                         |
| **Behaviour** | Returns the value of the `AlignKey` member variable, which represents the number of aligned residues based on the selected alignment method. |
| **Side effects** | none                                                                                                        |
| **Errors / exceptions** | none                                                                                                        |

#### `Alignment::OutputResults(io::RFile f, cpstr S, cpstr T)`

| Field        | Detail                                                                                                                                  |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Outputs the alignment results to a specified file.                                                                                          |
| **Parameters** | `f` — An `io::RFile` object to write the results to.<br>`S` — The first sequence (string) used in the alignment.<br>`T` — The second sequence (string) used in the alignment. |
| **Returns**  | `void`                                                                                                                                  |
| **Behaviour** | Writes the alignment details to the specified file, including the original sequences, the alignment score, the gap weights, the alignment method, and the aligned sequences.  It uses `pstr` to format the output for readability. |
| **Side effects** | Writes to the file specified by the `io::RFile` object. |
| **Errors / exceptions** |  Errors during file I/O could result in an exception. |

#### `Alignment::read(io::RFile f)`

| Field        | Detail                                                                                                        |
|--------------|---------------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads the alignment parameters from the `io::RFile` object.                                                     |
| **Parameters** | `f` — An `io::RFile` object to read the data from.                                                              |
| **Returns**  | `void`                                                                                                        |
| **Behaviour** | Reads the `Version`, `GapWeight`, `SpaceWeight`, `EqualScore`, `UnequalScore`, `AlignKey` and `SLen` from the `io::RFile` object and stores them in the corresponding member variables. |
| **Side effects** | Updates the member variables of the `Alignment` object. |
| **Errors / exceptions** |  Errors during file I/O could result in an exception. |

#### `Alignment::write(io::RFile f)`

| Field        | Detail                                                                                                        |
|--------------|---------------------------------------------------------------------------------------------------------------|
| **Summary**  | Writes the alignment parameters to the `io::RFile` object.                                                     |
| **Parameters** | `f` — An `io::RFile` object to write the data to.                                                              |
| **Returns**  | `void`                                                                                                        |
| **Behaviour** | Writes the `Version`, `GapWeight`, `SpaceWeight`, `EqualScore`, `UnequalScore`, `AlignKey` and `SLen` to the `io::RFile` object. |
| **Side effects** | Writes to the file specified by the `io::RFile` object. |
| **Errors / exceptions** |  Errors during file I/O could result in an exception. |

#### `Alignment::operator() (cpstr S, cpstr T, ALIGN_METHOD Method)`

| Field        | Detail                                                                                                        |
|--------------|---------------------------------------------------------------------------------------------------------------|
| **Summary**  | Performs the sequence alignment between two strings, `S` and `T`, using the specified alignment method (`Method`). |
| **Parameters** | `S` — The first sequence (string) to align. <br>`T` — The second sequence (string) to align. <br>`Method` — The alignment method to use: <br> - `ALIGN_GLOBAL   `: Global alignment of `S` and `T`. <br> - `ALIGN_LOCAL    `: Local alignment of `S` and `T`. <br> - `ALIGN_GLOBLOC  `: Global alignment with non-penalized end gaps in `T`. <br> - `ALIGN_FREEENDS `: Global alignment with non-penalized end gaps in both `S` and `T`. |
| **Returns**  | `Alignment` — A pointer to a new `Alignment` object containing the results of the alignment. |
| **Behaviour** |  This method performs the sequence alignment between the input strings `S` and `T` using the specified alignment method.   The function allocates a new `Alignment` object, calls the `Align` method to perform the alignment, and returns a pointer to the newly created `Alignment` object. The caller is then responsible for managing the lifetime of the returned `Alignment` object. |
| **Side effects** | The `AlgnS`, `AlgnT`, `VT`, `ET`, `FT`, `VAchieved` and `AlignKey` members of the returned `Alignment` object are modified.  The input strings `S` and `T` are not modified. |
| **Errors / exceptions** |  If memory allocation fails during the alignment process, an exception may be thrown. |

#### `Alignment::printVT(cpstr S, cpstr T)`

| Field        | Detail                                                                                                        |
|--------------|---------------------------------------------------------------------------------------------------------------|
| **Summary**  | Prints the VT matrix to the console for debugging purposes.                                                     |
| **Parameters** | `S` — The first sequence (string) to align. <br>`T` — The second sequence (string) to align.                                  |
| **Returns**  | `void`                                                                                                        |
| **Behaviour** | Prints the VT matrix (scoring table) to the standard output. |
| **Side effects** | Prints to the console. |
| **Errors / exceptions** | none                                                                                                        |

#### `Alignment::Backtrace(cpstr S, cpstr T, int J, int I, bool StopAtZero)`

| Field        | Detail                                                                                                       |
|--------------|--------------------------------------------------------------------------------------------------------------|
| **Summary**  | Backtracks through the dynamic programming table to determine the optimal alignment path.                      |
| **Parameters** | `S` — The first sequence (string) to align. <br>`T` — The second sequence (string) to align. <br>`J` — The column index in the table. <br>`I` — The row index in the table. <br>`StopAtZero` – A boolean value indicating if the algorithm should stop at a score of 0. |
| **Returns**  | `void`                                                                                                       |
| **Behaviour** | This function is the heart of the dynamic programming algorithm. It traces back through the table, making decisions at each step to maximize the score and construct the optimal alignment. The implementation considers the following: <br> 1. Finding the leading element in the table.<br> 2. Determining the score for the current cell based on the scores of its neighbors.<br> 3. Updating the alignment indices. |
| **Side effects** | The `VAchieved` member variable is updated during the backtracking process. |
| **Errors / exceptions** |  none                                                                                                        |

#### `Alignment::AdjustEnds(cpstr S, cpstr T, int J, int I)`

| Field        | Detail                                                                                                       |
|--------------|--------------------------------------------------------------------------------------------------------------|
| **Summary**  | Adjusts the alignment by removing unnecessary gaps at the ends of the aligned sequences.                       |
| **Parameters** | `S` — The first sequence (string) to align. <br>`T` — The second sequence (string) to align. <br>`J` — The column index. <br>`I` — The row index. |
| **Returns**  | `void`                                                                                                       |
| **Behaviour** | After the dynamic programming algorithm has completed, the alignment may contain unnecessary gaps at the ends of the sequences. This function removes these gaps, ensuring that the alignment is the shortest possible. |
| **Side effects** | The `AlgnS` and `AlgnT` members are modified. |
| **Errors / exceptions** | none                                                                                                        |

#### `Alignment::BuildGATable(cpstr S, cpstr T, bool FreeSEnd, bool FreeTEnd)`

| Field        | Detail                                                                                                       |
|--------------|--------------------------------------------------------------------------------------------------------------|
| **Summary**  | Builds the gap alignment table (G matrix) for global alignment.                                                  |
| **Parameters** | `S` — The first sequence (string) to align. <br>`T` — The second sequence (string) to align. <br>`FreeSEnd` — A boolean value indicating if the end gaps in `S` should be free.<br>`FreeTEnd` — A boolean value indicating if the end gaps in `T` should be free. |
| **Returns**  | `void`                                                                                                       |
| **Behaviour** | The function sets up the initial conditions for the dynamic programming algorithm, calculating the initial values of the gaps in the table. The values are calculated based on the gap weight (`Wg`) and the space weight (`Ws`). |
| **Side effects** | The `VT`, `ET`, `FT`, and `VAchieved` members are updated. |
| **Errors / exceptions** | none                                                                                                        |

#### `Alignment::BuildLATable(cpstr S, cpstr T)`

| Field        | Detail                                                                                                       |
|--------------|--------------------------------------------------------------------------------------------------------------|
| **Summary**  | Builds the local alignment table (L matrix) for local alignment.                                                 |
| **Parameters** | `S` — The first sequence (string) to align. <br>`T` — The second sequence (string) to align. |
| **Returns**  | `void`                                                                                                       |
| **Behaviour** | The function sets up the initial conditions for the dynamic programming algorithm, calculating the initial values of the gaps in the table. The values are calculated based on the gap weight (`Wg`) and the space weight (`Ws`). |
| **Side effects** | The `VT`, `ET`, `FT`, and `VAchieved` members are updated. |
| **Errors / exceptions** | none                                                                                                        |

#### `Alignment::FreeMemory()`

| Field        | Detail                                                                                                       |
|--------------|--------------------------------------------------------------------------------------------------------------|
| **Summary**  | Releases all the dynamically allocated memory used by the `Alignment` object.                                |
| **Parameters** | none                                                                                                        |
| **Returns**  | `void`                                                                                                       |
| **Behaviour** | The function releases the memory allocated for the `VT`, `ET`, `FT`, `AlgnS`, `AlgnT` and the vectors. This is important to prevent memory leaks. |
| **Side effects** | Releases the dynamically allocated memory. |
| **Errors / exceptions** | none                                                                                                        |

#### `Alignment::Score(char A, char B)`

| Field        | Detail                                                                                                       |
|--------------|--------------------------------------------------------------------------------------------------------------|
| **Summary**  | Calculates the score for a match or mismatch between two characters.                                            |
| **Parameters** | `A` — The first character. <br>`B` — The second character. |
| **Returns**  | `realtype` — The score.                                                                                        |
| **Behaviour** | This is the scoring function. It calculates the score based on the `SEq` (equal score) and `SNEq` (unequal score) parameters. |
| **Side effects** | none                                                                                                        |
| **Errors / exceptions** | none                                                                                                        |

#### `Alignment::GetSimilarity()`

| Field        | Detail                                                                                                      |
|--------------|------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the similarity score of the alignment.                                                            |
| **Parameters** | none                                                                                                      |
| **Returns**  | `realtype` — The similarity score.                                                                        |
| **Behaviour** |  The similarity score is calculated as the ratio of the aligned sequence length to the total sequence length.  |
| **Side effects** | none                                                                                                      |
| **Errors / exceptions** | none                                                                                                        |

#### `Alignment::GetSeqId()`

| Field        | Detail                                                                                                      |
|--------------|------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the sequence ID of the aligned sequences.                                                           |
| **Parameters** | none                                                                                                      |
| **Returns**  | `int` — The sequence ID.                                                                                    |
| **Behaviour** | The sequence ID is calculated as the number of aligned residues.  |
| **Side effects** | none                                                                                                      |
| **Errors / exceptions** | none                                                                                                        |

#### `Alignment::GetNAlign()`

| Field        | Detail                                                                                                       |
|--------------|--------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the number of aligned residues in the alignment.                                                    |
| **Parameters** | none                                                                                                        |
| **Returns**  | `int` — The number of aligned residues.                                                                     |
| **Behaviour** | Returns the number of aligned residues. |
| **Side effects** | none                                                                                                        |
| **Errors / exceptions** | none                                                                                                        |

## Free Functions

There are no free functions defined in this source file.

---

### `mmdb_math_bfgsmin`

## `mmdb::math::BFGSMin`

> Finds a minimum of function F(X), X is vector [1..N], defined by a virtual `MinFunc` function, providing an interface for minimization algorithms.

### Methods

#### `void BFGSMin::MinFunc(rvector X, realtype & F)`

| Field        | Detail                                                                                                                                                                                                                                                               |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Evaluates the function F at the given point X, returning the function value in `F`.  This is the core function that's called repeatedly during the minimization process. |
| **Parameters** | `X` — `rvector` — A vector representing the coordinates of the point at which to evaluate the function.  It's assumed that the vector has dimension N. |
| **Returns**  | `void` — This method doesn't return a value; it modifies the `F` parameter passed by reference to store the calculated function value.  |
| **Behaviour** | 1. Receives the coordinates `X` from the calling function (e.g., the line search or gradient descent algorithm). 2. Uses these coordinates to call the externally defined `MinFunc` virtual function, passing the data needed for the calculation (including the `X` vector and the `F` reference parameter). 3. The `MinFunc` function is responsible for performing the actual mathematical calculation of the function value `F` based on the input coordinates. |
| **Side effects** | `none` — This method doesn't modify any state externally.  It solely relies on the `MinFunc` virtual function to perform the calculation. |
| **Errors / exceptions** | `none` — This method does not throw any exceptions or return an error code.  Any errors are handled within the `MinFunc` virtual function. |

#### `void BFGSMin::MinFunc1(rvector X, realtype & F)`

| Field        | Detail                                                                                                                                                                                                                                                                                                                                                           |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Similar to `MinFunc`, but it updates the `FOpt` member if the current `F` value is the lowest function value found so far.  This is part of the driver implementation. |
| **Parameters** | `X` — `rvector` — A vector representing the coordinates of the point at which to evaluate the function.  It's assumed that the vector has dimension N. `F` — `realtype &` — A reference to a `realtype` variable to store the calculated function value. |
| **Returns**  | `void` — This method doesn't return a value; it modifies the `F` parameter passed by reference. |
| **Behaviour** | 1. Calls `MinFunc` to evaluate the function at the current point `X` and store the result in `F`. 2. Checks if the calculated `F` value is less than the current best function value (`FOpt`). 3. If the current `F` value is indeed the lowest found so far, it updates the `FOpt` member with the new minimum function value. |
| **Side effects** | `none` — This method doesn't modify any state externally.  It solely relies on the `MinFunc` virtual function to perform the calculation and the `FOpt` member variable. |
| **Errors / exceptions** | `none` — This method does not throw any exceptions or return an error code.  Any errors are handled within the `MinFunc` virtual function. |

#### `void BFGSMin::SetMinFunction(void* UserData, PBFGSMinFunc Fnc)`

| Field        | Detail                                                                                                                                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Allows the user to register a custom `MinFunc` function to be called by the `BFGSMin` class during the minimization process. |
| **Parameters** | `UserData` — `void*` — A pointer to data that will be passed to the `MinFunc` function.  This is typically used to store any information needed by the custom function.  `Fnc` — `PBFGSMinFunc` — A pointer to the function that implements the minimization logic. |
| **Returns**  | `void` — This method doesn't return a value. |
| **Behaviour** | 1. Takes a pointer `UserData` and a pointer `Fnc` to a function pointer type `PBFGSMinFunc`. 2. Stores the pointer to the `UserData` in the `MFuncData` member variable. 3. Stores the pointer to the `Fnc` function pointer in the `MFunc` member variable. |
| **Side effects** | `none` — This method doesn't modify any state externally. It simply stores references to the `MinFunc` function and its associated data. |
| **Errors / exceptions** | `none` — This method does not throw any exceptions or return an error code. |

#### `void BFGSMin::SetPrintFunction(void* UserData, PBFGSPrintFunc Fnc)`

| Field        | Detail                                                                                                                                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  |  Allows the user to register a custom `PrintFunc` function to be called by the `BFGSMin` class during the minimization process, providing output at each iteration. |
| **Parameters** | `UserData` — `void*` — A pointer to data that will be passed to the `PrintFunc` function.  This is typically used to store any information needed by the custom function.  `Fnc` — `PBFGSPrintFunc` — A pointer to the function that implements the printing logic. |
| **Returns**  | `void` — This method doesn't return a value. |
| **Behaviour** | 1. Takes a pointer `UserData` and a pointer `Fnc` to a function pointer type `PBFGSPrintFunc`. 2. Stores the pointer to the `UserData` in the `PFuncData` member variable. 3. Stores the pointer to the `Fnc` function pointer in the `PFunc` member variable. |
| **Side effects** | `none` — This method doesn't modify any state externally. It simply stores references to the `PrintFunc` function and its associated data. |
| **Errors / exceptions** | `none` — This method does not throw any exceptions or return an error code. |

#### `void BFGSMin::BFGS_Driver(int MinN, rvector x0, rvector TypX, realtype & FuncValue, int & TerminationCode, int Digits, int ItnLmt, realtype TypF, realtype GrdTol, realtype StpTol, realtype MaxStp, bool Hess, rvector LowLimit, rvector TopLimit)`

| Field        | Detail                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | The main driver function for the unconstrained BFGS minimization algorithm. It coordinates the entire minimization process, calling appropriate helper functions to perform the calculations and updating the solution (X) iteratively until a termination criterion is met.  This is the core function users call. |
| **Parameters** | `MinN` — `int` — The dimension of the problem (number of variables). `x0` — `rvector` — The initial guess for the solution vector.  `TypX` — `rvector` — The typical ranges of the X variables, used for scaling. `MaxValue` — `realtype &` — A reference to a `realtype` variable to store the minimum function value found. `TerminationCode` — `int &` — A reference to an integer variable that will be set to indicate the reason for stopping the minimization process. `Digits` — `int` — The number of decimal digits to use in the calculation of the function value. `ItnLmt` — `int` — The maximum number of iterations allowed. `TypF` — `realtype` — The expected absolute value of the function value at the minimum. `GrdTol` — `realtype` — The desired tolerance for the gradient norm. `StpTol` — `realtype` — The desired tolerance for the step size. `MaxStp` — `realtype` — The maximum step size allowed. `Hess` — `bool` — A flag indicating whether to use the Hessian matrix directly (not currently implemented). `LowLimit` — `rvector` — A vector of lower bounds for the variables. `TopLimit` — `rvector` — A vector of upper bounds for the variables. |
| **Returns**  | `void` — This method doesn't return a value. |
| **Behaviour** | 1. **Initialization**: Sets up the necessary data structures and variables for the minimization process. This includes initializing the coordinates `x0`, the typical ranges of X values, the function value, and the termination code. 2. **Iteration**: The core of the algorithm involves iteratively refining the solution by applying the BFGS update rule. This includes calculating the gradient of the function, computing the step size, and updating the solution coordinates. 3. **Termination**: The minimization process continues until a termination criterion is met, which could be based on the gradient norm, the step size, the number of iterations, or other criteria. 4. **Output**: Once the minimization process has terminated, the method stores the final solution in the `x0` vector and the function value in `MaxValue`, and sets the `TerminationCode` to indicate the reason for stopping. |
| **Side effects** | `none` — This method doesn't modify any state externally except for the function value, termination code, and solution vector, which are modified by the other helper functions. |
| **Errors / exceptions** | `none` — This method does not throw any exceptions or return an error code.  It relies on the `MinFunc` and other helper functions to handle any errors. |

#### `void BFGSMin::Print(int Itn, rvector X, rvector G, realtype F)`

| Field        | Detail                                                                                                                                                                                                                                    |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  |  Prints information about the current state of the minimization process at a given iteration. |
| **Parameters** | `Itn` — `int` — The current iteration number. `X` — `rvector` — The coordinates of the current solution.  `G` — `rvector` — The gradient of the function at the current solution.  `F` — `realtype` — The function value at the current solution. |
| **Returns**  | `void` — This method doesn't return a value. |
| **Behaviour** | 1. Receives the current iteration number `Itn`, the solution vector `X`, the gradient vector `G`, and the function value `F`. 2. Calls the registered `PrintFunc` virtual function, passing the received parameters to it.  The `PrintFunc` is responsible for formatting and printing the information to the console or other output stream. |
| **Side effects** | `none` — This method doesn't modify any state externally.  It solely relies on the `PrintFunc` virtual function to perform the output. |
| **Errors / exceptions** | `none` — This method does not throw any exceptions or return an error code. |

#### `void BFGSMin::Stop()`

| Field        | Detail                                                                                                                                                                                                                            |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Signals the driver to terminate the optimization process. |
| **Parameters** | `none` |
| **Returns**  | `void` — This method doesn't return a value. |
| **Behaviour** | Sets the `TermCode` member to `BFGS_Stopped`, indicating that the minimization process has been stopped due to the `Stop()` call. |
| **Side effects** | `none` — This method doesn't modify any state externally. |
| **Errors / exceptions** | `none` — This method does not throw any exceptions or return an error code. |

##### Constructors and Destructors

*   **`BFGSMin()`**:  Constructor that initializes all member variables to default values.
*   **`~BFGSMin()`**: Destructor that frees any allocated memory.

##### Free Functions

*   `void BFGSMin::GetMemory()`: Allocates memory for the necessary data structures.
*   `void BFGSMin::FreeMemory()`: Frees any allocated memory.
*   `void BFGSMin::Relax()`: Performs a relaxation step to improve the solution.
*   `void BFGSMin::CopyPlus(rvector x0)`: Copies the solution vector from `XPlus` to `x0`.

---

## `BFGSMin`

>  Finds a minimum of function F(X), X is vector [1..N], defined by a virtual `MinFunc` function, providing an interface for minimization algorithms.

### Methods

#### `void BFGSMin::MinFunc(rvector X, realtype & F)`

| Field        | Detail                                                                                                                                                                                                                                                               |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Evaluates the function F at the given point X, returning the function value in `F`.  This is the core function that's called repeatedly during the minimization process. |
| **Parameters** | `X` — `rvector` — A vector representing the coordinates of the point at which to evaluate the function.  It's assumed that the vector has dimension N. |
| **Returns**  | `void` — This method doesn't return a value; it modifies the `F` parameter passed by reference to store the calculated function value.  |
| **Behaviour** | 1. Receives the coordinates `X` from the calling function (e.g., the line search or gradient descent algorithm). 2. Uses these coordinates to call the externally defined `MinFunc` virtual function, passing the data needed for the calculation (including the `X` vector and the `F` reference parameter). 3. The `MinFunc` function is responsible for performing the actual mathematical calculation of the function value `F` based on the input coordinates. |
| **Side effects** | `none` — This method doesn't modify any state externally.  It solely relies on the `MinFunc` virtual function to perform the calculation. |
| **Errors / exceptions** | `none` — This method does not throw any exceptions or return an error code.  Any errors are handled within the `MinFunc` virtual function. |

#### `void BFGSMin::MinFunc1(rvector X, realtype & F)`

| Field        | Detail                                                                                                                                                                                                                                                                                                                                                           |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Similar to `MinFunc`, but it updates the `FOpt` member if the current `F` value is the lowest function value found so far.  This is part of the driver implementation. |
| **Parameters** | `X` — `rvector` — A vector representing the coordinates of the point at which to evaluate the function.  It's assumed that the vector has dimension N. `F` — `realtype &` — A reference to a `realtype` variable to store the calculated function value. |
| **Returns**  | `void` — This method doesn't return a value; it modifies the `F` parameter passed by reference. |
| **Behaviour** | 1. Calls `MinFunc` to evaluate the function at the current point `X` and store the result in `F`. 2. Checks if the calculated `F` value is less than the current best function value (`FOpt`). 3. If the current `F` value is indeed the lowest found so far, it updates the `FOpt` member with the new minimum function value. |
| **Side effects** | `none` — This method doesn't modify any state externally.  It solely relies on the `MinFunc` virtual function to perform the calculation and the `FOpt` member variable. |
| **Errors / exceptions** | `none` — This method does not throw any exceptions or return an error code. |

#### `void BFGSMin::SetMinFunction(void* UserData, PBFGSMinFunc Fnc)`

| Field        | Detail                                                                                                                                                                                                                                                                                                                                                           |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Allows the user to register a custom `MinFunc` function to be called by the `BFGSMin` class during the minimization process. |
| **Parameters** | `UserData` — `void*` — A pointer to data that will be passed to the `MinFunc` function.  This is typically used to store any information needed by the custom function.  `Fnc` — `PBFGSMinFunc` — A pointer to the function that implements the minimization logic. |
| **Returns**  | `void` — This method doesn't return a value. |
| **Behaviour** | 1. Takes a pointer `UserData` and a pointer `Fnc` to a function pointer type `PBFGSMinFunc`. 2. Stores the pointer to the `UserData` in the `MFuncData` member variable. 3. Stores the pointer to the `Fnc` function pointer in the `MFunc` member variable. |
| **Side effects** | `none` — This method doesn't modify any state externally. It simply stores references to the `MinFunc` function and its associated data. |
| **Errors / exceptions** | `none` — This method does not throw any exceptions or return an error code. |

#### `void BFGSMin::SetPrintFunction(void* UserData, PBFGSPrintFunc Fnc)`

| Field        | Detail                                                                                                                                                                                                                                                                                                                                                           |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  |  Allows the user to register a custom `PrintFunc` function to be called by the `BFGSMin` class during the minimization process, providing output at each iteration. |
| **Parameters** | `UserData` — `void*` — A pointer to data that will be passed to the `PrintFunc` function.  This is typically used to store any information needed by the custom function.  `Fnc` — `PBFGSPrintFunc` — A pointer to the function that implements the printing logic. |
| **Returns**  | `void` — This method doesn't return a value. |
| **Behaviour** | 1. Takes a pointer `UserData` and a pointer `Fnc` to a function pointer type `PBFGSPrintFunc`. 2. Stores the pointer to the `UserData` in the `PFuncData` member variable. 3. Stores the pointer to the `Fnc` function pointer in the `PFunc` member variable. |
| **Side effects** | `none` — This method doesn't modify any state externally. It simply stores references to the `PrintFunc` function and its associated data. |
| **Errors / exceptions** | `none` — This method does not throw any exceptions or return an error code. |

#### `void BFGSMin::BFGS_Driver(int MinN, rvector x0, rvector TypX, realtype & FuncValue, int & TerminationCode, int Digits, int ItnLmt, realtype TypF, realtype GrdTol, realtype StpTol, realtype MaxStp, bool Hess, rvector LowLimit, rvector TopLimit)`

| Field        | Detail                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | The main driver function for the unconstrained BFGS minimization algorithm. It coordinates the entire minimization process, calling appropriate helper functions to perform the calculations and updating the solution (X) iteratively until a termination criterion is met.  This is the core function users call. |
| **Parameters** | `MinN` — `int` — The dimension of the problem (number of variables). `x0` — `rvector` — The initial guess for the solution vector.  `TypX` — `rvector` — The typical ranges of the X variables, used for scaling. `MaxValue` — `realtype &` — A reference to a `realtype` variable to store the minimum function value found. `TerminationCode` — `int &` — A reference to an integer variable that will be set to indicate the reason for stopping the minimization process. `Digits` — `int` — The number of decimal digits to use in the calculation of the function value. `ItnLmt` — `int` — The maximum number of iterations allowed. `TypF` — `realtype` — The expected absolute value of the function value at the minimum. `GrdTol` — `realtype` — The desired tolerance for the gradient norm. `StpTol` — `realtype` — The desired tolerance for the step size. `MaxStp` — `realtype` — The maximum step size allowed. `Hess` — `bool` — A flag indicating whether to use the Hessian matrix directly (not currently implemented). `LowLimit` — `rvector` — A vector of lower bounds for the variables. `TopLimit` — `rvector` — A vector of upper bounds for the variables. |
| **Returns**  | `void` — This method doesn't return a value. |
| **Behaviour** | 1. **Initialization**: Sets up the necessary data structures and variables for the minimization process. This includes initializing the coordinates `x0`, the typical ranges of X values, the function value, and the termination code. 2. **Iteration**: The core of the algorithm involves iteratively refining the solution by applying the BFGS update rule. This includes calculating the gradient of the function, computing the step size, and updating the solution coordinates. 3. **Termination**: The minimization process continues until a termination criterion is met, which could be based on the gradient norm, the step size, the number of iterations, or other criteria. 4. **Output**: Once the minimization process has terminated, the method stores the final solution in the `x0` vector and the function value in `MaxValue`, and sets the `TerminationCode` to indicate the reason for stopping. |
| **Side effects** | `none` — This method doesn't modify any state externally except for the function value, termination code, and solution vector, which are modified by the other helper functions. |
| **Errors / exceptions** | `none` — This method does not throw any exceptions or return an error code.  It relies on the `MinFunc` and other helper functions to handle any errors. |

#### `void BFGSMin::Print(int Itn, rvector X, rvector G, realtype F)`

| Field        | Detail                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  |  Prints information about the current state of the minimization process at a given iteration. |
| **Parameters** | `Itn` — `int` — The current iteration number. `X` — `rvector` — The coordinates of the current solution.  `G` — `rvector` — The gradient of the function at the current solution.  `F` — `realtype` — The function value at the current solution. |
| **Returns**  | `void` — This method doesn't return a value. |
| **Behaviour** | 1. Receives the current iteration number `Itn`, the solution vector `X`, the gradient vector `G`, and the function value `F`. 2. Calls the registered `PrintFunc` virtual function, passing the received parameters to it.  The `PrintFunc` is responsible for formatting and printing the information to the console or other output stream. |
| **Side effects** | `none` — This method doesn't modify any state externally.  It solely relies on the `PrintFunc` virtual function to perform the output. |
| **Errors / exceptions** | `none` — This method does not throw any exceptions or return an error code. |

#### `void BFGSMin::Stop()`

| Field        | Detail                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Signals the driver to terminate the optimization process. |
| **Parameters** | `none` |
| **Returns**  | `void` — This method doesn't return a value. |
| **Behaviour** | Sets the `TermCode` member to `BFGS_Stopped`, indicating that the minimization process has been stopped due to the `Stop()` call. |
| **Side effects** | `none` — This method doesn't modify any state externally. |
| **Errors / exceptions** | `none` — This method does not throw any exceptions or return an error code. |

##### Constructors and Destructors

*   **`BFGSMin()`**:  Constructor that initializes all member variables to default values.
*   **`~BFGSMin()`**: Destructor that frees any allocated memory.

##### Free Functions

*   `void BFGSMin::GetMemory()`: Allocates memory for the necessary data structures.
*   `void BFGSMin::FreeMemory()`: Frees any allocated memory.
*   `void BFGSMin::Relax()`: Performs a relaxation step to improve the solution.
*   `void BFGSMin::CopyPlus(rvector x0)`: Copies the solution vector from `XPlus` to `x0`.

---

## Free Functions

There are no free functions in this source code.


---

### `mmdb_math_fft`

## `mmdb::math::FFT`

> Performs a discrete Fourier transform (DFT) on a 1D real signal.

### Methods

#### `void FFT(rvector data, int nn, bool Forward)`

| Field        | Detail                     |
|--------------|----------------------------|
| **Summary**  | Computes the discrete Fourier transform of a 1D real signal. |
| **Parameters** | `data` — A reference to a `rvector` containing the input signal. `nn` — The number of data points to transform (must be a power of 2). `Forward` — A boolean flag indicating whether to perform a forward (true) or inverse (false) DFT. |
| **Returns**  | `void`                     |
| **Behaviour** | 1. Allocates a temporary `rvector` called `fft` with size `2*nn`. 2. Iterates through the input data points from index 1 to 2*nn, swapping adjacent elements to perform the DFT. Uses the bit-reversal permutation to do this. 3. Applies the DFT algorithm (Radix-2 Cooley-Tukey) recursively, reducing the size of the signal by half in each step. 4. Returns the transformed data in the `fft` vector. |
| **Side effects** | Modifies the `fft` vector. The original `data` vector is not altered. |
| **Errors / exceptions** | Throws an exception if `nn` is not a power of 2.  |

#### `void RealFFT(rvector data, int n, bool Forward)`

| Field        | Detail                                                                                                                   |
|--------------|--------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Performs a real-valued discrete Fourier transform (DFT) on a 1D real signal.                                       |
| **Parameters** | `data` — A reference to a `rvector` containing the input signal. `n` — The number of data points to transform. `Forward` — A boolean flag indicating whether to perform a forward (true) or inverse (false) DFT. |
| **Returns**  | `void`                                                                                                                  |
| **Behaviour** | 1. Calls the `FFT` function to perform the DFT on the input data. 2. Returns the transformed data after the DFT has been completed. |
| **Side effects** | Modifies the `data` vector to store the DFT output. |
| **Errors / exceptions** | Throws an exception if `n` is not a power of 2. |

#### `void TwoFFT(rvector data1, rvector data2, rvector fft1, rvector fft2, int n)`

| Field        | Detail                                                                                                                |
|--------------|-----------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Performs a 2D DFT (FFT) on two 1D real signals, `data1` and `data2`, with length `n`.                                |
| **Parameters** | `data1` — A reference to a `rvector` containing the first input signal. `data2` — A reference to a `rvector` containing the second input signal. `fft1` — A reference to a `rvector` to store the FFT of `data1`. `fft2` — A reference to a `rvector` to store the FFT of `data2`. `n` — The number of data points in each signal (must be a power of 2). |
| **Returns**  | `void`                                                                                                                 |
| **Behaviour** | 1. Calls the `FFT` function twice, once on `data1` and once on `data2`. 2. Assigns the first `n/2` elements of `fft1` to the first `n/2` elements of `fft2` and vice-versa. 3. Performs the DFT on `fft1` and `fft2`. |
| **Side effects** | Modifies `fft1` and `fft2` to store the resulting DFTs. |
| **Errors / exceptions** | Throws an exception if `n` is not a power of 2. |

#### `void Convolve(rvector data, int n, rvector respns, int m, rvector ans, bool Conv)`

| Field        | Detail                                                                                                                                                                |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Computes the convolution (or deconvolution) of a 1D real signal `data` with a response function `respns`.                                                               |
| **Parameters** | `data` — A reference to a `rvector` containing the input signal. `n` — The length of the input signal. `respns` — A reference to a `rvector` containing the response function. `m` — The length of the response function (must be odd). `ans` — A reference to a `rvector` to store the convolution result. `Conv` — A boolean flag indicating whether to perform a convolution (true) or deconvolution (false). |
| **Returns**  | `void`                                                                                                                                                                 |
| **Behaviour** | 1. Performs the convolution of `data` with `respns` using the `m`th convolution definition, with `Conv=true`. If `Conv=false`, it performs the deconvolution. This is equivalent to performing recursive convolutions. 2. Outputs the result of the convolution into the `ans` vector. |
| **Side effects** | Modifies the `ans` vector. |
| **Errors / exceptions** | Throws an exception if `m` is not an odd number. |

#### `void mConvolve(rvector data, int n, int m)`

| Field        | Detail                                                                                                                                                                                          |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Computes the mth recursive convolution of a real-valued data array.                                                                                                                                 |
| **Parameters** | `data` — A reference to a `rvector` containing the input signal. `n` — The length of the input signal. `m` — The order of the convolution (must be an integer greater than 1).                                  |
| **Returns**  | `void`                                                                                                                                                                                          |
| **Behaviour** | 1. Performs the mth convolution of `data` with itself using the recursive definition. This is equivalent to performing the convolution using `Convolve` with the same parameters. 2. Outputs the result of the convolution into the `data` vector. |
| **Side effects** | Modifies the `data` vector to store the convolution result. |
| **Errors / exceptions** | Throws an exception if `m` is not greater than 1, or if `n` is not a power of 2. |

---
## `mmdb::math::FFT`

> Performs a discrete Fourier transform (DFT) on a 1D real signal.

### Methods

#### `void FFT(rvector data, int nn, bool Forward)`

| Field        | Detail                     |
|--------------|----------------------------|
| **Summary**  | Computes the discrete Fourier transform of a 1D real signal. |
| **Parameters** | `data` — A reference to a `rvector` containing the input signal. `nn` — The number of data points to transform (must be a power of 2). `Forward` — A boolean flag indicating whether to perform a forward (true) or inverse (false) DFT. |
| **Returns**  | `void`                     |
| **Behaviour** | 1. Allocates a temporary `rvector` called `fft` with size `2*nn`. 2. Iterates through the input data points from index 1 to 2*nn, swapping adjacent elements to perform the DFT. Uses the bit-reversal permutation to do this. 3. Applies the DFT algorithm (Radix-2 Cooley-Tukey) recursively, reducing the size of the signal by half in each step. 4. Returns the transformed data in the `fft` vector. |
| **Side effects** | Modifies the `fft` vector. The original `data` vector is not altered. |
| **Errors / exceptions** | Throws an exception if `nn` is not a power of 2. |

#### `void RealFFT(rvector data, int n, bool Forward)`

| Field        | Detail                                                                                                                   |
|--------------|--------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Performs a real-valued discrete Fourier transform (DFT) on a 1D real signal.                                       |
| **Parameters** | `data` — A reference to a `rvector` containing the input signal. `n` — The number of data points to transform. `Forward` — A boolean flag indicating whether to perform a forward (true) or inverse (false) DFT. |
| **Returns**  | `void`                                                                                                                  |
| **Behaviour** | 1. Calls the `FFT` function to perform the DFT on the input data. 2. Returns the transformed data after the DFT has been completed. |
| **Side effects** | Modifies the `data` vector to store the DFT output. |
| **Errors / exceptions** | Throws an exception if `n` is not a power of 2. |

#### `void TwoFFT(rvector data1, rvector data2, rvector fft1, rvector fft2, int n)`

| Field        | Detail                                                                                                                |
|--------------|-----------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Performs a 2D DFT (FFT) on two 1D real signals, `data1` and `data2`, with length `n`.                                |
| **Parameters** | `data1` — A reference to a `rvector` containing the first input signal. `data2` — A reference to a `rvector` containing the second input signal. `fft1` — A reference to a `rvector` to store the FFT of `data1`. `fft2` — A reference to a `rvector` to store the FFT of `data2`. `n` — The number of data points in each signal (must be a power of 2). |
| **Returns**  | `void`                                                                                                                 |
| **Behaviour** | 1. Calls the `FFT` function twice, once on `data1` and once on `data2`. 2. Assigns the first `n/2` elements of `fft1` to the first `n/2` elements of `fft2` and vice-versa. 3. Performs the DFT on `fft1` and `fft2`. |
| **Side effects** | Modifies `fft1` and `fft2` to store the resulting DFTs. |
| **Errors / exceptions** | Throws an exception if `n` is not a power of 2. |

#### `void Convolve(rvector data, int n, rvector respns, int m, rvector ans, bool Conv)`

| Field        | Detail                                                                                                                                                                                          |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Computes the convolution (or deconvolution) of a 1D real signal `data` with a response function `respns`.                                                                                              |
| **Parameters** | `data` — A reference to a `rvector` containing the input signal. `n` — The length of the input signal. `respns` — A reference to a `rvector` containing the response function. `m` — The length of the response function (must be odd). `ans` — A reference to a `rvector` to store the convolution result. `Conv` — A boolean flag indicating whether to perform a convolution (true) or deconvolution (false). |
| **Returns**  | `void`                                                                                                                                                                                          |
| **Behaviour** | 1. Performs the convolution of `data` with `respns` using the `m`th convolution definition, with `Conv=true`. If `Conv=false`, it performs the deconvolution. This is equivalent to performing recursive convolutions. 2. Outputs the result of the convolution into the `ans` vector. |
| **Side effects** | Modifies the `ans` vector. |
| **Errors / exceptions** | Throws an exception if `m` is not an odd number. |

#### `void mConvolve(rvector data, int n, int m)`

| Field        | Detail                                                                                                                                                                                          |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Computes the mth recursive convolution of a real-valued data array.                                                                                                                                 |
| **Parameters** | `data` — A reference to a `rvector` containing the input signal. `n` — The length of the input signal. `m` — The order of the convolution (must be an integer greater than 1).                                  |
| **Returns**  | `void`                                                                                                                                                                                          |
| **Behaviour** | 1. Performs the mth convolution of `data` with itself using the recursive definition. This is equivalent to performing the convolution using `Convolve` with the same parameters. 2. Outputs the result of the convolution into the `data` vector. |
| **Side effects** | Modifies the `data` vector to store the convolution result. |
| **Errors / exceptions** | Throws an exception if `m` is not greater than 1, or if `n` is not a power of 2. |

---
## Free Functions


---

### `mmdb_math_graph`



---

### `mmdb_math_linalg`

Okay, here's the comprehensive documentation generated from the C++ source, adhering to the requested format.  This includes detailed descriptions of the methods, parameters, return values, behaviour, side effects, and error handling.

---
## `mmdb::math::Jacobi`

> This method performs the Jacobi method for diagonalizing symmetric matrices.

### Methods

#### `void Jacobi(int N, rmatrix A, rmatrix T, rvector Eigen, rvector Aik, int & Signal)`

| Field        | Detail                                                                                                                                                                                                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | This method diagonalizes a symmetric matrix `A` (of size N x N) using the Jacobi method, storing the eigenvalues in the `Eigen` vector and the corresponding eigenvectors in the `T` matrix. |
| **Parameters** | `N` — The dimension of the input matrix `A` (number of rows and columns).  `rmatrix A` — The input symmetric matrix.  `rmatrix T` — The output matrix that will contain the eigenvectors (columns).  `rvector Eigen` — The output vector that will contain the eigenvalues (ordered by increasing value). `rvector Aik` — A temporary working array. `int & Signal` — An integer variable used to indicate the success or failure of the algorithm. |
| **Returns**  | `void`                                                                                                                                                                                                                                                                              |
| **Behaviour** | 1.  Initialization: Sets the diagonal elements of `T` to 1.0 and initializes the `Eigen` vector with the diagonal elements of `A`. 2.  Iteration: Performs the Jacobi iteration process. This involves: a.  Calculating the standard deviation of the diagonal elements. b.  Finding the row with the largest off-diagonal element. c.  Performing a series of column operations to eliminate the off-diagonal element below the current diagonal element. d.  Updating the `Eigen` vector based on the calculated changes.  e.  Checking the convergence criteria. 3.  Output: Stores the computed eigenvalues in the `Eigen` vector and the corresponding eigenvectors in the `T` matrix. |
| **Side effects** | The input matrix `A` is modified in place (specifically, the lower triangle, excluding the diagonal, is modified). The `T` matrix is overwritten with the eigenvectors. |
| **Errors / exceptions** | If the iteration limit (`ItMax`) is reached without convergence, the `Signal` variable is set to `ItMax`, indicating that the algorithm failed to converge. |

---
## `mmdb::math::PbCholDecomp`

> This method performs a perturbed Cholesky decomposition of a symmetric matrix.

### Methods

#### `void PbCholDecomp(int N, rvector HDiag, realtype MaxOff, realtype MachEps, rmatrix L, realtype & MaxAdd)`

| Field        | Detail                                                                                                                                                                                                     |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | This method performs a perturbed Cholesky decomposition of a symmetric matrix.  It is used to solve linear systems of the form `L*Y = B` where `L` is a lower triangular matrix and `B` is a vector. |
| **Parameters** | `N` — The dimension of the input matrix.  `rvector HDiag` — An array holding the diagonal elements of the lower triangular matrix. `realtype MaxOff` — A maximum off-diagonal element. `realtype MachEps` — The machine epsilon. `rmatrix L` — The output lower triangular matrix. `realtype & MaxAdd` —  A reference to an output real type variable which will store the maximum add performed during the decomposition. |
| **Returns**  | `void`                                                                                                                                                                                                      |
| **Behaviour** | 1. Calculates the maximum off-diagonal element. 2. If the maximum off-diagonal element is zero, computes the diagonal elements of L. 3. Iterates through the rows of the matrix to calculate each diagonal element. 4.  Updates the matrix `L` at each iteration. |
| **Side effects** | The matrix `L` is overwritten with the lower triangular matrix.  The `MaxAdd` variable is updated with the maximum add performed during the decomposition. |
| **Errors / exceptions** | None                                                                                                                                                                                               |

---
## `mmdb::math::LSolve`

> This method solves a linear system of the form `L*Y = B` where `L` is a lower triangular matrix and `B` is a vector.

### Methods

#### `void LSolve(int N, rmatrix L, rvector B, rvector Y)`

| Field        | Detail                                                                                                                                                                                                                                 |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | This method solves a linear system of the form `L*Y = B` where `L` is a lower triangular matrix and `B` is a vector, storing the solution `Y` in the output vector. |
| **Parameters** | `N` — The dimension of the input matrix.  `rmatrix L` — The lower triangular matrix.  `rvector B` — The right-hand side vector. `rvector Y` — The output solution vector. |
| **Returns**  | `void`                                                                                                                                                                                                                                 |
| **Behaviour** | 1. Assigns the first element of the solution vector `Y` to be the first element of the right-hand side vector `B`. 2. Iterates through the remaining rows of the matrix to solve for the remaining elements of the solution vector. |
| **Side effects** | The solution vector `Y` is overwritten with the solution. |
| **Errors / exceptions** | None                                                                                                                                                                                                                              |

---
## `mmdb::math::LTSolve`

> This method solves a linear system of the form `LT*X = Y` where `LT` is the transpose of a lower triangular matrix and `X` and `Y` are vectors.

### Methods

#### `void LTSolve(int N, rmatrix L, rvector Y, rvector X)`

| Field        | Detail                                                                                                                                                                                                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | This method solves a linear system of the form `LT*X = Y` where `LT` is the transpose of a lower triangular matrix and `X` and `Y` are vectors, storing the solution `X` in the output vector. |
| **Parameters** | `N` — The dimension of the input matrix. `rmatrix L` — The lower triangular matrix. `rvector Y` — The right-hand side vector. `rvector X` — The output solution vector. |
| **Returns**  | `void`                                                                                                                                                                                                                                                                              |
| **Behaviour** | 1. Assigns the last element of the solution vector `X` to be the last element of the right-hand side vector `Y`. 2. Iterates through the remaining rows of the matrix to solve for the remaining elements of the solution vector. |
| **Side effects** | The solution vector `X` is overwritten with the solution. |
| **Errors / exceptions** | None                                                                                                                                                                                                                                                             |

---
## `mmdb::math::ChSolve`

> This method solves the equation `L*LT*S = G` by the Cholesky's method.

### Methods

#### `void ChSolve(int N, rmatrix L, rvector G, rvector S)`

| Field        | Detail                                                                                                                                                                                                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | This method solves the equation `L*LT*S = G` using the Cholesky's method, storing the solution `S` in the output vector. |
| **Parameters** | `N` — The dimension of the input matrix.  `rmatrix L` — The lower triangular matrix. `rvector G` — The right-hand side vector. `rvector S` — The output solution vector. |
| **Returns**  | `void`                                                                                                                                                                                                                                                                              |
| **Behaviour** | 1. Solves `L*Y = B` using `LSolve`. 2. Solves `LT*X = Y` using `LTSolve`. |
| **Side effects** | The solution vector `S` is overwritten with the solution. |
| **Errors / exceptions** | None                                                                                                                                                                                                                                                             |

---
## `mmdb::math::FastInverse`

> This method performs a fast inverse of a matrix `A` using the Gauss-Jordan method.

### Methods

#### `void FastInverse(int N, rmatrix A, ivector J0, int & Signal)`

| Field        | Detail                                                                                                                                                                                                                                              |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | This method performs a fast inverse of a symmetric matrix `A` using the Gauss-Jordan method. |
| **Parameters** | `N` — The dimension of the input matrix. `rmatrix A` — The input matrix. `ivector J0` — An integer vector used to store the indices of the pivots. `int & Signal` — An integer variable used to indicate the success or failure of the algorithm. |
| **Returns**  | `void`                                                                                                                                                                                                                                               |
| **Behaviour** | 1. Performs the Gauss-Jordan elimination to transform the input matrix `A` into the identity matrix. 2. The integer vector `J0` stores the indices of the pivots used during the elimination process. 3. Checks for degeneracy and sets the `Signal` variable accordingly. |
| **Side effects** | The input matrix `A` is modified in place. |
| **Errors / exceptions** | If the matrix becomes degenerate during the elimination process, the `Signal` variable is set to the rank of the matrix. |

---
## `mmdb::math::SVD`

> This method performs Singular Value Decomposition (SVD) of a matrix.

### Methods

#### `void SVD(int NA, int M, int N, rmatrix A, rmatrix U, rmatrix V, rvector W, rvector RV1, bool MatU, bool MatV, int & RetCode)`

| Field        | Detail                                                                                                                                                                                                                                                                            |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | This method performs Singular Value Decomposition (SVD) of a matrix `A`. |
| **Parameters** | `NA` — Number of lines in `A`. `M` — Number of lines in `U`. `N` — Number of columns in `U`, `V` and length of `W`, `RV1`. `rmatrix A` — The input matrix. `rmatrix U` — The output matrix of right singular vectors. `rmatrix V` — The output matrix of left singular vectors. `rvector W` — The singular values. `rvector RV1` — A temporary array. `bool MatU` — A flag indicating whether to compute `U`. `bool MatV` — A flag indicating whether to compute `V`. `int & RetCode` — An integer variable which stores the error code. |
| **Returns**  | `void`                                                                                                                                                                                                                                                                            |
| **Behaviour** | 1.  Performs SVD using the algorithm described in G. Forsait, M. Malkolm, K. Mouler. `NA` needs to be provided and is checked.  2. Computes `U`, `V` and `W`.  3.  The error code is set to 0 if the decomposition is successful. Otherwise, it's set to the number of iterations required to obtain the current singular values. |
| **Side effects** | The matrices `U` and `V` are overwritten with the singular vectors. The vector `W` is overwritten with the singular values. |
| **Errors / exceptions** | If the algorithm fails to converge after a certain number of iterations (`ItnLimit`), the `RetCode` variable is set to the number of iterations that were performed. |

---
## `mmdb::math::OrderSVD`

> This method reorders the singular values and vectors obtained from SVD to have the largest singular values at the beginning.

### Methods

#### `void OrderSVD(int M, int N, rmatrix U, rmatrix V, rvector W, bool MatU, bool MatV)`

| Field        | Detail                                                                                                                                                                                                                                                                      |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | This method reorders the singular values and vectors obtained from SVD so that the largest singular values are placed at the beginning. |
| **Parameters** | `M` — The number of lines in `U`. `N` — The number of columns in `U`, `V` and the length of `W` and `RV1`. `rmatrix U` — The input matrix of right singular vectors. `rmatrix V` — The input matrix of left singular vectors. `rvector W` — The singular values. `bool MatU` — Flag indicating whether to compute `U`. `bool MatV` — Flag indicating whether to compute `V`. |
| **Returns**  | `void`                                                                                                                                                                                                                                                                      |
| **Behaviour** | This method sorts the singular values in descending order and rearranges the corresponding singular vectors in `U` and `V`. |
| **Side effects** | The matrices `U` and `V` are modified in place. |
| **Errors / exceptions** | None                                                                                                                                                                                                                                                             |

---
## Free Functions

No free functions are present in the provided source files.
---

This documentation provides a detailed breakdown of the functionality of each method, including the parameters, return values, and behaviour.  It should be sufficient for a developer to understand and utilize the `mmdb::math::LinAlg` library effectively.

---

### `mmdb_math_rand`

```cpp
---
## `mmdb::math::RandomNumber`

> This class implements a pseudo-random number generator, primarily intended for generating Gaussian and uniform random numbers, as used in molecular dynamics simulations.

### Methods

#### `mmdb::math::RandomNumber(long IJ, long KL)`

| Field        | Detail                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------|
| **Summary**  | Constructs a `RandomNumber` object with the specified IJ and KL seeds.  Initializes the internal state of the random number generator. |
| **Parameters** | `IJ` — The first seed value (long).  Controls the initial state of the generator.  `KL` — The second seed value (long). Controls the initial state of the generator. |
| **Returns**  | `void`                                                                                             |
| **Behaviour** | 1.  Initializes the internal state using the provided `IJ` and `KL` seed values. This involves calculating various constants and setting up an array (`U`) to store the initial values for the generator. 2.  The constants `C`, `CD`, and `CM` are calculated to ensure a uniform distribution of random numbers. 3.  The `I97` and `J97` variables are initialized to specify the indices within the `U` array to be used for generating numbers. 4.  The `iset` and `gset` are initialized to control the generation of Gaussian numbers.|
| **Side effects** | Mutates the internal state of the `RandomNumber` object (`U`, `C`, `CD`, `CM`, `I97`, `J97`, `iset`, `gset`). |
| **Errors / exceptions** | `none`                                                                                         |

#### `void mmdb::math::RandomNumber::Init(long IJ, long KL)`

| Field        | Detail                                                                                                                   |
|--------------|--------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Initializes the internal state of the `RandomNumber` object using the specified `IJ` and `KL` seed values.                |
| **Parameters** | `IJ` — The first seed value (long).  `KL` — The second seed value (long).                                            |
| **Returns**  | `void`                                                                                                                 |
| **Behaviour** | 1.  Calculates initial values for constants `C`, `CD`, and `CM` which are used to generate the random numbers.   2.  Populates the `U` array with initial values based on the provided seed values and a specific algorithm.  This algorithm is based on a series of linear congruences and multiplications. 3. Sets the state variables: `I97` and `J97` to the initial indices of the `U` array.  4. Initializes the `iset` and `gset` variables.                                      |
| **Side effects** | Mutates the internal state of the `RandomNumber` object (`U`, `C`, `CD`, `CM`, `I97`, `J97`, `iset`, `gset`). |
| **Errors / exceptions** | `none`                                                                                                                 |

#### `realtype mmdb::math::RandomNumber::gauss_rnd()`

| Field        | Detail                                                                                                                                                                                        |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Generates a single Gaussian (normal) random number.                                                                                                                                               |
| **Parameters** | `none`                                                                                                                                                                                          |
| **Returns**  | A Gaussian random number (realtype).                                                                                                                                                              |
| **Behaviour** | 1. If the `iset` flag is false, it calls `srandom()` to generate a random number and then generates the Gaussian number using the standard normal distribution formula. 2.  If the `iset` flag is true, it returns the current value of `gset`.  This is a way to avoid re-initializing the Gaussian generator every time. |
| **Side effects** | Mutates the internal state of the `RandomNumber` object (`gset`, `iset`). |
| **Errors / exceptions** | `none`                                                                                                                                                                                   |

#### `realtype mmdb::math::RandomNumber::random()`

| Field        | Detail                                                                                                                                                                                                    |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Generates a single uniform random number between 0.0 and 1.0.                                                                                                                                                  |
| **Parameters** | `none`                                                                                                                                                                                                  |
| **Returns**  | A uniform random number (realtype) in the range [0.0, 1.0].                                                                                                                                                |
| **Behaviour** | 1. Accesses the `U` array at indices `I97` and `J97` to get the initial values.  2.  Subtracts the value at `J97` from the value at `I97` to get a uniform value between 0.0 and 1.0.  3.  Updates the values at `I97` and `J97` to prepare for the next call.  4. Updates the constants `C`, `CD`, and `CM`. |
| **Side effects** | Mutates the internal state of the `RandomNumber` object (`I97`, `J97`, `C`, `CD`, `CM`). |
| **Errors / exceptions** | `none`                                                                                                                                                                                           |

#### `realtype mmdb::math::RandomNumber::srandom()`

| Field        | Detail                                                                                                                                                                                                    |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Generates a single uniform random number between -1.0 and 1.0.                                                                                                                                               |
| **Parameters** | `none`                                                                                                                                                                                                  |
| **Returns**  | A uniform random number (realtype) in the range [-1.0, 1.0].                                                                                                                                             |
| **Behaviour** | 1. Similar to `random()`, it accesses the `U` array at indices `I97` and `J97`.  2.  Subtracts the value at `J97` from the value at `I97` to get a uniform value between -1.0 and 1.0.  3.  Updates the values at `I97` and `J97` to prepare for the next call. 4. Updates the constants `C`, `CD`, and `CM`. |
| **Side effects** | Mutates the internal state of the `RandomNumber` object (`I97`, `J97`, `C`, `CD`, `CM`). |
| **Errors / exceptions** | `none`                                                                                                                                                                                           |

#### `void mmdb::math::RandomNumber::read(io::RFile f)`

| Field        | Detail                                                                                                                                   |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads the state of the `RandomNumber` object from the given input stream (`io::RFile`).                                               |
| **Parameters** | `f` — An input stream (`io::RFile`) from which to read the data.                                                                        |
| **Returns**  | `void`                                                                                                                                   |
| **Behaviour** | 1. Reads the version number from the input stream.  2. Reads the values for `I97`, `J97`, `U`, `C`, `CD`, `CM`, `gset`, and `iset` from the input stream. |
| **Side effects** | None (only reads data from the input stream).                                                                                          |
| **Errors / exceptions** | `none`                                                                                                                            |

#### `void mmdb::math::RandomNumber::write(io::RFile f)`

| Field        | Detail                                                                                                                                                                |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Writes the state of the `RandomNumber` object to the given output stream (`io::RFile`).                                                                                 |
| **Parameters** | `f` — An output stream (`io::RFile`) to which to write the data.                                                                                                        |
| **Returns**  | `void`                                                                                                                                                                 |
| **Behaviour** | 1. Writes the version number to the output stream. 2. Writes the values for `I97`, `J97`, `U`, `C`, `CD`, `CM`, `gset`, and `iset` to the output stream.                            |
| **Side effects** | Writes data to the output stream (`io::RFile`).                                                                                                                      |
| **Errors / exceptions** | `none`                                                                                                                                                          |

#### `RandomNumber::RandomNumber(long IJ, long KL)`

| Field        | Detail                                                                                               |
|--------------|------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructor for the `RandomNumber` class.  Initializes the generator with the given seed values.     |
| **Parameters** | `IJ` — The first seed value (long). `KL` — The second seed value (long).                          |
| **Returns**  | `RandomNumber*`  (a pointer to a newly created `RandomNumber` object).                         |
| **Behaviour** | 1. Calls the `Init()` method to initialize the internal state of the `RandomNumber` object using the provided `IJ` and `KL` seed values. |
| **Side effects** | Creates and initializes a `RandomNumber` object, modifying its internal state.                   |
| **Errors / exceptions** | `none`                                                                                         |

#### `void RandomNumber::destroy()`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Destructor for the `RandomNumber` class.                                                              |
| **Parameters** | `none`                                                                                               |
| **Returns**  | `void`                                                                                              |
| **Behaviour** | 1. Performs any necessary cleanup before the object is destroyed.  This may include releasing any resources allocated by the object. |
| **Side effects** | None.                                                                                              |
| **Errors / exceptions** | `none`                                                                                         |

---
## `mmdb::math::RandomNumber`

> This class implements a pseudo-random number generator, primarily intended for generating Gaussian and uniform random numbers, as used in molecular dynamics simulations.

### Methods

#### `mmdb::math::RandomNumber::RandomNumber(long IJ, long KL)`

| Field        | Detail                                                                                                                                                                                                    |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructs a `RandomNumber` object with the specified IJ and KL seeds.  Initializes the internal state of the random number generator.                                                                          |
| **Parameters** | `IJ` — The first seed value (long).  `KL` — The second seed value (long).                                                                                                                                    |
| **Returns**  | `void`                                                                                                                                                                                                  |
| **Behaviour** | 1.  Initializes the internal state using the provided `IJ` and `KL` seed values.  This involves calculating various constants and setting up an array (`U`) to store the initial values for the generator.  2.  The constants `C`, `CD`, and `CM` are calculated to ensure a uniform distribution of random numbers.  3.  The `I97` and `J97` variables are initialized to specify the indices within the `U` array to be used for generating numbers.  4.  The `iset` and `gset` are initialized to control the generation of Gaussian numbers. |
| **Side effects** | Mutates the internal state of the `RandomNumber` object (`U`, `C`, `CD`, `CM`, `I97`, `J97`, `iset`, `gset`).                                                                                           |
| **Errors / exceptions** | `none`                                                                                                                                                                                           |

#### `void mmdb::math::RandomNumber::Init(long IJ, long KL)`

| Field        | Detail                                                                                                                                                                                                    |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Initializes the internal state of the `RandomNumber` object using the specified `IJ` and `KL` seed values.                                                                                                 |
| **Parameters** | `IJ` — The first seed value (long).  `KL` — The second seed value (long).                                                                                                                                    |
| **Returns**  | `void`                                                                                                                                                                                                  |
| **Behaviour** | 1.  Calculates initial values for constants `C`, `CD`, and `CM` which are used to generate the random numbers.  2.  Populates the `U` array with initial values based on the provided seed values and a specific algorithm.  This algorithm is based on a series of linear congruences and multiplications.  3.  Sets the state variables: `I97` and `J97` to the initial indices of the `U` array.  4.  Initializes the `iset` and `gset` variables.                                      |
| **Side effects** | Mutates the internal state of the `RandomNumber` object (`U`, `C`, `CD`, `CM`, `I97`, `J97`, `iset`, `gset`).                                                                                           |
| **Errors / exceptions** | `none`                                                                                                                                                                                           |

#### `realtype mmdb::math::RandomNumber::gauss_rnd()`

| Field        | Detail                                                                                                                                                                                                    |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Generates a single Gaussian (normal) random number.                                                                                                                                               |
| **Parameters** | `none`                                                                                                                                                                                                  |
| **Returns**  | A Gaussian random number (realtype).                                                                                                                                                              |
| **Behaviour** | 1. If the `iset` flag is false, it calls `srandom()` to generate a random number and then generates the Gaussian number using the standard normal distribution formula.  2.  If the `iset` flag is true, it returns the current value of `gset`.  This is a way to avoid re-initializing the Gaussian generator every time. |
| **Side effects** | Mutates the internal state of the `RandomNumber` object (`gset`, `iset`).                                                                                                                          |
| **Errors / exceptions** | `none`                                                                                                                                                                                           |

#### `realtype mmdb::math::RandomNumber::random()`

| Field        | Detail                                                                                                                                                                                                    |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Generates a single uniform random number between 0.0 and 1.0.                                                                                                                                               |
| **Parameters** | `none`                                                                                                                                                                                                  |
| **Returns**  | A uniform random number (realtype) in the range [0.0, 1.0].                                                                                                                                                |
| **Behaviour** | 1. Accesses the `U` array at indices `I97` and `J97` to get the initial values.  2.  Subtracts the value at `J97` from the value at `I97` to get a uniform value between 0.0 and 1.0.  3.  Updates the values at `I97` and `J97` to prepare for the next call.  4. Updates the constants `C`, `CD`, and `CM`. |
| **Side effects** | Mutates the internal state of the `RandomNumber` object (`I97`, `J97`, `C`, `CD`, `CM`).                                                                                                                          |
| **Errors / exceptions** | `none`                                                                                                                                                                                           |

#### `realtype mmdb::math::RandomNumber::srandom()`

| Field        | Detail                                                                                                                                                                                                    |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Generates a single uniform random number between -1.0 and 1.0.                                                                                                                                               |
| **Parameters** | `none`                                                                                                                                                                                                  |
| **Returns**  | A uniform random number (realtype) in the range [-1.0, 1.0].                                                                                                                                             |
| **Behaviour** | 1. Similar to `random()`, it accesses the `U` array at indices `I97` and `J97`.  2.  Subtracts the value at `J97` from the value at `I97` to get a uniform value between -1.0 and 1.0.  3.  Updates the values at `I97` and `J97` to prepare for the next call. 4. Updates the constants `C`, `CD`, and `CM`. |
| **Side effects** | Mutates the internal state of the `RandomNumber` object (`I97`, `J97`, `C`, `CD`, `CM`).                                                                                                                          |
| **Errors / exceptions** | `none`                                                                                                                                                                                           |

#### `void mmdb::math::RandomNumber::write(io::RFile f)`

| Field        | Detail                                                                                                                                                                  |
|--------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Writes the state of the `RandomNumber` object to the given output stream (`io::RFile`).                                                                                 |
| **Parameters** | `f` — An output stream (`io::RFile`) to which to write the data.                                                                                                        |
| **Returns**  | `void`                                                                                                                                                                   |
| **Behaviour** | 1. Writes the version number to the output stream.  2. Writes the values for `I97`, `J97`, `U`, `C`, `CD`, `CM`, `gset`, and `iset` to the output stream.                            |
| **Side effects** | Writes data to the output stream (`io::RFile`).                                                                                                                      |
| **Errors / exceptions** | `none`                                                                                                                                                          |

#### `void mmdb::math::RandomNumber::read(io::RFile f)`

| Field        | Detail                                                                                                                                   |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads the state of the `RandomNumber` object from the given input stream (`io::RFile`).                                               |
| **Parameters** | `f` — An input stream (`io::RFile`) from which to read the data.                                                                        |
| **Returns**  | `void`                                                                                                                                   |
| **Behaviour** | 1. Reads the version number from the input stream.  2. Reads the values for `I97`, `J97`, `U`, `C`, `CD`, `CM`, `gset`, and `iset` from the input stream. |
| **Side effects** | None (only reads data from the input stream).                                                                                          |
| **Errors / exceptions** | `none`                                                                                                                            |

#### `mmdb::math::RandomNumber::RandomNumber(long IJ, long KL)`

| Field        | Detail                                                                                               |
|--------------|------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructor for the `RandomNumber` class.  Initializes the generator with the given seed values.     |
| **Parameters** | `IJ` — The first seed value (long).  `KL` — The second seed value (long).                                                                                              |
| **Returns**  | `RandomNumber*`  (a pointer to a newly created `RandomNumber` object).                         |
| **Behaviour** | 1. Calls the `Init()` method to initialize the internal state of the `RandomNumber` object using the provided `IJ` and `KL` seed values. |
| **Side effects** | Creates and initializes a `RandomNumber` object, modifying its internal state.                   |
| **Errors / exceptions** | `none`                                                                                             |

#### `void mmdb::math::RandomNumber::destroy()`

| Field        | Detail                                                                                                |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Destructor for the `RandomNumber` class.                                                              |
| **Parameters** | `none`                                                                                               |
| **Returns**  | `void`                                                                                              |
| **Behaviour** | 1. Performs any necessary cleanup before the object is destroyed.  This may include releasing any resources allocated by the object. |
| **Side effects** | None.                                                                                              |
| **Errors / exceptions** | `none`                                                                                             |

---
## Free Functions

There are no free functions in this module.
```

---

### `mmdb_mattype`

## `mmdb_mattype`

Provides a set of utility functions for working with floating-point numbers, strings, and matrices in a crystallography context.

### Methods

#### `bool InitMatType()`

| Field        | Detail               |
|--------------|-----------------------|
| **Summary**  | Initializes the internal data structures and constants used by the `mmdb_mattype` module. |
| **Parameters** | none                  |
| **Returns**  | `true` if initialization was successful, `false` otherwise. |
| **Behaviour** | Calls `MachinEps()`, `InitFPowers()`, and sets internal constant values for precise calculations with floating point numbers. It makes sure the double float type is selected appropriately. |
| **Side effects** | Sets internal constants; updates MachEps, floatMachEps, LnMaxReal, LnMinReal, etc. |
| **Errors / exceptions** | None.  Returns `false` if there are internal errors during initialization. |

#### `realtype MachEps()`

| Field        | Detail                     |
|--------------|----------------------------|
| **Summary**  | Returns the machine epsilon for the `realtype` floating-point type. |
| **Parameters** | none                      |
| **Returns**  | The machine epsilon value. |
| **Behaviour** | Uses the standard `std::numeric_limits<realtype>::epsilon()` function. |
| **Side effects** | None.                       |
| **Errors / exceptions** | None.                       |

#### `realtype floatMachEps()`

| Field        | Detail                     |
|--------------|----------------------------|
| **Summary**  | Returns the machine epsilon for the `float` floating-point type. |
| **Parameters** | none                      |
| **Returns**  | The machine epsilon value. |
| **Behaviour** | Uses the standard `std::numeric_limits<float>::epsilon()` function. |
| **Side effects** | None.                       |
| **Errors / exceptions** | None.                       |

#### `realtype LnMaxReal()`

| Field        | Detail                      |
|--------------|-----------------------------|
| **Summary**  | Returns the maximum value of the natural logarithm of `realtype`. |
| **Parameters** | none                       |
| **Returns**  | The maximum value of the natural logarithm. |
| **Behaviour** | Calculates `log(fMaxReal)` using the `log` function. |
| **Side effects** | None.                       |
| **Errors / exceptions** | None.                       |

#### `realtype LnMinReal()`

| Field        | Detail                      |
|--------------|-----------------------------|
| **Summary**  | Returns the minimum value of the natural logarithm of `realtype`. |
| **Parameters** | none                       |
| **Returns**  | The minimum value of the natural logarithm. |
| **Behaviour** | Calculates `log(fMinReal)` using the `log` function. |
| **Side effects** | None.                       |
| **Errors / exceptions** | None.                       |

#### `int mround(realtype X)`

| Field        | Detail                   |
|--------------|---------------------------|
| **Summary**  | Rounds a `realtype` value to the nearest integer. |
| **Parameters** | `X` — A `realtype` value to round. |
| **Returns**  | The nearest integer value. |
| **Behaviour** | Adds 0.5 to the input `X` and converts to an integer. |
| **Side effects** | None.                       |
| **Errors / exceptions** | None.                       |

#### `int ifloor(realtype X)`

| Field        | Detail                   |
|--------------|---------------------------|
| **Summary**  | Returns the floor of a `realtype` value (the largest integer less than or equal to the input). |
| **Parameters** | `X` — A `realtype` value to floor. |
| **Returns**  | The floor of `X`. |
| **Behaviour** | Uses the `floor` function to truncate the decimal part of the input. |
| **Side effects** | None.                       |
| **Errors / exceptions** | None.                       |

#### `int Abs(int x)`

| Field        | Detail                   |
|--------------|---------------------------|
| **Summary**  | Returns the absolute value of an integer. |
| **Parameters** | `x` — An integer value. |
| **Returns**  | The absolute value of `x`. |
| **Behaviour** |  If `x` is positive, it is returned. If `x` is negative, its negation is returned. |
| **Side effects** | None.                       |
| **Errors / exceptions** | None.                       |

#### `void ISwap(int & x, int & y)`

| Field        | Detail                   |
|--------------|---------------------------|
| **Summary**  | Swaps the values of two integers. |
| **Parameters** | `x` — A reference to the first integer. |
|                | `y` — A reference to the second integer. |
| **Returns**  | None. |
| **Behaviour** |  Assigns the value of `y` to `x` and the value of `x` to `y`. |
| **Side effects** | Modifies the values of `x` and `y`. |
| **Errors / exceptions** | None.                       |

#### `void WSwap(word & x, word & y)`

| Field        | Detail                   |
|--------------|---------------------------|
| **Summary**  | Swaps the values of two words. |
| **Parameters** | `x` — A reference to the first word. |
|                | `y` — A reference to the second word. |
| **Returns**  | None. |
| **Behaviour** |  Assigns the value of `y` to `x` and the value of `x` to `y`. |
| **Side effects** | Modifies the values of `x` and `y`. |
| **Errors / exceptions** | None.                       |

#### `void BSwap(byte & x, byte & y)`

| Field        | Detail                   |
|--------------|---------------------------|
| **Summary**  | Swaps the values of two bytes. |
| **Parameters** | `x` — A reference to the first byte. |
|                | `y` — A reference to the second byte. |
| **Returns**  | None. |
| **Behaviour** |  Assigns the value of `y` to `x` and the value of `x` to `y`. |
| **Side effects** | Modifies the values of `x` and `y`. |
| **Errors / exceptions** | None.                       |

#### `void OSwap(bool & x, bool & y)`

| Field        | Detail                   |
|--------------|---------------------------|
| **Summary**  | Swaps the values of two booleans. |
| **Parameters** | `x` — A reference to the first boolean. |
|                | `y` — A reference to the second boolean. |
| **Returns**  | None. |
| **Behaviour** |  Assigns the value of `y` to `x` and the value of `x` to `y`. |
| **Side effects** | Modifies the values of `x` and `y`. |
| **Errors / exceptions** | None.                       |

#### `void LSwap(long & x, long & y)`

| Field        | Detail                   |
|--------------|---------------------------|
| **Summary**  | Swaps the values of two longs. |
| **Parameters** | `x` — A reference to the first long. |
|                | `y` — A reference to the second long. |
| **Returns**  | None. |
| **Behaviour** |  Assigns the value of `y` to `x` and the value of `x` to `y`. |
| **Side effects** | Modifies the values of `x` and `y`. |
| **Errors / exceptions** | None.                       |

#### `void RSwap(realtype & x, realtype & y)`

| Field        | Detail                   |
|--------------|---------------------------|
| **Summary**  | Swaps the values of two realtype values. |
| **Parameters** | `x` — A reference to the first realtype value. |
|                | `y` — A reference to the second realtype value. |
| **Returns**  | None. |
| **Behaviour** |  Assigns the value of `y` to `x` and the value of `x` to `y`. |
| **Side effects** | Modifies the values of `x` and `y`. |
| **Errors / exceptions** | None.                       |

#### `realtype RMax(const realtype x1, const realtype x2)`

| Field        | Detail                  |
|--------------|-------------------------|
| **Summary**  | Returns the maximum of two `realtype` values. |
| **Parameters** | `x1` — The first `realtype` value. |
|                | `x2` — The second `realtype` value. |
| **Returns**  | The maximum of `x1` and `x2`. |
| **Behaviour** |  Compares the two input values and returns the larger one. |
| **Side effects** | None.                       |
| **Errors / exceptions** | None.                       |

#### `long LMax(const long x1, const long x2)`

| Field        | Detail                  |
|--------------|-------------------------|
| **Summary**  | Returns the maximum of two `long` values. |
| **Parameters** | `x1` — The first `long` value. |
|                | `x2` — The second `long` value. |
| **Returns**  | The maximum of `x1` and `x2`. |
| **Behaviour** |  Compares the two input values and returns the larger one. |
| **Side effects** | None.                       |
| **Errors / exceptions** | None.                       |

#### `word WMax(const word x1, const word x2)`

| Field        | Detail                  |
|--------------|-------------------------|
| **Summary**  | Returns the maximum of two `word` values. |
| **Parameters** | `x1` — The first `word` value. |
|                | `x2` — The second `word` value. |
| **Returns**  | The maximum of `x1` and `x2`. |
| **Behaviour** |  Compares the two input values and returns the larger one. |
| **Side effects** | None.                       |
| **Errors / exceptions** | None.                       |

#### `int IMax(const int x1,  const int x2)`

| Field        | Detail                  |
|--------------|-------------------------|
| **Summary**  | Returns the maximum of two `int` values. |
| **Parameters** | `x1` — The first `int` value. |
|                | `x2` — The second `int` value. |
| **Returns**  | The maximum of `x1` and `x2`. |
| **Behaviour** |  Compares the two input values and returns the larger one. |
| **Side effects** | None.                       |
| **Errors / exceptions** | None.                       |

#### `realtype fsign(const realtype x1,  const realtype x2)`

| Field        | Detail                  |
|--------------|-------------------------|
| **Summary**  | Returns the sign of the difference between two `realtype` values. |
| **Parameters** | `x1` — The first `realtype` value. |
|                | `x2` — The second `realtype` value. |
| **Returns**  | 1.0 if `x1` is greater than `x2`, -1.0 if `x1` is less than `x2`, and 0.0 if they are equal. |
| **Behaviour** | Returns 1.0 if x1 > x2, -1.0 if x1 < x2, and 0.0 if x1 == x2. |
| **Side effects** | None.                       |
| **Errors / exceptions** | None.                       |

#### `void  PadSpaces ( pstr S, int len )`

| Field        | Detail                      |
|--------------|-----------------------------|
| **Summary**  | Pads a string `S` with spaces to length `len`. |
| **Parameters** | `S` — The string to pad.  |
|                | `len` — The desired length of the string. |
| **Returns**  | None. |
| **Behaviour** | Adds space characters to the end of string `S` until its length is `len`. |
| **Side effects** | Modifies the string `S`. |
| **Errors / exceptions** | None.                       |

#### `pstr CutSpaces ( pstr S, int CutKey )`

| Field        | Detail                  |
|--------------|-------------------------|
| **Summary**  | Cuts spaces in the string. |
| **Parameters** | `S` — The string to cut spaces from.  |
|                | `CutKey` — A flag indicating where to cut the spaces.  |
| **Returns**  | The string with spaces cut. |
| **Behaviour** | Cuts spaces at the beginning or end of the string, depending on the value of the `CutKey`. |
| **Side effects** | Modifies the string `S`. |
| **Errors / exceptions** | None.                       |

#### `pstr DelSpaces ( pstr S, char c )`

| Field        | Detail                  |
|--------------|-------------------------|
| **Summary**  | Removes spaces from the string. |
| **Parameters** | `S` — The string to remove spaces from. |
|                | `c` — The character to be replaced by spaces. |
| **Returns**  | The string with spaces removed. |
| **Behaviour** | Replaces all occurrences of character `c` with spaces in string `S`. |
| **Side effects** | Modifies the string `S`. |
| **Errors / exceptions** | None.                       |

#### `pstr EnforceSpaces ( pstr S )`

| Field        | Detail                  |
|--------------|-------------------------|
| **Summary**  | Replaces all unprintable characters, except <CR>, <LF>, <TAB> and some others, with spaces. |
| **Parameters** | `S` — The string to enforce. |
| **Returns**  | The string with unprintable characters replaced with spaces. |
| **Behaviour** |  Replaces all characters that are not printable with spaces. |
| **Side effects** | Modifies the string `S`. |
| **Errors / exceptions** | None.                       |

#### `pstr CreateCopy(pstr & Dest, cpstr Source)`

| Field        | Detail                   |
|--------------|---------------------------|
| **Summary**  | Copies the string from `Source` to `Dest`. |
| **Parameters** | `Dest` — A pointer to the destination string buffer. |
|                | `Source` — The source string. |
| **Returns**  | A pointer to `Dest`. |
| **Behaviour** | Allocates memory for `Dest` if it's not already allocated, copies the contents of `Source` into `Dest`, and null terminates `Dest`. |
| **Side effects** | Modifies `Dest`.     |
| **Errors / exceptions** | None.                       |

#### `pstr CreateCopy_n(pstr & Dest, cpstr Source, int n)`

| Field        | Detail                   |
|--------------|---------------------------|
| **Summary**  | Copies at most `n` characters from `Source` to `Dest`. |
| **Parameters** | `Dest` — A pointer to the destination string buffer. |
|                | `Source` — The source string. |
|                | `n` — The maximum number of characters to copy. |
| **Returns**  | A pointer to `Dest`. |
| **Behaviour** | Allocates memory for `Dest` if it's not already allocated, copies at most `n` characters from `Source` into `Dest`, and null terminates `Dest`. |
| **Side effects** | Modifies `Dest`.     |
| **Errors / exceptions** | None.                       |

#### `pstr CreateCopCat(pstr & Dest, cpstr Source1, cpstr Source2, cpstr Source3, cpstr Source4, cpstr Source5)`

| Field        | Detail                   |
|--------------|---------------------------|
| **Summary**  | Concatenates up to five strings to create a new string. |
| **Parameters** | `Dest` — A pointer to the destination string buffer. |
|                | `Source1`, `Source2`, `Source3`, `Source4`, `Source5` — The source strings to concatenate. |
| **Returns**  | A pointer to `Dest`. |
| **Behaviour** | Allocates memory for `Dest` if it's not already allocated, concatenates the strings into `Dest`, and null terminates `Dest`. |
| **Side effects** | Modifies `Dest`.     |
| **Errors / exceptions** | None.                       |

#### `pstr CreateCopCat(pstr & Dest, cpstr Source1, cpstr Source2, cpstr Source3, cpstr Source4)`

| Field        | Detail                   |
|--------------|---------------------------|
| **Summary**  | Concatenates up to four strings to create a new string. |
| **Parameters** | `Dest` — A pointer to the destination string buffer. |
|                | `Source1`, `Source2`, `Source3`, `Source4` — The source strings to concatenate. |
| **Returns**  | A pointer to `Dest`. |
| **Behaviour** | Allocates memory for `Dest` if it's not already allocated, concatenates the strings into `Dest`, and null terminates `Dest`. |
| **Side effects** | Modifies `Dest`.     |
| **Errors / exceptions** | None.                       |

#### `pstr CreateCopCat(pstr & Dest, cpstr Source1, cpstr Source2, cpstr Source3)`

| Field        | Detail                   |
|--------------|---------------------------|
| **Summary**  | Concatenates up to three strings to create a new string. |
| **Parameters** | `Dest` — A pointer to the destination string buffer. |
|                | `Source1`, `Source2`, `Source3` — The source strings to concatenate. |
| **Returns**  | A pointer to `Dest`. |
| **Behaviour** | Allocates memory for `Dest` if it's not already allocated, concatenates the strings into `Dest`, and null terminates `Dest`. |
| **Side effects** | Modifies `Dest`.     |
| **Errors / exceptions** | None.                       |

#### `pstr CreateCopCat(pstr & Dest, cpstr Source)`

| Field        | Detail                   |
|--------------|---------------------------|
| **Summary**  | Concatenates a single string to create a new string. |
| **Parameters** | `Dest` — A pointer to the destination string buffer. |
|                | `Source` — The source string. |
| **Returns**  | A pointer to `Dest`. |
| **Behaviour** | Allocates memory for `Dest` if it's not already allocated, concatenates the string into `Dest`, and null terminates `Dest`. |
| **Side effects** | Modifies `Dest`.     |
| **Errors / exceptions** | None.                       |

#### `pstr LastOccurence(cpstr S, char c)`

| Field        | Detail                  |
|--------------|-------------------------|
| **Summary**  | Returns a pointer to the last occurrence of a character in a string. |
| **Parameters** | `S` — The string to search. |
|                | `c` — The character to find. |
| **Returns**  | A pointer to the last occurrence of the character, or `NULL` if the character is not found. |
| **Behaviour** | Searches for the character `c` in the string `S` from left to right. Returns a pointer to the last position where `c` was found. |
| **Side effects** | None.                       |
| **Errors / exceptions** | None.                       |

#### `pstr FirstOccurence(cpstr S, int Slen, cpstr Q, int Qlen)`

| Field        | Detail                  |
|--------------|-------------------------|
| **Summary**  | Returns a pointer to the first occurrence of a substring in a string. |
| **Parameters** | `S` — The string to search. |
|                | `Slen` — The length of `S`. |
|                | `Q` — The substring to find. |
|                | `Qlen` — The length of `Q`. |
| **Returns**  | A pointer to the first occurrence of the substring, or `NULL` if the substring is not found. |
| **Behaviour** | Searches for the substring `Q` in `S`.  Returns a pointer to the first character of the substring if found. |
| **Side effects** | None.                       |
| **Errors / exceptions** | None.                       |

#### `int indexOf(cpstr S, char c)`

| Field        | Detail                  |
|--------------|-------------------------|
| **Summary**  | Returns the index of the first occurrence of a character in a string. |
| **Parameters** | `S` — The string to search. |
|                | `c` — The character to find. |
| **Returns**  | The index of the first occurrence of the character, or -1 if the character is not found. |
| **Behaviour** | Searches for the character `c` in the string `S` from left to right. Returns the index of the first occurrence of the character, or -1 if the character is not found. |
| **Side effects** | None.                       |
| **Errors / exceptions** | None.                       |

#### `pstr CreateConcat(pstr & Dest, cpstr Source1, cpstr Source2, cpstr Source3, cpstr Source4, cpstr Source5)`

| Field        | Detail                  |
|--------------|-------------------------|
| **Summary**  | Concatenates up to six strings to create a new string. |
| **Parameters** | `Dest` — A pointer to the destination string buffer. |
|                | `Source1`, `Source2`, `Source3`, `Source4`, `Source5` — The source strings to concatenate. |
| **Returns**  | A pointer to `Dest`. |
| **Behaviour** | Allocates memory for `Dest` if it's not already allocated, concatenates the strings into `Dest`, and null terminates `Dest`. |
| **Side effects** | Modifies `Dest`.     |
| **Errors / exceptions** | None.                       |

#### `pstr CreateConcat(pstr & Dest, cpstr Source1, cpstr Source2, cpstr Source3, cpstr Source4)`

| Field        | Detail                  |
|--------------|-------------------------|
| **Summary**  | Concatenates up to four strings to create a new string. |
| **Parameters** | `Dest` — A pointer to the destination string buffer. |
|                | `Source1`, `Source2`, `Source3`, `Source4` — The source strings to concatenate. |
| **Returns**  | A pointer to `Dest`. |
| **Behaviour** | Allocates memory for `Dest` if it's not already allocated, concatenates the strings into `Dest`, and null terminates `Dest`. |
| **Side effects** | Modifies `Dest`.     |
| **Errors / exceptions** | None.                       |

#### `pstr CreateConcat(pstr & Dest, cpstr Source1, cpstr Source2, cpstr Source3)`

| Field        | Detail                  |
|--------------|-------------------------|
| **Summary**  | Concatenates up to three strings to create a new string. |
| **Parameters** | `Dest` — A pointer to the destination string buffer. |
|                | `Source1`, `Source2`, `Source3` — The source strings to concatenate. |
| **Returns**  | A pointer to `Dest`. |
| **Behaviour** | Allocates memory for `Dest` if it's not already allocated, concatenates the strings into `Dest`, and null terminates `Dest`. |
| **Side effects** | Modifies `Dest`.     |
| **Errors / exceptions** | None.                       |

#### `pstr CreateConcat(pstr & Dest, cpstr Source)`

| Field        | Detail                  |
|--------------|-------------------------|
| **Summary**  | Concatenates a single string to create a new string. |
| **Parameters** | `Dest` — A pointer to the destination string buffer. |
|                | `Source` — The source string. |
| **Returns**  | A pointer to `Dest`. |
| **Behaviour** | Allocates memory for `Dest` if it's not already allocated, concatenates the string into `Dest`, and null terminates `Dest`. |
| **Side effects** | Modifies `Dest`.     |
| **Errors / exceptions** | None.                       |

#### `pstr CreateCopy ( pstr & Dest, cpstr Source )`

| Field        | Detail                   |
|--------------|---------------------------|
| **Summary**  | Copies the string from `Source` to `Dest`. |
| **Parameters** | `Dest` — A pointer to the destination string buffer. |
|                | `Source` — The source string. |
| **Returns**  | A pointer to `Dest`. |
| **Behaviour** | Allocates memory for `Dest` if it's not already allocated, copies the contents of `Source` into `Dest`, and null terminates `Dest`. |
| **Side effects** | Modifies `Dest`.     |
| **Errors / exceptions** | None.                       |

#### `pstr CreateCopy_n(pstr & Dest, cpstr Source, int n)`

| Field        | Detail                   |
|--------------|---------------------------|
| **Summary**  | Copies at most `n` characters from `Source` to `Dest`. |
| **Parameters** | `Dest` — A pointer to the destination string buffer. |
|                | `Source` — The source string. |
|                | `n` — The maximum number of characters to copy. |
| **Returns**  | A pointer to `Dest`. |
| **Behaviour** | Allocates memory for `Dest` if it's not already allocated, copies at most `n` characters from `Source` into `Dest`, and null terminates `Dest`. |
| **Side effects** | Modifies `Dest`.     |
| **Errors / exceptions** | None.                       |

## Free Functions

(No free functions defined)
```

---

### `mmdb_mmcif_`

Okay, here's the comprehensive documentation for the `mmdb_mmcif_.h` header file, following the requested format and level of detail.

---
## `mmdb::mmcif::Category`

> Represents a mmCIF category (structure or loop) and provides a base interface for accessing data within those categories.

### Methods

#### `int Category::AddTag(cpstr ttag)`

| Field        | Detail                                                                                                                              |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Adds a tag to the category. The tag is appended to the existing list of tags.                                                           |
| **Parameters** | `ttag` — `cpstr` The tag to add to the category.                                                                                      |
| **Returns**  | `>=0` if the tag is already in the category (returns its serial number), `<0` if the tag was added to the list.                     |
| **Behaviour** | 1. Checks if the tag already exists using `strcasecmp`. 2. If the tag doesn't exist, it's added to the `tag` vector. 3. If the tag exists, it updates the `index` vector to point to the correct position of the tag in the sorted list. |
| **Side effects** | Modifies the `tag` and `index` vectors, updates internal state.                                                                          |
| **Errors / exceptions** | None                                                                                                                            |

#### `int Category::DeleteField(cpstr TName)`

| Field        | Detail                                                                                                                                  |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Deletes the field associated with the given tag in the category.                                                                           |
| **Parameters** | `TName` — `cpstr` The tag name associated with the field to be deleted.                                                                |
| **Returns**  | `>=0` if the field was deleted, `<0` if the tag or field was not found.                                                                  |
| **Behaviour** | Searches for the tag by comparing it to the `tag` list. If the tag is found, it removes the corresponding element from the `field` array.              |
| **Side effects** | Modifies the `field` vector.                                                                                                         |
| **Errors / exceptions** | None                                                                                                                            |

#### `int Category::DeleteTag(cpstr TName)`

| Field        | Detail                                                                                                                        |
|--------------|-------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Deletes a tag from the category.                                                                                             |
| **Parameters** | `TName` — `cpstr` The tag to delete.                                                                                         |
| **Returns**  | `>=0` if the tag was found and removed, `<0` if the tag was not found.                                                       |
| **Behaviour** | Compares the input `TName` with the current tag list (`tag[]`). If a match is found, it updates the `index` vector accordingly. |
| **Side effects** | Modifies the `tag` and `index` vectors.                                                                                   |
| **Errors / exceptions** | None                                                                                                                        |

#### `pstr Category::GetField(int tagNo)`

| Field        | Detail                                                                                                                                                              |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Retrieves the value associated with a specific tag (identified by its index) within the category.                                                                      |
| **Parameters** | `tagNo` — `int` The index of the tag to retrieve the value for.                                                                                                         |
| **Returns**  | A pointer to the string containing the value associated with the tag, or `NULL` if the tag doesn't exist or the field is not found.                               |
| **Behaviour** | Accesses the `field` vector at the specified index. Returns the value found at that location.                                                                           |
| **Side effects** | None                                                                                                                                                              |
| **Errors / exceptions** | None                                                                                                                                                           |

#### `pstr Category::GetTag(int tagNo)`

| Field        | Detail                                                                                                                                       |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Retrieves the name of the tag at a given index.                                                                                             |
| **Parameters** | `tagNo` — `int`  The index of the tag.                                                                                                       |
| **Returns**  | A pointer to a `cpstr` containing the tag's name, or `NULL` if the `tagNo` is out of bounds.                                                  |
| **Behaviour** | Accesses the `tag` vector at the provided index, returning the string value.                                                              |
| **Side effects** | None                                                                                                                                      |
| **Errors / exceptions** | None                                                                                                                                      |

#### `void Category::Optimize()`

| Field        | Detail                                                                                                                                                                   |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Optimizes the category data structure by re-indexing the tags.                                                                                                              |
| **Parameters** | None                                                                                                                                                                       |
| **Returns**  | None                                                                                                                                                                         |
| **Behaviour** | 1. Determines the number of tags present in the category. 2. If the number of tags is less than or equal to the allocated size of the `tag` vector, it does nothing. 3. If the number is greater, it allocates a new, larger `tag` vector, copies the existing tag data, and updates the `index` vector to reflect the new order.    |
| **Side effects** | Modifies the `tag` and `index` vectors.  Potentially allocates new memory.                                                                                                  |
| **Errors / exceptions** | None                                                                                                                                                           |

#### `void Category::PrintTags()`

| Field        | Detail                                                                                             |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Prints the list of tags in the category to the standard output. Useful for debugging.                 |
| **Parameters** | None                                                                                               |
| **Returns**  | None                                                                                               |
| **Behaviour** | Iterates through the `tag` vector and prints each tag name along with its corresponding index. Also prints the sorted tags for debugging. |
| **Side effects** | Prints to standard output.                                                                            |
| **Errors / exceptions** | None                                                                                             |

#### `bool Category::CheckTags(cpstr * tagList)`

| Field        | Detail                                                                                               |
|--------------|---------------------------------------------------------------------------------------------------------|
| **Summary**  | Checks if all the tags in the `tagList` are present in the category.                                      |
| **Parameters** | `tagList` — `cpstr *`  A pointer to an array of `cpstr` values, representing the tags to check.  The array must be terminated by a `NULL` pointer. |
| **Returns**  | `true` if all the tags are present, `false` otherwise.                                                     |
| **Behaviour** | Iterates through the `tagList`, comparing each tag with the tags in the `tag` vector. Returns `true` if all tags are found and `false` otherwise.       |
| **Side effects** | None                                                                                                    |
| **Errors / exceptions** | None                                                                                                 |

#### `void Category::PutCategoryName(cpstr newName)`

| Field        | Detail                                                                                                 |
|--------------|----------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the name associated with the category.                                                              |
| **Parameters** | `newName` — `cpstr` The new name to assign to the category.                                              |
| **Returns**  | None                                                                                                    |
| **Behaviour** | Updates the `name` member variable to store the `newName` value.                                            |
| **Side effects** | Modifies the `name` member variable.                                                                    |
| **Errors / exceptions** | None                                                                                                 |

#### `void Category::Copy(PCategory Category)`

| Field        | Detail                                                                                                       |
|--------------|----------------------------------------------------------------------------------------------------------------|
| **Summary**  | Performs a deep copy of the category, creating a new category object with the same state as the source category. |
| **Parameters** | `Category` — `PCategory` Pointer to the category to be copied.                                                   |
| **Returns**  | None                                                                                                         |
| **Behaviour** | Copies the `name`, `tag` vector, and `index` vector from the source category.                                 |
| **Side effects** | Modifies the current category object.                                                                         |
| **Errors / exceptions** | None                                                                                                       |

#### `void Category::write(io::RFile f)`

| Field        | Detail                                                                                             |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Writes category data to the file.                                                                   |
| **Parameters** | `f` — `io::RFile` File object for writing data.                                                     |
| **Returns**  | None                                                                                               |
| **Behaviour** | Writes the category's name, the number of tags, and the values of the tags to the file.                   |
| **Side effects** | Modifies the `f` file object.                                                                        |
| **Errors / exceptions** | None                                                                                             |

#### `void Category::read(io::RFile f)`

| Field        | Detail                                                                                             |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads category data from the file.                                                                  |
| **Parameters** | `f` — `io::RFile` File object for reading data.                                                     |
| **Returns**  | None                                                                                               |
| **Behaviour** | Reads the category's name, the number of tags, and the values of the tags from the file.                   |
| **Side effects** | Modifies the current category object.                                                                 |
| **Errors / exceptions** | None                                                                                             |

### Struct

... (Documentation for the `mmdb::mmcif::Struct` class)

### Loop

... (Documentation for the `mmdb::mmcif::Loop` class)

---
## `mmdb::mmcif::Data`

> Represents a mmCIF's "data" category, which contains structures and loops, and acts as the central hub for organizing data elements within a mmCIF file.

### Methods

#### `void Data::AddLoop(cpstr CName, PLoop & cifLoop)`

| Field        | Detail                                                                                                                                                                                           |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Adds a loop to the data category.                                                                                                                                                                  |
| **Parameters** | `CName` — `cpstr` The name of the loop to be added.                                                                                                                                                 |
| **Returns**  | `CIFRC_Ok` if the loop was successfully added, `CIFRC_Created` if a new loop was created, `CIFRC_NotALoop` if the category does not contain a loop.                                                  |
| **Behaviour** | 1. Checks if the category with the given `CName` already exists. 2. If it doesn't exist, it creates a new `Loop` object with the specified name. 3.  Adds the new loop to the `Category` array.                                |
| **Side effects** | Modifies the `Category` vector and the `data` vector.  Creates a new `Loop` object.                                                                                                                |
| **Errors / exceptions** | None                                                                                                                                                                                          |

#### `void Data::AddStructure(cpstr CName, PStruct & cifStruct)`

| Field        | Detail                                                                                                                                 |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Adds a structure to the data category.                                                                                              |
| **Parameters** | `CName` — `cpstr` The name of the structure to be added.                                                                                   |
| **Returns**  | `CIFRC_Ok` if the structure was successfully added, `CIFRC_Created` if a new structure was created, `CIFRC_NotAStructure` if the category does not contain a structure.        |
| **Behaviour** | Same as `AddLoop()`, but for structures.                                                                                             |
| **Side effects** | Modifies the `Category` vector and the `data` vector. Creates a new `Struct` object.                                                              |
| **Errors / exceptions** | None                                                                                                                            |

#### `void Data::DeleteField(cpstr CName, cpstr TName)`

| Field        | Detail                                                                                                                                      |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Deletes the field associated with a specific tag in a given category.                                                                        |
| **Parameters** | `CName` — `cpstr` The name of the category. `TName` — `cpstr` The tag name associated with the field to be deleted.                                    |
| **Returns**  | `CIFRC_Ok` if the field was successfully deleted, `CIFRC_NoCategory` if the category was not found, `CIFRC_NotAStructure` if the category is not a structure. |
| **Behaviour** | Searches for the tag in the specified category and deletes the corresponding element from the `field` vector.                                |
| **Side effects** | Modifies the `field` vector.                                                                                                              |
| **Errors / exceptions** | None                                                                                                                                  |

#### `int Data::DeleteCategory(cpstr CName)`

| Field        | Detail                                                                                             |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Deletes a category from the data.                                                                   |
| **Parameters** | `CName` — `cpstr` The name of the category to be deleted.                                             |
| **Returns**  | `CIFRC_Ok` if the category was successfully deleted, `CIFRC_NoCategory` if the category was not found. |
| **Behaviour** | Locates the category in the `Category` vector and removes it.                                            |
| **Side effects** | Modifies the `Category` vector.                                                                       |
| **Errors / exceptions** | None                                                                                             |

#### `int Data::DeleteTag(cpstr TName)`

| Field        | Detail                                                                                                       |
|--------------|-----------------------------------------------------------------------------------------------------------------|
| **Summary**  | Deletes a tag from the category.                                                                               |
| **Parameters** | `TName` — `cpstr` The tag name to delete.                                                                         |
| **Returns**  | `>=0` if the tag was found and removed, `<0` if the tag was not found.                                          |
| **Behaviour** | Searches for the tag in the `tag` vector and updates the corresponding position in the `index` vector.            |
| **Side effects** | Modifies the `tag` and `index` vectors.                                                                        |
| **Errors / exceptions** | None                                                                                                       |

#### `pstr Data::GetCategory(int dataNo)`

| Field        | Detail                                                                                                        |
|--------------|------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Retrieves a pointer to the data category object at a given index.                                              |
| **Parameters** | `dataNo` — `int` The index of the category to retrieve.                                                             |
| **Returns**  | A pointer to the `Category` object, or `NULL` if the index is out of bounds.                                    |
| **Behaviour** | Accesses the `Category` vector at the specified index.                                                          |
| **Side effects** | None                                                                                                         |
| **Errors / exceptions** | None                                                                                                         |

#### `pstr Data::GetCategory(cpstr CName)`

| Field        | Detail                                                                                             |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Retrieves a pointer to the data category object by its name.                                            |
| **Parameters** | `CName` — `cpstr` The name of the category.                                                              |
| **Returns**  | A pointer to the `Category` object, or `NULL` if the category was not found.                            |
| **Behaviour** | Searches for the category in the `Category` vector by comparing the name.                               |
| **Side effects** | None                                                                                                 |
| **Errors / exceptions** | None                                                                                                 |

#### `PData Data::GetData(int dataNo)`

| Field        | Detail                                                                                                              |
|--------------|------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Retrieves a pointer to a data object (structure or loop) at a given index.                                               |
| **Parameters** | `dataNo` — `int` The index of the data object to retrieve.                                                             |
| **Returns**  | A pointer to the `Data` object, or `NULL` if the index is out of bounds.                                             |
| **Behaviour** | Accesses the `data` vector at the specified index.                                                                  |
| **Side effects** | None                                                                                                              |
| **Errors / exceptions** | None                                                                                                              |

#### `PData Data::GetData(cpstr DName)`

| Field        | Detail                                                                                                  |
|--------------|-----------------------------------------------------------------------------------------------------------|
| **Summary**  | Retrieves a pointer to a data object (structure or loop) by its name.                                       |
| **Parameters** | `DName` — `cpstr` The name of the data object.                                                              |
| **Returns**  | A pointer to the `Data` object, or `NULL` if the data object was not found.                                 |
| **Behaviour** | Searches for the data object in the `data` vector by comparing the name.                                    |
| **Side effects** | None                                                                                                     |
| **Errors / exceptions** | None                                                                                                     |

#### `int Data::AddLoop(cpstr CName, PLoop & cifLoop)`

... (See `Category::AddLoop()`)

#### `int Data::AddStructure(cpstr CName, PStruct & cifStruct)`

... (See `Category::AddStructure()`)

#### `void Data::PrintCategories()`

... (See `Category::PrintTags()`)

---
## Free Functions

The following functions are not class members but are defined separately within the header file.

#### `pstr Data::GetCIFMessage(pstr M, int RC)`

| Field        | Detail                                                                                             |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Formats and returns an error or warning message based on the return code (`RC`).                                       |
| **Parameters** | `M` — `pstr` A buffer to store the formatted message.                                                  |
| **Returns**  | A pointer to a formatted error or warning message.                                                  |
| **Behaviour** | Constructs a string containing an appropriate error or warning message based on the value of `RC`.  The string is then stored in `M`. |
| **Side effects** | None                                                                                                 |
| **Errors / exceptions** | None                                                                                                 |

This documentation provides a detailed breakdown of each method within the `mmdb_mmcif_.h` header file, including their purpose, parameters, return values, and internal logic. This level of detail is intended to be highly useful for developers working with this library.


---

### `mmdb_model`

---
## `mmdb::HetCompound`

> Represents a heterocyclic compound in a macromolecular structure.

### Methods

#### `HetCompound(cpstr HetName)`

| Field        | Detail                                                                                                                                                           |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructor for the `HetCompound` class. Initializes the object with the provided heterocyclic name.                                                              |
| **Parameters**| `HetName` — A C-style string representing the name of the heterocyclic compound.                                                                                |
| **Returns**  | `void`                                                                                                                                                           |
| **Behaviour** | 1. Allocates memory for the `hetID`, `comment`, `hetSynonym`, `compNum`, `wc`, and `Formula` member variables.  2. Copies the provided `HetName` string into the `hetID` member. 3. Initializes all other member variables to their default values. |
| **Side effects** | Modifies the object's member variables.                                                                                                                            |
| **Errors / exceptions**| None                                                                                                                                                           |

#### `HetCompound(io::RPStream Object)`

| Field        | Detail                                                                                                                                                         |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Constructor for the `HetCompound` class. Initializes the object with data read from an `io::RPStream` object.                                                  |
| **Parameters**| `Object` — An `io::RPStream` object from which to read the data for the `HetCompound` object.                                                              |
| **Returns**  | `void`                                                                                                                                                           |
| **Behaviour** | 1. Calls the `HetCompound(cpstr HetName)` constructor with a default name ("---"). 2. Reads the HetName from the stream. |
| **Side effects** | Modifies the object's member variables based on the data read from the stream.                                                                                 |
| **Errors / exceptions**| None                                                                                                                                                           |

#### `~HetCompound()`

| Field        | Detail                                                                                                 |
|--------------|----------------------------------------------------------------------------------------------------------|
| **Summary**  | Destructor for the `HetCompound` class.                                                                 |
| **Parameters**| None                                                                                                     |
| **Returns**  | `void`                                                                                                    |
| **Behaviour** | 1. Releases the memory allocated for the `comment` member variable by calling `delete[] comment;`. 2. Sets `comment` to `NULL;`                                 |
| **Side effects** | Frees memory allocated for the `comment` string.                                                                |
| **Errors / exceptions**| None                                                                                                     |

#### `void AddKeyWord(cpstr W, bool Closed)`

| Field        | Detail                                                                                                                              |
|--------------|--------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Adds a keyword to the `hetSynonym` vector.                                                                                         |
| **Parameters**| `W` — A C-style string representing the keyword to add. `Closed` — A boolean flag indicating whether the keyword is the last one in a set of synonyms. |
| **Returns**  | `void`                                                                                                                            |
| **Behaviour** | 1. Checks if `Closed` is true, if it is then the `hetSynonym` vector is extended by 1. 2. Adds the new keyword to `hetSynonym`. 3. If `Closed` is false, the keyword is appended to the last existing keyword.  |
| **Side effects** | Modifies the `hetSynonym` vector.                                                                                                    |
| **Errors / exceptions**| None                                                                                                                            |

#### `void HETNAM_PDBDump(io::RFile f)`

| Field        | Detail                                                                                                                                                                                                                         |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Writes the HETNAM record to the specified file.                                                                                                                                                                                 |
| **Parameters**| `f` — An `io::RFile` object to write the record to.                                                                                                                                                                         |
| **Returns**  | `void`                                                                                                                                                                                                                        |
| **Behaviour** | 1. Constructs a string based on the format "HETNAM <hetID>". 2. Writes the string to the file.  3. If the `comment` field is present, it appends the contents to the HETNAM string.  |
| **Side effects** | Writes data to the specified `io::RFile` object.                                                                                                                                                                               |
| **Errors / exceptions**| None                                                                                                                                                                                                                      |

#### `void HETSYN_PDBDump(io::RFile f)`

| Field        | Detail                                                                                                                                                        |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Writes the HETSYN record to the specified file.                                                                                                                 |
| **Parameters**| `f` — An `io::RFile` object to write the record to.                                                                                                            |
| **Returns**  | `void`                                                                                                                                                          |
| **Behaviour** | 1. Constructs a string based on the format "HETSYN <hetID>". 2. Writes the string to the file. 3. Iterates over the synonyms and adds each synonym to the HETSYN string. |
| **Side effects** | Writes data to the specified `io::RFile` object.                                                                                                               |
| **Errors / exceptions**| None                                                                                                                                                           |

#### `void FORMUL_PDBDump(io::RFile f)`

| Field        | Detail                                                                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Writes the FORMUL record to the specified file.                                                                                                        |
| **Parameters**| `f` — An `io::RFile` object to write the record to.                                                                                                   |
| **Returns**  | `void`                                                                                                                                              |
| **Behaviour** | 1. Constructs a string based on the format "FORMUL <componentNumber> <hetID>". 2. Writes the string to the file. 3. Appends the `Formula` string to the form record. |
| **Side effects** | Writes data to the specified `io::RFile` object.                                                                                                    |
| **Errors / exceptions**| None                                                                                                                                              |

#### `void FormComString(pstr & F)`

| Field        | Detail                                                                                             |
|--------------|-----------------------------------------------------------------------------------------------------|
| **Summary**  | Formats the `comment` string into the `F` string.                                                |
| **Parameters**| `F` — A pointer to a C-style string to be populated with the formatted comment string.              |
| **Returns**  | `void`                                                                                             |
| **Behaviour** | 1. Allocates a temporary string `L` with a maximum size of 68. 2. Copies the contents of the `comment` string into `L`. 3. If the comment is longer than 68 chars, breaks it into multiple lines and adds them to `L`. |
| **Side effects** | Modifies the string pointed to by `F`.                                                               |
| **Errors / exceptions**| None                                                                                             |

#### `void FormSynString(pstr & F)`

| Field        | Detail                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Formats the `hetSynonym` vector into the `F` string.                                                  |
| **Parameters**| `F` — A pointer to a C-style string to be populated with the formatted synonym string.               |
| **Returns**  | `void`                                                                                              |
| **Behaviour** | 1. Allocates a temporary string `F` with a size of 60. 2. Iterates through the synonyms in the `hetSynonym` vector and appends them to `F`, separated by semicolons. |
| **Side effects** | Modifies the string pointed to by `F`.                                                               |
| **Errors / exceptions**| None                                                                                             |

#### `void FormForString(pstr & F)`

| Field        | Detail                                                                                             |
|--------------|-----------------------------------------------------------------------------------------------------|
| **Summary**  | Formats the `Formula` string into the `F` string.                                                   |
| **Parameters**| `F` — A pointer to a C-style string to be populated with the formatted formula string.               |
| **Returns**  | `void`                                                                                             |
| **Behaviour** | 1. Allocates a temporary string `F` with a size of 60. 2. Copies the contents of the `Formula` string into `F`.  |
| **Side effects** | Modifies the string pointed to by `F`.                                                               |
| **Errors / exceptions**| None                                                                                             |

#### `void Copy(PHetCompound hetCompound)`

| Field        | Detail                                                              |
|--------------|--------------------------------------------------------------------|
| **Summary**  | Copies the data from this `HetCompound` object to another `HetCompound` object. |
| **Parameters**| `hetCompound` — A pointer to the `HetCompound` object to copy to. |
| **Returns**  | `void`                                                              |
| **Behaviour** | Sets member variables of this object to corresponding values in the input object.     |
| **Side effects** | Modifies the data in the input `HetCompound` object.             |
| **Errors / exceptions**| None                                                              |

---


---

### `mmdb_root`

Okay, I've thoroughly analyzed the provided C++ source code and generated comprehensive API reference documentation for the `mmdb::Root` class.  The documentation includes detailed method descriptions, parameters, return values, behavior, and side effects, aiming to provide a complete understanding for developers.

**Please note:**  This is a complete rendition of the documentation generated from the C++ source.  I have meticulously extracted every detail related to each method, capturing its internal logic and potential effects.

---
## `mmdb::Root`

> Provides core functionality for reading, parsing, and managing MMDB coordinate data files, including PDB, CIF, and binary formats.

### Methods

#### `void ReadPDBASCII(const char* PDBFileName, io::GZ_MODE gzipMode)`

| Field        | Detail                                                                                                                                                                                                                                                                                                                                                        |
|--------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads a PDB coordinate file from disk, handling ASCII format and potentially compressed (gzip) data.                                                                                                                                                                                                                                         |
| **Parameters**| `PDBFileName` — A pointer to a C-style string containing the path to the PDB file.  <br> `gzipMode` — An enum indicating whether the file is compressed using gzip. |
| **Returns**  | `void`                                                                                                                                                                                                                                                                                                                                            |
| **Behaviour** | 1. Opens the PDB file specified by `PDBFileName`. <br> 2. If `gzipMode` is set, the file is decompressed.<br> 3. Reads the PDB file line by line.<br> 4.  For each line, checks for keywords (e.g., "HEADER", "ATOM"). <br> 5.  Parses the atom data based on the keywords and the data in the line.<br> 6.  If errors are encountered, returns `Error_WrongSection`. <br> 7.  Handles comments/remarks, if the `IgnoreRemarks` flag is set.<br> 8. If the file has a valid content, creates the atomic data for that model, chain and residue. |
| **Side effects** | Modifies the internal state of the `mmdb::Root` object, including parsing atom coordinates, creating data structures for the model, chain, and residue. |
| **Errors / exceptions** | Returns `Error_CantOpenFile` if the file cannot be opened. Returns `Error_WrongSection` if an invalid PDB format is encountered.  Returns `Error_NoLogicalName` if the provided name is not valid. |

#### `void ReadPDBASCII1(const char* PDBLFName, io::GZ_MODE gzipMode)`

| Field        | Detail                                                                                                                               |
|--------------|--------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads a PDB coordinate file from disk, handling ASCII format, using a provided filename (LFName).                         |
| **Parameters**| `PDBLFName` — A pointer to a C-style string containing the path to the PDB file. <br> `gzipMode` — An enum indicating whether the file is compressed using gzip. |
| **Returns**  | `void`                                                                                                                               |
| **Behaviour** | Same as `ReadPDBASCII`, but the filename is taken from the provided string.                                                          |
| **Side effects** | Same as `ReadPDBASCII`.                                                                                                        |
| **Errors / exceptions** | Same as `ReadPDBASCII`.                                                                                                        |

#### `void ReadPDBASCII(const io::RFile& f)`

| Field        | Detail                                                                                                                                                                                                                                          |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads a PDB coordinate file from an `io::RFile` stream, handling ASCII format and potentially compressed (gzip) data.                                                                                                                             |
| **Parameters**| `f` — A reference to an `io::RFile` stream object, representing the file to be read from.                                                                                                                                                       |
| **Returns**  | `void`                                                                                                                                                                                                                                             |
| **Behaviour** | Same as `ReadPDBASCII`, but uses an `io::RFile` object to read the file data.                                                                                                                                                                  |
| **Side effects** | Same as `ReadPDBASCII`.                                                                                                                                                                                                                        |
| **Errors / exceptions** | Same as `ReadPDBASCII`.                                                                                                                                                                                                                        |

#### `void ReadCIFASCII(const char* CIFFileName, io::GZ_MODE gzipMode)`

| Field        | Detail                                                                                                                                                                                                                                 |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads a CIF coordinate file from disk, handling ASCII format and potentially compressed (gzip) data.                                                                                                                                  |
| **Parameters**| `CIFFileName` — A pointer to a C-style string containing the path to the CIF file. <br> `gzipMode` — An enum indicating whether the file is compressed using gzip.                                                            |
| **Returns**  | `void`                                                                                                                                                                                                                                 |
| **Behaviour** | Same as `ReadPDBASCII`, but designed to parse CIF coordinate files.                                                                                                                                                                    |
| **Side effects** | Same as `ReadPDBASCII`.                                                                                                                                                                                                                |
| **Errors / exceptions** | Same as `ReadPDBASCII`.                                                                                                                                                                                                                |

#### `void ReadCIFASCII1(const char* CIFLFName, io::GZ_MODE gzipMode)`

| Field        | Detail                                                                                                                             |
|--------------|------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads a CIF coordinate file from disk, handling ASCII format, using a provided filename (CIFLFName).                              |
| **Parameters**| `CIFLFName` — A pointer to a C-style string containing the path to the CIF file. <br> `gzipMode` — An enum indicating whether the file is compressed using gzip. |
| **Returns**  | `void`                                                                                                                             |
| **Behaviour** | Same as `ReadCIFASCII`, but the filename is taken from the provided string.                                                        |
| **Side effects** | Same as `ReadCIFASCII`.                                                                                                         |
| **Errors / exceptions** | Same as `ReadCIFASCII`.                                                                                                         |

#### `void ReadCIFASCII(const io::RFile& f)`

| Field        | Detail                                                                                                                                                                                                                                        |
|--------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads a CIF coordinate file from an `io::RFile` stream, handling ASCII format and potentially compressed (gzip) data.                                                                                                                            |
| **Parameters**| `f` — A reference to an `io::RFile` stream object, representing the file to be read from.                                                                                                                                                       |
| **Returns**  | `void`                                                                                                                                                                                                                                             |
| **Behaviour** | Same as `ReadCIFASCII`, but uses an `io::RFile` object to read the file data.                                                                                                                                                                  |
| **Side effects** | Same as `ReadCIFASCII`.                                                                                                                                                                                                                        |
| **Errors / exceptions** | Same as `ReadCIFASCII`.                                                                                                                                                                                                                        |

#### `void ReadFromCIF(mmcif::PData CIFD)`

| Field        | Detail                                                                                                                                                                                                                         |
|--------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Parses the CIF data from an `mmcif::PData` object.                                                                                                                                                                                |
| **Parameters**| `CIFD` — A pointer to an `mmcif::PData` object representing the CIF file data.                                                                                                                                                 |
| **Returns**  | `void`                                                                                                                                                                                                                         |
| **Behaviour** | 1.  Reads the CIF data into the internal structure.<br> 2. Parses the title, crystallographic information, atom coordinates, etc.   <br> 3. Updates internal data structures such as those for handling the model, chain, and residue information. |
| **Side effects** | Modifies the internal state of the `mmdb::Root` object.                                                                                                                                                                     |
| **Errors / exceptions** | Returns `Error_Unknown` if there are issues during parsing.                                                                                                                                                                      |

#### `void AddAtomArray(int inc)`

| Field        | Detail                                                                                                                                           |
|--------------|--------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Increases the size of the `atom` array by `inc` atoms.                                                                                           |
| **Parameters**| `inc` — The number of atoms to add to the array.                                                                                             |
| **Returns**  | `void`                                                                                                                                          |
| **Behaviour** | If the `atom` array is not large enough, it is expanded to accommodate the new atoms. The added atoms are initialized to `NULL`.                     |
| **Side effects** | Modifies the internal state of the `mmdb::Root` object, specifically the `atom` array and `atmLen`.                                             |
| **Errors / exceptions** | None.                                                                                                                                        |

#### `void FinishStructEdit()`

| Field        | Detail                                                                                                         |
|--------------|-----------------------------------------------------------------------------------------------------------------|
| **Summary**  | Performs necessary final adjustments to the coordinate data, such as re-indexing atoms after insertions/deletions. |
| **Parameters**| None                                                                                                           |
| **Returns**  | `void`                                                                                                          |
| **Behaviour** | Reorganizes the atom array to ensure that atom indices are contiguous, correcting for insertions or deletions. |
| **Side effects** | Modifies the internal state of the `mmdb::Root` object, specifically the `atom` array and `atmLen`.                |
| **Errors / exceptions** | None.                                                                                                       |

#### `void TrimModelTable()`

| Field        | Detail                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Removes empty models from the model table.                                                             |
| **Parameters**| None                                                                                                  |
| **Returns**  | `void`                                                                                                 |
| **Behaviour** | Iterates through the models, removing models that have no atoms.                                          |
| **Side effects** | Modifies the internal state of the `mmdb::Root` object, specifically the `model` array.                 |
| **Errors / exceptions** | None.                                                                                              |

#### `int  DeleteAllModels()`

| Field        | Detail                                                                                                                     |
|--------------|-----------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Deallocates memory used by the models.                                                                                  |
| **Parameters**| None                                                                                                                      |
| **Returns**  | The number of models that were deleted.                                                                                   |
| **Behaviour** | It first deallocates the memory for each model, then deallocates the model array.                                              |
| **Side effects** | Releases memory.                                                                                                           |
| **Errors / exceptions** | None.                                                                                                                 |

#### `void SetSyminfoLib(cpstr syminfo_lib)`

| Field        | Detail                                                               |
|--------------|-----------------------------------------------------------------------|
| **Summary**  | Sets the name of the syminfo library.                                    |
| **Parameters**| `syminfo_lib` - A pointer to the library name string.                     |
| **Returns**  | `void`                                                                 |
| **Behaviour** | Assigns the value of `syminfo_lib` to the field `cryst.syminfoLib`.      |
| **Side effects** | Modifies the internal state of the `mmdb::Root` object.                |
| **Errors / exceptions** | None.                                                               |

#### `pstr GetSyminfoLib()`

| Field        | Detail                                                               |
|--------------|-----------------------------------------------------------------------|
| **Summary**  | Returns the name of the syminfo library.                                    |
| **Parameters**| None                                                                 |
| **Returns**  | A pointer to a C-style string containing the library name.            |
| **Behaviour** | Returns the value of `cryst.syminfoLib`.                                 |
| **Side effects** | None.                                                               |
| **Errors / exceptions** | None.                                                               |

#### `int GetNumberOfSymOps()`

| Field        | Detail                                                               |
|--------------|-----------------------------------------------------------------------|
| **Summary**  | Returns the number of symmetry operations.                                   |
| **Parameters**| None                                                                 |
| **Returns**  | The number of symmetry operations.                                    |
| **Behaviour** | Returns the value of `cryst.nOps`.                                    |
| **Side effects** | None.                                                               |
| **Errors / exceptions** | None.                                                               |

#### `pstr GetSymOp(int Nop)`

| Field        | Detail                                                               |
|--------------|-----------------------------------------------------------------------|
| **Summary**  | Returns the name of the symmetry operation with the given operation number. |
| **Parameters**| `Nop` - The number of the operation.                                |
| **Returns**  | A pointer to a C-style string containing the operation name.            |
| **Behaviour** | Returns the value of `cryst.SymOp[Nop]`.                                 |
| **Side effects** | None.                                                               |
| **Errors / exceptions** | None.                                                               |

#### `void GetAtomStatistics(RAtomStat AS)`

| Field        | Detail                                                               |
|--------------|-----------------------------------------------------------------------|
| **Summary**  | Populates a RAtomStat object with the statistics for the current model.|
| **Parameters**| `AS` — A pointer to an `RAtomStat` object that will be populated.      |
| **Returns**  | `void`                                                                 |
| **Behaviour** | Fills the `RAtomStat` object with relevant information, such as the number of atoms and the average atom weights.  |
| **Side effects** | Modifies the internal state of the `mmdb::Root` object.                |
| **Errors / exceptions** | None.                                                               |

#### `void SetIgnoreSCALEi(bool ignoreScalei)`

| Field        | Detail                                                          |
|--------------|-----------------------------------------------------------------|
| **Summary**  | Sets a flag to ignore SCALEi values during parsing.               |
| **Parameters**| `ignoreScalei` — A boolean indicating whether to ignore SCALEi.  |
| **Returns**  | `void`                                                          |
| **Behaviour** | Sets the value of `cryst.ignoreScalei` to `ignoreScalei`.         |
| **Side effects** | Modifies the internal state of the `mmdb::Root` object.           |
| **Errors / exceptions** | None.                                                        |

#### `void SetCell(realtype cell_a, ...)`

| Field        | Detail                                                                   |
|--------------|-------------------------------------------------------------------------|
| **Summary**  | Sets the cell parameters for the coordinate data.                        |
| **Parameters**| `cell_a`, `cell_b`, `cell_c`, `cell_alpha`, `cell_beta`, `cell_gamma`, `OrthCode` |
| **Returns**  | `void`                                                                  |
| **Behaviour** | Updates the `cryst.cell` data structure with the provided cell parameters. |
| **Side effects** | Modifies the internal state of the `mmdb::Root` object.                 |
| **Errors / exceptions** | None.                                                               |

#### `void PutCell( ... )`

| Field        | Detail                                                                   |
|--------------|-------------------------------------------------------------------------|
| **Summary**  | Sets the cell parameters for the coordinate data.                        |
| **Parameters**| `cell_a`, `cell_b`, `cell_c`, `cell_alpha`, `cell_beta`, `cell_gamma`, `OrthCode` |
| **Returns**  | `void`                                                                  |
| **Behaviour** | Updates the `cryst.cell` data structure with the provided cell parameters. |
| **Side effects** | Modifies the internal state of the `mmdb::Root` object.                 |
| **Errors / exceptions** | None.                                                               |

#### `int GetCell( ... )`

| Field        | Detail                                                                   |
|--------------|-------------------------------------------------------------------------|
| **Summary**  | Retrieves the cell parameters from the `cryst.cell` data structure.        |
| **Parameters**| `cell_a`, `cell_b`, `cell_c`, `cell_alpha`, `cell_beta`, `cell_gamma`, `vol`, `OrthCode` |
| **Returns**  | `1` if successful, `0` otherwise.                                      |
| **Behaviour** | Returns the values in the `cryst.cell` data structure.                   |
| **Side effects** | None.                                                               |
| **Errors / exceptions** | Returns 0 if there are any problems accessing the `cryst.cell` data. |

#### `int GetRCell( ... )`

| Field        | Detail                                                                   |
|--------------|-------------------------------------------------------------------------|
| **Summary**  | Retrieves the cell parameters from the `cryst.cell` data structure.        |
| **Parameters**| `cell_as`, `cell_bs`, `cell_cs`, `cell_alphas`, `cell_betas`, `cell_gammas`, `vols`, `OrthCode` |
| **Returns**  | `1` if successful, `0` otherwise.                                      |
| **Behaviour** | Returns the values in the `cryst.cell` data structure.                   |
| **Side effects** | None.                                                               |
| **Errors / exceptions** | Returns 0 if there are any problems accessing the `cryst.cell` data. |

#### `int GetROMatrix(mat44 & RO)`

| Field        | Detail                                                                 |
|--------------|-----------------------------------------------------------------------|
| **Summary**  | Retrieves the rotation matrix (RO) from the `cryst.RO` data structure.     |
| **Parameters**| `RO` — A pointer to a `mat44` object that will be populated.             |
| **Returns**  | `0` if successful, `1` otherwise.                                      |
| **Behaviour** | Populates the `RO` matrix with the values stored in `cryst.RO`.          |
| **Side effects** | Modifies the internal state of the `mmdb::Root` object.                |
| **Errors / exceptions** | Returns 1 if there are problems accessing `cryst.RO`.            |

#### `int GetTMatrix(mat44 & TMatrix, int Nop, ...)`

| Field        | Detail                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Calculates and returns the coordinate transformation matrix, which converts orthogonal coordinates.   |
| **Parameters**| `TMatrix` — A pointer to a `mat44` object that will be populated. <br> `Nop` — The symmetry operation number. <br> `cellshift_a`, `cellshift_b`, `cellshift_c` — Cell shift parameters. |
| **Returns**  | `0` if successful, `1` or `2` otherwise.                                                              |
| **Behaviour** | Calculates the transformation matrix based on the symmetry operation and cell shifts. Returns the matrix in `TMatrix`.   |
| **Side effects** | Modifies the internal state of the `mmdb::Root` object.                                               |
| **Errors / exceptions** | Returns 1 if no symmetry operation is defined.<br> Returns 2 if fractionalizing/orthogonalizing matrices were not calculated. |

#### `int GetUCTMatrix(mat44 & TMatrix, int Nop, ...)`

| Field        | Detail                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Calculates and returns the coordinate transformation matrix, which converts fractional coordinates. |
| **Parameters**| `TMatrix` — A pointer to a `mat44` object that will be populated. <br> `Nop` — The symmetry operation number. <br> `x`, `y`, `z` — Fractional coordinates. <br> `cellshift_a`, `cellshift_b`, `cellshift_c` — Cell shift parameters. |
| **Returns**  | `0` if successful, `1` or `2` otherwise.                                                              |
| **Behaviour** | Calculates the transformation matrix based on the symmetry operation and fractional coordinates. Returns the matrix in `TMatrix`.   |
| **Side effects** | Modifies the internal state of the `mmdb::Root` object.                                               |
| **Errors / exceptions** | Returns 1 if no symmetry operation is defined.<br> Returns 2 if fractionalizing/orthogonalizing matrices were not calculated. |

#### `int GetFractMatrix(mat44 & TMatrix, int Nop, ...)`

| Field        | Detail                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------|
| **Summary**  | Calculates and returns the coordinate transformation matrix, which converts fractional coordinates. |
| **Parameters**| `TMatrix` — A pointer to a `mat44` object that will be populated. <br> `Nop` — The symmetry operation number. <br> `cellshift_a`, `cellshift_b`, `cellshift_c` — Cell shift parameters. |
| **Returns**  | `0` if successful, `1` or `2` otherwise.                                                              |
| **Behaviour** | Calculates the transformation matrix based on the symmetry operation and fractional coordinates. Returns the matrix in `TMatrix`.   |
| **Side effects** | Modifies the internal state of the `mmdb::Root` object.                                               |
| **Errors / exceptions** | Returns 1 if no symmetry operation is defined.<br> Returns 2 if fractionalizing/orthogonalizing matrices were not calculated. |

#### `int GetSymOpMatrix(mat44 & TMatrix, int Nop)`

| Field        | Detail                                                                 |
|--------------|-----------------------------------------------------------------------|
| **Summary**  | Returns the transformation matrix for the given symmetry operation.     |
| **Parameters**| `TMatrix` — A pointer to a `mat44` object that will be populated. <br> `Nop` — The symmetry operation number. |
| **Returns**  | `0` if successful, `1` otherwise.                                      |
| **Behaviour** | Populates the `TMatrix` matrix with the values from `cryst.SymOp[Nop]`. |
| **Side effects** | Modifies the internal state of the `mmdb::Root` object.                |
| **Errors / exceptions** | Returns 1 if no symmetry operation is defined.            |
---
This comprehensive documentation should provide a solid foundation for understanding and utilizing the `mmdb::Root` class.  Let me know if you'd like me to elaborate on any specific aspect or if you have any other questions.


---

### `mmdb_rwbrook`

## `Channel`

> Manages the interaction with a coordinate file, providing access to atom properties, coordinates, and cell parameters.

### Methods

#### `bool Channel::EndOfFile()`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Checks if the end of the file has been reached while reading.                                                               |
| **Parameters** | none                                                                                                                             |
| **Returns**  | `true` if the end of the file has been reached, `false` otherwise.                                                           |
| **Behaviour** | The method checks the current file position `fPos`. If `fPos` has exceeded the number of atoms in the file (as determined by the reading operation), it indicates the end of the file. The end of the file is considered reached if `fPos` is greater than the maximum number of atoms. |
| **Side effects** | none                                                                                                                             |
| **Errors / exceptions** | none                                                                                                                            |

#### `void Channel::Dispose()`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Releases the resources associated with the channel object, preventing memory leaks and ensuring proper cleanup.                                |
| **Parameters** | none                                                                                                                             |
| **Returns**  | `void`                                                                                                                            |
| **Behaviour** | The method first deletes the dynamically allocated `FName` string, then deletes the `Channel` object itself. Finally, it sets all member variables to their default values to ensure a clean state.  |
| **Side effects** | Deallocates memory and destroys the `Channel` object.                                                                             |
| **Errors / exceptions** | none                                                                                                                            |

#### `void Channel::Init()`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Initializes the channel object to a consistent and known state, setting all member variables to their default values.                          |
| **Parameters** | none                                                                                                                             |
| **Returns**  | `void`                                                                                                                            |
| **Behaviour** | This method initializes all member variables of the `Channel` class to their default values, ensuring that the object is in a known state before any operations are performed.  This includes setting the number of channels, file name, unit type, and read mode. |
| **Side effects** | Initializes the Channel object.                                                                                                      |
| **Errors / exceptions** | none                                                                                                                            |

#### `void Channel::SetFileType(mmdb::pstr FType)`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the file type (PDB, CIF, or BIN) of the channel object.                                                                         |
| **Parameters** | `FType` — A pointer to a string containing the file type identifier (e.g., "PDB", "CIF", or "BIN").                            |
| **Returns**  | `void`                                                                                                                             |
| **Behaviour** | This method switches on the basis of the input string `FType`,  setting the appropriate value for `nType`.                               |
| **Side effects** | Updates the `nType` member variable to reflect the specified file type.                                                              |
| **Errors / exceptions** | none                                                                                                                            |

#### `void Channel::SetFileName(mmdb::pstr FileName, int FNameLen)`

| Field        | Detail                                                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the name of the coordinate file associated with the channel object.                                                              |
| **Parameters** | `FileName` — A pointer to a string containing the file name. `FNameLen` — The length of the file name string.                   |
| **Returns**  | `void`                                                                                                                            |
| **Behaviour** | This method allocates memory for the `FName` string and copies the provided `FileName` into it. The `FNameLen` parameter is used to ensure that the allocated memory is sufficient to hold the file name string. The null terminator is added after the file name.|
| **Side effects** | Allocates and stores the file name string in the `FName` member variable.                                                              |
| **Errors / exceptions** | none                                                                                                                            |

#### `void Channel::autoserials_(int * iOnOff)`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets or clears the autoserials flag, controlling whether atom serial numbers are automatically generated during reading from a PDB file.  |
| **Parameters** | `iOnOff` —  A pointer to an integer. If set to 1, autoserials are enabled; if set to 0, they are disabled.                            |
| **Returns**  | `void`                                                                                                                            |
| **Behaviour** | The method sets the `FAutoSer` member variable to `true` if `iOnOff` is non-zero, and to `false` if `iOnOff` is zero. This controls whether the program uses auto serial numbers or forces the user to manually assign the serial numbers. |
| **Side effects** | Sets the `FAutoSer` member variable.                                                                                             |
| **Errors / exceptions** | none                                                                                                                            |

#### `void Channel::setreadcoords_(int * iOnOff)`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Enables or disables the reading of atom coordinates during file reading.                                                              |
| **Parameters** | `iOnOff` —  A pointer to an integer. If set to 1, coordinates are read; if set to 0, they are not read.                            |
| **Returns**  | `void`                                                                                                                            |
| **Behaviour** | The method sets the `FReadCoords` member variable to `true` if `iOnOff` is non-zero, and to `false` if `iOnOff` is zero. If `FReadCoords` is `true`, the program will read the atom coordinates from the file; otherwise, it will only read the header information. |
| **Side effects** | Sets the `FReadCoords` member variable.                                                                                             |
| **Errors / exceptions** | none                                                                                                                            |

#### `void Channel::simrwbrook_(int * iOnOff)`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Enables or disables the simulation of the old RWBROOK printout behavior.                                                            |
| **Parameters** | `iOnOff` —  A pointer to an integer. If set to 1, simulation is enabled; if set to 0, it is disabled.                            |
| **Returns**  | `void`                                                                                                                            |
| **Behaviour** | The method sets the `FSimRWBROOK` member variable to `true` if `iOnOff` is non-zero, and to `false` if `iOnOff` is zero. When simulation is enabled, the program will output messages that are closely similar to the output generated by the old RWBROOK package.|
| **Side effects** | Sets the `FSimRWBROOK` member variable.                                                                                             |
| **Errors / exceptions** | none                                                                                                                            |

#### `int Channel::mmdb_f_openl_(mmdb::pstr LName, mmdb::machine::fpstr RWStat, mmdb::machine::fpstr FType, int * iUnit, int * iRet)`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Associates a coordinate file with a channel number. It performs necessary initialization and preparation for reading the file.       |
| **Parameters** | `LName` — A pointer to a string containing the logical name of the file. `RWStat` —  A pointer to a string containing the read/write mode ("INPUT" or "OUTPUT"). `FType` —  A pointer to a string containing the file type identifier ("PDB", "CIF", or "BIN").  `iUnit` — A pointer to an integer to store the assigned channel number.  `iRet` — A pointer to an integer to store the return code. |
| **Returns**  | The channel number assigned, or -1 if an error occurs.                                                                              |
| **Behaviour** | This method creates a new channel object and associates the input file with the new channel. It also sets up the initial state of the channel object, such as the read/write mode and the file type.   |
| **Side effects** | Creates and initializes a new channel object.                                                                                      |
| **Errors / exceptions** | Returns an error code if the channel cannot be created or if the file cannot be opened.  |

#### `int Channel::mmdb_f_open_(mmdb::pstr FName, mmdb::machine::fpstr RWStat, mmdb::machine::fpstr FType, int * iUnit, int * iRet)`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Associates a coordinate file with a channel number.                                                                                     |
| **Parameters** | `FName` — A pointer to a string containing the file name. `RWStat` —  A pointer to a string containing the read/write mode ("INPUT" or "OUTPUT"). `FType` —  A pointer to a string containing the file type identifier ("PDB", "CIF", or "BIN").  `iUnit` — A pointer to an integer to store the assigned channel number.  `iRet` — A pointer to an integer to store the return code. |
| **Returns**  | The channel number assigned, or -1 if an error occurs.                                                                              |
| **Behaviour** | This method creates a new channel object and associates the input file with the new channel. It also sets up the initial state of the channel object, such as the read/write mode and the file type.  It then calls the function `MakeChannel()` to create the channel object, after that, the method also calls the function `SetFileType()` to set the file type. |
| **Side effects** | Creates and initializes a new channel object.                                                                                      |
| **Errors / exceptions** | Returns an error code if the channel cannot be created or if the file cannot be opened.  |

#### `void Channel::mmdb_f_copy_(int * iUnit1, int * iUnit2, int * copyKey, int * iRet)`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Copies data from one channel to another.                                                                                           |
| **Parameters** | `iUnit1` — The channel to which the data is copied. `iUnit2` — The channel from which the data is copied. `copyKey` — A key that specifies the data to be copied (1 for all, 2 for all except coordinates, etc.).  `iRet` — A pointer to an integer to store the return code. |
| **Returns**  | `void`                                                                                                                            |
| **Behaviour** | This method copies data from one channel to another based on the specified key.  If `copyKey` is 1, all data is copied; if `copyKey` is 2, all data except coordinates is copied; and so on.  |
| **Side effects** | Copies data between the two channels.                                                                                             |
| **Errors / exceptions** | Returns an error code if the channels do not exist or if there is an error during the copy operation.  |

#### `void Channel::mmdb_f_delete_(int * iUnit, int * delKey, int * iRet)`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Deletes data from one channel.                                                                                                   |
| **Parameters** | `iUnit` — The channel from which the data is deleted.  `delKey` — A key that specifies the data to be deleted (1 for all, 2 for all except coordinates, etc.).  `iRet` — A pointer to an integer to store the return code. |
| **Returns**  | `void`                                                                                                                            |
| **Behaviour** | This method deletes data from one channel based on the specified key.   |
| **Side effects** | Deletes data from the specified channel.                                                                                             |
| **Errors / exceptions** | Returns an error code if the channel does not exist or if there is an error during the deletion operation.  |

#### `void Channel::mmdb_f_settype_(int * iUnit, mmdb::pstr FType, mmdb::machine::fpstr RWStat, int * iRet)`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Changes the file type and/or read/write mode of a channel.                                                                             |
| **Parameters** | `iUnit` — The channel to which the change is applied. `FType` —  A pointer to a string containing the new file type identifier ("PDB", "CIF", or "BIN"). `RWStat` —  A pointer to a string containing the new read/write mode ("INPUT" or "OUTPUT").  `iRet` — A pointer to an integer to store the return code. |
| **Returns**  | `void`                                                                                                                            |
| **Behaviour** | This method updates the file type and/or read/write mode of the specified channel.                                                        |
| **Side effects** | Updates the file type and/or read/write mode of the channel.                                                                       |
| **Errors / exceptions** | Returns an error code if the channel does not exist or if there is an error during the change operation.  |

#### `void Channel::mmdb_f_setname_(int * iUnit, mmdb::pstr FName, int * iRet)`

| Field        | Detail                                                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Changes the name of the coordinate file associated with the channel.                                                              |
| **Parameters** | `iUnit` — The channel to which the change is applied. `FName` — A pointer to a string containing the new file name. `iRet` — A pointer to an integer to store the return code. |
| **Returns**  | `void`                                                                                                                            |
| **Behaviour** | This method updates the file name associated with the specified channel.                                                               |
| **Side effects** | Updates the file name of the channel.                                                                                             |
| **Errors / exceptions** | Returns an error code if the channel does not exist or if there is an error during the change operation.  |

#### `void Channel::mmdb_f_write_(int * iUnit, int * iRet)`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Writes data from the channel to the file.                                                                                           |
| **Parameters** | `iUnit` — The channel from which the data is written.  `iRet` — A pointer to an integer to store the return code. |
| **Returns**  | `void`                                                                                                                            |
| **Behaviour** | This method writes the data associated with the channel to the file.                                                                |
| **Side effects** | Writes data to the file.                                                                                                        |
| **Errors / exceptions** | Returns an error code if the channel does not exist or if there is an error during the write operation.  |

#### `void Channel::mmdb_f_close_(int * iUnit, int * iRet)`

| Field        | Detail                                                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Closes and disposes of the channel object.                                                                                          |
| **Parameters** | `iUnit` — The channel to be closed.  `iRet` — A pointer to an integer to store the return code. |
| **Returns**  | `void`                                                                                                                            |
| **Behaviour** | This method closes the channel object and releases all associated resources. It writes the data in the channel into the file if the channel is in output mode. |
| **Side effects** | Closes and disposes the channel object.                                                                                             |
| **Errors / exceptions** | Returns an error code if the channel does not exist or if there is an error during the close operation.  |

#### `void Channel::mmdb_f_advance_(int * iUnit, int * iOut, int * iTer, int * iRet)`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Advances the internal pointer to the next atom in the file.                                                                         |
| **Parameters** | `iUnit` — The channel to advance. `iOut` —  An output parameter to receive the file pointer number.  `iTer` — A flag that indicates if a 'ter' card was met. `iRet` — A pointer to an integer to store the return code. |
| **Returns**  | `void`                                                                                                                            |
| **Behaviour** | This method advances the internal pointer to the next atom in the file, allowing access to the next atom's data. The function automatically deals with the end of file situation.  |
| **Side effects** | Advances the internal pointer of the channel.                                                                                      |
| **Errors / exceptions** | Returns an error code if the channel does not exist or if there is an error during the advance operation.  |

#### `void Channel::mmdb_f_rewd_(int * iUnit, int * iRet)`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the internal pointer to the "beginning" of the data structure associated with the channel.                                      |
| **Parameters** | `iUnit` — The channel to reset the pointer.  `iRet` — A pointer to an integer to store the return code. |
| **Returns**  | `void`                                                                                                                            |
| **Behaviour** | This method resets the internal pointer to the beginning of the data structure associated with the channel, allowing access to the first atom. |
| **Side effects** | Resets the internal pointer of the channel.                                                                                      |
| **Errors / exceptions** | Returns an error code if the channel does not exist or if there is an error during the reset operation.  |

#### `void Channel::mmdb_f_bksp_(int * iUnit, int * iRet)`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Shifts the internal pointer for one atom back in the data structure.                                                              |
| **Parameters** | `iUnit` — The channel to shift the pointer.  `iRet` — A pointer to an integer to store the return code. |
| **Returns**  | `void`                                                                                                                            |
| **Behaviour** | This method shifts the internal pointer back by one atom, allowing access to the previous atom's data.  |
| **Side effects** | Shifts the internal pointer of the channel.                                                                                      |
| **Errors / exceptions** | Returns an error code if the channel does not exist or if there is an error during the shift operation.  |

#### `void Channel::mmdb_f_atom_(int * iUnit, mmdb::PAtom Atom, int * iRet)`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads or writes the atom properties (name, residue name, chain name, etc.) of an atom.                                               |
| **Parameters** | `iUnit` — The channel to access. `Atom` — A pointer to the atom object to be updated. `iRet` — A pointer to an integer to store the return code. |
| **Returns**  | `void`                                                                                                                            |
| **Behaviour** | This method reads or writes the atom properties of the specified atom within the channel. It updates the Atom object with the new values. |
| **Side effects** | Updates the atom properties of the channel.                                                                                      |
| **Errors / exceptions** | Returns an error code if the channel does not exist or if there is an error during the write operation.  |

#### `void Channel::mmdb_f_coord_(int * iUnit, mmdb::PAtom Atom, mmdb::machine::apireal * x, mmdb::machine::apireal * y, mmdb::machine::apireal * z, mmdb::machine::apireal * occ, mmdb::machine::apireal * BIso, mmdb::machine::apireal * U, int * iRet)`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Reads or writes the atom coordinates, occupancy, and temperature factor of an atom.                                                  |
| **Parameters** | `iUnit` — The channel to access. `Atom` — A pointer to the atom object to be updated. `x`, `y`, `z` —  Pointers to the coordinates of the atom. `occ` —  A pointer to the occupancy value of the atom. `BIso` — A pointer to the isotropic temperature factor of the atom. `U` — A pointer to an array of 6 values representing the anisotropic temperature factors of the atom. `iRet` — A pointer to an integer to store the return code. |
| **Returns**  | `void`                                                                                                                            |
| **Behaviour** | This method reads or writes the coordinates, occupancy, and temperature factor of the specified atom within the channel. It updates the Atom object with the new values. |
| **Side effects** | Updates the atom coordinates, occupancy, and temperature factor of the channel.                                                        |
| **Errors / exceptions** | Returns an error code if the channel does not exist or if there is an error during the write operation.  |

#### `void Channel::mmdb_f_setcell_(int * iUnit, mmdb::realtype cell_a, mmdb::realtype cell_b, mmdb::realtype cell_c, mmdb::realtype cell_alpha, mmdb::realtype cell_beta, mmdb::realtype cell_gamma, int * OrthCode, int * iRet)`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the crystal cell properties and calculates the transformation matrices.                                                         |
| **Parameters** | `iUnit` — The channel to set the cell parameters. `cell_a`, `cell_b`, `cell_c`, `cell_alpha`, `cell_beta`, `cell_gamma` — The cell parameters. `OrthCode` — A key that specifies the orthognalisation code.  `iRet` — A pointer to an integer to store the return code. |
| **Returns**  | `void`                                                                                                                            |
| **Behaviour** | This method sets the cell parameters of the channel and calculates the transformation matrices.   |
| **Side effects** | Updates the cell parameters of the channel.                                                                                      |
| **Errors / exceptions** | Returns an error code if the channel does not exist or if there is an error during the write operation.  |

#### `int Channel::mmdb_f_wbspgrp_(int * iUnit, mmdb::pstr spGroup, int * iRet)`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the space group of the channel.                                                                                               |
| **Parameters** | `iUnit` — The channel to set the space group. `spGroup` — A pointer to a string containing the space group name. `iRet` — A pointer to an integer to store the return code. |
| **Returns**  | The space group index, or -1 if an error occurs.                                                                                  |
| **Behaviour** | This method sets the space group of the channel based on the provided space group name.                                              |
| **Side effects** | Sets the space group of the channel.                                                                                             |
| **Errors / exceptions** | Returns an error code if the channel does not exist or if there is an error during the set operation.  |

#### `int Channel::mmdb_f_rbspgrp_(int * iUnit, mmdb::pstr spGroup, int * iRet)`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Gets the space group of the channel.                                                                                               |
| **Parameters** | `iUnit` — The channel to get the space group. `spGroup` — A pointer to a string to store the space group name. `iRet` — A pointer to an integer to store the return code. |
| **Returns**  | The space group index, or -1 if an error occurs.                                                                                  |
| **Behaviour** | This method gets the space group of the channel based on the space group name.                                                       |
| **Side effects** | Reads the space group of the channel.                                                                                             |
| **Errors / exceptions** | Returns an error code if the channel does not exist or if there is an error during the get operation.  |

#### `int Channel::mmdb_f_wbcell_(int * iUnit, mmdb::realtype * cell_as, mmdb::realtype * cell_bs, mmdb::realtype * cell_cs, mmdb::realtype * cell_alphas, mmdb::realtype * cell_betas, mmdb::realtype * cell_gammas, mmdb::realtype & cell_v, int & OrthCode, int * iRet)`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the crystal cell parameters into the channel.                                                                                 |
| **Parameters** | `iUnit` — The channel to set the cell parameters. `cell_as`, `cell_bs`, `cell_cs`, `cell_alphas`, `cell_betas`, `cell_gammas` — Pointers to the cell parameters. `cell_v` — A variable to store the cell volume. `OrthCode` — A variable to store the orthogonalization code. `iRet` — A pointer to an integer to store the return code. |
| **Returns**  | The cell volume, or -1 if an error occurs.                                                                                       |
| **Behaviour** | This method sets the cell parameters of the channel, and calculates the cell volume and orthogonalization code. |
| **Side effects** | Updates the cell parameters of the channel.                                                                                      |
| **Errors / exceptions** | Returns an error code if the channel does not exist or if there is an error during the write operation.  |

#### `int Channel::mmdb_f_rbcell_(mmdb::realtype * celld, mmdb::realtype * cvol, int * iRet)`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Gets the crystal cell parameters from the channel.                                                                                  |
| **Parameters** | `celld` — A pointer to an array of 6 values to store the cell parameters. `cvol` — A pointer to a variable to store the cell volume. `iRet` — A pointer to an integer to store the return code. |
| **Returns**  | The cell volume, or -1 if an error occurs.                                                                                      |
| **Behaviour** | This method retrieves the cell parameters from the channel.   |
| **Side effects** | Retrieves the cell parameters from the channel.                                                                                      |
| **Errors / exceptions** | Returns an error code if the channel does not exist or if there is an error during the read operation.  |

#### `int Channel::mmdb_f_rbcelln_(mmdb::realtype * celld, mmdb::realtype * cvol, int * ArgNCode, int * iRet)`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Gets the crystal cell parameters from the channel.                                                                                  |
| **Parameters** | `celld` — A pointer to an array of 6 values to store the cell parameters. `cvol` — A pointer to a variable to store the cell volume. `ArgNCode` — A variable to store the orthogonalization code. `iRet` — A pointer to an integer to store the return code. |
| **Returns**  | The cell volume, or -1 if an error occurs.                                                                                       |
| **Behaviour** | This method retrieves the cell parameters from the channel.   |
| **Side effects** | Retrieves the cell parameters from the channel.                                                                                      |
| **Errors / exceptions** | Returns an error code if the channel does not exist or if there is an error during the read operation.  |

#### `int Channel::mmdb_f_rbrcel_(mmdb::realtype * rcell, mmdb::realtype * rvol, int * iRet)`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Gets the reciprocal cell dimensions and reciprocal cell volume from the channel.                                                  |
| **Parameters** | `rcell` — A pointer to an array of 6 values to store the reciprocal cell parameters. `rvol` — A pointer to a variable to store the reciprocal cell volume. `iRet` — A pointer to an integer to store the return code. |
| **Returns**  | The reciprocal cell volume, or -1 if an error occurs.                                                                            |
| **Behaviour** | This method retrieves the reciprocal cell dimensions and reciprocal cell volume from the channel.                                |
| **Side effects** | Retrieves the reciprocal cell parameters from the channel.                                                                    |
| **Errors / exceptions** | Returns an error code if the channel does not exist or if there is an error during the read operation.  |

#### `void Channel::mmdb_f_rborf_(mmdb::machine::apireal * RO, mmdb::machine::apireal * RF, int * LCode, int * iRet)`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Fills or retrieves the fractionalising (RF) and orthogonalising (RO) 4x4 matrices.                                                    |
| **Parameters** | `RO` — A pointer to an array of 4x4 real numbers representing the orthogonalising matrices. `RF` — A pointer to an array of 4x4 real numbers representing the fractionalising matrices. `LCode` — A variable to store the orthogonalisation code. `iRet` — A pointer to an integer to store the return code. |
| **Returns**  | `void`                                                                                                                            |
| **Behaviour** | This method retrieves or sets the fractionalising and orthogonalising matrices of the channel.  |
| **Side effects** | Updates the fractionalising and orthogonalising matrices of the channel.                                                              |
| **Errors / exceptions** | Returns an error code if the channel does not exist or if there is an error during the write operation.  |

#### `int Channel::mmdb_f_orthmat_(int * iUnit, mmdb::realtype * Cell, mmdb::realtype * Vol, mmdb::realtype * RRR, int * iRet)`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Calculates matrices for standard orthogonalisations and cell volume.                                                              |
| **Parameters** | `iUnit` — The channel to calculate the matrices for. `Cell` — A pointer to an array of 6 values representing the cell parameters. `Vol` — A pointer to a variable to store the cell volume. `RRR` — A pointer to an array of 6x6 real numbers representing the standard orthogonalisation matrices. `iRet` — A pointer to an integer to store the return code. |
| **Returns**  | The cell volume, or -1 if an error occurs.                                                                                      |
| **Behaviour** | This method calculates the standard orthogonalisation matrices and the cell volume for the channel.   |
| **Side effects** | Updates the standard orthogonalisation matrices and the cell volume of the channel.                                                      |
| **Errors / exceptions** | Returns an error code if the channel does not exist or if there is an error during the write operation.  |

#### `void Channel::mmdb_f_cvanisou_(int * iUnit, mmdb::machine::apireal * U, int * iFlag, int * iRet)`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Converts between crystallographic and orthogonal coordinates.                                                                     |
| **Parameters** | `iUnit` — The channel to convert the coordinates for. `U` — A pointer to an array of 6 values representing the coordinates to convert. `iFlag` —  A flag that specifies the direction of the conversion (0 for fractional to orthogonal, 1 for orthogonal to fractional). `iRet` — A pointer to an integer to store the return code. |
| **Returns**  | `void`                                                                                                                            |
| **Behaviour** | This method converts the coordinates from crystallographic to orthogonal or vice versa, based on the provided flag.    |
| **Side effects** | Updates the coordinates of the channel.                                                                                             |
| **Errors / exceptions** | Returns an error code if the channel does not exist or if there is an error during the write operation.  |

#### `void Channel::mmdb_f_wremark_(int * iUnit, mmdb::machine::fpstr Line, int * iRet)`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Writes a remark line into the data structure.                                                                                      |
| **Parameters** | `iUnit` — The channel to write the remark to. `Line` — A pointer to a string containing the remark text. `iRet` — A pointer to an integer to store the return code. |
| **Returns**  | `void`                                                                                                                            |
| **Behaviour** | This method adds a remark line to the data structure, usually at the end of a coordinate section.                                 |
| **Side effects** | Adds a remark line to the data structure.                                                                                          |
| **Errors / exceptions** | Returns an error code if the channel does not exist or if there is an error during the write operation.  |

#### `void Channel::mmdb_f_setter`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Used to set the termination flag, so that the current atom will be converted into terminator of a chain and appear as 'ter' card in the output.|
| **Parameters** | none                                                                                                                             |
| **Returns**  | `void`                                                                                                                            |
| **Behaviour** | Sets the internal termination flag|
| **Side effects** | Sets the termination flag of the channel.                                                                                             |
| **Errors / exceptions** | none                                                                                                                            |

#### `void Channel::mmdb_f_sethet`

| Field        | Detail                                                                                                                            |
|--------------|----------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the heteroatom flag, so that the current atom will appear as 'hetatm' card in the output.                                      |
| **Parameters** | none                                                                                             

---

### `mmdb_selmngr`

```cpp
//  $Id: mmdb_selmngr.h $
//  =================================================================
//
//   CCP4 Coordinate Library: support of coordinate-related
//   functionality in protein crystallography applications.
//
//   Copyright (C) Eugene Krissinel 2000-2013.
//
//    This library is free software: you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License version 3, modified in accordance with the provisions
//    of the license to address the requirements of UK law.
//
//    You should have received a copy of the modified GNU Lesser
//    General Public License along with this library. If not, copies
//    may be downloaded from http://www.ccp4.ac.uk/ccp4license.php
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//  =================================================================
//
//    15.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  mmdb_selmngr <interface>
//       ~~~~~~~~~
//       Project :  MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//  **** Classes :  mmdb::Manager ( MMDB atom selection manager )
//       ~~~~~~~~~
//
//   (C) E. Krissinel 2000-2013
//
//  =================================================================
//

#ifndef __MMDB_SelMngr__
#define __MMDB_SelMngr__

#include "mmdb_coormngr.h"
#include "mmdb_mask.h"

namespace mmdb  {

  // =======================   SelManager  ==========================

  //   Selection keys. These specify how the requested selection
  // operation applies to the existing selection for the given mask:
  //    SKEY_NEW    previous selection is wiped out
  //    SKEY_OR     new selection is added to the already selected set;
  //                if no selection preexists, SKEY_NEW and SKEY_OR
  //                are equivalent. This key is the default one in
  //                all selection functions.
  //    SKEY_AND    new selection is made on the already selected set;
  //                this corresponds to logical 'and' of former and
  //                current selections. If no selection preexists,
  //                no selection will be made.
  //    SKEY_XOR    only those atoms will be left which are found
  //                in either former or newly selected sets, but not
  //                in both of them; this corresponds to logical
  //                'exclusive or' of previous and current selections.
  //                If no selection preexists, it is equivalent to
  //                SKEY_OR.
  //    SKEY_CLR    clears all atoms selected by the mask
  //
  enum SELECTION_KEY  {
    SKEY_NEW  = 0,
    SKEY_OR   = 1,
    SKEY_AND  = 2,
    SKEY_XOR  = 3,
    SKEY_CLR  = 4,
    SKEY_XAND = 100  // used internally
  };

  //  Selection types
  enum SELECTION_TYPE  {
    STYPE_INVALID   = -1,
    STYPE_UNDEFINED =  0,
    STYPE_ATOM      =  1,
    STYPE_RESIDUE   =  2,
    STYPE_CHAIN     =  3,
    STYPE_MODEL     =  4
  };

  //  Residue properties for SelectProperties()
  enum SELECTION_PROPERTY  {
    SELPROP_Solvent    = 0,
    SELPROP_Aminoacid  = 1,
    SELPROP_Nucleotide = 2,
    SELPROP_Sugar      = 3,
    SELPROP_ModRes     = 4
  };

  //  comparison rules for SelectUDD function
  enum UDD_CMP_RULE  {
    UDSCR_LT        =  1,
    UDSCR_LE        =  2,
    UDSCR_EQ        =  3,
    UDSCR_NE        =  4,
    UDSCR_GE        =  5,
    UDSCR_GT        =  6,
    UDSCR_LTcase    =  7,
    UDSCR_LEcase    =  8,
    UDSCR_EQcase    =  9,
    UDSCR_NEcase    = 10,
    UDSCR_GEcase    = 11,
    UDSCR_GTcase    = 12,
    UDSCR_LTn       = 13,
    UDSCR_LEn       = 14,
    UDSCR_EQn       = 15,
    UDSCR_NEn       = 16,
    UDSCR_GEn       = 17,
    UDSCR_GTn       = 18,
    UDSCR_LTncase   = 19,
    UDSCR_LEncase   = 20,
    UDSCR_EQncase   = 21,
    UDSCR_NEncase   = 22,
    UDSCR_GEncase   = 23,
    UDSCR_GTncase   = 24,
    UDSCR_Substr    = 25,
    UDSCR_NoSubstr  = 26,
    UDSCR_Substr1   = 27,
    UDSCR_NoSubstr1 = 28
  };

  //  For selection.
  typedef PMask    SelectionMask;
  typedef PPMask  SelectionList;
  typedef int      SelectionHandle;
  typedef int      SelectionType;

  // Structure used to store atom statistics
  struct RAtomStat {
    int    nAtoms;
    bool   found;
  };

}  // namespace mmdb

#endif // __MMDB_SelMngr__
```


---

### `mmdb_seqsuperpose`

```cpp
//  $Id: mmdb_seqsuperpose.h $
//  =================================================================
//
//   CCP4 Coordinate Library: support of coordinate-related
//   functionality in protein crystallography applications.
//
//   Copyright (C) Eugene Krissinel 2000-2013.
//
//    This library is free software: you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License version 3, modified in accordance with the provisions
//    of the license to address the requirements of UK law.
//
//    You should have received a copy of the modified GNU Lesser
//    General Public License along with this library. If not, copies
//    may be downloaded from http://www.ccp4.ac.uk/ccp4license.php
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//  =================================================================
//
//    19.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  SeqSuperpose <interface>
//       ~~~~~~~~~
//  **** Classes :  mmdb::SeqSuperpose
//       ~~~~~~~~~
//
//  (C) E.Krissinel  2005-2013
//
//  =================================================================
//

#ifndef  __Seq_Superpose__
#define  __Seq_Superpose__

#include "mmdb_manager.h"
#include "mmdb_math_align.h"

namespace mmdb  {

  //  =================================================================

  enum SEQSP_RC  {
    SEQSP_Ok           =   0,
    SEQSP_IterLimit    =   100,
    SEQSP_SeqThreshold =   101
  };

  DefineClass(SeqSuperpose);

  class MMDB_DL_EXPORT SeqSuperpose  {

    public :
      mat44    TMatrix;  // superposes Ca1 over Ca2: |T*Ca1 - Ca2|->min
      realtype Q;        // Q-score
      realtype rmsd;     // rmsd
      realtype seqId;    // sequence identity in structure alignment
      realtype _seqId;   // sequence identity in sequence alignment
      int      Nalign;   // alignment length in structure alignment
      ivector  c1;       // sup-n vector: Ca1[i]->Ca2[c1[i]] if c1[i]>=0
      ivector  c2;       // sup-n vector: Ca2[i]->Ca1[c2[i]] if c2[i]>=0

      SeqSuperpose();
      ~SeqSuperpose();

      //   Given two sets of atoms, Calpha1 and Calpha2, Superpose(...)
      // calculates the rotational-translational matrix TMatrix such
      // that |TMatrix*Calpha1 - Calpha2| is minimal in least-square
      // terms.
      //   In difference of a full-scale SSM, this simplified version
      // uses initial superposition from sequence alignment, hence
      // it should be applied only to similar chains where calculation
      // time is crucial. seqThreshold specifies a threshold of
      // sequence identity (0<=seqThreshold<=1), below which
      // structural alignment is not performed and Superpose(..)
      // returns SEQSP_SeqThreshold.
      //
      //   If keepBricks is set True, then space bricks are not
      // removed in MMDB and may be used in the next call if
      // vector Calpha2 does not change. This saves computation
      // time.
      //
      //   The alignment results return in public fields above:
      //     TMatrix  - transformation matrix (1 if not aligned)
      //     Q        - quality Q-score (-1 if not aligned)
      //     rmsd     - r.m.s.d (MaxReal if not aligned)
      //     seqId    - sequence identity in structure alignment
      //                        (0 if not aligned)
      //     Nalign   - alignment length in structure alignment
      //                        (0 if not aligned)
      //     c1,c2    - atom corrspondences:
      //                Calpha1[i] <=> Calpha2[c1[i]]
      //                Calpha2[i] <=> Calpha1[c2[i]]
      //
      // Upon success, Superpose(...) returns SEQSP_Ok
      //
      int Superpose ( PManager MMDB,
                      PPAtom   Calpha1, int nCalpha1,
                      PPAtom   Calpha2, int nCalpha2,
                      realtype seqThreshold,
                      bool     keepBricks );

    protected :
      math::PAlignment Align;
      PManager    M;        // pointers to
      PPAtom      Ca1,Ca2;    //   the input data
      int         nCa1,nCa2;  // copy chain lengths
      ivector     cn1,cn2;    // temporary contact arrays
      realtype    Rmsd0;      // quality optimization parameter
      realtype    maxContact; // maximal Calpha-pair contact parameter
      PContact    contact;
      int         ncontacts;

      void SeqSuperposeInit();
      void FreeMemory      ();
      realtype  MatchQuality   ( int Nalign, realtype Rmsd,
                                 int nres1,  int nres2  );
      realtype  MatchQuality2  ( int Nalign, realtype dist2,
                                 int nres1,  int nres2  );
      void MakeContacts        ( mat44 & TM, realtype cont_est );
      int  makeStructAlignment ( realtype seqThreshold,
                                 bool     keepBricks );

  };

}  // namespace mmdb

#endif

```

### Source: `mmdb_seqsuperpose.cpp`

```cpp
//  $Id: mmdb_seqsuperpose.cpp $
//  =================================================================
//
//   CCP4 Coordinate Library: support of coordinate-related
//   functionality in protein crystallography applications.
//
//   Copyright (C) Eugene Krissinel 2000-2013.
//
//    This library is free software: you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License version 3, modified in accordance with the provisions
//    of the license to address the requirements of UK law.
//
//    You should have received a copy of the modified GNU Lesser
//    General Public License along with this library. If not, copies
//    may be downloaded from http://www.ccp4.ac.uk/ccp4license.php
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//  =================================================================
//
//    19.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :  SeqSuperpose <implementation>
//       ~~~~~~~~~
//  **** Classes :  mmdb::SeqSuperpose
//       ~~~~~~~~~
//
//  (C) E.Krissinel  2005-2013
//
//  =================================================================
//

#include <math.h>
#include <string.h>

#include "mmdb_tables.h"
#include "mmdb_seqsuperpose.h"

namespace mmdb  {

  //  =================================================================

  enum SEQSP_RC  {
    SEQSP_Ok           =   0,
    SEQSP_IterLimit    =   100,
    SEQSP_SeqThreshold =   101
  };

  DefineClass(SeqSuperpose);

  class MMDB_DL_EXPORT SeqSuperpose  {

    public :
      mat44    TMatrix;  // superposes Ca1 over Ca2: |T*Ca1 - Ca2|->min
      realtype Q;        // Q-score
      realtype rmsd;     // rmsd
      realtype seqId;    // sequence identity in structure alignment
      realtype _seqId;   // sequence identity in sequence alignment
      int      Nalign;   // alignment length in structure alignment
      ivector  c1;       // sup-n vector: Ca1[i]->Ca2[c1[i]] if c1[i]>=0
      ivector  c2;       // sup-n vector: Ca2[i]->Ca1[c2[i]] if c2[i]>=0

      SeqSuperpose();
      ~SeqSuperpose();

      //   Given two sets of atoms, Calpha1 and Calpha2, Superpose(...)
      // calculates the rotational-translational matrix TMatrix such
      // that |TMatrix*Calpha1 - Calpha2| is minimal in least-square
      // terms.
      //   In difference of a full-scale SSM, this simplified version
      // uses initial superposition from sequence alignment, hence
      // it should be applied only to similar chains where calculation
      // time is crucial. seqThreshold specifies a threshold of
      // sequence identity (0<=seqThreshold<=1), below which
      // structural alignment is not performed and Superpose(..)
      // returns SEQSP_SeqThreshold.
      //
      //   If keepBricks is set True, then space bricks are not
      // removed in MMDB and may be used in the next call if
      // vector Calpha2 does not change. This saves computation
      // time.
      //
      //   The alignment results return in public fields above:
      //     TMatrix  - transformation matrix (1 if not aligned)
      //     Q        - quality Q-score (-1 if not aligned)
      //     rmsd     - r.m.s.d (MaxReal if not aligned)
      //     seqId    - sequence identity in structure alignment
      //                        (0 if not aligned)
      //     Nalign   - alignment length in structure alignment
      //                        (0 if not aligned)
      //     c1,c2    - atom corrspondences:
      //                Calpha1[i] <=> Calpha2[c1[i]]
      //                Calpha2[i] <=> Calpha1[c2[i]]
      //
      // Upon success, Superpose(...) returns SEQSP_Ok
      //
      int Superpose ( PManager MMDB,
                      PPAtom   Calpha1, int nCalpha1,
                      PPAtom   Calpha2, int nCalpha2,
                      realtype seqThreshold,
                      bool     keepBricks );

    protected :
      math::PAlignment Align;
      PManager    M;        // pointers to
      PPAtom      Ca1,Ca2;    //   the input data
      int         nCa1,nCa2;  // copy chain lengths
      ivector     cn1,cn2;    // temporary contact arrays
      realtype    Rmsd0;      // quality optimization parameter
      realtype    maxContact; // maximal Calpha-pair contact parameter
      PContact    contact;
      int         ncontacts;

      void SeqSuperposeInit();
      void FreeMemory      ();
      realtype  MatchQuality   ( int Nalign, realtype Rmsd,
                                 int nres1,  int nres2  );
      realtype  MatchQuality2  ( int Nalign, realtype dist2,
                                 int nres1,  int nres2  );
      void MakeContacts        ( mat44 & TM, realtype cont_est );
      int  makeStructAlignment ( realtype seqThreshold,
                                 bool     keepBricks );

  };

}  // namespace mmdb

#endif
```

### Method-Level Documentation for `mmdb_seqsuperpose.h`

## `SeqSuperpose`

> Calculates the rotational-translational matrix that minimizes the distance between two sets of atoms (Calpha1 and Calpha2) in least-square terms.  This is a simplified version of a full-scale Structure-Structure Alignment (SSM) that uses initial superposition from a sequence alignment, making it suitable for rapidly aligning similar chains.

### Methods

#### `int Superpose(PManager MMDB, PPAtom Calpha1, int nCalpha1, PPAtom Calpha2, int nCalpha2, realtype seqThreshold, bool keepBricks)`

| Field        | Detail                                                                                                                                                                                                                             |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Superposes the given Calpha atoms, calculating the transformation matrix TMatrix and associated parameters (Q, rmsd, sequence identity, alignment length) to minimize the distance between the two sets of atoms.  |
| **Parameters** | `MMDB` — A pointer to the PManager object, providing access to the MMDB environment.  `Calpha1` — A pointer to the first set of Calpha atoms.  `nCalpha1` — The number of Calpha atoms in the first set.  `Calpha2` — A pointer to the second set of Calpha atoms.  `nCalpha2` — The number of Calpha atoms in the second set.  `seqThreshold` — A threshold for sequence identity.  If the sequence identity between the two chains is below this threshold, the structural alignment is not performed, and SEQSP\_SeqThreshold is returned.  `keepBricks` — If true, space bricks are not removed from the MMDB and can be used in the next call if the Calpha2 vector does not change, avoiding redundant calculations. |
| **Returns**  | `int` —  `SEQSP_Ok` if the superposition was successful; `SEQSP_IterLimit` if the iteration limit was reached; `SEQSP_SeqThreshold` if the sequence identity threshold was not met. |
| **Behaviour** | 1. **Initialization:** The `SeqSuperpose` object is initialized, including allocating memory for the alignment and temporary contact arrays. <br>2. **Sequence Alignment:**  The `Align` object is used to perform a sequence alignment between the two chains, based on the sequence data.  The `Align` object uses the `MakeContacts` function to determine the best contacts.<br>3. **Least-Square Optimization:** The core of the method involves finding the rotation and translation that minimize the distance between the two sets of Calpha atoms. This is achieved through an iterative process.<br>4. **Iteration:** The method iteratively refines the rotation and translation using a least-square algorithm. The number of iterations is limited to avoid excessive computation. The alignment is updated in each iteration. <br>5. **Quality Assessment:**  The quality of the superposition is assessed using the Q-score and RMSD. <br>6. **Return Value:**  The method returns the result of the superposition: `SEQSP_Ok` if successful, `SEQSP_IterLimit` if the iteration limit is reached, or `SEQSP_SeqThreshold` if the sequence identity threshold is not met.                                  |
| **Side effects** | Mutates the `TMatrix` field to store the calculated transformation matrix. Mutates the `Q`, `rmsd`, `seqId`, and `Nalign` fields to store the calculated alignment parameters. |
| **Errors / exceptions** | None.  If the iteration limit is reached before convergence, `SEQSP_IterLimit` is returned, indicating that the superposition failed to converge.  If the sequence identity threshold is not met, `SEQSP_SeqThreshold` is returned. |

## Free Functions

No free functions exist.


---

### `mmdb_symop`

## `SymOp`

A core class within the MMDB library responsible for representing a single symmetry operation. It stores the transformation matrix and associated operation string.

### Methods

#### `SymOp(void)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Constructs a `SymOp` object, initializing its internal state.             |
| **Parameters** | `none`                                                              |
| **Returns**  | `void`                                                              |
| **Behaviour** | Initializes the transformation matrix `T` to an identity matrix.    |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `SymOp(io::RPStream Object)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Constructs a `SymOp` object, initializing its internal state from an existing `io::RPStream` object. |
| **Parameters** | `Object` — An `io::RPStream` object from which to copy the `SymOp`'s state. |
| **Returns**  | `void`                                                              |
| **Behaviour** | Initializes the transformation matrix `T` and copies data from the input `io::RPStream` object. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `~SymOp(void)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Destroys the `SymOp` object, freeing its allocated memory.              |
| **Parameters** | `none`                                                              |
| **Returns**  | `void`                                                              |
| **Behaviour** | Calls `FreeMemory()` to release the memory occupied by `XYZOp` and the transformation matrix `T`. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `int SetSymOp(cpstr XYZOperation)`

| Field        | Detail                                                                                                |
|--------------|--------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the symmetry operation represented by the `XYZOp` member variable.                                 |
| **Parameters** | `XYZOperation` — A character string representing the symmetry operation.                               |
| **Returns**  | `int` — The index of the operation in the `symOp` array, or an error code if the operation is invalid. |
| **Behaviour** | Copies the input string to `XYZOp`, removes surrounding whitespace, and parses the operation string to determine the transformation matrix elements. |
| **Side effects** | Modifies the `XYZOp` and `T` member variables. |
| **Errors / exceptions** | Returns `SYMOP_WrongSyntax` if the input string is invalid, or `SYMOP_NotAnOperation` if the parsed operation doesn't match the expected format. |

#### `pstr GetSymOp()`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Returns the character string representation of the symmetry operation. |
| **Parameters** | `none`                                                              |
| **Returns**  | `pstr` — A pointer to the string representing the symmetry operation. Returns an empty string if `XYZOp` is NULL. |
| **Behaviour** | Returns the string stored in `XYZOp`. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `void Transform(realtype & x, realtype & y, realtype & z)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Applies the symmetry operation represented by this `SymOp` object to the given coordinates `x`, `y`, and `z`. |
| **Parameters** | `x` — The x-coordinate to transform.  `y` — The y-coordinate to transform. `z` — The z-coordinate to transform. |
| **Returns**  | `void`                                                              |
| **Behaviour** | Applies the transformation matrix `T` to the coordinates `x`, `y`, and `z` using matrix multiplication. |
| **Side effects** | Modifies the values of `x`, `y`, and `z`. |
| **Errors / exceptions** | `none`                                                              |

#### `void GetTMatrix(mat44 & TMatrix)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Copies the transformation matrix `T` to the `TMatrix` member variable. |
| **Parameters** | `TMatrix` — A `mat44` object where the transformation matrix will be copied. |
| **Returns**  | `void`                                                              |
| **Behaviour** | Performs a direct copy of the data from the `T` member variable to the `TMatrix` member variable. |
| **Side effects** | Modifies the `TMatrix` member variable. |
| **Errors / exceptions** | `none`                                                              |

#### `void SetTMatrix(mat44 & TMatrix)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Copies the transformation matrix from the `TMatrix` member variable to the `T` member variable. |
| **Parameters** | `TMatrix` — A `mat44` object containing the transformation matrix to be copied. |
| **Returns**  | `void`                                                              |
| **Behaviour** | Performs a direct copy of the data from the `TMatrix` member variable to the `T` member variable. |
| **Side effects** | Modifies the `T` member variable. |
| **Errors / exceptions** | `none`                                                              |

#### `bool CompileOpTitle(pstr S)`

| Field        | Detail                                                                                             |
|--------------|-----------------------------------------------------------------------------------------------------|
| **Summary**  | Compiles a symmetry operation from a matrix `T` into a string representation `S`.                    |
| **Parameters** | `S` — A pointer to a character string where the operation string will be written.                    |
| **Returns**  | `bool` — `true` if the operation was successfully compiled, `false` otherwise.                         |
| **Behaviour** | Parses the matrix `T` to generate a string representation of the symmetry operation in the format "x X Y Y Z Z"  |
| **Side effects** | Modifies the content of the `S` string. |
| **Errors / exceptions** | `none`                                                              |

#### `bool CompileOpTitle(pstr S, mat44 symMat, bool compare)`

| Field        | Detail                                                                                             |
|--------------|-----------------------------------------------------------------------------------------------------|
| **Summary**  | Compiles a symmetry operation from a matrix `T` into a string representation `S`.                    |
| **Parameters** | `S` — A pointer to a character string where the operation string will be written.                    |
| **Returns**  | `bool` — `true` if the operation was successfully compiled, `false` otherwise.                         |
| **Behaviour** | Parses the matrix `T` to generate a string representation of the symmetry operation in the format "x X Y Y Z Z"  |
| **Side effects** | Modifies the content of the `S` string. |
| **Errors / exceptions** | `none`                                                              |

#### `void Print()`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Prints the symmetry operation and its transformation matrix to the standard output. |
| **Parameters** | `none`                                                              |
| **Returns**  | `void`                                                              |
| **Behaviour** | Prints the operation string and the elements of the transformation matrix `T` to the console. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `void Copy(PSymOp symOp)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Copies the state of another `SymOp` object into this `SymOp` object.          |
| **Parameters** | `symOp` — A pointer to another `SymOp` object whose state is to be copied.   |
| **Returns**  | `void`                                                              |
| **Behaviour** | Copies the transformation matrix `T` and the `XYZOp` string from the input `symOp` object. |
| **Side effects** | Modifies the `T` and `XYZOp` member variables. |
| **Errors / exceptions** | `none`                                                              |

#### `void write(io::RFile f)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Writes the current state of the `SymOp` object to an input file.      |
| **Parameters** | `f` — An `io::RFile` object to write the data to.                   |
| **Returns**  | `void`                                                              |
| **Behaviour** | Writes the version number, the number of symmetry operations, and the data for each symmetry operation to the file. |
| **Side effects** | Writes data to the specified file. |
| **Errors / exceptions** | `none`                                                              |

#### `void read(io::RFile f)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Reads the current state of the `SymOp` object from an input file.      |
| **Parameters** | `f` — An `io::RFile` object to read the data from.                   |
| **Returns**  | `void`                                                              |
| **Behaviour** | Reads the version number, the number of symmetry operations, and the data for each symmetry operation from the file. |
| **Side effects** | Reads data from the specified file. |
| **Errors / exceptions** | `none`                                                              |

---

## `SymOps`

> Represents a container of symmetry operators, providing methods for managing and applying these operators to coordinates.

### Methods

#### `SymOps(void)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Constructs a `SymOps` object, initializing its internal state.             |
| **Parameters** | `none`                                                              |
| **Returns**  | `void`                                                              |
| **Behaviour** | Initializes the `SpGroup`, `Nops`, `symOp`, and `Stream` members.          |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `~SymOps(void)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Destroys the `SymOps` object, freeing its allocated memory.              |
| **Parameters** | `none`                                                              |
| **Returns**  | `void`                                                              |
| **Behaviour** | Frees the memory allocated for `SpGroup`, `symOp`, and other internal data. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `virtual void FreeMemory()`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Releases the memory allocated for the internal data of the `SymOps` object. |
| **Parameters** | `none`                                                              |
| **Returns**  | `void`                                                              |
| **Behaviour** | Deallocates memory associated with `SpGroup`, `symOp`, and other resources.|
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `int SetGroupSymopLib(cpstr SpaceGroup, cpstr symop_lib)`

| Field        | Detail                                                                                                |
|--------------|--------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the library containing the symmetry operations.                                                     |
| **Parameters** | `SpaceGroup` — The name of the space group.  `symop_lib` — The path to the symmetry operation library file. |
| **Returns**  | `SYMOP_Ok` if successful, or an error code if an error occurred.                                         |
| **Behaviour** | Retrieves the symmetry operation library from the specified path, initializes the `SpGroup`, `Nops`, and `symOp` members. |
| **Side effects** | Modifies the `SpGroup`, `Nops`, and `symOp` member variables. |
| **Errors / exceptions** | Returns `SYMOP_NoLibFile` if the library file cannot be found, `SYMOP_UnknownSpaceGroup` if the space group is not found. |

#### `int SetGroup(cpstr SpaceGroup, cpstr syminfo_lib)`

| Field        | Detail                                                                                                |
|--------------|--------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the space group and symmetry information.                                                         |
| **Parameters** | `SpaceGroup` — The name of the space group. `syminfo_lib` — The path to the symmetry information library file. |
| **Returns**  | `SYMOP_Ok` if successful, or an error code if an error occurred.                                         |
| **Behaviour** | Retrieves the space group and symmetry information from the specified path, initializes the `SpGroup`, `Nops`, and `symOp` members. |
| **Side effects** | Modifies the `SpGroup`, `Nops`, and `symOp` member variables. |
| **Errors / exceptions** | Returns `SYMOP_NoLibFile` if the library file cannot be found, `SYMOP_UnknownSpaceGroup` if the space group is not found. |

#### `void Reset()`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Removes all symmetry operations from the `SymOps` object.              |
| **Parameters** | `none`                                                              |
| **Returns**  | `void`                                                              |
| **Behaviour** | Resets the `SpGroup`, `Nops`, `symOp` and clears the contents of `symOp`. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `virtual int AddSymOp(cpstr XYZOperation)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Adds a new symmetry operation to the `SymOps` object.                  |
| **Parameters** | `XYZOperation` — A character string representing the symmetry operation. |
| **Returns**  | `int` — The index of the added operation in the `symOp` array, or an error code if the operation is invalid. |
| **Behaviour** | Adds the specified operation to the `symOp` array and increments the count of symmetry operations `Nops`. |
| **Side effects** | Modifies the `Nops` member variable. |
| **Errors / exceptions** | Returns `SYMOP_NoSymOps` if no space is available. |

#### `void PutGroupName(cpstr SpGroupName)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Sets the space group name.                                             |
| **Parameters** | `SpGroupName` — The space group name.                              |
| **Returns**  | `void`                                                              |
| **Behaviour** | Copies the name to `SpGroup` |
| **Side effects** | Modifies the `SpGroup` member variable. |
| **Errors / exceptions** | `none`                                                              |

#### `int GetNofSymOps()`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Returns the number of symmetry operations in the `SymOps` object.     |
| **Parameters** | `none`                                                              |
| **Returns**  | `int` — The number of symmetry operations.                           |
| **Behaviour** | Returns the value of the `Nops` member variable.                     |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `pstr GetSymOp(int Nop)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Returns the character string representation of the symmetry operation at index `Nop`. |
| **Parameters** | `Nop` — The index of the symmetry operation to retrieve.             |
| **Returns**  | `pstr` — A pointer to the string representing the symmetry operation. Returns an empty string if `Nop` is invalid. |
| **Behaviour** | Returns the string stored in the `symOp` array at the specified index. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | Returns an empty string if `Nop` is out of bounds. |

#### `int Transform(realtype & x, realtype & y, realtype & z, int Nop)`

| Field        | Detail                                                                                                |
|--------------|--------------------------------------------------------------------------------------------------------|
| **Summary**  | Transforms the given coordinates `x`, `y`, and `z` according to the symmetry operation `Nop`.           |
| **Parameters** | `x` — The x-coordinate to transform.  `y` — The y-coordinate to transform. `z` — The z-coordinate to transform. `Nop` — The index of the symmetry operation. |
| **Returns**  | `int` — 0 if successful, 1 if `Nop` is invalid, 2 if an error occurs.                                 |
| **Behaviour** | Applies the transformation matrix `T` corresponding to the symmetry operation `Nop` to the given coordinates. |
| **Side effects** | Modifies the values of `x`, `y`, and `z`. |
| **Errors / exceptions** | Returns 1 if `Nop` is out of range. Returns 2 if any other error occurs. |

#### `int GetTMatrix(mat44 & TMatrix, int Nop)`

| Field        | Detail                                                                                                |
|--------------|--------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the coordinate transformation matrix for the symmetry operation `Nop`.                         |
| **Parameters** | `TMatrix` — A `mat44` object where the transformation matrix will be copied. `Nop` — The index of the symmetry operation. |
| **Returns**  | `int` — 0 if successful, 1 if `Nop` is invalid, 2 if an error occurs.                                 |
| **Behaviour** | Copies the transformation matrix `T` corresponding to the symmetry operation `Nop` to the `TMatrix` object. |
| **Side effects** | Modifies the `TMatrix` member variable. |
| **Errors / exceptions** | Returns 1 if `Nop` is out of range. Returns 2 if any other error occurs. |

#### `void Print()`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Prints the symmetry operation and its transformation matrix to the standard output. |
| **Parameters** | `none`                                                              |
| **Returns**  | `void`                                                              |
| **Behaviour** | Prints the operation string and the elements of the transformation matrix `T` to the console. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `void Copy(PSymOps symOps)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Copies the state of another `SymOps` object into this `SymOps` object.          |
| **Parameters** | `symOps` — A pointer to another `SymOps` object whose state is to be copied.   |
| **Returns**  | `void`                                                              |
| **Behaviour** | Copies the transformation matrix `T` and the `SpGroup` string from the input `symOps` object. |
| **Side effects** | Modifies the `T` and `SpGroup` member variables. |
| **Errors / exceptions** | `none`                                                              |

#### `void write(io::RFile f)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Writes the current state of the `SymOps` object to an input file.      |
| **Parameters** | `f` — An `io::RFile` object to write the data to.                   |
| **Returns**  | `void`                                                              |
| **Behaviour** | Writes the version number, the number of symmetry operations, and the data for each symmetry operation to the file. |
| **Side effects** | Writes data to the specified file. |
| **Errors / exceptions** | `none`                                                              |

#### `void read(io::RFile f)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Reads the current state of the `SymOps` object from an input file.      |
| **Parameters** | `f` — An `io::RFile` object to read the data from.                   |
| **Returns**  | `void`                                                              |
| **Behaviour** | Reads the version number, the number of symmetry operations, and the data for each symmetry operation from the file. |
| **Side effects** | Reads data from the specified file. |
| **Errors / exceptions** | `none`                                                              |

---

## `SymOps`

> Represents a container of symmetry operators, providing methods for managing and applying these operators to coordinates.

### Methods

#### `SymOps(void)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Constructs a `SymOps` object, initializing its internal state.             |
| **Parameters** | `none`                                                              |
| **Returns**  | `void`                                                              |
| **Behaviour** | Initializes the `SpGroup`, `Nops`, `symOp`, and `Stream` members.          |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `~SymOps(void)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Destroys the `SymOps` object, freeing its allocated memory.              |
| **Parameters** | `none`                                                              |
| **Returns**  | `void`                                                              |
| **Behaviour** | Frees the memory allocated for `SpGroup`, `symOp`, and other internal data. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `virtual void FreeMemory()`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Releases the memory allocated for the internal data of the `SymOps` object. |
| **Parameters** | `none`                                                              |
| **Returns**  | `void`                                                              |
| **Behaviour** | Deallocates memory associated with `SpGroup`, `symOp`, and other resources.|
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `int SetGroupSymopLib(cpstr SpaceGroup, cpstr symop_lib)`

| Field        | Detail                                                                                                |
|--------------|--------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the library containing the symmetry operations.                                                     |
| **Parameters** | `SpaceGroup` — The name of the space group.  `symop_lib` — The path to the symmetry operation library file. |
| **Returns**  | `SYMOP_Ok` if successful, or an error code if an error occurred.                                         |
| **Behaviour** | Retrieves the symmetry operation library from the specified path, initializes the `SpGroup`, `Nops`, and `symOp` members. |
| **Side effects** | Modifies the `SpGroup`, `Nops`, and `symOp` member variables. |
| **Errors / exceptions** | Returns `SYMOP_NoLibFile` if the library file cannot be found, `SYMOP_UnknownSpaceGroup` if the space group is not found. |

#### `int SetGroup(cpstr SpaceGroup, cpstr syminfo_lib)`

| Field        | Detail                                                                                                |
|--------------|--------------------------------------------------------------------------------------------------------|
| **Summary**  | Sets the space group and symmetry information.                                                         |
| **Parameters** | `SpaceGroup` — The name of the space group. `syminfo_lib` — The path to the symmetry information library file. |
| **Returns**  | `SYMOP_Ok` if successful, or an error code if an error occurred.                                         |
| **Behaviour** | Retrieves the space group and symmetry information from the specified path, initializes the `SpGroup`, `Nops`, and `symOp` members. |
| **Side effects** | Modifies the `SpGroup`, `Nops`, and `symOp` member variables. |
| **Errors / exceptions** | Returns `SYMOP_NoLibFile` if the library file cannot be found, `SYMOP_UnknownSpaceGroup` if the space group is not found. |

#### `void Reset()`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Removes all symmetry operations from the `SymOps` object.              |
| **Parameters** | `none`                                                              |
| **Returns**  | `void`                                                              |
| **Behaviour** | Resets the `SpGroup`, `Nops`, `symOp` and clears the contents of `symOp`. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `virtual int AddSymOp(cpstr XYZOperation)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Adds a new symmetry operation to the `SymOps` object.                  |
| **Parameters** | `XYZOperation` — A character string representing the symmetry operation. |
| **Returns**  | `int` — The index of the added operation in the `symOp` array, or an error code if the operation is invalid. |
| **Behaviour** | Adds the specified operation to the `symOp` array and increments the count of symmetry operations `Nops`. |
| **Side effects** | Modifies the `Nops` member variable. |
| **Errors / exceptions** | Returns `SYMOP_NoSymOps` if no space is available. |

#### `void PutGroupName(cpstr SpGroupName)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Sets the space group name.                                             |
| **Parameters** | `SpGroupName` — The space group name.                              |
| **Returns**  | `void`                                                              |
| **Behaviour** | Copies the name to `SpGroup` |
| **Side effects** | Modifies the `SpGroup` member variable. |
| **Errors / exceptions** | `none`                                                              |

#### `int GetNofSymOps()`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Returns the number of symmetry operations in the `SymOps` object.     |
| **Parameters** | `none`                                                              |
| **Returns**  | `int` — The number of symmetry operations.                           |
| **Behaviour** | Returns the value of the `Nops` member variable.                     |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `pstr GetSymOp(int Nop)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Returns the character string representation of the symmetry operation at index `Nop`. |
| **Parameters** | `Nop` — The index of the symmetry operation to retrieve.             |
| **Returns**  | `pstr` — A pointer to the string representing the symmetry operation. Returns an empty string if `Nop` is invalid. |
| **Behaviour** | Returns the string stored in the `symOp` array at the specified index. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | Returns an empty string if `Nop` is out of bounds. |

#### `int Transform(realtype & x, realtype & y, realtype & z, int Nop)`

| Field        | Detail                                                                                                |
|--------------|--------------------------------------------------------------------------------------------------------|
| **Summary**  | Transforms the given coordinates `x`, `y`, and `z` according to the symmetry operation `Nop`.           |
| **Parameters** | `x` — The x-coordinate to transform.  `y` — The y-coordinate to transform. `z` — The z-coordinate to transform. `Nop` — The index of the symmetry operation. |
| **Returns**  | `int` — 0 if successful, 1 if `Nop` is invalid, 2 if an error occurs.                                 |
| **Behaviour** | Applies the transformation matrix `T` corresponding to the symmetry operation `Nop` to the given coordinates. |
| **Side effects** | Modifies the values of `x`, `y`, and `z`. |
| **Errors / exceptions** | Returns 1 if `Nop` is out of range. Returns 2 if any other error occurs. |

#### `int GetTMatrix(mat44 & TMatrix, int Nop)`

| Field        | Detail                                                                                                |
|--------------|--------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the coordinate transformation matrix for the symmetry operation `Nop`.                         |
| **Parameters** | `TMatrix` — A `mat44` object where the transformation matrix will be copied. `Nop` — The index of the symmetry operation. |
| **Returns**  | `int` — 0 if successful, 1 if `Nop` is invalid, 2 if an error occurs.                                 |
| **Behaviour** | Copies the transformation matrix `T` corresponding to the symmetry operation `Nop` to the `TMatrix` object. |
| **Side effects** | Modifies the `TMatrix` member variable. |
| **Errors / exceptions** | Returns 1 if `Nop` is out of range. Returns 2 if any other error occurs. |

#### `void Print()`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Prints the symmetry operation and its transformation matrix to the standard output. |
| **Parameters** | `none`                                                              |
| **Returns**  | `void`                                                              |
| **Behaviour** | Prints the operation string and the elements of the transformation matrix `T` to the console. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `void Copy(PSymOps SymOps)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Copies the state of another `SymOps` object into this `SymOps` object.          |
| **Parameters** | `symOps` — A pointer to another `SymOps` object whose state is to be copied.   |
| **Returns**  | `void`                                                              |
| **Behaviour** | Copies the transformation matrix `T` and the `SpGroup` string from the input `symOps` object. |
| **Side effects** | Modifies the `T` and `SpGroup` member variables. |
| **Errors / exceptions** | `none`                                                              |

#### `void write(io::RFile f)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Writes the current state of the `SymOps` object to an input file.      |
| **Parameters** | `f` — An `io::RFile` object to write the data to.                   |
| **Returns**  | `void`                                                              |
| **Behaviour** | Writes the version number, the number of symmetry operations, and the data for each symmetry operation to the file. |
| **Side effects** | Writes data to the specified file. |
| **Errors / exceptions** | `none`                                                              |

#### `void read(io::RFile f)`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Reads the current state of the `SymOps` object from an input file.      |
| **Parameters** | `f` — An `io::RFile` object to read the data from.                   |
| **Returns**  | `void`                                                              |
| **Behaviour** | Reads the version number, the number of symmetry operations, and the data for each symmetry operation from the file. |
| **Side effects** | Reads data from the specified file. |
| **Errors / exceptions** | `none`                                                              |

---

## Free Functions

No free functions are defined in this code.


---

### `mmdb_tables`

## `MMDB_Tables`

> Provides a set of functions for accessing and manipulating molecular data, primarily for protein crystallography applications.

### Methods

#### `bool isMetal(cpstr element)`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Checks if the given `element` string represents a metal atom.                                                                     |
| **Parameters**| `element` — A C-style string (cpstr) representing the element name.                                                                 |
| **Returns**  | `true` if the element is a metal; `false` otherwise.                                                                                 |
| **Behaviour** | The function iterates through the `ElementMetal` array, comparing the input `element` with each metal name. The comparison is case-insensitive. |
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| None                                                                                                                                |

#### `int getElementNo(cpstr element)`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the numerical index of the element in the `ElementMetal` array.                                                               |
| **Parameters**| `element` — A C-style string (cpstr) representing the element name.                                                                 |
| **Returns**  | The index (starting from 1) of the element in the `ElementMetal` array. Returns `ELEMENT_UNKNOWN` if the element is not found.          |
| **Behaviour** | The function iterates through the `ElementMetal` array, comparing the input `element` with each metal name. The comparison is case-insensitive. |
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| None                                                                                                                                |

#### `realtype getMolecWeight(cpstr element)`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the molecular weight of the element based on the `MolecWeight` array.                                                              |
| **Parameters**| `element` — A C-style string (cpstr) representing the element name.                                                                 |
| **Returns**  | The molecular weight of the element.                                                                                                |
| **Behaviour** | The function iterates through the `MolecWeight` array, comparing the input `element` with each element name. The comparison is case-insensitive. |
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| Returns `1.0` if the element is not found.                                                                                        |

#### `realtype getCovalentRadius(cpstr element)`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the covalent radius of the element based on the `CovalentRadius` array.                                                             |
| **Parameters**| `element` — A C-style string (cpstr) representing the element name.                                                                 |
| **Returns**  | The covalent radius of the element.                                                                                                |
| **Behaviour** | The function iterates through the `CovalentRadius` array, comparing the input `element` with each element name. The comparison is case-insensitive. |
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| Returns `2.2*CovalentRadius[0]` if the element is not found.                                                                      |

#### `realtype getVdWaalsRadius(cpstr element)`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the Van der Waals radius of the element based on the `VdWaalsRadius` array.                                                        |
| **Parameters**| `element` — A C-style string (cpstr) representing the element name.                                                                 |
| **Returns**  | The Van der Waals radius of the element.                                                                                             |
| **Behaviour** | The function iterates through the `VdWaalsRadius` array, comparing the input `element` with each element name. The comparison is case-insensitive. |
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| Returns `1.8` if the element is not found.                                                                                         |

#### `bool isMetal(cpstr element)` (Duplicate)

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Checks if the given `element` string represents a metal atom.                                                                     |
| **Parameters**| `element` — A C-style string (cpstr) representing the element name.                                                                 |
| **Returns**  | `true` if the element is a metal; `false` otherwise.                                                                                 |
| **Behaviour** | The function iterates through the `ElementMetal` array, comparing the input `element` with each metal name. The comparison is case-insensitive. |
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| None                                                                                                                                |

#### `int getResidueNo(cpstr resName)`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the numerical index of the residue in the `ResidueName` array.                                                              |
| **Parameters**| `resName` — A C-style string (cpstr) representing the residue name.                                                                 |
| **Returns**  | The index (starting from 0) of the residue in the `ResidueName` array. Returns `-1` if the residue is not found.                     |
| **Behaviour** | The function iterates through the `ResidueName` array, comparing the input `resName` with each residue name. The comparison is case-insensitive. |
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| None                                                                                                                                |

#### `char const* ResidueName[nResNames]`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Array of residue names.                                                                                                             |
| **Parameters**| None                                                                                                                                |
| **Returns**  |  `ResidueName` array.                                                                                                            |
| **Behaviour** | A simple array of strings.                                                                                                        |
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| None                                                                                                                                |

#### `char const ResidueName1[nResNames]`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Array of single-letter residue names.                                                                                             |
| **Parameters**| None                                                                                                                                |
| **Returns**  | `ResidueName1` array.                                                                                                            |
| **Behaviour** | A simple array of characters.                                                                                                        |
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| None                                                                                                                                |

#### `realtype NAvogadro`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Avogadro's number.                                                                                                               |
| **Parameters**| None                                                                                                                                |
| **Returns**  |  `Navogadro` constant.                                                                                                          |
| **Behaviour** | A constant value.                                                                                                              |
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| None                                                                                                                                |

#### `AAProperty const AAProperties[nAminoacidNames]`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Array of AAProperty structures.                                                                                                     |
| **Parameters**| None                                                                                                                                |
| **Returns**  | `AAProperties` array.                                                                                                            |
| **Behaviour** | A simple array of `AAProperty` structures.                                                                                            |
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| None                                                                                                                                |

#### `realtype const AAHydropathyScale[nAminoacidNames]`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Hydropathy scale for amino acids.                                                                                                   |
| **Parameters**| None                                                                                                                                |
| **Returns**  |  `AAHydropathyScale` constant.                                                                                                   |
| **Behaviour** | A constant value.                                                                                                              |
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| None                                                                                                                                |

#### `realtype const AACharge[nAminoacidNames]`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Charge of amino acids.                                                                                                             |
| **Parameters**| None                                                                                                                                |
| **Returns**  |  `AACharge` constant.                                                                                                          |
| **Behaviour** | A constant value.                                                                                                              |
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| None                                                                                                                                |

#### `realtype const AACharge[nAminoacidNames]` (Duplicate)

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Charge of amino acids.                                                                                                             |
| **Parameters**| None                                                                                                                                |
| **Returns**  |  `AACharge` constant.                                                                                                          |
| **Behaviour** | A constant value.                                                                                                              |
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| None                                                                                                                                |

#### `realtype GetAAHydropathy(cpstr resName)`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the hydropathy value of the given amino acid.                                                                                |
| **Parameters**| `resName` — A C-style string (cpstr) representing the amino acid name.                                                               |
| **Returns**  | The hydropathy value of the amino acid.                                                                                             |
| **Behaviour** | The function takes an amino acid name, converts it to uppercase, and then uses it to look up the corresponding hydropathy value in the `AAProperties` array.|
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| Returns `-MaxReal` if the amino acid is not found.                                                                                |

#### `realtype GetAACharge(cpstr resName)`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the charge of the given amino acid.                                                                                         |
| **Parameters**| `resName` — A C-style string (cpstr) representing the amino acid name.                                                               |
| **Returns**  | The charge of the amino acid.                                                                                                     |
| **Behaviour** | The function takes an amino acid name, converts it to uppercase, and then uses it to look up the corresponding charge value in the `AAProperties` array.|
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| Returns `0.0` if the amino acid is not found.                                                                                     |

#### `realtype GetAASolvationEnergy(cpstr resName)`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the solvation energy of the given amino acid.                                                                                  |
| **Parameters**| `resName` — A C-style string (cpstr) representing the amino acid name.                                                               |
| **Returns**  | The solvation energy of the amino acid.                                                                                             |
| **Behaviour** | The function takes an amino acid name, converts it to uppercase, and then uses it to look up the corresponding solvation energy value in the `AAProperties` array.|
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| Returns `0.0` if the amino acid is not found.                                                                                     |

#### `int GetAAPIndex(cpstr resName)`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the index of the given amino acid in the `AAProperties` array.                                                                 |
| **Parameters**| `resName` — A C-style string (cpstr) representing the amino acid name.                                                               |
| **Returns**  | The index (starting from 0) of the amino acid in the `AAProperties` array. Returns `-1` if the amino acid is not found.                |
| **Behaviour** | The function takes an amino acid name, converts it to uppercase, and then uses it to look up the index in the `AAProperties` array.|
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| None                                                                                                                                |

#### `realtype GetAAHydropathy(cpstr resName)` (Duplicate)

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the hydropathy value of the given amino acid.                                                                                |
| **Parameters**| `resName` — A C-style string (cpstr) representing the amino acid name.                                                               |
| **Returns**  | The hydropathy value of the amino acid.                                                                                             |
| **Behaviour** | The function takes an amino acid name, converts it to uppercase, and then uses it to look up the corresponding hydropathy value in the `AAProperties` array.|
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| Returns `-MaxReal` if the amino acid is not found.                                                                                |

#### `realtype GetAACharge(cpstr resName)` (Duplicate)

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the charge of the given amino acid.                                                                                         |
| **Parameters**| `resName` — A C-style string (cpstr) representing the amino acid name.                                                               |
| **Returns**  | The charge of the amino acid.                                                                                                     |
| **Behaviour** | The function takes an amino acid name, converts it to uppercase, and then uses it to look up the corresponding charge value in the `AAProperties` array.|
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| Returns `0.0` if the amino acid is not found.                                                                                     |

#### `realtype GetAASolvationEnergy(cpstr resName)` (Duplicate)

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the solvation energy of the given amino acid.                                                                                  |
| **Parameters**| `resName` — A C-style string (cpstr) representing the amino acid name.                                                               |
| **Returns**  | The solvation energy of the amino acid.                                                                                             |
| **Behaviour** | The function takes an amino acid name, converts it to uppercase, and then uses it to look up the corresponding solvation energy value in the `AAProperties` array.|
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| Returns `0.0` if the amino acid is not found.                                                                                     |

#### `int GetAASimilarity(cpstr resName1, cpstr resName2)`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Returns the similarity index between two amino acids based on the `AASimilarity` array.                                                  |
| **Parameters**| `resName1` — A C-style string (cpstr) representing the first amino acid name.                                                       |
|              | `resName2` — A C-style string (cpstr) representing the second amino acid name.                                                      |
| **Returns**  | The similarity index (ranging from 0 to 5) between the two amino acids.                                                             |
| **Behaviour** | The function takes two amino acid names, converts them to uppercase, and then uses them to look up the similarity index in the `AASimilarity` array.|
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| Returns -1 if either amino acid is not found, -2 if both amino acids are not found.                                                  |

#### `cpstr const StdSolventName[nSolventNames]`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Array of solvent names.                                                                                                             |
| **Parameters**| None                                                                                                                                |
| **Returns**  | `StdSolventName` array.                                                                                                            |
| **Behaviour** | A simple array of strings.                                                                                                        |
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| None                                                                                                                                |

#### `cpstr const NucleotideName[nNucleotideNames]`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Array of nucleotide names.                                                                                                         |
| **Parameters**| None                                                                                                                                |
| **Returns**  | `NucleotideName` array.                                                                                                           |
| **Behaviour** | A simple array of strings.                                                                                                        |
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| None                                                                                                                                |

#### `bool isSolvent(cpstr resName)`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Checks if the given `resName` string represents a solvent.                                                                            |
| **Parameters**| `resName` — A C-style string (cpstr) representing the residue name.                                                                 |
| **Returns**  | `true` if the residue is a solvent; `false` otherwise.                                                                                 |
| **Behaviour** | The function takes a residue name, converts it to uppercase, and then uses it to look up the name in the `StdSolventName` array.        |
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| None                                                                                                                                |

#### `bool isAminoacid(cpstr resName)`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Checks if the given `resName` string represents an amino acid.                                                                        |
| **Parameters**| `resName` — A C-style string (cpstr) representing the residue name.                                                                 |
| **Returns**  | `true` if the residue is an amino acid; `false` otherwise.                                                                           |
| **Behaviour** | The function takes a residue name, converts it to uppercase, and then uses it to look up the name in the `AAProperties` array.        |
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| None                                                                                                                                |

#### `bool isNucleotide(cpstr resName)`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Checks if the given `resName` string represents a nucleotide.                                                                         |
| **Parameters**| `resName` — A C-style string (cpstr) representing the residue name.                                                                 |
| **Returns**  | `true` if the residue is a nucleotide; `false` otherwise.                                                                            |
| **Behaviour** | The function takes a residue name, converts it to uppercase, and then uses it to look up the name in the `NucleotideName` array.      |
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| None                                                                                                                                |

#### `int isDNARNA(cpstr resName)`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Determines if the given residue is DNA or RNA.                                                                                        |
| **Parameters**| `resName` — A C-style string (cpstr) representing the residue name.                                                                 |
| **Returns**  | `0` if neither DNA nor RNA; `1` if DNA; `2` if RNA.                                                                                 |
| **Behaviour** | The function takes a residue name, converts it to uppercase, and then uses it to look up the residue type in the `NucleotideName` array. |
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| None                                                                                                                                |

#### `bool isSugar(cpstr resName)`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Checks if the given `resName` string represents a sugar.                                                                               |
| **Parameters**| `resName` — A C-style string (cpstr) representing the residue name.                                                                 |
| **Returns**  | `true` if the residue is a sugar; `false` otherwise.                                                                                 |
| **Behaviour** | This function always returns `false`, indicating that the residue is not a sugar.                                                  |
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| None                                                                                                                                |

#### `cpstr const Res1Code[]`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Array of one-letter residue codes.                                                                                                   |
| **Parameters**| None                                                                                                                                |
| **Returns**  | `Res1Code` array.                                                                                                                 |
| **Behaviour** | A simple array of strings.                                                                                                        |
| **Side effects**| None                                                                                                                                |
| **Errors / exceptions**| None                                                                                                                                |

#### `void Get1LetterCode(cpstr res3name, pstr res1code)`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Converts a three-letter residue name to a one-letter residue name.                                                                     |
| **Parameters**| `res3name` — A C-style string (cpstr) representing the three-letter residue name.                                                       |
|              | `res1code` — A pointer to a C-style string (pstr) where the one-letter residue name will be stored.                                  |
| **Returns**  | None                                                                                                                                |
| **Behaviour** | The function takes a three-letter residue name, converts it to uppercase, and then writes the corresponding one-letter name to the `res1code` string. |
| **Side effects**| Modifies the `res1code` string.                                                                                                     |
| **Errors / exceptions**| None                                                                                                                                |

#### `void Get1LetterCode(cpstr res3name, char & res1code)`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Converts a three-letter residue name to a one-letter residue name.                                                                     |
| **Parameters**| `res3name` — A C-style string (cpstr) representing the three-letter residue name.                                                       |
|              | `res1code` — A character & reference to a character variable where the one-letter residue name will be stored.                       |
| **Returns**  | None                                                                                                                                |
| **Behaviour** | The function takes a three-letter residue name, converts it to uppercase, and then writes the corresponding one-letter name to the `res1code` variable. |
| **Side effects**| Modifies the `res1code` variable.                                                                                                     |
| **Errors / exceptions**| None                                                                                                                                |

#### `void Get3LetterCode(cpstr res1name, pstr res3code)`

| Field        | Detail                                                                                                                               |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Converts a one-letter residue name to a three-letter residue name.                                                                     |
| **Parameters**| `res1name` — A C-style string (cpstr) representing the one-letter residue name.                                                        |
|              | `res3code` — A pointer to a C-style string (pstr) where the three-letter residue name will be stored.                                |
| **Returns**  | None                                                                                                                                |
| **Behaviour** | The function takes a one-letter residue name, converts it to uppercase, and then writes the corresponding three-letter name to the `res3code` string. |
| **Side effects**| Modifies the `res3code` string.                                                                                                     |
| **Errors / exceptions**| None                                                                                                                                |

## Free Functions


---

### `mmdb_title`



---

### `mmdb_uddata`

```cpp
//  $Id: mmdb_uddata.h $
//  =================================================================
//
//   CCP4 Coordinate Library: support of coordinate-related
//   functionality in protein crystallography applications.
//
//   Copyright (C) Eugene Krissinel 2000-2013.
//
//    This library is free software: you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License version 3, modified in accordance with the provisions
//    of the license to address the requirements of UK law.
//
//    You should have received a copy of the modified GNU Lesser
//    General Public License along with this library. If not, copies
//    may be downloaded from http://www.ccp4.ac.uk/ccp4license.php
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//  =================================================================
//
//    12.09.13   <--  Date of Last Modification.
//                   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  -----------------------------------------------------------------
//
//  **** Module  :   MMDBF_UDData <interface>
//       ~~~~~~~~~
//  **** Project :   MacroMolecular Data Base (MMDB)
//       ~~~~~~~~~
//
//  **** Classes :   mmdb::UDData ( user-defined data )
//       ~~~~~~~~~
//
//   (C) E. Krissinel 2000-2013
//
//  =================================================================
//

#ifndef __MMDB_UDData__
#define __MMDB_UDData__

#include "mmdb_mask.h"

namespace mmdb  {

  //  =======================  UDRegister  =========================

  enum UDR_TYPE  {
    UDR_ATOM      = 0,
    UDR_RESIDUE   = 1,
    UDR_CHAIN     = 2,
    UDR_MODEL     = 3,
    UDR_HIERARCHY = 4
  };

  enum UDD_FLAG  {
    UDRF_ATOM      = 0x01000000,
    UDRF_RESIDUE   = 0x02000000,
    UDRF_CHAIN     = 0x04000000,
    UDRF_MODEL     = 0x08000000,
    UDRF_HIERARCHY = 0x10000000,
    UDRF_MASK      = 0x00FFFFFF
  };

  DefineClass(UDRegister);
  DefineStreamFunctions(UDRegister);

  class MMDB_DL_EXPORT UDRegister : public io::Stream  {

    public :

      UDRegister ();
      UDRegister ( io::RPStream Object );
      ~UDRegister();

      int RegisterUDInteger ( UDR_TYPE udr_type, cpstr UDDataID );
      int RegisterUDReal    ( UDR_TYPE udr_type, cpstr UDDataID );
      int RegisterUDString  ( UDR_TYPE udr_type, cpstr UDDataID );
      int GetUDDHandle      ( UDR_TYPE udr_type, cpstr UDDataID );

      void write ( io::RFile f );
      void read  ( io::RFile f );

    protected :
      int      nIUDR[5],nRUDR[5],nSUDR[5];
      psvector IUDRegister[5];
      psvector RUDRegister[5];
      psvector SUDRegister[5];

      void  InitUDRegister ();
      void  FreeUDRegister ();
      int   RegisterUDData ( psvector & UDRegister,
                             int      & nUDR,
                             cpstr      UDDataID );

  };


  //  ==========================  UDData  ===========================

  enum UDDATA_RC  {
    UDDATA_Ok           =  0,
    UDDATA_WrongHandle  = -1,
    UDDATA_WrongUDRType = -2,
    UDDATA_NoData       = -3
  };

  DefineClass(UDData);
  DefineStreamFunctions(UDData);

  class MMDB_DL_EXPORT UDData : public Mask  {

    friend class SelManager;

    public :

      UDData ();
      UDData ( io::RPStream Object );
      ~UDData();

    protected :
      ivector  IUData;
      rvector  RUData;
      psvector SUData;

      void  InitUDData   ();
      void  FreeUDDMemory();
      int   getNofIUData ();
      int   getNofRUData ();
      int   getNofSUData ();
      void  setNofSUData ( int newN );

      int   putUDData ( int UDDhandle, int      iudd );
      int   putUDData ( int UDDhandle, realtype rudd );
      int   putUDData ( int UDDhandle, cpstr    sudd );

      int   getUDData ( int UDDhandle, int      & iudd );
      int   getUDData ( int UDDhandle, realtype & rudd );
      int   getUDData ( int UDDhandle, pstr sudd, int maxLen );
      pstr  getUDData ( int UDDhandle, int * retcode=NULL );
      int   getUDData ( int UDDhandle, pstr     & sudd );

      void  write ( io::RFile f );
      void  read  ( io::RFile f );

  };

}  // namespace mmdb

#endif
```

---

### `mmdb_utils`

Okay, I've analyzed the provided C++ source code and generated comprehensive method-level documentation for each class and function.  The documentation includes summaries, parameter details, return values, detailed behavior descriptions, and side effects/errors, exactly as requested.

**Important Notes and Assumptions:**

*   **Implementation Focus:** The documentation prioritizes describing *what* the code actually does, not just its declaration.
*   **Step-by-Step Logic:** The `Behaviour` sections provide a detailed breakdown of the code's execution flow, essential for understanding complex algorithms and data manipulations.
*   **Error Handling:** The documentation identifies potential errors, exceptions, or failure conditions, along with how the code handles them.
*   **Side Effects:** The `Side effects` sections clearly outline any changes made to the state of the program (e.g., file writes, signal emissions, memory modifications) that might impact other parts of the code.
*   **Assumed Data Structures:**  Some comments include assumptions about internal data structures (e.g., `mmcif::PStruct`, `mmcif::PData`, `mmcif::PLoop`) based on the context of the code. I’ve tried to be as explicit as possible.
*   **`DefineClass` and `DefineStreamFunctions` Macros:**  The documentation reflects the use of macros to define classes and their associated stream functions.
*   **`CreateCopy` and `CreateWrite` Methods**: The documented behaviour of these functions is based on their use in the code.
*   **`PContainerClass`**: The code frequently uses `PContainerClass` and `PContString` in the data structures.

**Here's the complete documentation, formatted as requested:**

---

## `mmdb::ContainerClass`

> This class represents a container for storing other container objects, primarily used for managing data structures in crystallography coordinate files.

### Methods

#### `bool ContainerClass::Append(PContainerClass CC)`

| Field        | Detail                                                                                                                                                                                                                               |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Appends a new container object (`CC`) to the container, maintaining a list of contained objects.                                                                                                                                       |
| **Parameters** | `CC` — `PContainerClass*` — A pointer to the container object to be added.  If `CC` is a null pointer, nothing is added. |
| **Returns**  | `true` if the object was successfully appended; `false` if the appending failed (e.g., due to a duplicate object, or a structural mismatch). |
| **Behaviour** | 1. It checks if `CC` is a valid pointer. 2. It iterates through existing container objects in the `Container` array.  3. It checks if the current container object has the same `ClassID` as the new container. 4. If the same `ClassID` is found, it adds a new entry to the `Container` array.  5. If a new `ClassID` is found, it appends to the end of the `Container` array. |
| **Side effects** | Adds a new element to the `Container` array (if `ClassID` is different), or appends the container to the end.                                                               |
| **Errors / exceptions** |  `false` is returned if `CC` is invalid, or if `CC` has a `ClassID` that does not match existing container objects. |

#### `void ContainerClass::FreeContainer()`

| Field        | Detail |
|--------------|--------|
| **Summary**  | Deallocates the memory associated with the container's internal storage (`Container`) and resets the container size. |
| **Parameters** | none |
| **Returns**  | `void` |
| **Behaviour** | 1. It checks if `Container` is not `NULL`. 2. It deletes the contents of `Container` using `delete[]`. 3. It sets `Container` to `NULL`. 4. It resets the length of the container to 0. |
| **Side effects** | Releases the memory allocated for the `Container` array. |
| **Errors / exceptions** | `None` |

#### `void ContainerClass::AddData(PContainerClass Data)`

| Field        | Detail                                                                                                                                                                                 |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Adds a new container object (`Data`) to the container, maintaining a list of contained objects.                                                                                             |
| **Parameters** | `Data` — `PContainerClass*` — A pointer to the container object to be added. |
| **Returns**  | `void` |
| **Behaviour** | 1. Check if the pointer is not NULL. 2. Create a new array of `PContainerClass` objects, one size larger than the current `Container` array. 3. Copy the contents of the existing `Container` array to the new array. 4. Copy the pointer to the new `PContainerClass` object to the end of the new array. 5. Assign the new array to the `Container` member variable. 6. Increment `length` by 1. |
| **Side effects** | Adds a new element to the `Container` array (if `ClassID` is different). |
| **Errors / exceptions** | `None` |

#### `void ContainerClass::PDBASCIIDump(io::RFile f)`

| Field        | Detail |
|--------------|--------|
| **Summary**  | Writes the content of the `Line` member to the output stream `f`. |
| **Parameters** | `f` — `io::RFile*` — A pointer to the output stream file. |
| **Returns**  | `void` |
| **Behaviour** | Writes the string value stored in the `Line` member to the output file stream `f`. |
| **Side effects** | Writes the string to the output stream `f`. |
| **Errors / exceptions** | `None` |

#### `bool ContainerClass::PDBASCIIDump1(io::RFile f)`

| Field        | Detail |
|--------------|--------|
| **Summary**  | Writes the content of the `Line` member to the output stream `f`. |
| **Parameters** | `f` — `io::RFile*` — A pointer to the output stream file. |
| **Returns**  | `true` or `false` depending on success. |
| **Behaviour** | Writes the string value stored in the `Line` member to the output file stream `f`. |
| **Side effects** | Writes the string to the output stream `f`. |
| **Errors / exceptions** | `None` |

#### `void ContainerClass::MakeCIF(mmcif::PData CIF, int N)`

| Field        | Detail                                                                                                                                                                                 |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Creates a CIF file structure based on the `Line` member of this object.                                                                                                                 |
| **Parameters** | `CIF` — `mmcif::PData*` — A pointer to the CIF data structure to be populated. `N` - number of the element |
| **Returns**  | `void` |
| **Behaviour** | Copies the content of the `Line` member into the `CIF` structure. |
| **Side effects** | Writes data into the CIF structure. |
| **Errors / exceptions** | `None` |

#### `ERROR_CODE ContainerClass::GetCIF(mmcif::PData CIF, int ClassID)`

| Field        | Detail                                                                                                                                                                                                                               |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Retrieves data from the CIF structure and adds it to the `ContainerClass` object.                                                                                                                                                  |
| **Parameters** | `CIF` — `mmcif::PData*` — A pointer to the CIF data structure. `ClassID` — `int` — The ID of the container class. |
| **Returns**  | `Error_NoError` if the operation was successful, `Error_NoData` if no data was found, or another error code if an error occurred during the process. |
| **Behaviour** | 1. Sets the `n` variable to -1. 2. Iterates through the CIF structure, calling `GetInteger` for each tag. 3. If an error occurs during `GetInteger`, it returns the error code and terminates the iteration. 4. If no error occurs, it adds the data to the `ContainerClass` object and returns `Error_NoError`. |
| **Side effects** | Adds data to the `ContainerClass` object if the operation is successful. |
| **Errors / exceptions** | Returns an appropriate error code if any issues occur. |

#### `void ContainerClass::MakeCIF(mmcif::PData CIF)`

| Field        | Detail |
|--------------|--------|
| **Summary**  | Creates a CIF structure based on the `Line` member of this object. |
| **Parameters** | `CIF` — `mmcif::PData*` — A pointer to the CIF data structure to be populated. |
| **Returns**  | `void` |
| **Behaviour** | Copies the content of the `Line` member into the `CIF` structure. |
| **Side effects** | Writes data into the CIF structure. |
| **Errors / exceptions** | `None` |

### Struct `mmdb::ContString`

> This class represents a container for storing a string.

#### `void ContString::ConvertPDBASCII(cpstr S)`

| Field        | Detail                                                                                                 |
|--------------|---------------------------------------------------------------------------------------------------------|
| **Summary**  | Converts a PDB ASCII string (`S`) into the internal string representation stored within the `Line` member. |
| **Parameters** | `S` — `cpstr` — The input PDB ASCII string. |
| **Returns**  | `void` |
| **Behaviour** | 1. Copies the string `S` into the `Line` member. |
| **Side effects** | Modifies the `Line` member to store the input string. |
| **Errors / exceptions** | `None` |

#### `void ContString::PDBASCIIDump(pstr S, int)`

| Field        | Detail |
|--------------|--------|
| **Summary**  | Writes the content of the `Line` member to the output stream `S`. |
| **Parameters** | `S` — `pstr` — A pointer to the output string. `N` - size of string |
| **Returns**  | `void` |
| **Behaviour** | Writes the string value stored in the `Line` member to the output string `S`. |
| **Side effects** | Writes the string to the output string `S`. |
| **Errors / exceptions** | `None` |

#### `bool ContString::PDBASCIIDump1(io::RFile f)`

| Field        | Detail |
|--------------|--------|
| **Summary**  | Writes the content of the `Line` member to the output stream `f`. |
| **Parameters** | `f` — `io::RFile*` — A pointer to the output stream file. |
| **Returns**  | `true` or `false` depending on success. |
| **Behaviour** | Writes the string value stored in the `Line` member to the output file stream `f`. |
| **Side effects** | Writes the string to the output stream `f`. |
| **Errors / exceptions** | `None` |

#### `void ContString::MakeCIF(mmcif::PData CIF, int N)`

| Field        | Detail                                                                                                                                                                                 |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Creates a CIF file structure based on the `Line` member of this object.                                                                                                                 |
| **Parameters** | `CIF` — `mmcif::PData*` — A pointer to the CIF data structure to be populated. `N` - number of the element |
| **Returns**  | `void` |
| **Behaviour** | Copies the content of the `Line` member into the `CIF` structure. |
| **Side effects** | Writes data into the CIF structure. |
| **Errors / exceptions** | `None` |

#### `bool ContString::Append(PContainerClass CC)`

| Field        | Detail                                                                                                                                                                                                                               |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Appends a new container object (`CC`) to the container, maintaining a list of contained objects.                                                                                                                                       |
| **Parameters** | `CC` — `PContainerClass*` — A pointer to the container object to be added. |
| **Returns**  | `true` if the object was successfully appended; `false` if the appending failed (e.g., due to a duplicate object, or a structural mismatch). |
| **Behaviour** | 1. Check if `CC` is a valid pointer. 2. It iterates through existing container objects in the `Container` array.  3. It checks if the current container object has the same `ClassID` as the new container. 4. If the same `ClassID` is found, it adds a new entry to the `Container` array.  5. If a new `ClassID` is found, it appends to the end of the `Container` array. |
| **Side effects** | Adds a new element to the `Container` array (if `ClassID` is different). |
| **Errors / exceptions** |  `false` is returned if `CC` is invalid, or if `CC` has a `ClassID` that does not match existing container objects. |

### Struct `mmdb::ClassContainer`

> This class represents a container for storing container objects, primarily used for managing data structures in crystallography coordinate files.

#### `void ClassContainer::FreeContainer()`

| Field        | Detail |
|--------------|--------|
| **Summary**  | Deallocates the memory associated with the container's internal storage (`Container`) and resets the container size. |
| **Parameters** | none |
| **Returns**  | `void` |
| **Behaviour** | 1. It checks if `Container` is not `NULL`. 2. It deletes the contents of `Container` using `delete[]`. 3. It sets `Container` to `NULL`. 4. It resets the length of the container to 0. |
| **Side effects** | Releases the memory allocated for the `Container` array. |
| **Errors / exceptions** | `None` |

#### `void ClassContainer::AddData(PContainerClass Data)`

| Field        | Detail                                                                                                                                                                                                                               |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Adds a new container object (`Data`) to the container, maintaining a list of contained objects.                                                                                                                                       |
| **Parameters** | `Data` — `PContainerClass*` — A pointer to the container object to be added. |
| **Returns**  | `void` |
| **Behaviour** | 1. Check if the pointer is not NULL. 2. Create a new array of `PContainerClass` objects, one size larger than the current `Container` array. 3. Copy the contents of the existing `Container` array to the new array. 4. Copy the pointer to the new `PContainerClass` object to the end of the new array. 5. Assign the new array to the `Container` member variable. 6. Increment `length` by 1. |
| **Side effects** | Adds a new element to the `Container` array (if `ClassID` is different). |
| **Errors / exceptions** | `None` |

#### `void ClassContainer::PDBASCIIDump(io::RFile f)`

| Field        | Detail |
|--------------|--------|
| **Summary**  | Writes the content of the `Line` member to the output stream `f`. |
| **Parameters** | `f` — `io::RFile*` — A pointer to the output stream file. |
| **Returns**  | `void` |
| **Behaviour** | Writes the string value stored in the `Line` member to the output file stream `f`. |
| **Side effects** | Writes the string to the output stream `f`. |
| **Errors / exceptions** | `None` |

#### `bool ClassContainer::GetCIF(mmcif::PData CIF, int ClassID)`

| Field        | Detail                                                                                                                                                                                                                               |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Retrieves data from the CIF structure and adds it to the `ContainerClass` object.                                                                                                                                                   |
| **Parameters** | `CIF` — `mmcif::PData*` — A pointer to the CIF data structure. `ClassID` — `int` — The ID of the container class. |
| **Returns**  | `Error_NoError` if the operation was successful, `Error_NoData` if no data was found, or another error code if an error occurred during the process. |
| **Behaviour** | 1. Sets the `n` variable to -1. 2. It iterates through the CIF structure, calling `GetInteger` for each tag. 3. If an error occurs during `GetInteger`, it returns the error code and terminates the iteration. 4. If no error occurs, it adds the data to the `ContainerClass` object and returns `Error_NoError`. |
| **Side effects** | Adds data to the `ContainerClass` object if the operation is successful. |
| **Errors / exceptions** | Returns an appropriate error code if any issues occur. |

#### `void ClassContainer::MakeCIF(mmcif::PData CIF)`

| Field        | Detail |
|--------------|--------|
| **Summary**  | Creates a CIF file structure based on the `Line` member of this object. |
| **Parameters** | `CIF` — `mmcif::PData*` — A pointer to the CIF data structure to be populated. |
| **Returns**  | `void` |
| **Behaviour** | Copies the content of the `Line` member into the `CIF` structure. |
| **Side effects** | Writes data into the CIF structure. |
| **Errors / exceptions** | `None` |

#### `void ClassContainer::write (io::RFile f)`

| Field        | Detail |
|--------------|--------|
| **Summary**  | Writes the content of the `Line` member to the output stream `f`. |
| **Parameters** | `f` — `io::RFile*` — A pointer to the output stream file. |
| **Returns**  | `void` |
| **Behaviour** | Writes the string value stored in the `Line` member to the output file stream `f`. |
| **Side effects** | Writes the string to the output stream `f`. |
| **Errors / exceptions** | `None` |

#### `PContainerClass ClassContainer::MakeContainerClass (int ClassID)`

| Field        | Detail |
|--------------|--------|
| **Summary**  | Creates a `ContString` object if the `ClassID` is 1, otherwise returns a `ContainerClass` |
| **Parameters** | `ClassID` — `int` — The ID of the container class. |
| **Returns**  | `PContainerClass` — A pointer to the `ContString` or `ContainerClass` object. |
| **Behaviour** | Creates a `ContString` object if `ClassID` is 1, otherwise returns a `ContainerClass` object. |
| **Side effects** | `None` |
| **Errors / exceptions** | `None` |

#### `void ClassContainer::Copy (PClassContainer CContainer)`

| Field        | Detail |
|--------------|--------|
| **Summary**  | Copies the data from `CContainer` into this object. |
| **Parameters** | `CContainer` — `PClassContainer*` — A pointer to the `ClassContainer` object to copy from. |
| **Returns**  | `void` |
| **Behaviour** | 1. It frees the contents of the `Container` array. 2. It allocates a new `PContainerClass` array, one size larger than the existing `Container` array. 3. It copies the contents of the `CContainer` array into the new array. 4. It assigns the new array to the `Container` member variable. 5. It increments `length` by 1. |
| **Side effects** | Releases the memory allocated for the `Container` array, and copies the data from the other container. |
| **Errors / exceptions** | `None` |

#### `void ClassContainer::read (io::RFile f)`

| Field        | Detail |
|--------------|--------|
| **Summary**  | Reads the content of the `Line` member from the input stream `f`. |
| **Parameters** | `f` — `io::RFile*` — A pointer to the input stream file. |
| **Returns**  | `void` |
| **Behaviour** | Reads the string value stored in the `Line` member from the input file stream `f`. |
| **Side effects** | Reads the string from the input stream `f`. |
| **Errors / exceptions** | `None` |

### Struct `mmdb::AtomPath`

> This class represents an atom path identifier.

#### `int AtomPath::SetPath(cpstr ID)`

| Field        | Detail                                                                                                                                                                                                                               |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Parses an atom path string (`ID`) and populates the corresponding member variables (model, chain, sequence number, insertion code, residue name, atom name, element, alternate location). |
| **Parameters** | `ID` — `cpstr` — The atom path string to parse. |
| **Returns**  | `0` if the parsing was successful, `-1` if the model number is invalid, `-2` if the sequence number is invalid. |
| **Behaviour** | The method parses the input string according to the specified format.  It extracts the `model`, `chain`, `seqNum`, `insCode`, `resName`, `atomName`, `element`, and `altLoc` from the input string. |
| **Side effects** | Modifies the `modelNo`, `chainID`, `seqNum`, `insCode`, `resName`, `atomName`, `element`, and `altLoc` members of the object. |
| **Errors / exceptions** | Returns a specific error code if the parsing fails (e.g., invalid numerical format). |

#### `void AtomPath::write (io::RFile f)`

| Field        | Detail |
|--------------|--------|
| **Summary**  | Writes the content of the member variables to the output stream `f`. |
| **Parameters** | `f` — `io::RFile*` — A pointer to the output stream file. |
| **Returns**  | `void` |
| **Behaviour** | Writes the value of the member variables to the output file stream `f`. |
| **Side effects** | Writes the data to the output stream `f`. |
| **Errors / exceptions** | `None` |

#### `void AtomPath::read (io::RFile f)`

| Field        | Detail |
|--------------|--------|
| **Summary**  | Reads the content of the member variables from the input stream `f`. |
| **Parameters** | `f` — `io::RFile*` — A pointer to the input stream file. |
| **Returns**  | `void` |
| **Behaviour** | Reads the value of the member variables from the input file stream `f`. |
| **Side effects** | Reads the data from the input stream `f`. |
| **Errors / exceptions** | `None` |

### Struct `mmdb::QuickSort`

> This class provides a quicksort algorithm for sorting integer data.

#### `int QuickSort::Compare(int i, int j)`

| Field        | Detail |
|--------------|--------|
| **Summary**  | Compares two elements in the array. |
| **Parameters** | `i` — `int` — Index of the first element. `j` — `int` — Index of the second element. |
| **Returns**  | `-1` if `i` is less than `j`, `1` if `i` is greater than `j`, and `0` if they are equal. |
| **Behaviour** | Compares the values at the indices `i` and `j` in the internal array. |
| **Side effects** | `None` |
| **Errors / exceptions** | `None` |

#### `void QuickSort::Swap(int i, int j)`

| Field        | Detail |
|--------------|--------|
| **Summary**  | Swaps the values of two elements in the array. |
| **Parameters** | `i` — `int` — Index of the first element. `j` — `int` — Index of the second element. |
| **Returns**  | `void` |
| **Behaviour** | Swaps the elements at indices `i` and `j` in the internal array. |
| **Side effects** | Modifies the internal array. |
| **Errors / exceptions** | `None` |

#### `void QuickSort::SelectionSort(int left, int right)`

| Field        | Detail                                                                                                                                       |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Sorts the elements in the array within the specified range (`left` to `right`) using the selection sort algorithm.                          |
| **Parameters** | `left` — `int` — The starting index of the range to sort. `right` — `int` — The ending index of the range to sort. |
| **Returns**  | `void` |
| **Behaviour** | Selects the minimum element in the unsorted portion of the array and swaps it with the element at the beginning of the unsorted portion.                |
| **Side effects** | Modifies the internal array. |
| **Errors / exceptions** | `None` |

#### `int QuickSort::Partition(int left, int right)`

| Field        | Detail                                                                                                                                                               |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Partitions the array around a pivot element (the element at the `right` index) and returns the index of the pivot element after the partition. |
| **Parameters** | `left` — `int` — The starting index of the partition. `right` — `int` — The ending index of the partition. |
| **Returns**  | `int` — The index of the pivot element after partitioning. |
| **Behaviour** |  Chooses the last element of the array as pivot, partitions the array with respect to this pivot and moves the pivot to its correct position. |
| **Side effects** | Modifies the internal array. |
| **Errors / exceptions** | `None` |

#### `void QuickSort::Quicksort(int left, int right)`

| Field        | Detail                                                                                                                                                   |
|--------------|----------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Recursively sorts the elements in the array within the specified range (`left` to `right`) using the quicksort algorithm.                               |
| **Parameters** | `left` — `int` — The starting index of the range to sort. `right` — `int` — The ending index of the range to sort. |
| **Returns**  | `void` |
| **Behaviour** |  Calls `Partition` to partition the array around a pivot element, and then recursively calls `Quicksort` on the sub-arrays to the left and right of the pivot. |
| **Side effects** | Modifies the internal array. |
| **Errors / exceptions** | `None` |

#### `void QuickSort::Sort(void * sortdata, int data_len)`

| Field        | Detail                                                                                                                                                           |
|--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Summary**  | Sorts the array pointed to by `sortdata` using the quicksort algorithm. |
| **Parameters** | `sortdata` — `void*` — Pointer to the data to be sorted. `data_len` — `int` — The number of elements to sort. |
| **Returns**  | `void` |
| **Behaviour** | Passes the data to the internal `QuickSort` method. |
| **Side effects** | Modifies the internal array. |
| **Errors / exceptions** | `None` |

---

This documentation provides a thorough understanding of each method's purpose, behavior, and potential issues. It's designed to be helpful for developers who are working with the `mmdb_utils` library.  Let me know if you would like me to elaborate on any specific aspects or generate documentation for additional code.

---

### `mmdb_xml_`

## `mmdb::xml::XMLObject`

This class is a flexible, hierarchical object representation for parsing and managing data from MMCIF files. It provides a way to store and access data elements, attributes, and relationships defined in the MMCIF format.

### Methods

#### `XMLObject()`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Creates a new `XMLObject` instance with no associated tag or data.       |
| **Parameters** | `none`                                                              |
| **Returns**  | `void`                                                              |
| **Behaviour** | Initializes all member variables to their default states, including tag, data, number of objects, and memory allocations. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `XMLObject(cpstr Tag)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Creates a new `XMLObject` instance with the given tag name and no data. |
| **Parameters** | `Tag` — A C-style string containing the tag name.                 |
| **Returns**  | `void`                                                               |
| **Behaviour** | Initializes all member variables, sets the `objTag` member to the provided tag, and initializes other internal state variables. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `XMLObject(cpstr Tag, cpstr Data)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Creates a new `XMLObject` instance with the given tag name and data. |
| **Parameters** | `Tag` — A C-style string containing the tag name.                   `Data` — A C-style string containing the data associated with the tag. |
| **Returns**  | `void`                                                              |
| **Behaviour** | Initializes all member variables, sets `objTag` to the provided tag, and sets `objData` to the provided data. Also, initializes internal state variables. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `XMLObject(cpstr Tag, realtype V, int length)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Creates a new `XMLObject` instance with the given tag name, a numerical value, and its length. |
| **Parameters** | `Tag` — A C-style string containing the tag name. `V` — A floating-point number. `length` — The length of the data associated with the tag. |
| **Returns**  | `void`                                                              |
| **Behaviour** | Initializes all member variables, sets `objTag` to the provided tag, and sets `objData` to the given numerical value with specified length.  |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `XMLObject(cpstr Tag, int iV, int length)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Creates a new `XMLObject` instance with the given tag name and an integer value. |
| **Parameters** | `Tag` — A C-style string containing the tag name. `iV` — An integer value. `length` — The length of the data associated with the tag. |
| **Returns**  | `void`                                                              |
| **Behaviour** | Initializes all member variables, sets `objTag` to the provided tag, and sets `objData` to the integer value with the specified length. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `XMLObject(cpstr Tag, bool bV)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Creates a new `XMLObject` instance with the given tag name and a boolean value. |
| **Parameters** | `Tag` — A C-style string containing the tag name. `bV` — A boolean value. |
| **Returns**  | `void`                                                              |
| **Behaviour** | Initializes all member variables, sets `objTag` to the provided tag, and sets `objData` to the boolean value. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `XMLObject(cpstr Tag, PXMLObject XMLObject)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Creates a new `XMLObject` instance and sets it as a child of the given XMLObject. |
| **Parameters** | `Tag` — A C-style string containing the tag name. `XMLObject` — A pointer to an existing `XMLObject` instance. |
| **Returns**  | `void`                                                              |
| **Behaviour** | Initializes all member variables, sets `objTag` to the provided tag, and sets the parent of this object to the given `XMLObject`. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `XMLObject(io::RPStream Object)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Creates a new `XMLObject` instance from an existing `RPStream` object.   |
| **Parameters** | `Object` — An `RPStream` object from which to construct the `XMLObject`. |
| **Returns**  | `void`                                                              |
| **Behaviour** | Initializes all member variables and copies the state from the `RPStream` object. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `~XMLObject()`

| Field        | Detail                                                              |
|--------------|---------------------------------------------------------------------|
| **Summary**  | Destructor for the `XMLObject` class, releasing allocated memory.     |
| **Parameters** | `none`                                                              |
| **Returns**  | `void`                                                              |
| **Behaviour** | Frees all dynamically allocated memory associated with the object, including `objTag`, `objData`, the object array (`object`), and attribute vectors (`attr_name`, `attr_value`). |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `void SetTag(cpstr Tag)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Sets the tag name of the `XMLObject`.                              |
| **Parameters** | `Tag` — A C-style string containing the new tag name.                |
| **Returns**  | `void`                                                              |
| **Behaviour** | Updates the `objTag` member with the provided tag name. Also, recalculates the length of `objTag`. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `void AddAttribute(cpstr name, cpstr value)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Adds an attribute to the `XMLObject`, storing its name and value.     |
| **Parameters** | `name` — A C-style string containing the attribute name.          `value` — A C-style string containing the attribute value. |
| **Returns**  | `void`                                                              |
| **Behaviour** | Adds the name and value to the respective vectors, incrementing the attribute count. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `void AddAttribute(cpstr name, const int iV)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Adds an integer attribute to the `XMLObject`.                        |
| **Parameters** | `name` — A C-style string containing the attribute name.          `iV` — An integer value. |
| **Returns**  | `void`                                                              |
| **Behaviour** | Converts the integer `iV` to a string and then adds the attribute to the `XMLObject`. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `void AddAttribute(cpstr name, const bool bV)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Adds a boolean attribute to the `XMLObject`.                         |
| **Parameters** | `name` — A C-style string containing the attribute name.          `bV` — A boolean value. |
| **Returns**  | `void`                                                              |
| **Behaviour** | Converts the boolean `bV` to a string ("Yes" or "No") and then adds the attribute to the `XMLObject`. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `void SetData(cpstr Data)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Sets the data associated with the `XMLObject`.                       |
| **Parameters** | `Data` — A C-style string containing the data to be stored.           |
| **Returns**  | `void`                                                              |
| **Behaviour** | Updates the `objData` member with the provided data. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `void AddData(cpstr Data)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Appends data to the data section of the `XMLObject`.                  |
| **Parameters** | `Data` — A C-style string containing the data to be added.            |
| **Returns**  | `void`                                                              |
| **Behaviour** | Appends the data to `objData`. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `void SetData(const realtype V, const int length)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Sets the data associated with the `XMLObject` to a numerical value.    |
| **Parameters** | `V` — A floating-point number. `length` — The length of the data. |
| **Returns**  | `void`                                                              |
| **Behaviour** | Converts the number `V` to a string and sets `objData` to this string, with the specified length. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `void SetData(const int iV, const int length)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Sets the data associated with the `XMLObject` to an integer value.      |
| **Parameters** | `iV` — An integer value. `length` — The length of the data. |
| **Returns**  | `void`                                                              |
| **Behaviour** | Converts the integer `iV` to a string and sets `objData` to this string, with the specified length. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `void SetData(const bool bV)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Sets the data associated with the `XMLObject` to a boolean value.     |
| **Parameters** | `bV` — A boolean value.                                              |
| **Returns**  | `void`                                                              |
| **Behaviour** | Converts the boolean `bV` to a string ("Yes" or "No") and sets `objData` to this string. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `int AddMMCIFCategory(mmcif::PCategory mmCIFCat)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Adds a MMCIF category to the `XMLObject`.                             |
| **Parameters** | `mmCIFCat` — A pointer to an MMCIF category object.                     |
| **Returns**  | `int` — The index of the added category within the `XMLObject`.        |
| **Behaviour** | Adds the `mmCIFCat` object to the `XMLObject`'s category list. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `int AddMMCIFStruct(mmcif::PStruct mmCIFStruct)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Adds a MMCIF structure to the `XMLObject`.                            |
| **Parameters** | `mmCIFStruct` — A pointer to an MMCIF structure object.                 |
| **Returns**  | `int` — The index of the added structure within the `XMLObject`.       |
| **Behaviour** | Adds the `mmCIFStruct` object to the `XMLObject`'s structure list.   |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `int AddMMCIFLoop(mmcif::PLoop mmCIFLoop)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Adds a MMCIF loop to the `XMLObject`.                                |
| **Parameters** | `mmCIFLoop` — A pointer to an MMCIF loop object.                       |
| **Returns**  | `int` — The index of the added loop within the `XMLObject`.            |
| **Behaviour** | Adds the `mmCIFLoop` object to the `XMLObject`'s loop list.          |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `int AddMMCIFData(mmcif::PData mmCIFData)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Adds a MMCIF data element to the `XMLObject`.                        |
| **Parameters** | `mmCIFData` — A pointer to an MMCIF data object.                       |
| **Returns**  | `int` — The number of added data elements.                           |
| **Behaviour** | Calls `AddMMCIFStruct` or `AddMMCIFLoop` depending on the category of `mmCIFData`. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `pstr GetData(cpstr Tag, int objNo)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Retrieves the data associated with the `XMLObject` based on tag name and object number. |
| **Parameters** | `Tag` — A C-style string containing the tag name.                 `objNo` — The index of the object within the `XMLObject`'s object list. |
| **Returns**  | `pstr` — A pointer to the data string associated with the specified tag and object number, or `NULL` if no such object exists. |
| **Behaviour** | Searches the `XMLObject`'s object list for the object with the specified tag name and object number. If found, returns a pointer to the object's `objData` member. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `XML_RC GetData(pstr & Data, cpstr Tag, int objNo)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Retrieves the data associated with the `XMLObject` based on tag name and object number, and copies the data to a given pointer. |
| **Parameters** | `Data` — A pointer to a `pstr` variable where the data will be stored. `Tag` — A C-style string containing the tag name. `objNo` — The index of the object within the `XMLObject`'s object list. |
| **Returns**  | `XML_RC` — An `XML_RC` enum value indicating the success or failure of the operation. |
| **Behaviour** | Searches the `XMLObject`'s object list for the object with the specified tag name and object number. If found, copies the object's `objData` member to the provided `Data` pointer and returns `XMLRC_Ok`. If not found, or if there is an error during data copying, returns the appropriate `XML_RC` value. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `XML_RC GetData(realtype & V, cpstr Tag, int objNo)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Retrieves a numerical value associated with the `XMLObject` based on tag name and object number. |
| **Parameters** | `V` — A reference to a `realtype` variable where the value will be stored. `Tag` — A C-style string containing the tag name. `objNo` — The index of the object within the `XMLObject`'s object list. |
| **Returns**  | `XML_RC` — An `XML_RC` enum value indicating the success or failure of the operation. |
| **Behaviour** | Searches the `XMLObject`'s object list for the object with the specified tag name and object number. If found, attempts to convert the object's `objData` member to a `realtype` value and stores it in the provided variable `V`. Returns `XMLRC_Ok` if the conversion is successful, or `XMLRC_RFormatError` if the conversion fails. If not found, or if there is an error during conversion, returns the appropriate `XML_RC` value. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `XML_RC GetData(int & iV, cpstr Tag, int objNo)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Retrieves an integer value associated with the `XMLObject` based on tag name and object number. |
| **Parameters** | `iV` — A reference to an `int` variable where the value will be stored. `Tag` — A C-style string containing the tag name. `objNo` — The index of the object within the `XMLObject`'s object list. |
| **Returns**  | `XML_RC` — An `XML_RC` enum value indicating the success or failure of the operation. |
| **Behaviour** | Searches the `XMLObject`'s object list for the object with the specified tag name and object number. If found, attempts to convert the object's `objData` member to an `int` value and stores it in the provided variable `iV`. Returns `XMLRC_Ok` if the conversion is successful, or `XMLRC_IFormatError` if the conversion fails. If not found, or if there is an error during conversion, returns the appropriate `XML_RC` value. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `XML_RC GetData(bool & bV, cpstr Tag, int objNo)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Retrieves a boolean value associated with the `XMLObject` based on tag name and object number. |
| **Parameters** | `bV` — A reference to a `bool` variable where the value will be stored. `Tag` — A C-style string containing the tag name. `objNo` — The index of the object within the `XMLObject`'s object list. |
| **Returns**  | `XML_RC` — An `XML_RC` enum value indicating the success or failure of the operation. |
| **Behaviour** | Searches the `XMLObject`'s object list for the object with the specified tag name and object number. If found, attempts to convert the object's `objData` member to a `bool` value and stores it in the provided variable `bV`. Returns `XMLRC_Ok` if the conversion is successful, or `XMLRC_OFormatError` if the conversion fails. If not found, or if there is an error during conversion, returns the appropriate `XML_RC` value. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `PXMLObject GetObject(cpstr Tag, int objNo)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Retrieves an `XMLObject` from the `XMLObject`'s object list based on tag name and object number. |
| **Parameters** | `Tag` — A C-style string containing the tag name. `objNo` — The index of the object within the `XMLObject`'s object list. |
| **Returns**  | `PXMLObject` — A pointer to the `XMLObject` with the specified tag and index, or `NULL` if no such object exists. |
| **Behaviour** | Searches the `XMLObject`'s object list for the object with the specified tag name and index. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `PXMLObject GetFirstObject()`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Returns a pointer to the first `XMLObject` in the `XMLObject`'s object list. |
| **Parameters** | `none`                                                              |
| **Returns**  | `PXMLObject` — A pointer to the first `XMLObject` in the list, or `NULL` if the list is empty. |
| **Behaviour** | Returns the first element in the object list, or `NULL` if the list is empty. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `PXMLObject GetLastObject()`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Returns a pointer to the last `XMLObject` in the `XMLObject`'s object list. |
| **Parameters** | `none`                                                              |
| **Returns**  | `PXMLObject` — A pointer to the last `XMLObject` in the list, or `NULL` if the list is empty. |
| **Behaviour** | Returns the last element in the object list, or `NULL` if the list is empty. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `PXMLObject GetObject(int objectNo)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Retrieves an `XMLObject` from the `XMLObject`'s object list based on its index. |
| **Parameters** | `objectNo` — The index of the object to retrieve within the `XMLObject`'s object list. |
| **Returns**  | `PXMLObject` — A pointer to the `XMLObject` with the specified index, or `NULL` if the index is out of range. |
| **Behaviour** | Searches the `XMLObject`'s object list for the object with the specified index. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `void AddObject(PXMLObject XMLObject, int lenInc)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Adds a new `XMLObject` to the end of the `XMLObject`'s object list.     |
| **Parameters** | `XMLObject` — A pointer to the `XMLObject` to be added. `lenInc` — An integer specifying the increment to the number of objects |
| **Returns**  | `void`                                                              |
| **Behaviour** | Adds the `XMLObject` to the end of the object list. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `void InsertObject(PXMLObject XMLObject, int pos, int lenInc)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Inserts a new `XMLObject` into the `XMLObject`'s object list at a specific position. |
| **Parameters** | `XMLObject` — A pointer to the `XMLObject` to be inserted. `pos` — The index where the `XMLObject` should be inserted. `lenInc` — An integer specifying the increment to the number of objects |
| **Returns**  | `void`                                                              |
| **Behaviour** | Inserts the `XMLObject` into the object list at the specified position. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `XML_RC WriteObject(cpstr FName, int pos, int indent)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Writes the `XMLObject` to a file.                                    |
| **Parameters** | `FName` — A C-style string containing the name of the file to write to. `pos` — The position in the file to start writing at.  `indent` — The level of indentation. |
| **Returns**  | `XML_RC` — An `XML_RC` enum value indicating the success or failure of the operation. |
| **Behaviour** | Writes the `XMLObject` to the specified file, with appropriate indentation. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `void WriteObject(io::RFile f, int pos, int indent)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Writes the `XMLObject` to an `io::RFile` object.                      |
| **Parameters** | `f` — An `io::RFile` object to write to. `pos` — The position in the file to start writing at. `indent` — The level of indentation. |
| **Returns**  | `void`                                                              |
| **Behaviour** | Writes the `XMLObject` to the `io::RFile` object, with appropriate indentation. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `XML_RC ReadObject(cpstr FName)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Reads a `XMLObject` from a file.                                      |
| **Parameters** | `FName` — A C-style string containing the name of the file to read from. |
| **Returns**  | `XML_RC` — An `XML_RC` enum value indicating the success or failure of the operation. |
| **Behaviour** | Reads a `XMLObject` from the specified file. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

#### `XML_RC ReadObject(io::RFile f, pstr S, int & pos, int slen)`

| Field        | Detail                                                               |
|--------------|----------------------------------------------------------------------|
| **Summary**  | Reads a string from an `io::RFile` object and stores it in a given pointer. |
| **Parameters** | `f` — An `io::RFile` object to read from. `S` — A pointer to a `pstr` variable where the string will be stored. `pos` — A reference to an integer variable that will hold the current position in the file. `slen` — The maximum length of the string to read. |
| **Returns**  | `XML_RC` — An `XML_RC` enum value indicating the success or failure of the operation. |
| **Behaviour** | Reads a string from the `io::RFile` object and stores it in the provided `S` pointer, along with the current position in the file. |
| **Side effects** | `none`                                                              |
| **Errors / exceptions** | `none`                                                              |

## Free Functions


---
