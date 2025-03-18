# T3Design

A collection of script for TRIP and TWIP transformations.

## `TRIP.py`

`TRIP_unitcell(V_bcc: float, V_hcp: float = None, coa: float = np.sqrt(8/3), delta_shape: float = 0, delta_shuffle: float = 0) -> Atoms`

Compute a unit cell structure corresponding to a TRIP (Transformation-Induced Plasticity) effect.

This function calculates a unit cell structure for a TRIP effect, which involves the transformation
from a BCC (Body-Centered Cubic) to an HCP (Hexagonal Close-Packed) structure.

*Parameters*

- `V_bcc` (float):
    The specific volume (per atom) of a BCC unit cell.
- `V_hcp` (float, optional):
    The specific volume (per atom) of an HCP unit cell. Defaults to None.
- `coa` (float):
    c-over-a ratio of the HCP structure. Defaults to sqrt(8/3).
- `delta_shape` (float):
    A shape parameter. Defaults to 0.
- `delta_shuffle` (float):
    A shuffle parameter. Defaults to 0.

*Returns*
- Atoms
    The resulting unit cell structure represented as an ASE Atoms object.

`TRIP_unitcell_bcc_omega(V_bcc: float, V_omega: float = None, coa: float = 0.619, delta_shape: float = 0, delta_shuffle: float = 0) -> Atoms`

Compute a unit cell structure corresponding to a TRIP (Transformation-Induced Plasticity) effect.

This function calculates a unit cell structure for a TRIP effect, which involves the transformation
from a BCC (Body-Centered Cubic) to a hexagonal omega structure.

*Parameters*
- `V_bcc` : float
    The specific volume (per atom) of a BCC unit cell.
- `V_omega` : float, optional
    The specific volume (per atom) of an omega unit cell. Defaults to None.
- `coa` : float
    c-over-a ratio of the omega structure. Defaults to 0.619.
- `delta_shape` : float
    A shape parameter. Defaults to 0.
- `delta_shuffle` : float
    A shuffle parameter. Defaults to 0.

*Returns*
- Atoms
    The resulting unit cell structure represented as an ASE Atoms object.


## `TWIP.py`

`TWIP_332(V_bcc: float, zeta: float) -> Atoms`

Generate a crystal structure representing a 332-oriented BCC (body-centered cubic) lattice with a twinning transformation.

*Parameters*
- `V_bcc` (float):
    Volume of the BCC unit cell.
- `zeta` (float):
    Parameter for the twinning transformation.

*Returns*
- Atoms
    ASE Atoms object representing the crystal structure after twinning transformation.


`TWIP_112(V_bcc: float, zeta: float) -> Atoms`

Generate a crystal structure representing a 112-oriented BCC (body-centered cubic) lattice with a twinning transformation.

*Parameters*
- `V_bcc` : float
    Volume of the BCC unit cell.
`zeta` : float
    Parameter for the twinning transformation.

*Returns*
- Atoms
    ASE Atoms object representing the crystal structure after twinning transformation.


## `utils.py`

`map_atoms(str1, str2, species=True)`

Maps atoms from one atomic structure to another based on their scaled positions.

This function attempts to match each atom in `str1` to the closest corresponding atom
in `str2`, considering periodic boundary conditions by using images of `str2`.
If `species` is `True`, the mapping returns atomic symbols instead of indices.

*Parameters*
- `str1` : ASE Atoms object
    The first atomic structure to be mapped.

- `str2` : ASE Atoms object
    The second atomic structure, which serves as the reference for mapping.

- `species` : bool, optional (default=`True`)
    If `True`, returns atomic symbols from `str2` corresponding to the mapped indices.
    If `False`, returns the indices of the mapped atoms in `str2`.

*Returns*
- list
    A list of atomic indices or atomic symbols (if `species=True`), representing
    the best matches from `str1` to `str2`.

*Raises*
- Exception
    If `str1` and `str2` have different numbers of atoms.

*Warning*: If the mapping contains duplicates, indicating possible mismatches.
