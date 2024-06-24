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