from ase import Atoms
import numpy as np

def TRIP_unitcell(V_bcc: float, V_hcp: float = None, delta_shape: float = 0, delta_shuffle: float = 0) -> Atoms:
    """
    Compute a unit cell structure corresponding to a TRIP (Transformation-Induced Plasticity) effect.

    This function calculates a unit cell structure for a TRIP effect, which involves the transformation
    from a BCC (Body-Centered Cubic) to an HCP (Hexagonal Close-Packed) structure.

    Args:
        V_bcc (float): The volume of a BCC unit cell.
        V_hcp (float, optional): The volume of an HCP unit cell. Defaults to None.
        delta_shape (float): A shape parameter. Defaults to 0.
        delta_shuffle (float): A shuffle parameter. Defaults to 0.

    Returns:
        ase.Atoms: The resulting unit cell structure represented as an ASE Atoms object.

    """
    
    # create a bcc cell, in the orthorhombic setting, as described in 
    #   N. Abdoshahi, M. Dehghani, A.V. Ruban, M. Friák, M. Šob, J. Spitaler, 
    #   D. Holec, On the energetics of the cubic-to-hexagonal transformations 
    #   in TiAl+Mo alloys, Acta Mater. 240 (2022) 118268.
    # lattice parameter of conventional bcc cell
    a_bcc = (2*V_bcc)**(1/3)
    struct_bcc = Atoms(
        'Ti4',
        scaled_positions = [
            [0, 0, 0],
            [0, 1/2, 1/2],
            [1/2, 1/2, 0],
            [1/2, 0, 1/2]],
        cell = a_bcc*np.array([np.sqrt(2), 1, np.sqrt(2)]),
        pbc = [True, True, True]
        )

    if V_hcp == None:
        # no volume of the hcp structure was provided, let's take the same
        # specific volume as for the bcc cell
        V_hcp = V_bcc

if __name__ == '__main__':
    TRIP_unitcell(3.2**3/2)