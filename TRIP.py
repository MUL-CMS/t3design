from ase import Atoms
import numpy as np

def TRIP_unitcell(V_bcc: float, V_hcp: float = None, coa: float = np.sqrt(8/3),
                  delta_shape: float = 0, delta_shuffle: float = 0) -> Atoms:
    """
    Compute a unit cell structure corresponding to a TRIP (Transformation-Induced Plasticity) effect.

    This function calculates a unit cell structure for a TRIP effect, which involves the transformation
    from a BCC (Body-Centered Cubic) to an HCP (Hexagonal Close-Packed) structure.

    Parameters
    ----------
    V_bcc : float
        The specific volume (per atom) of a BCC unit cell.
    V_hcp : float, optional
        The specific volume (per atom) of an HCP unit cell. Defaults to None.
    coa : float
        c-over-a ratio of the HCP structure. Defaults to sqrt(8/3).
    delta_shape : float
        A shape parameter. Defaults to 0.
    delta_shuffle : float
        A shuffle parameter. Defaults to 0.

    Returns
    -------
    Atoms
        The resulting unit cell structure represented as an ASE Atoms object.
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
            [1/2, 1/2, 0],
            [0, 1/2, 1/2],
            [1/2, 0, 1/2]],
        cell = a_bcc*np.array([np.sqrt(2), 1, np.sqrt(2)]),
        pbc = [True, True, True]
    )

    if V_hcp == None:
        # no volume of the hcp structure was provided, let's take the same
        # specific volume as for the bcc cell
        V_hcp = V_bcc
    # 2*V_hcp = Valpha = sqrt(3)/2*coa*a_hcp^3
    a_hcp = (4*V_hcp/(coa*3**(1/2)))**(1/3)
    c_hcp = coa*a_hcp
    # create a hcp cell (2 unit cells) suitable for the Burger's transformation
    struct_hcp = Atoms(
        'Ti4',
        scaled_positions = [
            [0, 0, 0],
            [1/2, 1/2, 0],
            [1/6, 1/2, 1/2],
            [4/6, 0, 1/2]],
        cell = [np.sqrt(3)*a_hcp, a_hcp, c_hcp],
        pbc = [True, True, True]
    )
    
    # generate and return the final structure based on the shape and shuffle parameters
    struct = Atoms(
        'Ti4',
        scaled_positions = (1-delta_shuffle)*struct_bcc.get_scaled_positions() + \
                           delta_shuffle*struct_hcp.get_scaled_positions(),
        cell = (1-delta_shape)*struct_bcc.cell + delta_shape*struct_hcp.cell,
        pbc = [True, True, True]
    )
    return struct


if __name__ == '__main__':

    # Test: create a set of POSCARs which represent the bcc->hcp transformation 
    # using the Burger's OR; the output is written as large supercells for
    # nicer visualisation

    from os import mkdir, chdir
    from os.path import exists
    from ase.io import write

    V_beta = 3.25**3/2
    N_images = 15

    shape   = np.linspace(0, 1, N_images).tolist() + np.linspace(1, 0, N_images).tolist() + \
              np.linspace(0, 0, N_images).tolist() + np.linspace(0, 0, N_images).tolist() + \
              np.linspace(0, 1, N_images).tolist()
    shuffle = np.linspace(0, 0, N_images).tolist() + np.linspace(0, 0, N_images).tolist() + \
              np.linspace(0, 1, N_images).tolist() + np.linspace(1, 0, N_images).tolist() + \
              np.linspace(0, 1, N_images).tolist()
    
    if not exists('TRIP_path'):
        mkdir('TRIP_path')
    chdir('TRIP_path')
    
    structs = [TRIP_unitcell(V_beta, delta_shape=c, delta_shuffle=s) for c, s in zip(shape, shuffle)]
    for i, (c, s) in enumerate(zip(shape, shuffle)):
        struct = TRIP_unitcell(V_beta, delta_shape=c, delta_shuffle=s)
        write(filename=f'POSCAR_{i:05d}.vasp', images=struct*(5,5,4), format='vasp', direct=True)
    chdir('..')