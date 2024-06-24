from ase import Atoms
import numpy as np


def get_332_cell(a_bcc: float) -> Atoms:
    """
    Generate a crystal structure representing a 332-oriented BCC (body-centered cubic) lattice of titanium.

    Parameters
    ----------
    a_bcc : float
        Lattice parameter for the BCC lattice.

    Returns
    -------
    Atoms
        ASE Atoms object representing the crystal structure.
    """
    # Generate positions for the 332-orientated BCC lattice
    positions = [np.array([i, 0, (22 - i * 4) % 22]) for i in range(11)] + \
                 [np.array([i, 1, (11 + 22 - i * 4) % 22]) for i in range(11)]

    # Normalize positions
    positions /= np.array([11, 2, 22])

    # Create and return the Atoms object
    return Atoms(
        cell=a_bcc * np.array([5.5 / (3 ** 2 + 1 + 1) ** 0.5, 2 / (1 + 1 + 0) ** 0.5, 22 / (2 ** 2 + 3 ** 2 + 3 ** 2) ** 0.5]),
        symbols='Ti22',
        scaled_positions=np.array(positions),
        pbc=True
    )



def TWIP_332(V_bcc: float, zeta: float) -> Atoms:
    """
    Generate a crystal structure representing a 332-oriented BCC (body-centered cubic) lattice with a twinning transformation.

    Parameters
    ----------
    V_bcc : float
        Volume of the BCC unit cell.
    zeta : float
        Parameter for the twinning transformation.

    Returns
    -------
    Atoms
        ASE Atoms object representing the crystal structure after twinning transformation.
    """
    # Based on following paper and its supplementary material:
    # P. Kwasniak, F. Sun, S. Mantri, R. Banerjee, F. Prima
    # Polymorphic nature of {332}〈113〉twinning mode in BCC alloys
    # Materials Research Letters. 10 (2022) 334–342
    
    # get bcc lattice parameter
    a_bcc = (V_bcc * 2) ** (1 / 3)

    # construct a prototype
    twin_cell = get_332_cell(a_bcc=a_bcc)
    # parameters for the twinning transformation
    bbeta = twin_cell.cell.cellpar()[1]
    vbeta = twin_cell.cell.cellpar()[0]
    wbeta = twin_cell.cell.cellpar()[2]
    # factor 5.5 is correcting the displacement delta4 for the supercell size
    delta4 = (3 * a_bcc ** 2 - bbeta ** 2) / (4 * vbeta) * 5.5

    # construct model
    twin = twin_cell.copy()
    pos = twin.get_positions()
    for i in range(len(pos)):
        # applying delta_4
        pos[i, 0] += pos[i, 2] * zeta * delta4 / wbeta
        if abs(pos[i, 1]) < 0.1:
            # applying delta_2
            pos[i, 2] += zeta * wbeta / 22
        else:
            # applying delta_3
            pos[i, 2] -= zeta * wbeta / 22
    twin.set_positions(pos)

    # new cell corresponding to the above displacements
    twin.set_cell([[vbeta, 0, 0], [0, bbeta, 0], [zeta*delta4, 0, wbeta]], scale_atoms=False)
    # apply additional shear, thereby recovering the PBC
    twin.set_cell([[vbeta, 0, 0], [0, bbeta, 0], [zeta*vbeta, 0, wbeta]], scale_atoms=True)

    return twin




def get_112_cell(a_bcc: float) -> Atoms:
    """
    Generate a crystal structure representing a 112-oriented BCC (body-centered cubic) lattice of titanium.

    Parameters
    ----------
    a_bcc : float
        Lattice parameter for the BCC lattice.

    Returns
    -------
    Atoms
        ASE Atoms object representing the crystal structure.
    """
    # Generate positions for the 112-orientated BCC lattice
    positions = [np.array([i, 0, (6-2*i)%6]) for i in range(3)] + [np.array([i, 1, (6-2*i+3)%6]) for i in range(3)]

    # Normalize positions
    positions /= np.array([3, 2, 6])

    # Create and return the Atoms object
    return Atoms(
        cell = a_bcc*np.array([(1+1+1)**0.5/2, (1+1+0)**0.5, (2**2+1+1)**0.5]),
        symbols = 'H6',
        scaled_positions=np.array(positions),
        pbc=True
    )



def TWIP_112(V_bcc: float, zeta: float) -> Atoms:
    """
    Generate a crystal structure representing a 112-oriented BCC (body-centered cubic) lattice with a twinning transformation.

    Parameters
    ----------
    V_bcc : float
        Volume of the BCC unit cell.
    zeta : float
        Parameter for the twinning transformation.

    Returns
    -------
    Atoms
        ASE Atoms object representing the crystal structure after twinning transformation.
    """
    # Based on following paper (Figure 4):
    # H. Tobe, H.Y. Kim, T. Inamura, H. Hosoda, S. Miyazaki
    # Origin of 332 twinning in metastable β-Ti alloys
    # Acta Mater. 64 (2014) 345–3550 (2022) 334–342
    
    # get bcc lattice parameter
    a_bcc = (V_bcc * 2) ** (1 / 3)

    # construct a prototype
    twin_cell = get_112_cell(a_bcc=a_bcc)
    # parameters for the twinning transformation
    bbeta = twin_cell.cell.cellpar()[1]
    vbeta = twin_cell.cell.cellpar()[0]
    wbeta = twin_cell.cell.cellpar()[2]
  
    # magnitudes of shuffling for a fully transformed twin  
    shifts = np.array([
        np.array([0, 0, 0]),
        np.array([0, 0, 0]),
        np.array([0, 0, 0]),
        np.array([vbeta/2, 0, 0]),
        np.array([-vbeta/2, 0, 0]),
        np.array([-vbeta/2, 0, 0])
    ])

    # construct model
    twin = twin_cell.copy()
    pos = twin.get_positions()
    # apply shear
    twin.set_cell([[vbeta, 0, 0], [0, bbeta, 0], [-zeta*vbeta, 0, wbeta]], scale_atoms=True)
    # apply additional shuffling
    twin.set_positions(twin.positions + zeta*shifts)

    return twin