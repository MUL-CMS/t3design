def map_atoms(str1, str2, species=True):
    """
    Maps atoms from one atomic structure to another based on their scaled positions.

    This function attempts to match each atom in `str1` to the closest corresponding atom
    in `str2`, considering periodic boundary conditions by using images of `str2`.
    If `species` is True, the mapping returns atomic symbols instead of indices.

    Parameters:
    -----------
    str1 : ASE Atoms object
        The first atomic structure to be mapped.

    str2 : ASE Atoms object
        The second atomic structure, which serves as the reference for mapping.

    species : bool, optional (default=True)
        If True, returns atomic symbols from `str2` corresponding to the mapped indices.
        If False, returns the indices of the mapped atoms in `str2`.

    Returns:
    --------
    list
        A list of atomic indices or atomic symbols (if `species=True`), representing
        the best matches from `str1` to `str2`.

    Raises:
    -------
    Exception
        If `str1` and `str2` have different numbers of atoms.

    Warning
        If the mapping contains duplicates, indicating possible mismatches.
    """

    from itertools import product
    from numpy import dot, concatenate
    from numpy.linalg import norm
    
    if len(str1) != len(str2):
        raise Exception(f'Different sizes of the two provided structures: {len(str1)} != {len(str2)}')
        return None
    Nat = len(str1)
    
    
    pos2 = []
    ind2 = []
    for im in product(range(-1, 2), range(-1, 2), range(-1, 2)):
        # shift = dot(str2.cell, im)
        if len(pos2) == 0:
            # pos2 = str2.positions+shift
            pos2 = str2.get_scaled_positions()+im
        else:
            # pos2 = concatenate([pos2, str2.positions+shift])
            pos2 = concatenate([pos2, str2.get_scaled_positions()+im])
        ind2 += range(Nat)
        
    # pos1 = str1.positions
    pos1 = str1.get_scaled_positions()
    
    mapping = []
    for p1 in pos1:
        dist = [norm(p2 - p1) for p2 in pos2]
        mapping.append(ind2[dist.index(min(dist))])
    
    if len(set(mapping)) != Nat:
        raise Warning(f'Careful, mapping may failed as it contains duplicates!')
    
    if species:
        mapping = [str2[i].symbol for i in mapping]
        
    return mapping
