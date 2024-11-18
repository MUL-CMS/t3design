def map_atoms(str1, str2, species=True):
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