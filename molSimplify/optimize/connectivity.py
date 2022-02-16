import numpy as np
import ase.data
import itertools


def find_connectivity(atoms, threshold=1.25, connect_fragments=True):
    """
    Follows the algorithm outlined in Section II A of
    Billeter et al. Phys. Chem. Chem. Phys., 2000, 2, 2177-2186
    """
    bonds = []
    N = len(atoms)
    xyzs = atoms.get_positions()
    types = atoms.get_atomic_numbers()

    r2 = np.zeros((N, N))
    np.fill_diagonal(r2, np.inf)
    neighbors_per_atom = [[] for _ in range(N)]

    for i in range(N):
        for j in range(i+1, N):
            r2[i, j] = r2[j, i] = np.sum((xyzs[i, :] - xyzs[j, :])**2)
            # The paper clearly states that squared distances are compared
            # without squaring the threshold:
            if (r2[i, j] < threshold * (
                    ase.data.covalent_radii[types[i]]
                    + ase.data.covalent_radii[types[j]])**2):
                bonds.append((i, j))
                neighbors_per_atom[i].append(j)
                neighbors_per_atom[j].append(i)

    if connect_fragments:
        # Check for disconnected fragments
        def depth_first_search(fragment, i, visited):
            visited[i] = True
            fragment.append(i)
            for j in neighbors_per_atom[i]:
                if not visited[j]:
                    fragment = depth_first_search(fragment, j, visited)
            return fragment

        visited = [False] * N
        fragments = []
        for i in range(N):
            if not visited[i]:
                fragment = []
                fragments.append(depth_first_search(fragment, i, visited))

        # If there are more than 1 fragment connect shortest distance atoms
        while len(fragments) > 1:
            # Merge first fragment with nearest fragment by adding a bond
            r_min = np.inf
            bond_to_add = ()
            for frag_ind, other_frag in enumerate(fragments[1:]):
                for i in fragments[0]:
                    for j in other_frag:
                        if r2[i, j] < r_min:
                            r_min = r2[i, j]
                            # Plus one because of iterating fragments[1:]
                            fragment_to_merge = frag_ind + 1
                            bond_to_add = tuple(sorted([i, j]))
            bonds.append(bond_to_add)
            fragments[0].extend(fragments.pop(fragment_to_merge))
    return bonds


def find_primitives(xyzs, bonds, linear_threshold=5., planar_threshold=0.95):
    """
    Finds primitive internals given a reference geometry and connectivity list.
    Follows the algorithm outlined in Section II A of
    Billeter et al. Phys. Chem. Chem. Phys., 2000, 2, 2177-2186.
    Throughout the code literal citations of this paper appear in quotes "".
    linear_threshold: threshold angle in degrees to determine when to replace
                      close to linear bends with a coplanar and perpendicular
                      bend coordinate.
    planar_threshold: Set larger 1.0 to turn off detection of planar
                      structures.
    """
    neighbors = [[] for _ in range(len(xyzs))]
    for b in bonds:
        neighbors[b[0]].append(b[1])
        neighbors[b[1]].append(b[0])
    neighbors = list(map(sorted, neighbors))
    bends = []
    linear_bends = []
    torsions = []
    planars = []
    # Transform the threshold angle to cos(angle) to avoid repeated
    # error prone calls to arccos.
    cos_lin_thresh = np.abs(np.cos(linear_threshold*np.pi/180.))

    def cos_angle(r1, r2):
        """
        Helper function that calculates the cosine of the angle between
        two vectors
        """
        return np.dot(r1, r2)/(np.linalg.norm(r1)*np.linalg.norm(r2))

    for a in range(len(xyzs)):
        for i, ai in enumerate(neighbors[a]):
            for aj in neighbors[a][i+1:]:
                r_ai = xyzs[ai, :] - xyzs[a, :]
                r_aj = xyzs[aj, :] - xyzs[a, :]
                cos_theta = np.abs(cos_angle(r_ai, r_aj))
                # "Almost linear bends are dropped." Implemented by comparsion
                # to a threshold angle.
                if (cos_theta < cos_lin_thresh):
                    bends.append((ai, a, aj))
                    for ak in neighbors[ai]:
                        if (ak != a and ak != aj
                                and (aj, a, ai, ak) not in torsions
                                and (ak, ai, a, aj) not in torsions):
                            # Check if (ak, ai, a) is linear
                            r_ik = xyzs[ak, :] - xyzs[ai, :]
                            cos_phi = cos_angle(r_ik, r_ai)
                            if np.abs(cos_phi) < 0.99:
                                torsions.append((ak, ai, a, aj))
                    for ak in neighbors[aj]:
                        if (ak != a and ak != ai
                                and (ak, aj, a, ai) not in torsions
                                and (ai, a, aj, ak) not in torsions):
                            # Check if (a, aj, ak) is linear
                            r_jk = xyzs[ak, :] - xyzs[aj, :]
                            cos_phi = cos_angle(r_jk, r_aj)
                            if np.abs(cos_phi) < 0.99:
                                torsions.append((ai, a, aj, ak))
                else:
                    # "If the centre of such an angle is connected to exactly
                    # two atoms, the dropped angle is replaced by..." a linear
                    # bend.
                    if len(neighbors[a]) == 2:
                        # Add linear bend
                        linear_bends.append((ai, a, aj))
                        # Do not add torsions on linear bends, instead try
                        # using ai and aj as center bond. TODO: This does not
                        # detect longer linear chains
                        for ak in neighbors[ai]:
                            for am in neighbors[aj]:
                                if (ak != a and am != a
                                        and (ak, ai, aj, am) not in torsions
                                        and (am, aj, ai, ak) not in torsions):
                                    r_ij = xyzs[aj, :] - xyzs[ai, :]
                                    r_ik = xyzs[ak, :] - xyzs[ai, :]
                                    r_jm = xyzs[am, :] - xyzs[aj, :]
                                    cos_t1 = np.abs(cos_angle(r_ij, r_ik))
                                    cos_t2 = np.abs(cos_angle(r_ij, r_jm))
                                    if cos_t1 < 0.99 and cos_t2 < 0.99:
                                        torsions.append((ak, ai, aj, am))
        # "If an atom is connected to more than two atoms" and the threshold
        # is small enough to be exceeded by the product of two unit vectors.
        if len(neighbors[a]) > 2 and planar_threshold < 1.0:
            # "and lies in the centre of a planar system" characterized by
            # a central atom 'a' and the plane spanned by three neighboring
            # atoms (ai, aj, ak)
            for (ai, aj, ak) in itertools.combinations(neighbors[a], 3):
                r_ai = xyzs[ai, :] - xyzs[a, :]
                r_aj = xyzs[aj, :] - xyzs[a, :]
                r_ak = xyzs[ak, :] - xyzs[a, :]
                n1 = np.cross(r_ai, r_aj)
                n1 /= np.linalg.norm(n1)
                n2 = np.cross(r_aj, r_ak)
                n2 /= np.linalg.norm(n2)
                n3 = np.cross(r_ak, r_ai)
                n3 /= np.linalg.norm(n3)
                if (np.abs(n1.dot(n2)) > planar_threshold
                        or np.abs(n2.dot(n3)) > planar_threshold
                        or np.abs(n3.dot(n1)) > planar_threshold):
                    # Remove bend
                    bends.remove((ai, a, aj))
                    # Try to find an improper (b, a, c, d)
                    # such neither the angle t1 between (b, a, c)
                    # nor t2 between (a, c, d) is close to linear
                    for (b, c, d) in itertools.permutations([ai, aj, ak], 3):
                        r_ab = xyzs[b, :] - xyzs[a, :]
                        r_ac = xyzs[c, :] - xyzs[a, :]
                        r_cd = xyzs[d, :] - xyzs[c, :]
                        cos_t1 = cos_angle(r_ab, r_ac)
                        cos_t2 = cos_angle(r_ac, r_cd)
                        if np.abs(cos_t1) < 0.95 and np.abs(cos_t2 < 0.95):
                            planars.append((a, b, c, d))
                            break
                    # Break after one improper has been added
                    break
    return bends, linear_bends, torsions, planars
