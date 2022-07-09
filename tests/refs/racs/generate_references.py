import pickle
from molSimplify.Classes.mol3D import mol3D
from molSimplify.Classes.atom3D import atom3D
from molSimplify.Classes.ligand import ligand
from molSimplify.Scripts.io import lig_load
from molSimplify.Informatics.RACassemble import assemble_connectivity_from_parts
from molSimplify.Informatics.lacRACAssemble import get_descriptor_vector


def generate_features(metal_str, eq_ligands, ax_ligands):
    metal = mol3D()
    metal.addAtom(atom3D(Sym=metal_str))

    eq_ligand_list = []
    for lig_str in eq_ligands:
        # Note: only works if connecting atom is 0
        lig = ligand(mol3D(), [0], 1)
        lig.mol = lig_load(lig_str)[0]
        eq_ligand_list.append(lig)

    ax_ligand_list = []
    for lig_str in ax_ligands:
        # Note: only works if connecting atom is 0
        lig = ligand(mol3D(), [0], 1)
        lig.mol = lig_load(lig_str)[0]
        ax_ligand_list.append(lig)

    print(ax_ligand_list)
    ligand_dict = {'eq_ligand_list': eq_ligand_list,
                   'ax_ligand_list': ax_ligand_list,
                   'eq_con_int_list': [h.mol.cat for h in eq_ligand_list],
                   'ax_con_int_list': [h.mol.cat for h in ax_ligand_list]
                   }
    mol = assemble_connectivity_from_parts(metal, ligand_dict)
    names, racs = get_descriptor_vector(mol, custom_ligand_dict=ligand_dict)
    features = dict(zip(names, racs))
    return features


def main():
    features = generate_features('Fe',
                                 ['carbonyl']*4,
                                 ['carbonyl']*2)
    with open('racs_Fe_(CO)_6.pickle', 'wb') as f:
        pickle.dump(features, f)

    features = generate_features('Mn',
                                 ['furan', 'water', 'ammonia', 'furan'],
                                 ['water', 'ammonia'])
    with open('racs_Mn_furan_water_ammonia_furan_water_ammonia.pickle', 'wb') as f:
        pickle.dump(features, f)


if __name__ == '__main__':
    main()
