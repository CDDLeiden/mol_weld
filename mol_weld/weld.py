from rdkit import Chem
from collections import defaultdict
from rdkit.Chem.rdchem import EditableMol

def weld_r_groups(core: Chem.rdchem.Mol, rgroups: Chem.rdchem.Mol) -> Chem.rdchem.Mol:
    """RDKit manipulations to put r-groups back on a core.

    Adapted from 
    https://sourceforge.net/p/rdkit/mailman/rdkit-discuss/thread/CANPfuvAqWxR%2BosdH6TUT-%2B1Fy85fUXh0poRddrEQDxXmguJJ7Q%40mail.gmail.com/
    AtomMapNum: number to keep track of what is connected to what
    Args:
        input_mol (Chem.rdchem.Mol): combined mol object of core and r-groups

    Returns:
        final_mol (Chem.rdchem.Mol|None): welded molecule
    """
    core = Chem.Mol(core)
    rgroups = Chem.Mol(rgroups)

    core_noatommap = Chem.Mol(core)
    for atom in core_noatommap.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num > 0:
            atom.SetAtomMapNum(0)

    rgroups = Chem.CombineMols(core_noatommap, rgroups)

    # loop over atoms and find the atoms with an AtomMapNum
    join_dict_core = defaultdict(
        list)  # list of atoms to connect for the core
    for atom in core.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num > 0:
            join_dict_core[map_num].append(atom)
    join_dict_rgroup = defaultdict(
        list)  # list of atoms to connect for each rgroup
    for atom in rgroups.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num > 0:
            join_dict_rgroup[map_num].append(atom)

    # loop over the bond between atom and dummy atom of r-group and save bond type
    bond_order_dict = defaultdict(list)
    chiral_tag_dict = {
    }  # dict to store chiral tag of the source of the branching off points
    source_bond_order_dict = {
    }  # dict to store bond order in the source of the branching off points

    # rgroups can have multiple attachment points
    for idx, atom_list in join_dict_rgroup.items():
        if len(atom_list) == 1:
            atm_1 = atom_list[0]
            bond_order_dict[idx].extend([x.GetBondType() for x in atm_1.GetBonds()])
        elif len(atom_list) == 2:
            atm_1, atm_2 = atom_list
            bond_order_dict[idx].extend([x.GetBondType() for x in atm_1.GetBonds()])
            bond_order_dict[idx].extend([x.GetBondType() for x in atm_2.GetBonds()])
        elif len(atom_list) == 3:
            atm_1, atm_2, atm_3 = atom_list
            bond_order_dict[idx].extend([x.GetBondType() for x in atm_1.GetBonds()])
            bond_order_dict[idx].extend([x.GetBondType() for x in atm_2.GetBonds()])
            bond_order_dict[idx].extend([x.GetBondType() for x in atm_3.GetBonds()])
        elif len(atom_list) == 4:
            atm_1, atm_2, atm_3, atm_4 = atom_list
            bond_order_dict[idx].extend([x.GetBondType() for x in atm_1.GetBonds()])
            bond_order_dict[idx].extend([x.GetBondType() for x in atm_2.GetBonds()])
            bond_order_dict[idx].extend([x.GetBondType() for x in atm_3.GetBonds()])
            bond_order_dict[idx].extend([x.GetBondType() for x in atm_4.GetBonds()])
    
    # transfer the atom maps to the neighbor atoms in core
    for idx, atom_list in join_dict_core.items():
        if len(atom_list) == 1:
            atm_1 = atom_list[0]
            nbr_1 = [x.GetOtherAtom(atm_1) for x in atm_1.GetBonds()]
            for nbr in nbr_1:
                nbr.SetAtomMapNum(idx)
        else:
            print('more atoms in list than expected')
    # transfer the atom maps to the neighbor atoms in rgroups
    for idx, atom_list in join_dict_rgroup.items():
        if len(atom_list) == 1:
            atm_1 = atom_list[0]
            nbr_1 = [x.GetOtherAtom(atm_1) for x in atm_1.GetBonds()]
            for nbr in nbr_1:
                nbr.SetAtomMapNum(idx)
        if len(atom_list) == 2:
            atm_1, atm_2 = atom_list
            nbr_1 = [x.GetOtherAtom(atm_1) for x in atm_1.GetBonds()]
            nbr_2 = [x.GetOtherAtom(atm_2) for x in atm_2.GetBonds()]
            for nbr in nbr_1 + nbr_2:
                nbr.SetAtomMapNum(idx)
        elif len(atom_list) == 3:
            atm_1, atm_2, atm_3 = atom_list
            nbr_1 = [x.GetOtherAtom(atm_1) for x in atm_1.GetBonds()]
            nbr_2 = [x.GetOtherAtom(atm_2) for x in atm_2.GetBonds()]
            nbr_3 = [x.GetOtherAtom(atm_3) for x in atm_3.GetBonds()]
            for nbr in nbr_1 + nbr_2 + nbr_3:
                nbr.SetAtomMapNum(idx)
        elif len(atom_list) == 4:
            atm_1, atm_2, atm_3, atm_4 = atom_list
            nbr_1 = [x.GetOtherAtom(atm_1) for x in atm_1.GetBonds()]
            nbr_2 = [x.GetOtherAtom(atm_2) for x in atm_2.GetBonds()]
            nbr_3 = [x.GetOtherAtom(atm_3) for x in atm_3.GetBonds()]
            nbr_4 = [x.GetOtherAtom(atm_4) for x in atm_4.GetBonds()]
            for nbr in nbr_1 + nbr_2 + nbr_3 + nbr_4:
                nbr.SetAtomMapNum(idx)
        chiral_tag_dict[idx] = Chem.rdchem.ChiralType.CHI_UNSPECIFIED
        source_bond_order_dict[idx] = None
    
    # remove the dummy atoms
    new_core = Chem.DeleteSubstructs(core, Chem.MolFromSmarts('[#0]'))

    # get the new atoms with AtomMapNum, these will be connected
    bond_join_dict_core = defaultdict(list)
    for atom in new_core.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num > 0:
            bond_join_dict_core[map_num].append(atom.GetIdx())

    # remove the dummy atoms
    new_rgroups = Chem.DeleteSubstructs(rgroups, Chem.MolFromSmarts('[#0]'))

    # get the new atoms with AtomMapNum, these will be connected
    bond_join_dict_rgroups = defaultdict(list)
    for atom in new_rgroups.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num > 0:
            bond_join_dict_rgroups[map_num].append(atom.GetIdx())

    # make an editable molecule and add bonds between atoms with correspoing AtomMapNum
    em = EditableMol(new_rgroups)
    for idx in bond_join_dict_core.keys():
        core_atms = bond_join_dict_core[idx]
        rgroup_atms = bond_join_dict_rgroups[idx]
        i = 0
        for core_atm in core_atms:
            for rgroup_atm in rgroup_atms:
                em.AddBond(core_atm, rgroup_atm, order=bond_order_dict[idx][i])
                i+=1

    combined_mol = em.GetMol()

    # postprocessing to fix number of hydrogens at attachment points
    for atom in combined_mol.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num == 0: continue
        if atom.GetAtomicNum(
        ) == 7 and not atom.IsInRing():  # for nitrogen atoms
            atom.SetNoImplicit(True)
            nbrs = list(atom.GetNeighbors())
            nonHs = [nbr.GetAtomicNum() != 1 for nbr in nbrs]
            bonds = list(atom.GetBonds())
            bondtypes = [bond.GetBondType() for bond in bonds]
            i = 0
            for bondtype in bondtypes:
                if bondtype == Chem.BondType.DOUBLE:
                    i += 1
                elif bondtype == Chem.BondType.TRIPLE:
                    i += 2
            numHs = 3 - len(nonHs) - i + atom.GetFormalCharge()
            if numHs < 0:
                return None
            atom.SetNumExplicitHs(numHs)
        if atom.GetAtomicNum(
        ) == 6 and not atom.IsInRing():  # for carbon atoms
            atom.SetNoImplicit(True)
            nbrs = list(atom.GetNeighbors())
            nonHs = [nbr.GetAtomicNum() != 1 for nbr in nbrs]
            bonds = list(atom.GetBonds())
            bondtypes = [bond.GetBondType() for bond in bonds]
            i = 0
            for bondtype in bondtypes:
                if bondtype == Chem.BondType.DOUBLE:
                    i += 1
                elif bondtype == Chem.BondType.TRIPLE:
                    i += 2
            numHs = 4 - len(nonHs) - i + atom.GetFormalCharge()
            if numHs < 0:
                return None
            atom.SetNumExplicitHs(numHs)
        if atom.GetAtomicNum(
        ) == 6 and atom.IsInRing():  # for carbon atoms
            atom.SetNoImplicit(True)
            nbrs = list(atom.GetNeighbors())
            nonHs = [nbr.GetAtomicNum() != 1 for nbr in nbrs]
            bonds = list(atom.GetBonds())
            bondtypes = [bond.GetBondType() for bond in bonds]
            i = 0
            for bondtype in bondtypes:
                if bondtype == Chem.BondType.DOUBLE:
                    i += 1
                elif bondtype == Chem.BondType.TRIPLE:
                    i += 2
            numHs = 4 - len(nonHs) - i + atom.GetFormalCharge()
            if numHs < 0:
                return None
            atom.SetNumExplicitHs(numHs)
    combined_mol_fixed = combined_mol

    # if molecule is invalid try replacing single bond tokens
    if Chem.MolFromSmiles(Chem.MolToSmiles(combined_mol)) is None:
        combined_mol_fixed = Chem.MolFromSmiles(
            Chem.MolToSmiles(combined_mol).replace('-', ''))
        if combined_mol_fixed is None:
            print(f'invalid molecule: {Chem.MolToSmiles(combined_mol)}')
        else:
            combined_mol = combined_mol_fixed

    # remove the AtomMapNum values
    for atom in combined_mol.GetAtoms():
        atom.SetAtomMapNum(0)

    # remove explicit Hs
    try:
        final_mol = Chem.RemoveHs(combined_mol)
    except Chem.rdchem.AtomValenceException:
        final_mol = Chem.RemoveHs(combined_mol, sanitize=False)
        print(final_mol)
    except Chem.rdchem.KekulizeException:
        return None

    #TODO maybe do somewhere else
    # restore bonds to aromatic type
    Chem.SanitizeMol(final_mol)

    return final_mol