from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem.MolStandardize import rdMolStandardize

enumerator = rdMolStandardize.TautomerEnumerator()


def get_product(reacts):
    # two step reaction, first nucleophilic addition by NCN to the aldehyde
    rxn1 = rdChemReactions.ReactionFromSmarts(
        '[#6:1](=[O]).[#7X3H1:6]~[#6X3:7](~[#7:8])~[#7X3H2:9]>>[#6:1][#7:8][#6X3:7](=[#7X3H0:6])~[#7X3H2:9]'
    )
    products = rxn1.RunReactants(reacts[:2])
    product1 = Chem.Mol(products[0][0])
    final = None
    # enumerate tautomers
    for taut in enumerator.Enumerate(product1):
        # second reaction, closing of ring by nucleophilic attack by primary amine onto carbonyl group
        rxn2 = rdChemReactions.ReactionFromSmarts(
            '[#7X3H2:1][#6X3:6]~[#7H0:5][#6:2].[#6][#6](=[#8])[#6][#6:3](=[#8])[*:4]>>[#6:2]1[#6]([#6:3](=[#8])[*:4])=[#6]([#6])[#7:1][#6X3:6]~[#7X3:5]1'
        )
        products2 = rxn2.RunReactants((taut, reacts[2]))
        # filter out reactions that do not yield a result or cannot be kekulized
        try:
            Chem.SanitizeMol(products2[0][0])
            print(products2[0][0])
            final = products2[0][0]
        except:
            pass

    return final

if __name__ == "__main__":
    reacts = (Chem.MolFromSmiles('O=Cc1ccco1'),Chem.MolFromSmiles('Nc1nnn[nH]1'),Chem.MolFromSmiles('CCOC(=O)CC(C)=O'))
    product = get_product(reacts)
    Chem.Kekulize(product, clearAromaticFlags=True)
    print(Chem.MolToSmiles(product))

    reacts = (Chem.MolFromSmiles('[H]C(C3=CC=CO3)=O'),Chem.MolFromSmiles('N=C(N)NC1=NC(C=CC=C2)=C2N1'),Chem.MolFromSmiles('O=C(C)CC(OCC)=O'))
    product = get_product(reacts)
    Chem.Kekulize(product, clearAromaticFlags=True)
    print(Chem.MolToSmiles(product))
