from rdkit import Chem

# Predefined dictionary of common covalent warheads and their SMARTS patterns
# Based on typical covalent docking reaction types in Schrödinger
WARHEAD_PATTERNS = {
    "Michael Addition": {
        "Michael Acceptor": "[C:1]=[C:2]-[C:3]=[O,S:4]",
        "Acrylate": "[C;H2]=C-C(=O)O",
        "Vinyl Sulfone": "C=CS(=O)=O",
        "Vinyl Sulfonamide": "C=CS(=O)(=O)N",
        "Propiolamide": "C#CC(=O)N"
    },
    "Nucleophilic Substitution": {
        "Haloacetamide": "[Cl,Br,I]CC(=O)N",
        "Halomethyl Ketone": "[Cl,Br,I]CC(=O)",
        "Epoxide": "C1OC1",
        "Azetidine": "C1CN1",
        "Sulfonate Ester": "COS(=O)(=O)C"
    },
    "Covalent Reversible": {
        "Nitrile": "C#N",
        "Boronic Acid": "B(O)O",
        "Aldehyde": "[C;H1]=O"
    }
}

class WarheadDetector:
    def __init__(self):
        self.patterns = {}
        for category, sub_patterns in WARHEAD_PATTERNS.items():
            for name, smarts in sub_patterns.items():
                mol_pattern = Chem.MolFromSmarts(smarts)
                if mol_pattern:
                    self.patterns[name] = {
                        "mol": mol_pattern,
                        "category": category,
                        "smarts": smarts
                    }

    def detect(self, smiles_list):
        """
        Scans a list of SMILES for warheads.
        Returns a summary dictionary.
        """
        results = {}
        detected_categories = set()
        
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(str(smiles))
            if not mol:
                continue
            
            for name, data in self.patterns.items():
                if mol.HasSubstructMatch(data["mol"]):
                    results[name] = results.get(name, 0) + 1
                    detected_categories.add(str(data["category"]))
        
        return results, detected_categories

    def print_summary(self, results, detected_categories):
        """Prints a user-friendly summary of detected warheads."""
        print("\n" + "="*40)
        print("COVALENT WARHEAD DETECTION SUMMARY")
        print("="*40)
        
        if not results:
            print("No known covalent warheads detected.")
            return

        for category, sub_patterns in WARHEAD_PATTERNS.items():
            found_names = [name for name in sub_patterns.keys() if name in results]
            if found_names:
                print(f"\n[{category}]")
                for name in found_names:
                    count = results[name]
                    print(f"  - {name}: {count} ligands")
        
        print("\n" + "="*40)
