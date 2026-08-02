import os,sys,glob,json,subprocess
import numpy as np
from Bio.PDB import PDBList, PDBParser, PDBIO, Select
from Bio.PDB.NeighborSearch import NeighborSearch

ligand_files = glob.glob("ligand/*.pdbqt")
if not ligand_files:
    print(" ❌ Error: No ligand .pdbqt file found in 'ligand/' directory.")
    sys.exit(1)
docking_ligand = os.path.basename(ligand_files[0]).replace(".pdbqt", "")

pdb_ids = ['8PQ9'] #<================================================================================================================= PDB ID

_BASE = os.path.join(os.path.dirname(__file__), "ADFR")
ADFR_BIN = os.path.join(_BASE, "bin")
ADFR_PYTHON = os.path.join(ADFR_BIN, "python")
PREPARE_RECEPTOR = os.path.join(ADFR_BIN, "prepare_receptor")


def find_prepare_flex_script():
    search_roots = []

    conda_prefix = os.environ.get("CONDA_PREFIX") or os.environ.get("MAMBA_ROOT_PREFIX")
    if conda_prefix:
        search_roots.append(conda_prefix)

    search_roots.append(sys.prefix)

    search_roots.append(_BASE)

    seen = set()
    for root in search_roots:
        if not root or root in seen or not os.path.isdir(root):
            continue
        seen.add(root)
        matches = glob.glob(os.path.join(root, "**", "prepare_flexreceptor4.py"), recursive=True)
        if matches:
            return matches[0]

    return None


PREPARE_FLEX_SCRIPT = find_prepare_flex_script()

STANDARD_RESIDUES = [
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS',
    'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO',
    'SER', 'THR', 'TRP', 'TYR', 'VAL'
]

MODIFIED_RESIDUES = ['CSD', 'CME', 'CSO', 'OCS', 'CSS', 'CSX', 'CYG', 'KCX', 'LLP', 'PCA', 'SEP', 'TPO', 'PTR', 'MSE'  ]

class RigidProteinSelect(Select):
    def __init__(self, keep_chains, flex_residues, ligand_chain_id=None):
        self.keep_chains = keep_chains
        self.flex_residues = flex_residues
        self.ligand_chain_id = ligand_chain_id 

    def accept_chain(self, chain):
        if self.ligand_chain_id and chain.get_id() == self.ligand_chain_id:
            return False
        return True

    def accept_residue(self, residue):
        res_name = residue.get_resname().strip().upper()
        return res_name in STANDARD_RESIDUES or res_name in MODIFIED_RESIDUES
    
    def accept_atom(self, atom):
        if atom.get_altloc() not in [' ', 'A']:
            return False
        atom.set_altloc(' ')
        return True


def prepare_receptor_adfr(rigid_pdb, rigid_pdbqt, flex_pdbqt, flex_res_objs):
    env = os.environ.copy()
    env["PATH"] = ADFR_BIN + os.pathsep + env.get("PATH", "")
    env["LD_LIBRARY_PATH"] = os.path.join(_BASE, "lib") + os.pathsep + env.get("LD_LIBRARY_PATH", "")

    if not PREPARE_FLEX_SCRIPT or not os.path.isfile(PREPARE_FLEX_SCRIPT):
        print("❌ Error: prepare_flexreceptor4.py could not be located automatically.")
        print("❌ Make sure the docking_env environment is active, then locate it manually with:")
        print('❌ find "$CONDA_PREFIX" -name "prepare_flexreceptor4.py" 2>/dev/null')
        print("❌ Once you have the path, hardcode it into PREPARE_FLEX_SCRIPT at the top of this file.")
        return False

    try:
        print("Running ADFR prepare_receptor...")
        subprocess.run(
            [PREPARE_RECEPTOR, "-r", rigid_pdb, "-o", rigid_pdbqt, "-A", "hydrogens", "-U", "waters"],
            check=True, env=env, stdout=subprocess.DEVNULL
        )   
        
        rigid_mol_name = os.path.basename(rigid_pdbqt).replace(".pdbqt", "")

        chain_groups = {}
        for r in flex_res_objs:
            chain_id = r.parent.id
            res_name_num = f"{r.get_resname().strip()}{r.id[1]}"
            chain_groups.setdefault(chain_id, []).append(res_name_num)

        flex_str = ",".join(
            f"{rigid_mol_name}:{cid}:{'_'.join(res_list)}"
            for cid, res_list in chain_groups.items()
        )

        print(f"DEBUG: Running ADFR prepare_flexreceptor for: {flex_str}")
        subprocess.run(
            [ADFR_PYTHON, PREPARE_FLEX_SCRIPT, "-r", rigid_pdbqt, "-s", flex_str,
             "-g", rigid_pdbqt, "-x", flex_pdbqt],
            check=True, env=env
        )
        remove_orphan_ha_and_fix_anchors(rigid_pdbqt, flex_pdbqt)
        return True

    except subprocess.CalledProcessError as e:
        print(f"❌ ADFR preparation failed for {rigid_pdb}: {e}")
        return False

def remove_orphan_ha_and_fix_anchors(rigid_pdbqt, flex_pdbqt):
    flex_ca_residues = set()
    root_atoms = []
    
    with open(flex_pdbqt, 'r') as f:
        in_root = False
        for line in f:
            if line.startswith("ROOT"):
                in_root = True
            elif line.startswith("ENDROOT"):
                in_root = False
            elif in_root and line.startswith(('ATOM', 'HETATM')):
                root_atoms.append(line)
                if line[12:16].strip() == 'CA':
                    flex_ca_residues.add((line[21], line[22:26].strip()))

    kept = []
    with open(rigid_pdbqt, 'r') as f:
        for line in f:
            if (line.startswith(('ATOM', 'HETATM'))
                    and line[12:16].strip() == 'HA'
                    and (line[21], line[22:26].strip()) in flex_ca_residues):
                continue
            kept.append(line)

    with open(rigid_pdbqt, 'w') as f:
        f.writelines(kept)
        f.writelines(root_atoms)

def process_targets(pdb_list, docking_ligand_name, cutoff_distance=5.0, padding=8.0, n_flex_residues=5):
    pdbl = PDBList()
    parser = PDBParser(QUIET=True)
    io = PDBIO()

    os.makedirs("raw_pdb", exist_ok=True)
    os.makedirs("clean_proteins_pdb", exist_ok=True)
    os.makedirs("clean_proteins", exist_ok=True)
    os.makedirs("configs", exist_ok=True)

    active_chains_map = {}

    for pdb_id in pdb_list:
        print(f"\n--- Processing {pdb_id} ---")
        try:
            pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir='raw_pdb', file_format='pdb')
            structure = parser.get_structure(pdb_id, pdb_file)
        except Exception as e:
            print(f"❌ Could not download/parse {pdb_id}: {e}. Skipping.")
            sys.exit(1)
        model = structure[0]
        ligand_atoms = []
        ligand_name = "Unknown"
        ligand_chain_id_for_removal = None
        peptide_chain = None

        for chain in model:
            res_count = sum(1 for r in chain if r.get_resname().strip().upper() in STANDARD_RESIDUES or r.get_resname().strip().upper() in MODIFIED_RESIDUES)
            
            if 3 <= res_count <= 35: 
                print(f"✅ Peptide ligand detected! Chain '{chain.id}' with {res_count} residues.")
                peptide_chain = chain
                break 

        if peptide_chain:
            ligand_atoms = list(peptide_chain.get_atoms())
            ligand_name = f"Peptide_Chain_{peptide_chain.id}"
            ligand_chain_id_for_removal = peptide_chain.id
            
        else:
            largest_ligand = None
            max_atoms = 0
            for residue in model.get_residues():
                res_name = residue.get_resname().strip().upper()
                if res_name not in STANDARD_RESIDUES and res_name not in MODIFIED_RESIDUES and res_name not in ['HOH', 'WAT']:
                    atom_count = len(list(residue.get_atoms()))
                    if atom_count > max_atoms and atom_count > 5:
                        max_atoms = atom_count
                        largest_ligand = residue

            if largest_ligand:
                ligand_name = largest_ligand.get_resname().strip()
                ligand_atoms = list(largest_ligand.get_atoms())
                print(f"✅ Small molecule ligand detected: {ligand_name} with {max_atoms} atoms.")
            else:
                print(f" ❌ No valid ligand (peptide or small molecule) found in {pdb_id}. Skipping.")
                continue

        coords = np.array([atom.coord for atom in ligand_atoms])

        all_atoms = list(model.get_atoms())
        ns = NeighborSearch(all_atoms)
        interacting_chains = set()
        residue_min_dist = {}

        for lig_atom in ligand_atoms:
            neighbors = ns.search(lig_atom.coord, cutoff_distance, level='R')
            for res in neighbors:
                if res.get_resname().strip().upper() not in STANDARD_RESIDUES:
                    continue
                
                if ligand_chain_id_for_removal and res.parent.id == ligand_chain_id_for_removal:
                    continue

                interacting_chains.add(res.parent.get_id())
                dist = min(
                    np.linalg.norm(atom.coord - lig_atom.coord) for atom in res.get_atoms()
                )
                if res not in residue_min_dist or dist < residue_min_dist[res]:
                    residue_min_dist[res] = dist

        interacting_residues = sorted(residue_min_dist, key=residue_min_dist.get)

        flex_residues = interacting_residues[:n_flex_residues]
        flex_coords = np.array([a.coord for r in flex_residues for a in r.get_atoms()])

        all_coords = np.vstack([coords, flex_coords])
        box_min = all_coords.min(axis=0)
        box_max = all_coords.max(axis=0)
        box_size = (box_max - box_min) + padding
        
        center = (box_max + box_min) / 2.0
        
        print(f"Grid box center: X={center[0]:.3f}, Y={center[1]:.3f}, Z={center[2]:.3f}")

        clean_pdb = f"clean_proteins_pdb/{pdb_id}_clean.pdb"

        io.set_structure(structure)
        io.save(clean_pdb, RigidProteinSelect(interacting_chains, flex_residues, ligand_chain_id_for_removal))

        rigid_pdbqt = f"clean_proteins/{pdb_id}_rigid.pdbqt"
        flex_pdbqt = f"clean_proteins/{pdb_id}_flex.pdbqt"

        success = prepare_receptor_adfr(clean_pdb, rigid_pdbqt, flex_pdbqt, flex_residues)
        if not success:
            print(f" ⚠️ Skipping {pdb_id} due to ADFR preparation failure.")
            continue
        config_path = f"configs/{pdb_id}_config.txt"
        with open(config_path, "w") as f:
            f.write(f"receptor = {rigid_pdbqt}\n")
            f.write(f"flex = {flex_pdbqt}\n")
            f.write(f"ligand = ligand/{docking_ligand_name}.pdbqt\n\n")
            f.write(f"center_x = {center[0]:.3f}\n")
            f.write(f"center_y = {center[1]:.3f}\n")
            f.write(f"center_z = {center[2]:.3f}\n\n")
            f.write(f"size_x = {box_size[0]:.3f}\n")
            f.write(f"size_y = {box_size[1]:.3f}\n")
            f.write(f"size_z = {box_size[2]:.3f}\n\n")
            f.write("cpu = 10\n")
            f.write("exhaustiveness = 64\n")
            f.write("num_modes = 9\n")
            f.write("energy_range = 3\n")
            f.write("seed = 12345\n")

        chains_string = ", ".join(sorted(interacting_chains))
        flex_res_names = [f"{r.get_resname()}_{r.id[1]}_{r.parent.id}" for r in flex_residues]

        active_chains_map[pdb_id] = {
            "Interacting_Chains": chains_string,
            "Flexible_Residues": ", ".join(flex_res_names)
        }
        print(f"✅ {pdb_id} done. Flex residues: {flex_res_names}")

    json_path = "configs/active_chains.json"
    with open(json_path, "w") as f:
        json.dump(active_chains_map, f, indent=4)
    print(f"\n✅ All active chains and flex residues exported to {json_path}")


if __name__ == "__main__":
    process_targets(pdb_ids, docking_ligand)