import os,glob,json,subprocess
import numpy as np
from Bio.PDB import PDBParser
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import PDBQTMolecule, RDKitMolCreate 

AFFINITY_THRESHOLD = -6.0
HBOND_DISTANCE_CUTOFF = 3.6  
MAX_STRAIN_PER_ATOM = 1.5  
POLAR_AD_TYPES = ['N', 'NA', 'O', 'OA', 'S', 'SA', 'F', 'CL', 'BR', 'I']

def extract_best_pose(pdbqt_file):
    affinity = None
    polar_coords = []
    heavy_atom_count = 0
    
    with open(pdbqt_file, 'r') as f:
        in_model_1 = False
        for line in f:
            if line.startswith("MODEL 1"):
                in_model_1 = True
            elif line.startswith("ENDMDL") and in_model_1:
                break
                
            elif line.startswith("BEGIN_RES"):
                break
            
            if in_model_1:
                if line.startswith("REMARK VINA RESULT:"):
                    affinity = float(line.split()[3])
                elif line.startswith(("ATOM", "HETATM")):
                    ad_type = line[77:].strip().upper()
                    
                    if ad_type not in ['H', 'HD']:
                        heavy_atom_count += 1
                        
                    if ad_type in POLAR_AD_TYPES:
                        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                        polar_coords.append(np.array([x, y, z]))
                        
    return affinity, polar_coords, heavy_atom_count

def get_protein_polar_atoms(pdb_file, target_residues_list):
    parser = PDBParser(QUIET=True)
    model = parser.get_structure("protein", pdb_file)[0]
    target_coords = []
    for res_info in target_residues_list:
        parts = res_info.split('_')
        res_num, chain_id = int(parts[1]), parts[2]
        try:
            for atom in model[chain_id][(' ', res_num, ' ')].get_atoms():
                if atom.element.upper() in ['N', 'O', 'S', 'F', 'CL']:
                    target_coords.append(atom.coord)
        except KeyError:
            pass
    return target_coords

def check_interactions(ligand_coords, protein_coords, cutoff):
    if not ligand_coords or not protein_coords: 
        return False, None
    
    min_dist = float('inf')
    for l_atom in ligand_coords:
        for p_atom in protein_coords:
            dist = np.linalg.norm(l_atom - p_atom)
            if dist < min_dist:
                min_dist = dist
                
    return (min_dist <= cutoff), min_dist

def calculate_ligand_strain(docked_pdbqt):
    try:
        pdbqt_mol = PDBQTMolecule.from_file(docked_pdbqt)
        
        rdkit_mols = None
        if hasattr(RDKitMolCreate, 'from_pdbqt_mol'):
            rdkit_mols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
        elif hasattr(RDKitMolCreate, 'from_pdbqt'):
            rdkit_mols = RDKitMolCreate.from_pdbqt(pdbqt_mol)
        else:
            creator = RDKitMolCreate()
            rdkit_mols = creator.from_pdbqt_mol(pdbqt_mol) if hasattr(creator, 'from_pdbqt_mol') else creator.from_pdbqt(pdbqt_mol)
            
        if not rdkit_mols: 
            return None
            
        rdkit_mol = rdkit_mols[0]
        Chem.SanitizeMol(rdkit_mol)
        mp = AllChem.MMFFGetMoleculeProperties(rdkit_mol, mmffVariant='MMFF94')
        ff = AllChem.MMFFGetMoleculeForceField(rdkit_mol, mp)
        
        if ff is None: return None
            
        e_docked = ff.CalcEnergy()
        ff.Minimize(maxIts=500)
        e_minimized = ff.CalcEnergy()
        
        return e_docked - e_minimized
    except Exception as e:
        print(f" ⚠️ Strain calculation failed: {e}")
        return None
def extract_and_convert_to_mol2(pdbqt_file, output_mol2_path):
    try:
        pdbqt_mol = PDBQTMolecule.from_file(pdbqt_file)
        
        rdkit_mols = None
        if hasattr(RDKitMolCreate, 'from_pdbqt_mol'):
            rdkit_mols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
        elif hasattr(RDKitMolCreate, 'from_pdbqt'):
            rdkit_mols = RDKitMolCreate.from_pdbqt(pdbqt_mol)
        else:
            creator = RDKitMolCreate()
            rdkit_mols = creator.from_pdbqt_mol(pdbqt_mol) if hasattr(creator, 'from_pdbqt_mol') else creator.from_pdbqt(pdbqt_mol)

        if not rdkit_mols:
            print(" ⚠️ Failed to parse ligand with Meeko.")
            return

        rdkit_mol = rdkit_mols[0]
        Chem.SanitizeMol(rdkit_mol)
        
        rdkit_mol = Chem.AddHs(rdkit_mol, addCoords=True)
        
        temp_sdf = pdbqt_file.replace('.pdbqt', '_temp.sdf')
        writer = Chem.SDWriter(temp_sdf)
        writer.write(rdkit_mol)
        writer.close()
        
        subprocess.run(
            ['obabel', '-i', 'sdf', temp_sdf, '-o', 'mol2', '-O', output_mol2_path],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
        
    except Exception as e:
        print(f" ⚠️ Error in extracting and converting pose: {e}")
    finally:
        if 'temp_sdf' in locals() and os.path.exists(temp_sdf):
            os.remove(temp_sdf)

def main():
    print("--- Starting Post-Docking Biophysical Filtering ---")

    os.makedirs("best_pose_candidat", exist_ok=True)

    with open("configs/active_chains.json", 'r') as f:
        active_chains_map = json.load(f)
        
    result_files = glob.glob("docking_results/*_out.pdbqt")
    passed_candidates = []
    
    for res_file in result_files:
        parts = os.path.basename(res_file).replace("_out.pdbqt", "").split("_")
        pdb_id, ligand_name = parts[0], "_".join(parts[1:])
        
        print(f"\n Analyzing: Target {pdb_id} | Ligand {ligand_name}")
        
        affinity, lig_polar_coords, heavy_atoms = extract_best_pose(res_file)
        if affinity is None or affinity > AFFINITY_THRESHOLD:
            print(f" ❌ REJECTED: Affinity too weak ({affinity}).")
            continue
            
        target_residues = [r.strip() for r in active_chains_map.get(pdb_id, {}).get("Flexible_Residues", "").split(",") if r]
        clean_pdb = f"clean_proteins_pdb/{pdb_id}_clean.pdb"
        
        prot_polar_coords = get_protein_polar_atoms(clean_pdb, target_residues)
        has_hbond, min_dist = check_interactions(lig_polar_coords, prot_polar_coords, HBOND_DISTANCE_CUTOFF)
        
        if not has_hbond:
            if min_dist:
                print(f" ❌ REJECTED: No close H-Bond (Closest distance: {min_dist:.2f} Å)")
            else:
                print(" ❌ REJECTED: No polar atoms available.")
            continue
            
        print(f"✅ Passed H-Bond test (Distance: {min_dist:.2f} Å). Checking Strain...")
        
        strain = calculate_ligand_strain(res_file)
        
        if strain is not None and heavy_atoms > 0:
            per_atom_strain = strain / heavy_atoms
            if per_atom_strain > MAX_STRAIN_PER_ATOM:
                print(f" ❌ REJECTED: Ligand is highly strained! (Per Atom: {per_atom_strain:.2f} > Limit: {MAX_STRAIN_PER_ATOM}) Total: {strain:.2f}")
                continue
        elif strain is not None and heavy_atoms == 0:
            continue
            
        if heavy_atoms > 0:
            sile = abs(affinity) / (heavy_atoms ** 0.3)
        else:
            sile = 0.0
            
        if strain is not None:
            print(f"✅ ACCEPTED: Low strain (Per Atom: {per_atom_strain:.2f}). Total: {strain:.2f}. SILE: {sile:.2f}")
        else:
            print(f"✅ ACCEPTED: Strain Unknown. SILE: {sile:.2f}")
            
        passed_candidates.append({
            "ligand": ligand_name,
            "target": pdb_id,
            "affinity": affinity,
            "strain": round(strain, 2) if strain is not None else "Unknown",
            "sile": round(sile, 2),
            "heavy_atoms": heavy_atoms
        })
        mol2_filename = f"best_pose_candidat/{pdb_id}_{ligand_name}_best.mol2"
        extract_and_convert_to_mol2(res_file, mol2_filename)
        print(f"✅ Best pose saved to: {mol2_filename}")

    passed_candidates.sort(key=lambda x: x["sile"], reverse=True)

    print("\n=========================================================================")
    print(f"--- Top Validated Candidates (Ranked by SILE) - {len(passed_candidates)} found ---")
    print(f"{'Rank':<5} | {'Ligand':<20} | {'Affinity':<10} | {'Strain':<8} | {'Heavy Atoms':<12} | {'SILE':<6}")
    print("-" * 73)
    
    for idx, cand in enumerate(passed_candidates, 1):
        print(f"{idx:<5} | {cand['ligand']:<20} | {cand['affinity']:<10} | {cand['strain']:<8} | {cand['heavy_atoms']:<12} | {cand['sile']:<6}")
    print("=========================================================================")

if __name__ == "__main__":
    main()