import os,glob,json
from Bio.PDB import PDBParser

MODIFIED_RESIDUES = [
    'CSD', 'CME', 'CSO', 'OCS', 'CSS', 'CSX', 'CYG', 
    'KCX', 'LLP', 'PCA', 'SEP', 'TPO', 'PTR', 'MSE'  
]

def update_flexible_residues_and_fix_gaps():
    print("--- Starting Protein Flexibility Updater & Gap Fixer ---")
    
    out_dir = "receptor_flex_updated"
    os.makedirs(out_dir, exist_ok=True)
    
    json_path = "configs/active_chains.json"
    if not os.path.exists(json_path):
        print(f" ❌ Error: {json_path} not found!")
        return
        
    with open(json_path, 'r') as f:
        active_chains_map = json.load(f)

    best_candidates = glob.glob("best_pose_candidat/*_best.mol2")
    if not best_candidates:
        print(" ⚠️ No validated candidates found.")
        return
    
    for cand_file in best_candidates:
        base_name = os.path.basename(cand_file).replace("_best.mol2", "")
        parts = base_name.split("_")
        pdb_id = parts[0]
        ligand_name = "_".join(parts[1:])
        
        dock_file = f"docking_results/{base_name}_out.pdbqt"
        rigid_pdb = f"clean_proteins_pdb/{pdb_id}_clean.pdb"
        
        if not os.path.exists(dock_file) or not os.path.exists(rigid_pdb):
            print(f" ❌ Missing files for {base_name}, skipping...")
            continue
            
        print(f"\nProcessing VIP Target: {pdb_id} | Ligand: {ligand_name}")

        target_info = active_chains_map.get(pdb_id, {})
        raw_chains = target_info.get("Interacting_Chains", "A")
        allowed_chains = {c.strip() for c in raw_chains.split(",")}
        print(f"✅ Keeping Active Chain(s) Only: {allowed_chains}")

        flex_coords = {}
        in_model_1 = False
        in_flex_res = False
        
        with open(dock_file, 'r') as f:
            for line in f:
                if line.startswith("MODEL 1"): in_model_1 = True
                elif line.startswith("ENDMDL") and in_model_1: break
                
                if in_model_1:
                    if line.startswith("BEGIN_RES"):
                        in_flex_res = True
                        parts_line = line.split()
                        current_chain = parts_line[2]
                        current_resnum = int(parts_line[3])
                    elif line.startswith("END_RES"):
                        in_flex_res = False
                    elif in_flex_res and line.startswith(("ATOM", "HETATM")):
                        atom_name = line[12:16].strip()
                        if atom_name.startswith('H'): continue
                        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                        flex_coords[(current_chain, current_resnum, atom_name)] = (x, y, z)
        
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", rigid_pdb)
        output_pdb = os.path.join(out_dir, f"{pdb_id}_{ligand_name}_rec.pdb")
        
        updated_atoms = 0

        atom_serial_counter = 1 
        
        with open(output_pdb, 'w') as out_f:
            for model in structure:
                for chain in model:
                    if chain.id not in allowed_chains:
                        continue 
                        
                    prev_res_num = None
                    
                    for res in chain:
                        res_name = res.get_resname().strip().upper()
                        
                        if res.id[0] != ' ' and res_name not in MODIFIED_RESIDUES:
                            continue
                            
                        curr_res_num = res.id[1]
                        
                        if prev_res_num is not None and (curr_res_num - prev_res_num > 1):
                            out_f.write("TER\n")
                            print(f"✅ Inserted TER gap-break in Chain {chain.id}: between {prev_res_num} and {curr_res_num}")
                            
                        prev_res_num = curr_res_num
                        
                        for atom in res:
                            key = (chain.id, curr_res_num, atom.name.strip())
                            if key in flex_coords:
                                atom.set_coord(flex_coords[key])
                                updated_atoms += 1
                            
                            atom_name = atom.name
                            atom_fmt = f" {atom_name:<3}" if len(atom_name) < 4 else f"{atom_name:<4}"
                            x, y, z = atom.coord
                            element = atom.element if atom.element else atom_name[0]
                            
                            record_type = "ATOM  "
                            
                            line = (f"{record_type}{atom_serial_counter:>5} {atom_fmt} {res.resname:>3} "
                                    f"{chain.id}{curr_res_num:>4}    "
                                    f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
                                    f"{atom.occupancy:>6.2f}{atom.bfactor:>6.2f}          "
                                    f"{element:>2}\n")
                            out_f.write(line)
                            atom_serial_counter += 1
                            
                    out_f.write("TER\n") 
            out_f.write("END\n")
            
        print(f"✅ Updated {updated_atoms} flex atoms & fixed gaps.")
        print(f"✅ Saved to: {output_pdb}")

if __name__ == "__main__":
    update_flexible_residues_and_fix_gaps()