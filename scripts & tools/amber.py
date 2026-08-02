import os,glob,subprocess

_AMBER_BIN = "/home/angelommd/miniconda3/envs/amber_env/bin"
AMBER_ENV = {
    **os.environ,
    "AMBERHOME": "/home/angelommd/miniconda3/envs/amber_env",
    "PATH": f"{_AMBER_BIN}:{os.environ.get('PATH', '')}",
}

def sanitize_receptor_for_amber(input_pdb, output_pdb):
    res_map = {
        'CSD': 'CYS', 'CME': 'CYS', 'CSO': 'CYS', 'OCS': 'CYS', 'CSS': 'CYS', 'CSX': 'CYS', 'CYG': 'CYS',
        'KCX': 'LYS', 'LLP': 'LYS', 'PCA': 'GLU', 'SEP': 'SER', 'TPO': 'THR', 'PTR': 'TYR', 'MSE': 'MET'
    }
    std_atoms = {
        'CYS': {'N', 'CA', 'C', 'O', 'CB', 'SG'},
        'LYS': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ'},
        'GLU': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2'},
        'SER': {'N', 'CA', 'C', 'O', 'CB', 'OG'},
        'THR': {'N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2'},
        'TYR': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'},
        'MET': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'SD', 'CE'}
    }

    with open(input_pdb, 'r') as f_in, open(output_pdb, 'w') as f_out:
        for line in f_in:
            if line.startswith(("ATOM", "HETATM")):
                res_name = line[17:20].strip()
                atom_name = line[12:16].strip()
                
                if res_name in res_map:
                    std_res = res_map[res_name]
                    if atom_name not in std_atoms[std_res]:
                        continue
                    line = line[:17] + f"{std_res:>3}" + line[20:]
                    if line.startswith("HETATM"):
                        line = "ATOM  " + line[6:]
            f_out.write(line)

def prepare_mmgbsa_topologies():
    print("--- Starting AmberTools Topology Builder  ---")
    
    ligand_dir = "best_pose_candidat"
    receptor_dir = "receptor_flex_updated"
    out_dir = "gbsa_topologies"
    os.makedirs(out_dir, exist_ok=True)
    
    best_ligands = glob.glob(f"{ligand_dir}/*_best.mol2")
    if not best_ligands:
        print(f" ⚠️ No ligand files found in '{ligand_dir}'.")
        return

    junk_files = ["sqm.in", "sqm.out", "sqm.pdb", "ATOMTYPE.INF", "leap.log"]

    for lig_file in best_ligands:
        base_name = os.path.basename(lig_file).replace("_best.mol2", "")
        rec_file = os.path.join(receptor_dir, f"{base_name}_rec.pdb")
        
        if not os.path.exists(rec_file):
            print(f" ❌ Missing receptor file for {base_name}. Skipping!")
            continue
            
        print(f"\n Processing Complex: {base_name}")
        complex_out_dir = os.path.join(out_dir, base_name)
        os.makedirs(complex_out_dir, exist_ok=True)
        
        expected_files = [
            "comp.prmtop", "comp.inpcrd", 
            "rec.prmtop", "rec.inpcrd", 
            "lig.prmtop", "lig.inpcrd"
        ]
        
        already_processed = True
        for ef in expected_files:
            if not os.path.exists(os.path.join(complex_out_dir, ef)):
                already_processed = False
                break
                
        if already_processed:
            print(f"✅ Topologies for {base_name} already exist. Skipping to next...")
            continue

        lig_bcc_mol2 = os.path.join(complex_out_dir, "lig_bcc.mol2")
        lig_frcmod = os.path.join(complex_out_dir, "lig.frcmod")
        tleap_in = os.path.join(complex_out_dir, "tleap.in")
        sanitized_rec = os.path.join(complex_out_dir, f"{base_name}_sanitized_rec.pdb")
        sanitize_receptor_for_amber(rec_file, sanitized_rec)

        try:
            print(" Calculating AM1-BCC charges (antechamber)...")
            subprocess.run([
                f"{_AMBER_BIN}/antechamber", "-i", lig_file, "-fi", "mol2",
                "-o", lig_bcc_mol2, "-fo", "mol2", 
                "-c", "bcc", "-s", "2", "-nc", "0"
            ], check=True, capture_output=True, text=True, env=AMBER_ENV)
            
            print("     Generating parameters (parmchk2)...")
            subprocess.run([
                f"{_AMBER_BIN}/parmchk2", "-i", lig_bcc_mol2, "-f", "mol2",
                "-o", lig_frcmod
            ], check=True, capture_output=True, text=True, env=AMBER_ENV)
            
            print(" Building Topologies (tleap)...")
            with open(tleap_in, "w") as f:
                f.write("source leaprc.protein.ff14SB\n")
                f.write("source leaprc.gaff\n")
                f.write(f"rec = loadpdb {os.path.basename(sanitized_rec)}\n")
                f.write(f"lig = loadmol2 lig_bcc.mol2\n")
                f.write(f"loadamberparams lig.frcmod\n")
                f.write("comp = combine {rec lig}\n")
                f.write("saveamberparm rec rec.prmtop rec.inpcrd\n")
                f.write("saveamberparm lig lig.prmtop lig.inpcrd\n")
                f.write("saveamberparm comp comp.prmtop comp.inpcrd\n")
                f.write("quit\n")
                
            subprocess.run([f"{_AMBER_BIN}/tleap", "-f", "tleap.in"], cwd=complex_out_dir, env=AMBER_ENV,
                check=True, capture_output=True, text=True)
                           
            print(f"  ✅ Successfully prepared topologies for {base_name}")
            
        except subprocess.CalledProcessError as e:
            print(f"  ❌ Error processing {base_name} at a specific stage.")
            print("  --- ⚠️ AMBERTOOLS ERROR LOG ---")
            error_msg = e.stderr if e.stderr else e.stdout
            print(error_msg)
            print("  -------------------------------")
            
        finally:
            for jf in junk_files:
                if os.path.exists(jf):
                    os.remove(jf)
                    
            for ant_file in glob.glob("ANTECHAMBER_*"):
                if os.path.exists(ant_file):
                    if os.path.isdir(ant_file):
                        import shutil
                        shutil.rmtree(ant_file)
                    else:
                        os.remove(ant_file)

if __name__ == "__main__":
    prepare_mmgbsa_topologies()