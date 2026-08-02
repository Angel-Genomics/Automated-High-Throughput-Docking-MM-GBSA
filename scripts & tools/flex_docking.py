import os,sys,glob,subprocess

ligand_files = glob.glob("ligand/*.pdbqt")
if not ligand_files:
    print(" ❌ Error: No ligand .pdbqt file found in the 'ligand/' directory.")
    sys.exit(1)

config_files = glob.glob("configs/*_config.txt")
if not config_files:
    print(" ❌ Error: No config files found in 'configs/' directory.")
    sys.exit(1)

output_dir = "docking_results"
os.makedirs(output_dir, exist_ok=True)

print(f"--- Starting GPU-Accelerated Docking Pipeline with Uni-Dock ---")
print(f"Found {len(ligand_files)} ligands and {len(config_files)} targets to process.\n")

for ligand_path in ligand_files:
    ligand_name = os.path.basename(ligand_path).replace(".pdbqt", "")
    print(f"\n Starting Batch for Ligand: {ligand_name}")
    
    for config in config_files:
        base_name = os.path.basename(config)
        pdb_id = base_name.replace("_config.txt", "")
        
        out_pdbqt = os.path.join(output_dir, f"{pdb_id}_{ligand_name}_out.pdbqt")
        log_file = os.path.join(output_dir, f"{pdb_id}_{ligand_name}_log.txt")
        
        print(f"  Target: {pdb_id} | Ligand: {ligand_name}")
        
        if os.path.exists(out_pdbqt) and os.path.getsize(out_pdbqt) > 0:
            print(f"✅ SKIPPED: Already docked! (File exists: {os.path.basename(out_pdbqt)})")
            continue  
        cmd_unidock = [
    "/home/angelommd/miniconda3/bin/unidock","--config", config,"--ligand", ligand_path,"--out", out_pdbqt]  
        
        print(f" Running Uni-Dock (GPU Mode)...")
        try:
            with open(log_file, "w") as f:
                process = subprocess.run(cmd_unidock, stdout=f, stderr=subprocess.STDOUT)
                
            if process.returncode == 0 and os.path.exists(out_pdbqt):
                print(f"✅ Docking Successful! Saved to {os.path.basename(out_pdbqt)}")
            else:
                print(f" ❌ Error: Docking failed for {pdb_id} with {ligand_name}.")
                print("  =" * 50)
                with open(log_file, "r") as err_log:
                    print(err_log.read())
                print("  =" * 50)
                continue 
                
        except Exception as e:
            print(f" ❌ System Error (Uni-Dock execution failed): {e}")
            sys.exit(1)
            
        print("  " + "-" * 40)

print("\n--- ✅ GPU Docking Pipeline Completed Successfully! ---")