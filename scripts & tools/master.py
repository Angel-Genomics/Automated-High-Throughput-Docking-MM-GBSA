import os,sys,glob,json,subprocess,time,shutil,re
import pandas as pd

STAGE_EXPECTATIONS = {
    "ligand.py":       {"path": "ligand", "is_dir": True},
    "clean_pro.py":    {"path": "clean_proteins_pdb", "is_dir": True},
    "flex_docking.py": {"path": "docking_results", "is_dir": True},
    "filter.py":       {"path": "best_pose_candidat", "is_dir": True},
    "rec_pdb.py":      {"path": "receptor_flex_updated", "is_dir": True},
    "amber.py":        {"path": "gbsa_topologies", "is_dir": True},
    "gbsa.py":         {"path": "Final_Outputs/MMGBSA_Final_Ranking.csv", "is_dir": False}
}

def verify_stage_output(script_name):
    if script_name not in STAGE_EXPECTATIONS:
        return True 
        
    expected = STAGE_EXPECTATIONS[script_name]
    target_path = expected["path"]
    
    if expected["is_dir"]:
        if os.path.isdir(target_path) and len(os.listdir(target_path)) > 0:
            return True
        return False
    else:
        if os.path.isfile(target_path) and os.path.getsize(target_path) > 0:
            return True
        return False

def run_script_and_capture(script_name):
    print(f"\n{'='*60}")
    print(f" STARTING STAGE: {script_name}")
    print(f"{'='*60}\n")
    
    try:
        process = subprocess.Popen(
            [sys.executable, script_name], 
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT, 
            text=True
        )
        
        for line in process.stdout:
            print(line, end='')  
            
        process.wait()
        
        if process.returncode != 0:
            print(f"\n ❌ CRITICAL ERROR: '{script_name}' crashed internally (Exit code {process.returncode}).")
            return False
            
        if not verify_stage_output(script_name):
            print(f"\n ❌ FAKE SUCCESS DETECTED: '{script_name}' finished without Python errors, BUT the expected output files are missing or empty in '{STAGE_EXPECTATIONS[script_name]['path']}'!")
            return False
            
        print(f"\n✅ STAGE VERIFIED: '{script_name}' successfully generated physical outputs.")
        return True
            
    except FileNotFoundError:
        print(f" ❌ ERROR: Could not find the file '{script_name}'.")
        return False

def parse_vina_log(log_file):
    try:
        with open(log_file, 'r') as f:
            lines = f.readlines()
            
        for i, line in enumerate(lines):
            if line.strip().startswith('----'):
                results_line = lines[i + 1]
                parts = results_line.split()
                if len(parts) >= 2:
                    return float(parts[1])
    except Exception:
        pass
    return None

def extract_protein_name(pdb_id, raw_pdb_dir="raw_pdb"):
    file_path = os.path.join(raw_pdb_dir, f"pdb{pdb_id.lower()}.ent")
    molecule_name = "Unknown"
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith("COMPND"):
                    if "MOLECULE:" in line:
                        molecule_name = line.split("MOLECULE:")[1].strip().replace(";", "")
                        break 
    except FileNotFoundError:
        pass
        
    return molecule_name

def extract_and_sort_results(active_chains_map, log_folder=".", raw_pdb_dir="raw_pdb"):
    print(f"\n{'='*60}")
    print(f" EXTRACTING AND UPDATING FLEX-DOCKING CSV")
    print(f"{'='*60}\n")

    output_dir = "Final_Outputs"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "flex_docking_summary.csv")

    existing_data = {}
    if os.path.exists(output_path):
        try:
            df_old = pd.read_csv(output_path)
            for _, row in df_old.iterrows():
                existing_data[str(row['PDB_ID'])] = row.to_dict()
        except Exception as e:
            print(f" ❌ Could not read old CSV: {e}")

    log_files = glob.glob(os.path.join(log_folder, "*_log.txt"))
    
    if not log_files:
        print(" ⚠️ No docking log files found to extract results.")
        return

    results = []
    for log_file in log_files:
        filename = os.path.basename(log_file)
        pdb_id = filename.split('_')[0].upper() 
        
        try:
            ligand_name = filename.split('_')[1]
        except IndexError:
            ligand_name = "Unknown"
        
        protein_name = extract_protein_name(pdb_id, raw_pdb_dir)
        
        target_info = active_chains_map.get(pdb_id, {})
        if isinstance(target_info, dict):
            active_chains = target_info.get("Interacting_Chains", "Unknown")
            flex_residues = target_info.get("Flexible_Residues", "None")
        else:
            active_chains = str(target_info)
            flex_residues = "None"
        
        if pdb_id in existing_data:
            old_row = existing_data[pdb_id]
            if protein_name == "Unknown" and str(old_row.get('Protein_Name', '')) not in ["Unknown", "nan", ""]:
                protein_name = old_row['Protein_Name']
            if active_chains == "Unknown" and str(old_row.get('Active_Chains', '')) not in ["Unknown", "nan", ""]:
                active_chains = old_row['Active_Chains']
            if flex_residues == "None" and str(old_row.get('Flexible_Residues', '')) not in ["None", "nan", ""]:
                flex_residues = old_row['Flexible_Residues']
        
        best_affinity = parse_vina_log(log_file)
        
        if best_affinity is not None:
            results.append({
                'PDB_ID': pdb_id,
                'Protein_Name': protein_name,
                'Ligand': ligand_name,
                'Active_Chains': active_chains,
                'Flexible_Residues': flex_residues,
                'Best_Affinity_kcal/mol': best_affinity
            })

    if not results:
        print(" ❌ Could not extract any valid results.")
        return

    df = pd.DataFrame(results)
    df_sorted = df.sort_values(by='Best_Affinity_kcal/mol', ascending=True)
    df_sorted.to_csv(output_path, index=False)
    print(f"✅ Docking Results successfully extracted and saved to: {output_path}")

def get_project_title(raw_pdb_dir="raw_pdb"):
    files = glob.glob(os.path.join(raw_pdb_dir, "*.ent")) + glob.glob(os.path.join(raw_pdb_dir, "*.pdb"))
    if not files:
        return None
        
    filepath = files[0]
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith("TITLE"):
                raw_title = line[6:].strip()
                safe_title = re.sub(r'[\\/*?:"<>|]', "", raw_title).strip()
                return safe_title
    return "UNKNOWN_PROJECT"
if __name__ == "__main__":
    print("--- STARTING FULL COMPUTATIONAL PIPELINE ---")
    start_time = time.time()
    
    scripts_to_run = [
        "ligand.py",
        "clean_pro.py",
        "flex_docking.py",
        "filter.py",
        "rec_pdb.py",
        "amber.py",
        "gbsa.py"
    ]
    
    project_title = "UNKNOWN_PROJECT" 

    for script in scripts_to_run:
        success = run_script_and_capture(script)
        if not success:
            print(f"\n PIPELINE HALTED: A critical stage failed its verification check. Fix the issue and restart.")
            sys.exit(1)
            
        if script == "clean_pro.py":
            project_title = get_project_title()
            if project_title:
                final_dest = os.path.join("final_results", project_title)
                
                if os.path.exists(final_dest):
                    print(f"\n STOPPING PIPELINE: The project '{project_title}' ALREADY EXISTS in final_results/")
                    print(" Skipping Docking and MM-GBSA to save computational resources.")
                    sys.exit(0) 

    active_chains_map = {}
    try:
        with open("configs/active_chains.json", "r") as f:
            active_chains_map = json.load(f)
    except FileNotFoundError:
        print(" ⚠️ Warning: active_chains.json not found. Chains will be marked as Unknown.")

    extract_and_sort_results(active_chains_map, log_folder="docking_results", raw_pdb_dir="raw_pdb")

    print(f"\n{'='*60}")
    print(f" ARCHIVING RESULTS AND CLEANING WORKSPACE")
    print(f"{'='*60}")

    final_dest = os.path.join("final_results", project_title)
    os.makedirs(final_dest, exist_ok=True)
    print(f" Created project archive at: {final_dest}/")

    folders_to_keep = [
        "best_pose_candidat", 
        "configs", 
        "docking_results", 
        "Final_Outputs", 
        "gbsa_topologies", 
        "receptor_flex_updated"
    ]
    
    for folder in folders_to_keep:
        if os.path.exists(folder):
            dest_folder = os.path.join(final_dest, folder)
            if os.path.exists(dest_folder):
                shutil.rmtree(dest_folder) 
                
            shutil.move(folder, dest_folder)
    folders_to_nuke = [
        "raw_pdb", 
        "clean_proteins_pdb", 
        "clean_proteins",     
        "ligand"
    ]
    
    print("\n✅ Nuking temporary files...")
    for folder in folders_to_nuke:
        if os.path.exists(folder):
            shutil.rmtree(folder)
            

    end_time = time.time()
    print(f"\n✅ ALL STAGES COMPLETED & ARCHIVED SUCCESSFULLY in {(end_time - start_time) / 60:.2f} minutes.")
    print(f"✅ Your final files are safely stored in: {final_dest}/")