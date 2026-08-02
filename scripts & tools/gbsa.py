import os,glob,subprocess
import pandas as pd

_AMBER_BIN = "/home/angelommd/miniconda3/envs/amber_env/bin"
AMBER_ENV = {
    **os.environ,
    "AMBERHOME": "/home/angelommd/miniconda3/envs/amber_env",
    "PATH": f"{_AMBER_BIN}:{os.environ.get('PATH', '')}",
}

def run_minimization_and_gbsa():
    print("--- Starting Pipeline: Minimization -> MM-GBSA -> Extraction ---")
    
    out_dir = "gbsa_topologies"
    mmpbsa_in = "mmpbsa.in"
    min_in = "min.in"
    
    if not os.path.exists(min_in):
        min_content = """Minimize complex to remove clashes
 &cntrl
  imin=1, maxcyc=500, ncyc=250,
  cut=999.0, rgbmax=999.0,
  igb=5, saltcon=0.150,
  ntpr=100,
 /
"""
        with open(min_in, "w") as f:
            f.write(min_content)

    if not os.path.exists(mmpbsa_in):
        mmpbsa_content = """Input file for running MM-GBSA
&general
   endframe=1, keep_files=0,
/
&gb
   igb=5, saltcon=0.150,
/
"""
        with open(mmpbsa_in, "w") as f:
            f.write(mmpbsa_content)

    complex_dirs = [d for d in glob.glob(f"{out_dir}/*") if os.path.isdir(d)]
    
    if not complex_dirs:
        print(f" ⚠️ No complex directories found in '{out_dir}'.")
        return False
    
    for complex_dir in complex_dirs:
        base_name = os.path.basename(complex_dir)
        print(f"\n Processing: {base_name}")
        
        cp = os.path.join(complex_dir, "comp.prmtop")
        y  = os.path.join(complex_dir, "comp.inpcrd")
        min_crd = os.path.join(complex_dir, "min.ncrst")
        out_dat = os.path.join(complex_dir, "FINAL_RESULTS_MMPBSA.dat")
        
        if os.path.exists(min_crd) and os.path.getsize(min_crd) > 0:
            print("✅ Already minimized (Valid output found).")
        else:
            print(" Relaxing structure (Minimization)...")
            try:
                cmd_min = [
                    f"{_AMBER_BIN}/sander", "-O", 
                    "-i", f"../../{min_in}", 
                    "-o", "min.out", 
                    "-p", "comp.prmtop", 
                    "-c", "comp.inpcrd", 
                    "-r", "min.ncrst"
                ]
                subprocess.run(cmd_min, cwd=complex_dir, check=True, env=AMBER_ENV)
            except subprocess.CalledProcessError:
                print(" ❌ Minimization failed! Check min.out.")
                if os.path.exists(min_crd): os.remove(min_crd) 
                continue

        if os.path.exists(out_dat) and os.path.getsize(out_dat) > 0:
            print("✅ Energy already calculated (Valid output found).")
        else:
            print(" Calculating MM-GBSA Energy...")
            try:
                cmd_gbsa = [
                    f"{_AMBER_BIN}/MMPBSA.py", "-O", "-i", f"../../{mmpbsa_in}", 
                    "-o", "FINAL_RESULTS_MMPBSA.dat",
                    "-cp", "comp.prmtop", "-rp", "rec.prmtop",
                    "-lp", "lig.prmtop", "-y", "min.ncrst"
                ]
                subprocess.run(cmd_gbsa, cwd=complex_dir, check=True, stdout=subprocess.DEVNULL, env=AMBER_ENV)
                print("✅ Calculation successful!")
            except subprocess.CalledProcessError:
                print(" ❌ MM-GBSA failed.")
                if os.path.exists(out_dat): os.remove(out_dat) 
                continue
            finally:
                for temp_file in glob.glob(f"{complex_dir}/_MMPBSA_*"):
                    if os.path.exists(temp_file):
                        os.remove(temp_file)

    return True

def generate_final_report():
    print("\n---  Generating Final Comprehensive MM-GBSA Report ---")
    out_dir = "gbsa_topologies"
    data = []
    complex_dirs = [d for d in glob.glob(f"{out_dir}/*") if os.path.isdir(d)]
    
    for complex_dir in complex_dirs:
        ligand_name = os.path.basename(complex_dir)
        pdb_id = ligand_name.split('_')[0] 
        
        res_file = os.path.join(complex_dir, "FINAL_RESULTS_MMPBSA.dat")
        decomp_file = os.path.join(complex_dir, "FINAL_DECOMP_MMPBSA.dat")
        
        delta_g = None
        if os.path.exists(res_file):
            with open(res_file, 'r') as f:
                for line in f:
                    if line.startswith("DELTA TOTAL"):
                        try:
                            delta_g = float(line.split()[2])
                        except (IndexError, ValueError):
                            pass
                        break
        
        res_mapping = build_residue_map(pdb_id)
        top_residues_str = get_top_5_residues_string(decomp_file, res_mapping)
        
        if delta_g is not None:
            data.append({
                "System": ligand_name, 
                "Delta_G_Total (kcal/mol)": delta_g,
                "Top_5_Interacting_Residues": top_residues_str
            })
            print(f"✅ {ligand_name:18s} | Energy: {delta_g:7.2f} | Top: {top_residues_str}")
        else:
            print(f" ⚠️ Could not find Energy for {ligand_name}")

    if data:
        df = pd.DataFrame(data)
        df = df.sort_values(by="Delta_G_Total (kcal/mol)", ascending=True)
        
        results_dir = "Final_Outputs"
        os.makedirs(results_dir, exist_ok=True)
        csv_filename = os.path.join(results_dir, "MMGBSA_Final_Ranking.csv")
        df.to_csv(csv_filename, index=False)
        print(f"\n✅ Final Report saved to: '{csv_filename}'")

def build_residue_map(pdb_id):
    clean_pdb_path = f"clean_proteins_pdb/{pdb_id}_clean.pdb"
    res_map = {}
    
    if not os.path.exists(clean_pdb_path):
        return res_map

    amber_index = 1
    last_res_num = None
    
    with open(clean_pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                res_num = line[22:26].strip() 
                if res_num != last_res_num:
                    res_map[amber_index] = res_num
                    last_res_num = res_num
                    amber_index += 1
                    
    return res_map

def run_decomposition():
    print("---  Starting Per-Residue Energy Decomposition ---")
    out_dir = "gbsa_topologies"
    mmpbsa_decomp_in = "mmpbsa_decomp.in"

    if not os.path.exists(mmpbsa_decomp_in):
        decomp_content = """Input file for running MM-GBSA Per-Residue Decomposition
&general
   endframe=1, keep_files=0,
/
&gb
   igb=5, saltcon=0.150,
/
&decomp
   idecomp=1,
   dec_verbose=0,
/
"""
        with open(mmpbsa_decomp_in, "w") as f:
            f.write(decomp_content)

    complex_dirs = [d for d in glob.glob(f"{out_dir}/*") if os.path.isdir(d)]
    
    for complex_dir in complex_dirs:
        base_name = os.path.basename(complex_dir)
        out_decomp = os.path.join(complex_dir, "FINAL_DECOMP_MMPBSA.dat")
        main_decomp_out = os.path.join(complex_dir, "FINAL_RESULTS_MMPBSA_DECOMP.dat")
        
        if os.path.exists(out_decomp) and os.path.getsize(out_decomp) > 0:
            print(f"\n✅ Decomposition already calculated for {base_name}.")
        else:
            print(f"\n✅ Calculating Decomposition for: {base_name} ")
            try:
                cmd_gbsa = [
                    f"{_AMBER_BIN}/MMPBSA.py", "-O", 
                    "-i", f"../../{mmpbsa_decomp_in}", 
                    "-o", "FINAL_RESULTS_MMPBSA_DECOMP.dat", 
                    "-do", "FINAL_DECOMP_MMPBSA.dat",
                    "-cp", "comp.prmtop", "-rp", "rec.prmtop",
                    "-lp", "lig.prmtop", "-y", "min.ncrst"
                ]
                subprocess.run(cmd_gbsa, cwd=complex_dir, check=True, stdout=subprocess.DEVNULL, env=AMBER_ENV)
                
            except subprocess.CalledProcessError:
                print(f" ❌ Decomposition failed for {base_name}. Skipping to next...")
                if os.path.exists(out_decomp): os.remove(out_decomp)
                if os.path.exists(main_decomp_out): os.remove(main_decomp_out)
                continue
            finally:
                for temp_file in glob.glob(f"{complex_dir}/_MMPBSA_*"):
                    if os.path.exists(temp_file):
                        os.remove(temp_file)

def get_top_5_residues_string(decomp_file, res_mapping):
    residues = []
    in_delta_section = False
    
    if not os.path.exists(decomp_file):
        return "No Decomposition Data"

    with open(decomp_file, 'r') as f:
        lines = f.readlines()
        
    for line in lines:
        if "DELTAS:" in line:
            in_delta_section = True
            continue
            
        if in_delta_section:
            if "Total Energy" in line or "Residue,Location" in line or "Avg.,Std" in line:
                continue
                
            parts = line.split(',')
            if len(parts) >= 18 and parts[1].startswith('R '):
                fake_res_str = parts[0].strip()
                res_name = fake_res_str[:3]
                fake_num = int(fake_res_str[3:].strip())
                
                real_num = res_mapping.get(fake_num, f"?({fake_num})")
                final_res_name = f"{res_name}{real_num}"
                
                try:
                    total_energy = float(parts[17])
                    residues.append((final_res_name, total_energy))
                except ValueError:
                    pass
            elif len(parts) > 1 and parts[1].startswith('L '):
                break

    if residues:
        residues.sort(key=lambda x: x[1]) 
        top_5 = residues[:5]
        return ", ".join([f"{res} ({energy:.2f})" for res, energy in top_5])
        
    return "None"

if __name__ == "__main__":
    if run_minimization_and_gbsa():
        run_decomposition()         
        generate_final_report()     
        
        for f in ["min.in", "mmpbsa.in", "mmpbsa_decomp.in"]:
            if os.path.exists(f):
                os.remove(f)