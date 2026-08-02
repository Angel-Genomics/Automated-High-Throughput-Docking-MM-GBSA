import os, time, warnings, requests, urllib.parse, re
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors  
from meeko import MoleculePreparation

warnings.filterwarnings("ignore", category=DeprecationWarning)

ligands_list = ["adamantane", "aspirin", "atorvastatin", "avapritinib", "caffeine", "cholesterol", "ibuprofen", "imatinib", "nilotinib", "sunitinib"]
output_dir = "ligand"
os.makedirs(output_dir, exist_ok=True)

failed_ligands = []
successful_count = 0

print("--- Starting High-Throughput Ligand Preparation ---")
print(f"Total ligands to process: {len(ligands_list)}\n")

for ligand_name in ligands_list:
    print(f"\n Processing: {ligand_name}")
    final_pdbqt = os.path.join(output_dir, f"{ligand_name}.pdbqt")
    
    if os.path.exists(final_pdbqt) and os.path.getsize(final_pdbqt) > 0:
        print(f"✅ SKIPPED: Valid {ligand_name}.pdbqt already exists.")
        successful_count += 1
        continue
    
    safe_name = urllib.parse.quote(ligand_name)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{safe_name}/SDF?record_type=3d"
    time.sleep(0.5) 
    
    try:
        print("  Downloading 3D structure from PubChem...")
        response = requests.get(url)
        
        if response.status_code != 200:
            print(f" ❌ ERROR: Download failed. (HTTP {response.status_code})")
            failed_ligands.append(f"{ligand_name} (Download Failed: HTTP {response.status_code})")
            continue  
            
        sdf_data = response.text
        mol = Chem.MolFromMolBlock(sdf_data, removeHs=False)
        
        if mol is None:
            print(" ❌ ERROR: RDKit could not read the 3D structure.")
            failed_ligands.append(f"{ligand_name} (Invalid 3D Structure)")
            continue  

        print("✅ Optimizing ligand geometry (Energy Minimization with MMFF94)...")
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=1000, nonBondedThresh=10.0)
            print("✅ Energy minimization successful.")
        except Exception as e:
            print(f" ⚠️ Warning: Minimization failed, continuing with raw structure. Error: {e}")
            
        print("✅ Converting to PDBQT with Meeko (Assigning charges & rotatable bonds)...")
        preparator = MoleculePreparation()
        preparator.prepare(mol)
        pdbqt_string = preparator.write_pdbqt_string()
        
        match = re.search(r"(\d+)\s+active torsion", pdbqt_string, re.IGNORECASE)
        if match:
            active_torsions = int(match.group(1))
            if active_torsions > 0:
                print(f"✅ Ligand Analysis: Found {active_torsions} active rotatable bond(s). Flexibility enabled.")
            else:
                print("✅ Ligand Analysis: 0 rotatable bonds found. Ligand is fully rigid.")
        else:
            rdkit_torsions = rdMolDescriptors.CalcNumRotatableBonds(mol)
            print(f"✅ Ligand Analysis (via RDKit): Found {rdkit_torsions} rotatable bond(s). Flexibility enabled.")
        
        with open(final_pdbqt, "w") as f:
            f.write(pdbqt_string)
            
        print(f"✅ Success! Saved to {final_pdbqt}")
        successful_count += 1
        
    except Exception as e:
        print(f" ❌ CRITICAL: Unexpected error while processing {ligand_name}: {e}")
        failed_ligands.append(f"{ligand_name} (Critical Error: {e})")

print("\n===========================================")
print("--- ✅ Batch Ligand Preparation Complete! ---")
print(f"Successfully processed: {successful_count} / {len(ligands_list)}")

if failed_ligands:
    print("\n--- ⚠️ Failed Ligands Report ---")
    for failure in failed_ligands:
        print(f" - {failure}")
        
    report_file = "rejected_ligands_report.txt"
    with open(report_file, "w") as f:
        f.write("Failed Ligands Report\n")
        f.write("======================\n")
        for failure in failed_ligands:
            f.write(f"{failure}\n")
    print(f"\n ⚠️ A detailed report of failed ligands has been saved to '{report_file}'")
print("===========================================")