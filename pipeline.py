import subprocess
import numpy as np
from pathlib import Path
from Bio.PDB import PDBParser, Select, PDBIO
from typing import List, Tuple
import os
import shutil
import warnings
import concurrent.futures
from functools import partial
import pandas as pd

warnings.filterwarnings("ignore")

def extract_chain(input_path: str, chain_id: str, output_path: str):
    """Extract a specific chain from a PDB file using Biopython"""
    class ChainSelector(Select):
        def accept_chain(self, chain):
            return chain.id == chain_id

    parser = PDBParser()
    structure = parser.get_structure('input', input_path)
    
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_path, ChainSelector())

def pdb_to_pdbqt(input_path: str, output_path: str):
    """Convert PDB to PDBQT using Open Babel with receptor parameters"""
    if not Path(input_path).exists():
        raise FileNotFoundError(f"Input PDB file {input_path} not found")
    
    result = subprocess.run(
        ["obabel", input_path, "-O", output_path, "-xr"],
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        raise RuntimeError(f"Open Babel failed: {result.stderr}")

def sdf_to_pdbqt(input_path: str, output_path: str):
    """Convert SDF to PDBQT using Open Babel with ligand parameters"""
    if not Path(input_path).exists():
        raise FileNotFoundError(f"Input SDF file {input_path} not found")
    
    result = subprocess.run(
        ["obabel", input_path, "-O", output_path, "-h", "--gen3d"],
        capture_output=True,
        text=True
    )
    
    if result.returncode != 0:
        raise RuntimeError(f"Open Babel failed: {result.stderr}")


def dock(protein_path: str, substrate_path: str, output_path: str, log_path, motif: str = None, timeout=600):
    """Run AutoDock Vina with optional motif-based or blind docking"""
    # Calculate binding site coordinates
    coords = []
    
    if motif:
        # Motif-based docking
        try:
            motif_residues = list(map(int, motif.strip().split(',')))
        except ValueError:
            raise ValueError("Invalid motif format. Use comma-separated residue numbers")

        # Collect coordinates from motif residues
        with open(protein_path) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        res_id = int(line[22:26].strip())
                        if res_id in motif_residues:
                            x = float(line[30:38])
                            y = float(line[38:46])
                            z = float(line[46:54])
                            coords.append([x, y, z])
                    except (ValueError, IndexError):
                        continue
        box_size = 20  # Smaller box for known binding site
    else:
        # Blind docking - use entire protein
        with open(protein_path) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append([x, y, z])
                    except (ValueError, IndexError):
                        continue
        box_size = 100  # Larger box for whole protein search

    if not coords:
        raise ValueError("No valid coordinates found in receptor")

    center = np.mean(coords, axis=0)
    
    vina_cmd = [
        "vina",
        "--receptor", protein_path,
        "--ligand", substrate_path,
        "--center_x", f"{center[0]:.3f}",
        "--center_y", f"{center[1]:.3f}",
        "--center_z", f"{center[2]:.3f}",
        "--size_x", str(box_size),
        "--size_y", str(box_size),
        "--size_z", str(box_size),
        "--exhaustiveness", "32" if motif else "64",  # Higher for blind docking
        "--num_modes", "10",
        "--energy_range", "5",  # Wider range for blind docking
        "--out", output_path,
    ]

    # result = subprocess.run(vina_cmd, capture_output=True, text=True)
    success = True

    try:
        result = subprocess.run(
            vina_cmd,
            capture_output=True,
            text=True,
            timeout=timeout  # None=no timeout, or set in seconds
        )
        with open(log_path, "a") as f:
            f.write(result.stdout)
            f.write(result.stderr)

        if result.returncode != 0:
            print(f"Vina docking failed: {result.stderr}")
            success = False
    except subprocess.TimeoutExpired as e:
        print(f"Docking timed out after {timeout} seconds")
        success = False

    if success:
        print(f"Docking completed! Results in {output_path}")
    
    return success

def pipeline(protein_file, substrate_file, motif, chain_id, out_dir):
    protein_name = os.path.basename(protein_file).split(".")[0]
    substrate_name = os.path.basename(substrate_file).split(".")[0]

    dirname = f"{protein_name}_{substrate_name}"
    output_dir = os.path.join(out_dir, dirname)

    os.makedirs(output_dir, exist_ok=True)
    
    chain_pdb_path = os.path.join(output_dir, "chain.pdb")
    chain_pdbqt_path = os.path.join(output_dir, "chain.pdbqt")
    substrate_pdbqt_path = os.path.join(output_dir, "substrate.pdbqt")
    dock_output = os.path.join(output_dir, "docking_results.pdbqt")
    dock_log_file = os.path.join(output_dir, "docking.log")

    # Extract chain A from PDB
    extract_chain(protein_file, chain_id, chain_pdb_path)
    
    # Convert to PDBQT
    pdb_to_pdbqt(chain_pdb_path, chain_pdbqt_path)
    sdf_to_pdbqt(substrate_file, substrate_pdbqt_path)
    
    # Run docking with motif residues
    success = dock(
        protein_path=chain_pdbqt_path,
        substrate_path=substrate_pdbqt_path,
        output_path=dock_output,
        log_path=dock_log_file,
        motif=motif
    )
    
    # if success:
    #     print(f"Visualize command: pymol {chain_pdb_path} {dock_output}")
    # else:
    #     print(f"Failure: {chain_pdb_path} {dock_output}")


def process_single_pair(args: Tuple, out_dir: str) -> str:
    """Wrapper function for parallel execution"""
    protein_file, substrate_file, motif, chain_id = args
    try:
        pipeline(protein_file, substrate_file, motif, chain_id, out_dir)
        return f"Success: {protein_file} + {substrate_file}"
    except Exception as e:
        return f"Failed: {protein_file} + {substrate_file} - {str(e)}"



def run_parallel_docking(task_list: List[Tuple], out_dir: str, workers: int = None):
    """
    Run multiple docking pipelines in parallel
    Args:
        task_list: List of tuples (protein_file, substrate_file, motif, chain_id)
        out_dir: Output directory for all results
        workers: Number of parallel processes (default: CPU count)
    """
    if not workers:
        workers = os.cpu_count() or 4

    print(f"Starting parallel docking with {workers} workers")
    print(f"Total tasks: {len(task_list)}")

    # Create process pool
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        # Create partial function with fixed out_dir
        process_fn = partial(process_single_pair, out_dir=out_dir)
        
        # Submit all tasks
        futures = [executor.submit(process_fn, task) for task in task_list]
        
        # Monitor progress
        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()
                print(result)
            except Exception as e:
                print(f"Task failed with exception: {str(e)}")

    print("All tasks completed")


# Example usage
if __name__ == "__main__":
    df = pd.read_csv("./crawler/samples.csv")

    TASKS = []
    for row in df.to_dict('records'):
        TASKS.append((
            f"./crawler/pdb_files/{row['pdb_entry']}.pdb",
            f"./crawler/chebi_sdf/{row['substrate_entry']}.sdf",
            None,
            row['pdb_chain']
        ))
    
    TASKS = TASKS[:3]

    OUTPUT_DIR = "docking_results"
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    run_parallel_docking(
        task_list=TASKS,
        out_dir=OUTPUT_DIR,
        workers=4  # Adjust based on your available cores
    )

# # Example usage
# if __name__ == "__main__":
    # pipeline(
    #     protein_file="./data/pdb/5b6m.pdb",
    #     substrate_file="./data/sdf/ChEBI_16240.sdf",
    #     motif="20,25,36,52,63,68,81,85,87,89,91,100,102,104,105,162,163,176,188,191,192,196,198,199,209,223,225,228,235,252,292,311,334",
    #     chain_id="A",
    #     out_dir="./output/bind_motif"
    # )

    # pipeline(
    #     protein_file="./data/pdb/5b6m.pdb",
    #     substrate_file="./data/sdf/ChEBI_16240.sdf",
    #     motif=None,
    #     chain_id="A",
    #     out_dir="./output/bind_no_motif"
    # )

    # # timeout
    # pipeline(
    #     protein_file="./data/pdb/4im4.pdb",
    #     substrate_file="./data/sdf/ChEBI_57540.sdf",
    #     motif="2,4,5,7,10,11,21,25,30,32,33,36,37,38,39,40,41,42,43,45,46,47,56,57,58,59,61,62,64,65,66,68,70,71,75,76,77,78,96,97,99,101,106,107,111,113,122,130,131,132,140,143,147,150,151,154,155,156,157,165,174,177",
    #     chain_id="A",
    #     out_dir="./output/no_bind_motif"
    # )

    # # timeout
    # pipeline(
    #     protein_file="./data/pdb/4im4.pdb",
    #     substrate_file="./data/sdf/ChEBI_57540.sdf",
    #     motif=None,
    #     chain_id="A",
    #     out_dir="./output/no_bind_no_motif"
    # )

    # pipeline(
    #     protein_file="./data/pdb/1v98.pdb",
    #     substrate_file="./data/sdf/ChEBI_57746.sdf",
    #     motif="2,5,7,8,10,12,13,14,18,21,32,36,37,38,40,41,46,48,57,58,60,61,65,70,72,73,79,84",
    #     chain_id="A",
    #     out_dir="./output/no_bind_motif"
    # )

    # # timeout
    # pipeline(
    #     protein_file="./data/pdb/1v98.pdb",
    #     substrate_file="./data/sdf/ChEBI_57746.sdf",
    #     motif=None,
    #     chain_id="A",
    #     out_dir="./output/no_bind_no_motif"
    # )

    # # timeout
    # pipeline(
    #     protein_file="./data/pdb/1cl0.pdb",
    #     substrate_file="./data/sdf/ChEBI_58349.sdf",
    #     motif="3,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23,24,25,26,28,31,32,33,35,39,40,45,48,49,50,52,53,56,58,60,67,69,73,77,78,83,84,85,86,87,90,92,93,95,97,98,101,102,103,105,107,108,109,110,111,112,113,114,115,118,119,120,121,123,128,130,131,132,133,134,135,137,140,143,144,147,149,150,151,152,153,157,158,159,160,161,163,166,167,168,169,170,173,174,175,176,177,180,181,183,184,187,188,189,196,197,198,200,202,203,208,209,211,212,216,218,230,235,236,237,239,241,242,243,244,245,246,248,249,250,251,260,261,262,263,264,268,271,273,277,279,280,281,282,283,284,285,286,288,290,292,293,294,295,296,298,299,301,304,305,306,307,308,309,312",
    #     chain_id="A",
    #     out_dir="./output/bind_motif"
    # )

    # # timeout
    # pipeline(
    #     protein_file="./data/pdb/1cl0.pdb",
    #     substrate_file="./data/sdf/ChEBI_58349.sdf",
    #     motif=None,
    #     chain_id="A",
    #     out_dir="./output/bind_no_motif"
    # )