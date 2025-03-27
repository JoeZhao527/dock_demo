import os
import requests
import time
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
import json

def download_pdb(pdb_id, save_dir="pdb_files", max_retries=5, initial_delay=1):
    """Download a PDB file given its ID with exponential backoff for rate limits."""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    retries = 0
    delay = initial_delay
    
    while retries < max_retries:
        response = requests.get(url)
        
        if response.status_code == 200:
            os.makedirs(save_dir, exist_ok=True)
            file_path = os.path.join(save_dir, f"{pdb_id}.pdb")
            
            with open(file_path, "wb") as file:
                file.write(response.content)
            print(f"Downloaded: {pdb_id}.pdb")
            return
        elif response.status_code == 429:
            print(f"Rate limit exceeded for {pdb_id}, retrying in {delay} seconds...")
            time.sleep(delay)
            delay *= 2  # Exponential backoff
        else:
            print(f"Failed to download: {pdb_id} (Status Code: {response.status_code})")
            return
        
        retries += 1
    print(f"Failed to download {pdb_id} after {max_retries} retries.")

def download_pdb_files(pdb_list, save_dir="pdb_files", max_workers=5):
    """Download multiple PDB files using multithreading."""
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        executor.map(lambda pdb_id: download_pdb(pdb_id, save_dir), pdb_list)

def download_chebi_sdf(chebi_id, save_dir="chebi_sdf", max_retries=5, initial_delay=1):
    """Download an SDF file from the ChEBI database with exponential backoff for rate limits."""
    url = f"https://www.ebi.ac.uk/chebi/saveStructure.do?defaultImage=true&chebiId={chebi_id}&format=sdf"
    retries = 0
    delay = initial_delay
    
    while retries < max_retries:
        response = requests.get(url)
        
        if response.status_code == 200:
            os.makedirs(save_dir, exist_ok=True)
            file_path = os.path.join(save_dir, f"{chebi_id}.sdf")
            
            with open(file_path, "wb") as file:
                file.write(response.content)
            print(f"Downloaded: {chebi_id}.sdf")
            return
        elif response.status_code == 429:
            print(f"Rate limit exceeded for {chebi_id}, retrying in {delay} seconds...")
            time.sleep(delay)
            delay *= 2  # Exponential backoff
        else:
            print(f"Failed to download: {chebi_id} (Status Code: {response.status_code})")
            return
        
        retries += 1
    print(f"Failed to download {chebi_id} after {max_retries} retries.")

def download_chebi_sdf_files(chebi_list, save_dir="chebi_sdf", max_workers=5):
    """Download multiple ChEBI SDF files using multithreading."""
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        executor.map(lambda chebi_id: download_chebi_sdf(chebi_id, save_dir), chebi_list)


with open("../enzygen_bench/enzyme_substrate_data_lucky_best.json", 'r') as f:
    data = json.load(f)

selected_samples = []

for ec_key in data:
    for trn_key in data[ec_key]:
        data_dict = data[ec_key][trn_key]
        ec_size = len(data_dict['pdb'])
        for i in range(min(5, ec_size)):
            selected_samples.append({
                k: data_dict[k][i] for k in data_dict
            })

samples = pd.DataFrame(selected_samples)
samples['pdb'] = samples['pdb'].apply(lambda x: x.replace('\n', '').upper())
samples['pdb_entry'] = samples['pdb'].apply(lambda x: x.split('.')[0])
samples['pdb_chain'] = samples['pdb'].apply(lambda x: x.split('.')[1])
samples['substrate_entry'] = samples['substrate'].apply(lambda x: x.split('_')[-1].split('.')[0])

pdb_entries = list(samples['pdb_entry'].unique())
download_pdb_files(pdb_entries)

chebi_entries = list(samples['substrate_entry'].unique())
download_chebi_sdf_files(chebi_entries)

samples.to_csv("./samples.csv", index=False)