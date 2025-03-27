import csv
import os
import sys

def parse_vina_log(log_path, dir_name):
    """Parse a Vina log file and return data with directory name as filename"""
    data = []
    try:
        with open(log_path, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        return data
    
    in_table = False
    for line in lines:
        stripped_line = line.strip()
        if in_table:
            if not stripped_line:  # End of table
                break
            parts = stripped_line.split()
            if len(parts) < 4 or not parts[0].isdigit():
                continue  # Skip non-data lines
            try:
                entry = {
                    'mode': parts[0],
                    'affinity': parts[1],
                    'rmsd_lb': parts[2],
                    'rmsd_ub': parts[3],
                    'filename': dir_name
                }
                data.append(entry)
            except IndexError:
                continue
        elif stripped_line.startswith('-----'):
            in_table = True  # Table starts after this line
    return data

def main():
    if len(sys.argv) != 2:
        print("Usage: python parse_vina.py <input_directory>")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    if not os.path.isdir(input_dir):
        print(f"Error: {input_dir} is not a valid directory")
        sys.exit(1)
    
    all_data = []
    # Iterate through all subdirectories in input directory
    for entry in os.listdir(input_dir):
        subdir_path = os.path.join(input_dir, entry)
        if os.path.isdir(subdir_path):
            log_path = os.path.join(subdir_path, "docking.log")
            if os.path.isfile(log_path):
                print(f"Processing: {entry}")
                all_data.extend(parse_vina_log(log_path, entry))
    
    # Write to CSV
    with open('docking_metrics.csv', 'w', newline='') as csvfile:
        fieldnames = ['mode', 'affinity', 'rmsd_lb', 'rmsd_ub', 'filename']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_data)
    
    print(f"Successfully processed {len(all_data)} entries")

if __name__ == '__main__':
    main()