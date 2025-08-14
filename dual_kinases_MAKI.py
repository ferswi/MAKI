import os
from importlib.metadata import files
from pathlib import Path

def clearing_duplicates(input_dir,verbose):
    input_dir = Path(input_dir)
    kinase_groups = {}
    for kinase_file in os.listdir(input_dir):
        file_path = input_dir / kinase_file
        if not file_path.is_file():
            continue
        kinase_name = kinase_file.split('_')[1]
        kinase_groups.setdefault(kinase_name, []).append(file_path)

    for kinase_name,files in kinase_groups.items() :
        if len(files)<=1:
            continue
        keep_file=None
        for file in files:
            if "PDB" in file.name:
                keep_file=file
                break
        if keep_file is None:
            keep_file=files[0]
        for file in files:
            if file != keep_file:
                verbose and print(f"removing file :{file}")
                os.remove(file)
    remaining_files = [f for f in input_dir.iterdir() if f.is_file()]
    print("Amount of files left ! : ",(len(remaining_files)))

if __name__ == '__main__':
    clearing_duplicates(input_dir=r"C:\Users\sophi\Documents\NRG\KINASES\Kincore_dual_PDBAF2_Active",verbose=True)
