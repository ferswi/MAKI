import os
import requests
import time
import json
import csv

def download_alphafold_structures(output_dir=r"C:\Users\sophi\Documents\NRG\fake_kinases",
                                  taxon_id=9606,
                                  min_len=200,
                                  max_len=300,
                                  exclude_kinase=True,
                                  limit=250):
    os.makedirs(output_dir, exist_ok=True)

    # Use the "stream" endpoint for robust querying
    base_url = "https://rest.uniprot.org/uniprotkb/stream"
    query = f"(organism_id:{taxon_id}) AND (length:[{min_len} TO {max_len}])"
    if exclude_kinase:
        query += " AND (NOT keyword:Kinase)"

    params = {
        "compressed": "false",
        "format": "json",
        "query": query
    }

    print("Sending query to UniProt stream endpoint...")
    r = requests.get(base_url, params=params)
    r.raise_for_status()
    results = r.json()

    print(f"Retrieved {len(results.get('results', []))} entries.")

    count = 0
    summary = []
    for entry in results.get("results", []):
        uniprot_id = entry["primaryAccession"]
        length = entry.get("sequence", {}).get("length")
        name = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value")

        pdb_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
        pdb_path = os.path.join(output_dir, f"{uniprot_id}.pdb")

        try:
            pdb_response = requests.get(pdb_url)
            if pdb_response.status_code == 200 and "ATOM" in pdb_response.text:
                with open(pdb_path, "w") as f:
                    f.write(pdb_response.text)
                print(f"[{count+1}] Downloaded {uniprot_id}")
                summary.append((uniprot_id, name, length, True))
                count += 1
            else:
                print(f"Skipping {uniprot_id}: no PDB file")
                summary.append((uniprot_id, name, length, False))
        except Exception as e:
            print(f"Error downloading {uniprot_id}: {e}")
            summary.append((uniprot_id, name, length, False))

        if count >= limit:
            break
        time.sleep(0.2)

    # Save a summary CSV
    csv_path = os.path.join(output_dir, "download_summary.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["UniProt ID", "Protein Name", "Length", "Downloaded"])
        writer.writerows(summary)

    print(f"\n✅ Completed! {count} PDB files downloaded.")
    print(f"Summary log: {csv_path}")

if __name__ == '__main__':
    download_alphafold_structures()