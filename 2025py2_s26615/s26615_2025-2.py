from Bio import Entrez, SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import time

email_address = input("Email: ")
api_key_input = input("API Key: ")
tax_id = input("TaxID: ")
try:
    min_seq_len = int(input("Min len: "))
    max_seq_len = int(input("Max len: "))
except ValueError:
    print("Invalid input.")
    exit()
Entrez.email = email_address
Entrez.api_key = api_key_input
Entrez.tool = "GenBankRetrieverScript"
try:
    Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml").close()
    search_handle = Entrez.esearch(db="nucleotide", term=f"txid{tax_id}[Organism]", usehistory="y")
    search_results = Entrez.read(search_handle)
    web_env, query_key, total_count = search_results["WebEnv"], search_results["QueryKey"], int(search_results["Count"])
except Exception as e:
    print(f"Error occurred: {e}")
    exit()
sequence_records = []
fetched_count = 0
batch_size = 200
max_sequences = 100
for start_index in range(0, total_count, batch_size):
    if fetched_count >= max_sequences:
        break
    try:
        fetch_handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text",
                                     retstart=start_index, retmax=batch_size, webenv=web_env, query_key=query_key)
        for record in SeqIO.parse(fetch_handle, "genbank"):
            sequence_length = len(record.seq)
            if min_seq_len <= sequence_length <= max_seq_len:
                sequence_records.append({"accession_number": record.id, "sequence_length": sequence_length, "description": record.description})
                fetched_count += 1
                if fetched_count >= max_sequences:
                    break
        time.sleep(0.4)
    except Exception as e:
        print(f"Error fetching data: {e}")

if not sequence_records:
    print(":(")
    exit()
csv_filename = f"filtered_taxid_{tax_id}.csv"
df_records = pd.DataFrame(sequence_records)
df_records.to_csv(csv_filename, index=False)
print(f"CSV file saved as {csv_filename}")
plot_filename = f"plot_taxid_{tax_id}.png"
df_sorted = df_records.sort_values(by='sequence_length', ascending=False)
plt.figure(figsize=(10, 5))
plt.plot(df_sorted["accession_number"], df_sorted["sequence_length"], marker='o')
plt.xticks(rotation=90)
plt.xlabel("Accession Number")
plt.ylabel("Seq Len")
plt.title("Seq Len Distribution")
plt.tight_layout()
plt.savefig(plot_filename)
plt.close()