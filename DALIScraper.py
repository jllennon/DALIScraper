import requests
from Bio.PDB import PDBList

class protein:
    def __init__(self, number, pdb_id, chain, l_ali, n_res, description):
        self._number = number
        self._pdb_id = pdb_id
        self._chain = chain
        self._l_ali = l_ali
        self._n_res = n_res
        self._description = description
        self._pdb_file = None

    def setPDBFile(self, pdb_file):
        self._pdb_file = pdb_file

def scrapePage(url, max_proteins):
    page = requests.get(url)

    split_text = page.text.splitlines()

    protein_list = []

    found = 0
    for i in range(len(split_text)):
        if found == 0 and "    No:  Chain   Z    rmsd lali nres  %id PDB  Description" == split_text[i]:
            found = 1
            continue

        elif found > 0:
            split_line = split_text[i].split()

            begin = split_line[14].find(">") + 1
            end = split_line[14].find("<")
            number = int(split_line[14][begin:end])

            pdb_id = split_line[15][:4]
            chain = split_line[15][-1:]
            l_ali = int(split_line[18])
            n_res = int(split_line[19])
            description = " ".join(split_line[24:])
            description = description[:-1]  # Remove the trailing ';'

            protein_list.append(protein(number, pdb_id, chain, l_ali, n_res, description))

            found += 1

        if found > max_proteins:
            break

    return protein_list

def getAndWritePDBFiles(proteins, l_ali_thresh, n_res_thresh, description_req, output_dir):
    pdb_names = []
    pdb = PDBList()

    for protein in proteins:
        if protein._l_ali >= l_ali_thresh or protein._n_res >= n_res_thresh and protein._pdb_id not in pdb_list:
            pdb.retrieve_pdb_file(protein._pdb_id, file_format="pdb", pdir=output_dir)
            pdb_names.append(protein._pdb_id)


proteins = scrapePage('http://ekhidna2.biocenter.helsinki.fi/barcosel/tmp//2v5wA/2v5wA.html', 5)

if len(proteins > 0):
    getAndWritePDBFiles(proteins, 360, 360, ("HISTONE", "HDAC"), "/Users/work/dali")
    print("Done.")
else:
    print("No proteins were found that match your search criteria.")



