import sys
import requests
from Bio.PDB import PDBList


class Protein:
    def __init__(self, number, pdb_id, chain, l_ali, n_res, description):
        self._number = number
        self._pdb_id = pdb_id
        self._chain = chain
        self._l_ali = l_ali
        self._n_res = n_res
        self._description = description
        self._pdb_file = None

    def get_number(self):
        return self._number

    def get_pdb_id(self):
        return self._pdb_id

    def get_l_ali(self):
        return self._l_ali

    def get_n_res(self):
        return self._n_res

    def get_chain(self):
        return self._chain

    def get_desc(self):
        return self._description

    def set_pdb_file(self, pdb_file):
        self._pdb_file = pdb_file

def getInputs(args):
    url = args.argv[1]
    output_dir = args.argv[2]
    nres = int(args.argv[3])
    l_ali = int(args.argv[4])
    max_num = int(args.argv[5])

    descs = []

    for i in range(len(args.argv[6:])):
        descs.append(args.argv[i])

    return url, output_dir, nres, l_ali, max_num, descs

def scrape_page(url, max_proteins):
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

            protein_list.append(Protein(number, pdb_id, chain, l_ali, n_res, description))

            found += 1

        if found > max_proteins:
            break

    return protein_list


def get_pdb_files(proteins, l_ali_thresh, n_res_thresh, description_req, output_dir):
    pdb_names = []
    structures = ''

    pdb = PDBList()

    for protein in proteins:
        if protein.get_l_ali() >= l_ali_thresh \
                and protein.get_n_res() >= n_res_thresh \
                and protein.get_pdb_id() not in pdb_names:

            desc_found = False
            for word in description_req:
                if word in description_req:
                    desc_found = True
                    break

            if not desc_found:
                continue

            pdb.retrieve_pdb_file(protein.get_pdb_id(), file_format="pdb", pdir=output_dir)
            pdb_names.append(protein.get_pdb_id())

            if len(structures) > 0:
                structures += '\n'

            structures += protein.get_pdb_id() + ' ' + protein.get_chain()

    return structures


def write_pdb_list(structures, output_dir):
    with open(output_dir + '/structures_list.txt', 'w+') as f:
        f.write(structures)


url, output_dir, l_ali, n_res, max_num, desc = getInputs(sys)
# proteins = scrape_page('http://ekhidna2.biocenter.helsinki.fi/barcosel/tmp//2v5wA/2v5wA.html', max_num)
proteins = scrape_page(url, max_num)

if len(proteins) > 0:
    structures = get_pdb_files(proteins, l_ali, n_res, desc, output_dir)
    write_pdb_list(structures, output_dir)

    print("Done.")
else:
    print("No proteins were found that match your search criteria.")



