import sys
import re
import gzip
import subprocess


def download_datasets(f):
    """Download cds from embl using embl ids fetched from .dat file"""

    # open up protein / nucleotide files which we will write to
    swiss_nuc_out = gzip.open('test_data/uniprot_nuc.fa.gz', 'wb')
    swiss_prot_out = gzip.open('test_data/uniprot_prot.fa.gz', 'wb')
    ignore_pattern = re.compile("(NOT_ANNOTATED|ALT|JOINED|mRNA)")
    reading_protein = False
    block = f.readlines(50000)
    term_list = []
    term = None
    sequence_counter = 0

    while block:
        # parse through the swissprot dat file in chunks
        for line in block:

            # start of a new uniprot entry
            if line[0:2] == "ID":

                # if we're at 400 terms, move onto fetching ebi data
                if len(term_list) >= 400:
                    sequence_counter += 400
                    if sequence_counter % 100000 == 0:
                        print("Sequences downloaded:" + str(sequence_counter))
                    fetch_ebi_data(term_list, swiss_nuc_out, swiss_prot_out)
                    del term_list[:]

                # otherwise, append to term list if it had an embl entry
                elif term and term['embl']:
                    term_list.append(term)
                reading_protein = False

                # create a new term and assign name; start reading sequences
                id = re.match("([^\s]*)", line[5:]).group(1)
                term = dict(name=id, embl="", sequence="")
            if not term:
                continue
            # parse EMBL entry (don't add if it's not annotated)
            elif line[5:9] == "EMBL" and not ignore_pattern.search(line):
                elements = line.split(";")
                #embl = re.sub("\..*","", elements[2].replace(" ",""))
                embl = elements[2][1:]
                term['embl'] = embl

            # only record sequence data if it has an embl entry
            elif term['embl'] and line[0:2] == "SQ":
                reading_protein = True
                continue

            if reading_protein and line[0] != "/":
                term['sequence'] += line.replace(" ","")

        block = f.readlines(50000)

    # close files
    swiss_prot_out.close()
    swiss_nuc_out.close()


def fetch_ebi_data(term_list, nuc_out, prot_out, retry=0):
    """Attempt to download ebi entries (will retry multiple times)"""
    # construct query for ebi
    base = "http://www.ebi.ac.uk/ena/data/view/"
    query = base + ",".join([term['embl'] for term in term_list])
    query += "&display=fasta"

    # send query to ebi
    proc = subprocess.Popen(["curl", "-s", query], stdout=subprocess.PIPE)
    (out, err) = proc.communicate()

    # if something went wrong, try a couple more times
    if retry >= 2:
        return
    if not out or err:
        fetch_ebi_data(term_list, nuc_out, prot_out, retry=(retry + 1))


    # construct dictionary with fasta header as key, sequence as value
    fasta_pattern = re.compile("(\>.*?)\n([^>]*)")
    acc_pattern = re.compile(".*\|([^\s]+)")
    fasta_match = fasta_pattern.findall(out)
    if not fasta_match:
        return
    fasta_dict = {}
    for (key, value) in fasta_match:
        acc_match = acc_pattern.match(key)
        if acc_match:
            fasta_dict[acc_match.group(1)] = value


    for term in term_list:

        # check to see if our terms were able to fetch cds
        embl_sequence = fasta_dict.get(term['embl'])
        if not embl_sequence:
            continue

        # write them to output if so
        nuc_out.write(">" + term['name'] + "\n" + embl_sequence)
        prot_out.write(">" + term['name'] + "\n" + term['sequence'])

if __name__ == '__main__':
    with gzip.open(sys.argv[1], 'rb') as swiss_data_file:
        download_datasets(swiss_data_file)

