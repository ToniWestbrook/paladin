import re
import dataset
import gzip
import sys

def create_swissprot_database(file):
    """Creates a new sqlite database from supplied uniprot .dat file"""
    mrna_pattern = re.compile("(NM_|XM_)")
    # db = dataset.connect('sqlite:///test_data/paladin.db')
    db = dataset.connect('sqlite:///paladin.db')
    swissprot_table = db.get_table('swissprot')
    ac_table = db.get_table('ac')
    entries = []
    ac_entries = []
    term = None
    million = 0
    linecounter = 0
    acs = ""

    # iterate over dat file
    lines = file.readlines(100000)
    while lines:
        for line in lines:
            header = line[0:2]
            if header == 'ID':
                if term:
                    linecounter += 1
                    term['go_terms'] = term['go_terms'].rstrip(',')
                    entries.append(term)

                    acs = acs.replace("\n","").split(";")
                    acs.pop()
                    for i in acs:
                        ac_entries.append(dict(ac=i, name=term['name']))

                    if linecounter % 100000 == 0:
                        print("Current number of entries added: {}".format(
                            linecounter))
                        db.begin()
                        swissprot_table.insert_many(entries, 10000)
                        ac_table.insert_many(ac_entries, 10000)
                        db.commit()
                        del entries[:]
                        del ac_entries[:]

                name = re.match("([^\s]*)", line[5:]).group(1)
                term = dict(name=name, refseq="", taxa="", go_terms="")
                acs = ""

            # parse out information to add to dictionary (and then to database)
            if header == 'AC':
                acs += line[5:].rstrip()
            elif header == 'OC':
                term['taxa'] += line[5:].rstrip()
            elif header == 'DR':
                if line[5:7] == 'GO':
                    m = re.search("GO:(\d+)", line)
                    term["go_terms"] += m.group(1) + ","
                elif line[5:11] == 'RefSeq' and not term['refseq']:
                    m = mrna_pattern.search(line)
                    if m:
                        term['refseq'] = line.rstrip()

        lines = file.readlines(100000)

    # insert any leftover entries
    db.begin()
    swissprot_table.insert_many(entries, 10000)
    ac_table.insert_many(ac_entries, 10000)
    db.commit()

    # index ac / name for faster queries
    ac_table.create_index(['ac'], "acname")
    swissprot_table.create_index(['name'], "acname")

if __name__ == "__main__":
    with gzip.open(sys.argv[1], 'rb') as f:
        create_swissprot_database(f)
