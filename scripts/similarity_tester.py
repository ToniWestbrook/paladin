#!/usr/bin/python2
import argparse
import re
import sys
import networkx as nx
import dataset
from subprocess import PIPE, Popen
from collections import namedtuple

# BIOLOGICAL_PROCESS = "0008150"
# CELLULAR_COMPONENT = "0005575"
# MOLECULAR_FUNCTION = "0003674"
GO = namedtuple('GO', ['id', 'desc'])

class Node(object):
    """Simple node class for representing parent/child relationships"""
    def __init__(self, id, parent, depth):
        self.parent = parent
        self.depth = depth
        self.children = []
        self.id = id


class GoTerm(object):
    """Container for GO term information extracted from go.obo"""
    def __init__(self, id):
        self.id = str(id)
        self.name = ""
        self.desc = bytearray()
        self.synonym = bytearray()
        self.is_a = []
        self.namespace = bytearray()
        self.children = []
        self.parents = []
        self.depth = 0


def create_godict(filename):
    """Parse .obo file for GO terms and stores them in a dict"""
    f = open(filename, "r+b")
    term = None
    go_dict = {}
    start_reading = False
    while 1:
        line = f.readline()
        if line == '' or line == '[Typedef]\n':
            break
        if line == '[Term]\n':
            start_reading = True

        if not start_reading or line == "\n" or line[0] == '[':
            continue

        header = re.match("(.*?):", line).group(1)
        if header == 'id':
            name = re.search("GO:(\w*)", line).group(1)
            term = GoTerm(name)
            go_dict[str(name)] = term
        elif header == 'name':
            term.name = line[6:]
        elif header == 'namespace':
            term.namespace = line[11:]
        elif header == 'def':
            term.desc = re.search('\"(.*?)\"', line).group(1)
        elif header == 'is_a':
            isa_id = re.search('GO:([0-9]*)', line).group(1)
            term.is_a.append(str(isa_id))
        elif header == 'synonym':
            term.synonym.extend(re.search('\"(.*?)\"', line).group(1))
        elif header == 'is_obsolete':
            go_dict.pop(term.id)
            start_reading = False
    f.close()
    return go_dict


def find_children(dict):
    """Adds parent/child information to goterm objects"""
    for i in dict:
        goterm = dict[i]
        if goterm.is_a:
            for parent_id in goterm.is_a:
                parent = dict[parent_id]
                if goterm.id not in parent.children:
                    parent.children.append(goterm)
                    goterm.parents.append(parent)


def parse_line(line, db, go_dict, go_tree, results):
    """Extracts GO terms from mapped read and updates similarity score"""

    # check to see which database we're using
    if line[0:2] == 'sp':
        uniprot_pattern = re.compile('sp\|.*?\|(\w+?),')
    else:
        uniprot_pattern = re.compile('(\w+?),')
    uniprot_id = uniprot_pattern.match(line)

    reference_id = re.search('<source>(.*?)</source>', line)
    ac_id = re.search('UniProtKB:(.*?);', line)

    # break out if we can't find a match
    if not (uniprot_id and reference_id and ac_id):
        return

    # extract groups
    ac_id = ac_id.group(1)
    uniprot_id = uniprot_id.group(1)
    reference_id = reference_id.group(1)
    tar_string = "SELECT * FROM uniprot WHERE name='" + str(uniprot_id + "'")
    ac_string = "SELECT * FROM ac WHERE ac='" + str(ac_id + "'")

    # try to retrieve uniprot entries from target and source
    try:
        target_entry = db.query(tar_string).next()
        ac_entry = db.query(ac_string).next()
        sor_string = "SELECT * FROM uniprot WHERE name='" \
                     + str(ac_entry['name'] + "'")
        source_entry = db.query(sor_string).next()
    except StopIteration:
        return

    initialize_source(results, reference_id)

    # start creating subgraphs
    source_nodes = []
    target_nodes = []
    for i in target_entry["go_terms"].split(","):
        if i in go_tree:
            get_subgraph_nodes(go_tree, target_nodes, i)
    for i in source_entry["go_terms"].split(","):
        if i in go_tree:
            get_subgraph_nodes(go_tree, source_nodes, i)

    # start creating trees
    source_tree = go_tree.subgraph(source_nodes)
    target_tree = go_tree.subgraph(target_nodes)

    # create sets from nodes in these trees
    source_set = set(source_tree.nodes())
    target_set = set(target_tree.nodes())

    # now find where they intersect
    intersection = source_set.intersection(target_set)
    union = source_set.union(target_set)
    source_c, source_b, source_f = 0, 0, 0
    target_b, target_c, target_f = 0, 0, 0
    intersection_c, intersection_b, intersection_f = 0, 0, 0
    union_c, union_b, union_f = 0, 0, 0

    # count number of cell/mol/bio go_terms from source
    for i in source_set:
        go_type = go_dict[i].namespace.rstrip()
        if go_type == "cellular_component":
            source_c += 1
        elif go_type == "molecular_function":
            source_b += 1
        else:
            source_f += 1

    # count number of cell/mol/bio go_terms from target
    for i in target_set:
        go_type = go_dict[i].namespace.rstrip()
        if go_type == "cellular_component":
            target_c += 1
        elif go_type == "molecular_function":
            target_b += 1
        else:
            target_f += 1

    # count cell/mol/bio go terms from intersection of target/source
    for i in intersection:
        go_type = go_dict[i].namespace.rstrip()
        if go_type == "cellular_component":
            intersection_c += 1
        elif go_type == "molecular_function":
            intersection_b += 1
        else:
            intersection_f += 1

    # count cell/mol/bio go terms from union of target/source
    for i in union:
        go_type = go_dict[i].namespace.rstrip()
        if go_type == "cellular_component":
            union_c += 1
        elif go_type == "molecular_function":
            union_b += 1
        else:
            union_f += 1

    # update cellular component similarity score
    if intersection_c:
        c_ratio = float(intersection_c) / max(1, union_c)
        c_dist = float(union_c - intersection_c) / max(1, union_c)
        results['c_count'] += 1
        results['c_similarity'] += c_ratio
        results['c_distance'] += c_dist

    # update biological process similarity score
    if intersection_b:
        b_ratio = float(intersection_b) / max(1, union_b)
        b_dist = float(union_b - intersection_b) / max(1, union_b)
        results['b_count'] += 1
        results['b_similarity'] += b_ratio
        results['b_distance'] += b_dist

    # update molecular function similarity score
    if intersection_f:
        f_ratio = float(intersection_f) / max(1, union_f)
        f_dist = float(union_f - intersection_f) / max(1, union_f)
        results['f_count'] += 1
        results['f_similarity'] += f_ratio
        results['f_distance'] += f_dist


def create_go_tree(go_dict):
    """Create tree containing go term parent/child relations"""
    go_tree = nx.DiGraph()
    for i in go_dict:
        goterm = go_dict[i]
        go_tree.add_node(goterm.id)
        for parent in goterm.is_a:
            go_tree.add_edge(parent, goterm.id)
    return go_tree


def remap_dict(dict):
    """create secondary dictionary where AC-ids are used as the keys instead"""
    rmdict = {}
    for key in dict:
        term = dict[key]
        acs = (term.ac.replace("\n","")).split(';')
        acs.pop()
        for ac in acs:
            rmdict[(ac.replace(' ', ''))] = term
    return rmdict


def initialize_source(dict, source):
    if dict.get(source):
        return

    dict[source] = \
              {'b_count':0, 'b_union':0, 'b_intersection':0, 'b_similarity':0.0,
               'c_count':0, 'c_union':0, 'c_intersection':0, 'c_similarity':0.0,
               'c_distance':0, 'b_distance':0, 'f_distance':0, 'read_counter':0,
               'f_count':0, 'f_union':0, 'f_intersection':0, 'f_similarity':0.0,
               'unmapped':0}
    dict['references'].append(source)


def get_subgraph_nodes(tree, list, node):
    """Creates a subgraph which includes target node and its ancestors"""
    list.append(node)
    parents = [i for i in tree.predecessors_iter(node)]
    while parents:
        p = parents.pop()
        if p in list:
            continue

        list.append(p)
        for i in tree.predecessors_iter(p):
            parents.append(i)


# override error class for argparse
class HelpParser(argparse.ArgumentParser):
    def error(self,message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


# parse through arguments
def create_parser():
    parser = HelpParser()
    parser.add_argument('mapper',
                        help='type of read mapper ("paladin" or "bwa")')
    parser.add_argument('sam_loc',
                        help='location of sam file')

    return parser.parse_args()

if __name__ == "__main__":

    parser = create_parser()

    # required files
    sam_type = parser.mapper
    sam_loc = parser.sam_loc
    base = '/home/genome/airjordan/paladin_data/similarity_references'
    go_file = "/home/genome/airjordan/paladin_data/similarity_references/go.obo"
    godict = create_godict(go_file)
    script_loc = "/home/genome/airjordan/paladin_data/" \
                 "similarity_references/listMappedCDS.py"

    # associate parent/child relationships with GoTerms
    find_children(godict)

    # load databases
    go_tree = create_go_tree(godict)
    db = dataset.connect('sqlite:///' + base + '/paladin.db')
    uniprot_table = db.get_table('uniprot')
    ac_table = db.get_table('ac')

    # results that we'll analyze later
    results = {'b_count':0, 'b_union':0, 'b_intersection':0, 'b_similarity':0.0,
               'c_count':0, 'c_union':0, 'c_intersection':0, 'c_similarity':0.0,
               'c_distance':0, 'b_distance':0, 'f_distance':0, 'references':[],
               'f_count':0, 'f_union':0, 'f_intersection':0, 'f_similarity':0.0}

    # run listMapped script and begin to pipe results into script
    counter = 0
    p = Popen([script_loc, sam_loc, "-r", sam_type], stdout=PIPE, bufsize=1)
    with p.stdout:
        for line in iter(p.stdout.readline, b''):
            parse_line(line, db, godict, go_tree, results)
    p.wait()

    # report average of results
    b_sim = str(float(results['b_similarity'] / results['b_count']))
    c_sim = str(float(results['c_similarity'] / results['c_count']))
    f_sim = str(float(results['f_similarity'] / results['f_count']))
    total_sim = str(float(results['b_similarity'] + results['c_similarity'] +
                   results['f_similarity']) /
                    (results['f_count'] + results['c_count'] +
                     results['b_count']))
    f_dist = str(float(results['f_distance'] / results['f_count']))
    c_dist = str(float(results['c_distance'] / results['c_count']))
    b_dist = str(float(results['b_distance'] / results['b_count']))
    total_dist = str(float(results['f_distance'] + results['c_distance'] +
                          results['b_distance']) /
                    (results['f_count'] + results['c_count'] +
                     results['b_count']))
    name = sam_loc.replace(".sam", "")
    base = re.sub(".*/", "", name)
    name = re.sub(".*/", "", name) + "_total"

    print("Name,b_sim,f_sim,c_sim,total_sim,b_dist,f_dist,c_dist,total_dist")
    print(name + "," + b_sim + "," + f_sim + "," + c_sim + "," + total_sim +
          "," + b_dist + "," + f_dist + "," + c_dist + "," + total_dist)

