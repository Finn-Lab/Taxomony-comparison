from anytree import Node, RenderTree
from anytree.exporter import DotExporter
#filename_mgrast = 'test_mgrast.tsv'
#filename_mgnify = 'test_mgnify.tsv'
filename_mgrast = 'amp_fq_test.700.annotation.lca.abundance'  #'mgrast.abundance'
filename_mgnify = 'ERR1255869_MERGED_FASTQ_SSU.fasta.mseq.tsv'  # 'mgnify.tsv'

def create_tree_mgrast(filename_mgrast):
    # MG-RAST
    prefixes = ['sk__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__', 'extra__']
    taxonomy = {}
    name_root = ""
    with open(filename_mgrast, 'r') as file_in:
        for line in file_in:
            line = line.strip().split('\t')
            groups = line[0].split(';')
            number = int(line[1])

            for index in range(0, len(groups)):
                if name_root == "":  # make root
                    name_root = prefixes[0] + groups[0]
                    taxonomy[name_root] = Node(name_root, edge=0)
                else:
                    pred_name = prefixes[index - 1] + groups[index - 1]
                    cur_name = prefixes[index] + groups[index]
                    if groups[index] == '-':
                        taxonomy[pred_name].edge = number
                        break
                    if cur_name not in taxonomy:
                        taxonomy[cur_name] = Node(cur_name, parent=taxonomy[pred_name], edge=0)
                if index == len(groups)-1:
                    taxonomy[cur_name].edge = number

            # set numbers to parents
            for index in range(len(groups)-1, -1, -1):
                if groups[index] == '-': continue
                parent_name = prefixes[index] + groups[index]
                if index != len(groups)-1:
                    if groups[index+1] == '-': continue
                else: continue
                taxonomy[parent_name].edge += number
    return taxonomy


def create_tree_mgnify(filename_mgnify):
    # MGNIFY
    name_root = ""
    taxonomy_mf = {}
    with open(filename_mgnify, 'r') as file_in:
        for line in file_in:
            line = line.strip().split('\t')
            groups = line[2].split(';')
            number = int(float(line[1]))

            for index in range(0, len(groups)):
                if name_root == "":  # make root
                    name_root = groups[0]
                    taxonomy_mf[name_root] = Node(name_root, edge=0)
                else:
                    pred_name = '.'.join(groups[0:index])
                    cur_name = '.'.join(groups[0:index+1])
                    if cur_name not in taxonomy_mf:
                        taxonomy_mf[cur_name] = Node(groups[index], parent=taxonomy_mf[pred_name], edge=0)
                if index == len(groups) - 1:
                    cur_name = '.'.join(groups[0:index+1])
                    taxonomy_mf[cur_name].edge = number
            # set numbers to parents
            for index in range(len(groups)-2, -1, -1):
                parent_name = '.'.join(groups[0:index+1])
                taxonomy_mf[parent_name].edge += number
    return taxonomy_mf

superkingdom = 'sk__Bacteria'

# MG-RAST
taxonomy = create_tree_mgrast(filename_mgrast)
print('MG-RAST')
DotExporter(taxonomy[superkingdom]).to_picture("mgrast.png")
#DotExporter(taxonomy['sk__Eukaryota']).to_dotfile("mgrast.dot")
for pre, fill, node in RenderTree(taxonomy[superkingdom]):
    print("%s%s [%s]" % (pre, node.name, str(node.edge)))

# MGNIFY
print('MGNIFY')
taxonomy_mf = create_tree_mgnify(filename_mgnify)
DotExporter(taxonomy_mf[superkingdom]).to_picture("mgnify.png")
for pre, fill, node in RenderTree(taxonomy_mf[superkingdom]):
    print("%s%s [%s]" % (pre, node.name, str(node.edge)))


def get_child(cur_tree):
    children = []
    for index_mn in range(len(cur_tree)):
        children += cur_tree[index_mn].children
    return children


tree_mgnify = taxonomy_mf[superkingdom]
tree_mgrast = taxonomy[superkingdom]

print('\t'.join(['NAME', 'MGNIFY', 'MG-RAST', 'common_name',
                 'common_MGNIFY', 'common_MG-RAST', 'different_MGNIFY', 'different_MG-RAST']))

# TODO ROOT comparison if > 1 sk
print('\t'.join([superkingdom, str(tree_mgnify.edge), str(tree_mgrast.edge)]))

# ====== skipping kingdom in MGNIFY =======

kingdoms = taxonomy_mf[superkingdom].children
tree_mgnify = get_child(kingdoms)
tree_mgrast = taxonomy[superkingdom].children

levels = ['sk', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'extra']
for level in range(2, 6, 1):  # TODO when > 1 super kingdom in taxonomy

    counts_mn = [tree_mgnify[index].edge for index in range(len(tree_mgnify))]
    counts_mr = [tree_mgrast[index].edge for index in range(len(tree_mgrast))]

    names_mgnify = [tree_mgnify[index].name for index in range(len(tree_mgnify))]
    names_mgrast = [tree_mgrast[index].name for index in range(len(tree_mgrast))]
    number_mgnify, number_mgrast = [0, 0]
    for name, index_mn in zip(names_mgnify, range(len(names_mgnify))):
        if name in names_mgrast:
            index_mg = [index for index in range(len(names_mgrast)) if name == names_mgrast[index]][0]
            number_mgnify += tree_mgnify[index_mn].edge
            number_mgrast += tree_mgrast[index_mg].edge
            print('\t'.join(['' for _ in range(3)] +
                            [name, str(tree_mgnify[index_mn].edge), str(tree_mgrast[index_mg].edge)]))

    print('\t'.join([levels[level], str(sum(counts_mn)), str(sum(counts_mr)), '', str(number_mgnify), str(number_mgrast),
                     str(sum(counts_mn) - number_mgnify), str(sum(counts_mr) - number_mgrast)]))

    children_mgnify_trees = get_child(tree_mgnify)
    children_mgrast_trees = get_child(tree_mgrast)
    tree_mgnify = children_mgnify_trees
    tree_mgrast = children_mgrast_trees
print(0)


