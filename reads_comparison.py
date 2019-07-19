import numpy as np
import distance

def mgnify_preparation(filename_mgnify):
    annotated_reads_mgnify = {}
    with open(filename_mgnify, 'r') as file_in:
        for line in file_in:
            line = line.strip().split('\t')
            if line[0][0] == '#': continue
            name = line[0].split('-')[0]
            if name not in annotated_reads_mgnify:
                annotated_reads_mgnify[name] = [line[len(line)-1]]
            else:
                if len(annotated_reads_mgnify[name][0]) < len(line[len(line)-1]):
                    annotated_reads_mgnify[name][0] = line[len(line)-1]

    ############## Make MGNIFY beautiful
    pair = 0
    for read in annotated_reads_mgnify:
        tax = annotated_reads_mgnify[read][pair].split(';')
        if len(tax[len(tax) - 1]) > 3: continue
        for index in range(len(tax) - 1, 0, -1):
            if len(tax[index]) != 3: break
        short_tax = ';'.join(tax[:index])
        annotated_reads_mgnify[read][pair] = short_tax
    return annotated_reads_mgnify


def mg_rast_preparation(filename_mgrast_annotation, filename_mgrast_clusters):
    classifications = ['sk__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__', 't__']
    annotations_mgrast = {}
    with open(filename_mgrast_annotation, 'r') as file_annotation:
        for line in file_annotation:
            line = line.strip().split('\t')
            annotations_mgrast[line[1]] = line[5]
    print('annotated clusters mgrast', len(annotations_mgrast))

    clustes = {}
    with open(filename_mgrast_clusters, 'r') as file_clustes:
        for line in file_clustes:
            line = line.strip().split('\t')
            cur_cluster = np.unique(line[1].split(','))
            clustes[line[0]] = list(cur_cluster)

    kol_no = 0
    annotated_reads_mgrast = {}
    with open('out_mgrast.txt', 'w') as file_out_mgrast:
        for cluster in clustes:
            items = clustes[cluster]
            if cluster in annotations_mgrast:
                taxonomy = annotations_mgrast[cluster].split(';')
                name = []
                for i, j in zip(taxonomy, range(len(taxonomy))):
                    if i == '-': break
                    if j == 1:
                        name.append(';'.join(['k__', classifications[j] + i,]))
                    else:
                        name.append(classifications[j] + i)
                cur_taxonomy = ';'.join(name)
                for read in items:
                    read_name = read.split('_')[0]
                    file_out_mgrast.write('\t'.join([read_name, cur_taxonomy]) + '\n')
                    annotated_reads_mgrast[read_name] = cur_taxonomy
            else:
                kol_no += 1
    print('non annotated clusters mgrast', kol_no)
    return annotated_reads_mgrast


def simple_statistics(annotated_reads_mgrast, annotated_reads_mgnify):
    print('Reads MGNIFY: ' + str(len(annotated_reads_mgnify)))
    print('Reads MG-RAST: ' + str(len(annotated_reads_mgrast)))
    common_reads = set(annotated_reads_mgrast).intersection(set(annotated_reads_mgnify))
    print('Common reads: ' + str(len(common_reads)))

    with open('comparison_common.txt', 'w') as file_comparison:
        file_comparison.write('\t'.join(['Read_name', 'MGNIFY', 'MG-RAST']) + '\n')
        for read in common_reads:
            file_comparison.write('\t'.join([read, annotated_reads_mgnify[read][0], annotated_reads_mgrast[read]]) + '\n')

    reads_only_mgnify = set(annotated_reads_mgnify).difference(set(annotated_reads_mgrast))
    print('only MGNIFY', len(reads_only_mgnify))
    with open('only_MGNIFY.txt', 'w') as file_mgnify:
        for read in list(reads_only_mgnify):
            file_mgnify.write('\t'.join([read, annotated_reads_mgnify[read][0]]) + '\n')

    reads_only_mgrast = set(annotated_reads_mgrast).difference(set(annotated_reads_mgnify))
    print('only MG-RAST', len(reads_only_mgrast))
    with open('only_MG-RAST.txt', 'w') as file_mgrast:
        for read in list(reads_only_mgrast):
            file_mgrast.write('\t'.join([read, annotated_reads_mgrast[read]]) + '\n')

filename_mgrast_clusters = 'amp_fq_test.440.cluster.rna97.mapping' # mgrast 'SRR6367227/amp_fq_test.440.cluster.rna97.mapping'
filename_mgrast_annotation = 'amp_fq_test.450.rna.expand.lca'
filename_mgnify = 'ERR1255869_MERGED_FASTQ_SSU.fasta.mseq'# 'SRR6367227/SRR6367227_MERGED_FASTQ_SSU.fasta.mseq'  #  # mgnify

annotated_reads_mgrast = mg_rast_preparation(filename_mgrast_annotation, filename_mgrast_clusters)
annotated_reads_mgnify = mgnify_preparation(filename_mgnify)


def check_unclassified(cur_unclassified, annotated_reads_mgrast, level):
    annotated, notannotated, uncultured = [[], [], []]
    for read in cur_unclassified:
        tax = annotated_reads_mgrast[read].split(';')
        while level < len(tax):
            if len(tax[level].split('unclassified')) == 1:  # NON unclassified = annotated
                annotated.append(read)
                break
            level += 1
        if level == len(tax):
            notannotated.append(read)

        if len(tax[len(tax)-1].split('uncultured bacterium')) > 1:
            uncultured.append(read)
    print('MG-Rast unclassified ', len(cur_unclassified))
    print('----> annotated', len(annotated))
    print('----> notannotated', len(notannotated))
    print('----------> uncultured', len(uncultured))
    return annotated, uncultured, notannotated


def print_dop_file(name, reads, level, num_read, annotated_reads_mgnify, annotated_reads_mgrast):
    with open(str(level)+name, 'w') as file_1:
        for read in reads:
            s = read + '\t' + annotated_reads_mgnify[read][num_read] + ' <---> ' + annotated_reads_mgrast[read]
            file_1.write(s + '\n')


def check_absolutely_different(absolutely_different, annotated_reads_mgrast, annotated_reads_mgnify, level, num_read):
    """
    Function detects and changes situations :
    "sk__Bacteria;k__;p__Actinobacteria;c__Actinobacteria <--> 
                                                       sk__Bacteria;k__;p__Actinobacteria;c__Actinobacteria (class)"
    :param absolutely_different: 
    :param annotated_reads_mgrast: 
    :param annotated_reads_mgnify: 
    :return: 
    """
    print_dop_file('_abs_different', absolutely_different, level, num_read, annotated_reads_mgnify, annotated_reads_mgrast)
    really_abs_dif = []
    for read in absolutely_different:
        name_mgnify = annotated_reads_mgnify[read][num_read].split(';')
        name_mgrast = annotated_reads_mgrast[read].split(';')
        if len(name_mgrast[level].split(name_mgnify[level])) <= 1:
            #name_mgrast[level] = name_mgnify[level]
            #annotated_reads_mgrast[read] = ';'.join(name_mgrast)
        #else: make in if > 1
            really_abs_dif.append(read)
    print('-----> Not involved to each other', len(really_abs_dif))
    print_dop_file('_0really_abs_dif', really_abs_dif, level, num_read, annotated_reads_mgnify,
                   annotated_reads_mgrast)

    class_difference = []
    dist_difference = []
    for read in really_abs_dif:  # check distance between words with (class)
        name_mgrast = annotated_reads_mgrast[read].split(';')[level]
        name_mgrast_spl = name_mgrast.split('(class)')
        name_mgnify = annotated_reads_mgnify[read][num_read].split(';')[level]
        if len(name_mgrast_spl) > 1:
            dist = distance.levenshtein(name_mgnify, name_mgrast_spl[0])
            if dist > 3:  # not similar
                class_difference.append(read)
        else:
            if distance.levenshtein(name_mgnify, name_mgrast) > 3:
                dist_difference.append(read)
    print('-----> Class difference', len(class_difference))
    print('-----> Distance difference', len(dist_difference))
    print_dop_file('_1class_diff', class_difference, level, num_read, annotated_reads_mgnify,
                   annotated_reads_mgrast)
    print_dop_file('_2distance_diff', dist_difference, level, num_read, annotated_reads_mgnify,
                   annotated_reads_mgrast)
    return class_difference + dist_difference


def compare_cur_union(cur_union, annotated_reads_mgrast, annotated_reads_mgnify, num_read):
    """
    Function decides which annotation deeper by number of consistent fields
    :return:
    """
    deeper_mgrast, deeper_mgnify, cur_common = [0 for _ in range(3)]
    cur_common_list = []
    for i in cur_union:
        name_mgnify = annotated_reads_mgnify[i][num_read].split(';')
        name_mgrast = annotated_reads_mgrast[i].split(';')
        if len(name_mgrast) > len(name_mgnify):
            deeper_mgrast += 1
        elif len(name_mgrast) < len(name_mgnify):
            deeper_mgnify += 1
        else:
            cur_common += 1
            cur_common_list.append(read)
    return deeper_mgrast, deeper_mgnify, cur_common, cur_common_list


def sk_stop_analysis(level_mgrast, level_mgnify, annotated_reads_mgrast, annotated_reads_mgnify):
    """
    Analysis of MG-RAST and MGnify in superkingdom level
    :param level_mgrast:
    :param level_mgnify:
    :param annotated_reads_mgrast:
    :param annotated_reads_mgnify:
    :return:
    """

    cur_common = set(level_mgnify).intersection(set(level_mgrast))
    level_mgrast = list(set(level_mgrast) - set(cur_common))
    level_mgnify = list(set(level_mgnify) - set(cur_common))
    passed_mgnify, passed_mgrast, stop_mgnify, stop_mgrast = [[] for _ in range(4)]
    with open('Stopped_MGnify.txt', 'w') as file_sm:
        for read in level_mgnify:
            file_sm.write(read + '\n')
    types = ['phylum', 'class', 'order', 'family', 'genus', 'species', 'extra1', 'extra2']
    type_levels_mgnify, type_levels_mgrast = [[0 for _ in range(8)] for _ in range(2)]
    # 0 - phylum, 1 - class, 2 - order, 3 - family, 4 - genus, 5 - species, ...

    passed_mgnify.append(len(level_mgrast))
    stop_mgnify.append(0)
    print('List from MG-RAST. What happened with MGnify', len(level_mgrast))
    for read in level_mgrast:
        type_levels_mgnify[len(annotated_reads_mgnify[read][0].split(';')) - 3] += 1
    for i in range(len(type_levels_mgnify)):
        print(types[i] + ' :' + str(type_levels_mgnify[i]))
        stop_mgnify.append(type_levels_mgnify[i])
    passed_mgnify += list(len(level_mgrast) - np.cumsum(type_levels_mgnify[:-1]))

    passed_mgrast.append(len(level_mgnify))
    stop_mgrast.append(0)
    print('List from MGnify. What happened with MG-RAST', len(level_mgnify))
    for read in level_mgnify:
        type_levels_mgrast[len(annotated_reads_mgrast[read].split(';')) - 3] += 1
    for i in range(len(type_levels_mgrast)):
        print(types[i] + ' :' + str(type_levels_mgrast[i]))
        stop_mgrast.append(type_levels_mgrast[i])
    passed_mgrast += list(len(level_mgnify) - np.cumsum(type_levels_mgrast[:-1]))
    #a = [annotated_reads_mgrast[read] for read in level_mgnify if len(annotated_reads_mgrast[read].split(';')) - 3 == 4]
    #with open('bacteria.txt', 'w') as file:
    #    file.write('\n'.join(a))
    return passed_mgnify, passed_mgrast, stop_mgnify, stop_mgrast



def report(report_left_reads, report_stop_mgnify, report_stop_mgrast, report_stop_both, report_abs_dif,
           report_uncultured, report_notannotated, passed_mgnify_sk, passed_mgrast_sk, stop_mgnify_sk, stop_mgrast_sk):

    # Main part
    Fields = ['All reads', 'Superkingdom', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    colours = {'left_reads': '#00FF00', 'stop_mgnify': '#FFFF00', 'stop_mgrast': '#0000FF', 'stop_both': '#8B008B',
               'different': '#FF1493', 'unclassified': '#FF0000', 'uncultured': '#FF0000', 'mg-rast_passed': "#FFFF00",
               "mgnify_passed": "#00CED1"}
    for index in range(len(report_left_reads)):
        if report_left_reads[index] == 0: break
        # left reads
        print(' '.join(
            [Fields[index], '[' + str(report_left_reads[index]) + ']', Fields[index + 1], colours['left_reads']]))
        # stop MGnify
        if report_stop_mgnify[index] - report_stop_both[index]:
            print(' '.join([Fields[index], '[' + str(report_stop_mgnify[index] - report_stop_both[index]) + ']',
                            '{' + Fields[index + 1][:3] + '}', 'Stop_MGnify', colours['stop_mgnify']]))
        # stop MG-RAST
        if report_stop_mgrast[index] - report_stop_both[index]:
            print(' '.join([Fields[index], '[' + str(report_stop_mgrast[index] - report_stop_both[index]) + ']',
                            '{' + Fields[index + 1][:3] + '}', 'Stop_MG-RAST', colours['stop_mgrast']]))
        # stop both
        if report_stop_both[index]:
            print(' '.join([Fields[index], '[' + str(report_stop_both[index]) + ']',
                            '{' + Fields[index + 1][:3] + '}', 'Stop_both', colours['stop_both']]))
        # Different
        if report_abs_dif[index]:
            print(' '.join([Fields[index], '[' + str(report_abs_dif[index]) + ']',
                            '{' + Fields[index + 1][:3] + '}', 'Different', colours['different']]))
        if report_uncultured[index]:
            print(' '.join([Fields[index], '[' + str(report_uncultured[index]) + ']',
                            '{' + Fields[index + 1][:3] + '}', 'Uncultured', colours['uncultured']]))
        if report_notannotated[index]:
            print(' '.join([Fields[index], '[' + str(report_notannotated[index]) + ']',
                            '{' + Fields[index + 1][:3] + '}', 'Unclassified', colours['unclassified']]))
        print('\n')

    # Not annotated by MGnify
    print(" ".join(["{" + Fields[1][:3] + "}", "Stop_MGnify", "["+str(passed_mgrast_sk[0])+"]", "{"+Fields[2][:3]+"}", "MG-RAST_passed_"+Fields[2], colours["mg-rast_passed"]]))
    for index in range(1, len(passed_mgrast_sk)-2):
        if passed_mgrast_sk[index]:
            print(" ".join(["{"+Fields[index+1][:3]+"}", "MG-RAST_passed_"+Fields[index+1], "["+str(passed_mgrast_sk[index])+"]", "{"+Fields[index+2][:3]+"}",
                            "MG-RAST_passed_" + Fields[index+2], colours['mg-rast_passed']]))
        if stop_mgrast_sk[index]:
            print(" ".join(["{" + Fields[index+1][:3] + "}", "MG-RAST_passed_" + Fields[index+1],
                            "[" + str(stop_mgrast_sk[index]) + "]", "{" + Fields[index+1][:3] + "}",
                            "MG-RAST_stop_" + Fields[index+2], "#FF0000"]))
    print('\n')
    # Not annotated by MG-RAST
    print(" ".join(["{" + Fields[1][:3] + "}", "Stop_MG-RAST", "["+str(passed_mgnify_sk[0])+"]", "{"+Fields[2][:3]+"}", "MGnify_passed_"+Fields[2], colours["mgnify_passed"]]))
    for index in range(1, len(passed_mgnify_sk)-2):
        if passed_mgnify_sk[index]:
            print(" ".join(["{"+Fields[index+1][:3]+"}", "MGnify_passed_"+Fields[index+1], "["+str(passed_mgnify_sk[index])+"]", "{"+Fields[index+2][:3]+"}",
                            "MGnify_passed_" + Fields[index+2], colours['mgnify_passed']]))
        if stop_mgnify_sk[index]:
            print(" ".join(["{" + Fields[index+1][:3] + "}", "MGnify_passed_" + Fields[index+1],
                            "[" + str(stop_mgnify_sk[index]) + "]", "{" + Fields[index+1][:3] + "}",
                            "MGnify_stop_" + Fields[index+2], "#FF0000"]))


# statistics by reads
simple_statistics(annotated_reads_mgrast, annotated_reads_mgnify)
common_reads = set(annotated_reads_mgrast).intersection(set(annotated_reads_mgnify))
print('========== Comparison ==========')
print('---> SUPER KINGDOM: ' + str(len(common_reads)) + ' <----> ' + str(len(common_reads)))
num_read = 0
levels = ['      KINGDOM', '    PHYLUM', '     CLASS', '    ORDER', '    FAMILY', '    GENUS', '     SPECIES', '    T']
exclude_set = []
report_left_reads, report_stop_mgnify, report_stop_mgrast, report_stop_both = [[] for _ in range(4)]
report_abs_dif, report_notannotated, report_uncultured = [[0] for _ in range(3)]
for level, num_level in zip(levels, range(len(levels))):
    an_mr, an_mn = [0, 0]
    number_of_comparisons, number_equal = [0, 0]
    num_different_annotations_on_level, num_unclassified_annotations_on_level = [0, 0]
    level_mgrast, level_mgnify = [[], []]
    absolutely_different, cur_unclassified = [[], []]
    for read in common_reads:
        if read in exclude_set: continue
        flag_mr, flag_mn = [False, False]
        tax_mgrast = annotated_reads_mgrast[read].split(';')
        tax_mgnify = annotated_reads_mgnify[read][num_read].split(';')
        if len(tax_mgrast) > num_level + 1:  # MG-RAST annotated to level
            an_mr += 1
            cur_mr = tax_mgrast[:num_level+2]
            flag_mr = True
        else:
            level_mgrast.append(read)  # stopped on previous level

        if len(tax_mgnify) > num_level + 1:  # MGNIFY annotated to level
            an_mn += 1
            cur_mn = tax_mgnify[:num_level+2]
            flag_mn = True
        else:
            level_mgnify.append(read)  # stopped on previous level

        if flag_mr and flag_mn:  # if read annotated to level
            number_of_comparisons += 1
            if ';'.join(cur_mr) == ';'.join(cur_mn):  # equal
                number_equal += 1
            else:  # not equal
                dop_mgrast = ';'.join(cur_mr)
                if len(dop_mgrast.split('unclassified')) > 1:  # because of "unclassified"
                    num_unclassified_annotations_on_level += 1
                    cur_unclassified.append(read)
                else:  # because of absolutely different trees
                    absolutely_different.append(read)
                    num_different_annotations_on_level += 1
    print('prev Level mgnify ', level_mgnify[:10])
    print('prev Level mgrast ', level_mgrast[:10])

    cur_union = set(level_mgrast).union(set(level_mgnify))
    # calc number of previous level
    deeper_mgrast, deeper_mgnify, cur_common, cur_common_list = compare_cur_union(cur_union, annotated_reads_mgrast,
                                                                                  annotated_reads_mgnify, num_read)
    # extra info about SK
    if num_level == 0:
        passed_mgnify_sk, passed_mgrast_sk, stop_mgnify_sk, stop_mgrast_sk = sk_stop_analysis(level_mgrast,
                                                        level_mgnify, annotated_reads_mgrast, annotated_reads_mgnify)

    print('---> ' + level + ': ' + str(an_mn) + '  <----> ' + str(an_mr))
    print('left on previous level: ' + str(len(level_mgnify)) + ' <--> ' + str(
        len(level_mgrast)) + ' from them common ' + str(cur_common))
    print('----> deeper MGNIFY: ' + str(deeper_mgnify))
    print('----> deeper MG-RAST: ' + str(deeper_mgrast))

    annotated, uncultured, notannotated = check_unclassified(cur_unclassified, annotated_reads_mgrast, num_level+1)
    report_notannotated.append(len(notannotated))
    report_uncultured.append(len(uncultured))
    new_abs_dif = check_absolutely_different(absolutely_different, annotated_reads_mgrast, annotated_reads_mgnify,
                                             num_level+1, num_read)
    report_abs_dif.append(len(new_abs_dif))
    print('+ num absolutely different ' + str(len(new_abs_dif)) +
          '(by name:' + str(num_different_annotations_on_level) + ')')
    old_size_excluded = len(exclude_set)
    exclude_set += list(cur_union)  # exclude
    exclude_set += new_abs_dif  # exclude absolutely different annotations
    exclude_set += list(notannotated)  # exclude not annotated
    exclude_set += uncultured  # exclude uncultured
    print('Size of excluded ' + str(old_size_excluded) + ' + ' + str(len(exclude_set)-old_size_excluded) +
          ' = ' + str(len(exclude_set)))

    print('number equal: ' + str(number_equal) + ' from ' + str(number_of_comparisons))
    with open(str(num_level) + '.txt', 'w') as file_dop:
        for read in new_abs_dif:
            file_dop.write(read + '\t' + annotated_reads_mgnify[read][num_read] +
                           ' <--> ' + annotated_reads_mgrast[read] + '\n')
    # Report
    report_left_reads.append(number_of_comparisons)
    report_stop_mgnify.append(len(level_mgnify))
    report_stop_mgrast.append(len(level_mgrast))
    report_stop_both.append(cur_common)

print('Copy all to http://sankeymatic.com/build/')
report(report_left_reads, report_stop_mgnify, report_stop_mgrast, report_stop_both, report_abs_dif,
       report_uncultured, report_notannotated, passed_mgnify_sk, passed_mgrast_sk, stop_mgnify_sk, stop_mgrast_sk)
