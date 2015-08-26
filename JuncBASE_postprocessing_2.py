#!usr/bin/env python

#
# this script will pull "junction sets" and their normalized read counts from _grouped_AS_exclusion_inclusion_counts_lenNorm.txt
#

import sys
import optparse 
import os
import pdb
import re
import shutil

#############
# CONSTANTS #
#############

COLORS_SMALL = ["0,0,0", "0,0,255", "0,100,0", "0,206,209", "255,0,0", "255,140,0", "148,0,211", "176,48,96", "205,92,92"] 

#################
# END CONSTANTS #
#################


class OptionParser(optparse.OptionParser):
    """
    Adding a method for required arguments.
    Taken from:
    http://www.python.org/doc/2.3/lib/optparse-extending-examples.html
    """
    def check_required(self, opt):
        option = self.get_option(opt)

        # Assumes the option's 'default' is set to None!
        if getattr(self.values, option.dest) is None:
            print "%s option not supplied" % option
            self.print_help()
            sys.exit(1)

def main():
	
    opt_parser = OptionParser()
   
    opt_parser.add_option("-i",
                          dest="input",
                          type="string",
                          help="input file",
                          default=None)
    opt_parser.add_option("-e",
                          dest="entropy",
                          type="string",
                          help="entropy_scores.txt with K/N",
                          default=None)
    opt_parser.add_option("--e_alt",
                          dest="e_alt",
                          action="store_true",
                          help="-e input is list of known introns (\t)",
                          default=False)
    opt_parser.add_option("-o",
                          dest="output",
                          type="string",
                          help="output file",
                          default=None)
    opt_parser.add_option("-b",
                          dest="bed",
                          type="string",
                          help="output bed file",
                          default=None)
    opt_parser.add_option("-c",
                          dest="complex",
                          type="string",
                          help="output file of complex events",
                          default=None)
    opt_parser.add_option("-n",
                          dest="name",
                          type="string",
                          help="track name (opt)",
                          default=False)
    opt_parser.add_option("--jcn_seq_len",
                          dest="jcn_seq_len",
                          type="int",
                          help="jcn_seq_len",
                          default=20)
    
    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-i")
    opt_parser.check_required("-e")
    opt_parser.check_required("-o")
    opt_parser.check_required("-c")
    opt_parser.check_required("-b")

    jcn_len = options.jcn_seq_len / 2
    name = "JuncBASE events"
    if options.name:
        name = options.name
        
    # First read in K/N annotation for events

    juncs_Known = []                     # (known juncs only)

    e_alt = options.e_alt
    e = open(options.entropy)
    elines = e.readlines()
    for eline in elines:
        eline = formatLine(eline).split("\t")
        if e_alt:
            j = "%s:%s-%s" % (fixChr(eline[0], "+"), eline[1], eline[2]);
            juncs_Known.append(j)
        else:
            if eline[2] == "K":
                junc = eline[0].split("_")
                j = "%s:%s-%s" % ("_".join(junc[0:-2]), junc[-2], junc[-1]);
                juncs_Known.append(j);
    e.close()
        

    # Next read in and parse all juncEvents

    groups = {}                     # groups[ASclass][groupID][eventID][gene/chr/str/exclJ/inclJ/exclE/inclE/data(list)]
    allLines = {}                   # allLines[eventID] = line
    
    f = open(options.input)
    header = f.readline()
    header = formatLine(header)
    samples = "\t".join(header.split("\t")[11:-2])
    lines = f.readlines()
    for line in lines:
        line = formatLine(line).split("\t")
        groups = addJuncEvent2groups(line, groups)
        allLines[line[-1]] = "\t".join(line)
    f.close
    
    # Next find junction sets

    juncSets, complexLines = getJuncSets(groups, juncs_Known)        # juncSets[gene][ASclass][groupID][juncSetID][chrom/strand/KN/junctions/exclusionJuncs/readCounts(list)]

    # Now print output and make bed file

    o = open(options.output, "w")
    b = open(options.bed, "w")
    c = open(options.complex, "w")

    print >>o, "gene\tASclass\tgroupID\tjuncSetID\tchr\tstrand\tknown_or_novel\tjunctions\texclusion_junctions\t%s" % samples
    print >>b, "track name=\'%s\' itemRgb=\'On\'" % name
    print >>c, header
    for gn in juncSets:
        for asclass in juncSets[gn]:
            for g in juncSets[gn][asclass]:
                i = 0
                for j in juncSets[gn][asclass][g]:
                    j_list = j.split("-")
                    gg = g + "_" + j_list[0]
                    jj = j_list[1]
                    print >>o, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (gn, asclass, gg, jj, juncSets[gn][asclass][g][j]["chrom"], juncSets[gn][asclass][g][j]["strand"],
                                                                           juncSets[gn][asclass][g][j]["KN"], juncSets[gn][asclass][g][j]["junctions"], 
                                                                           juncSets[gn][asclass][g][j]["exclusionJuncs"], "\t".join(juncSets[gn][asclass][g][j]["data"]))
                    print >>b, makeJuncSetBedLines(juncSets[gn][asclass][g][j], asclass + "-" + gg + "-" + jj, asclass, i, jcn_len),
                    i += 1
                    if i > 8:
                        i = 0

    for ev in complexLines:
        print >>c, allLines[ev]
        
    o.close()
    b.close()

############
# END_MAIN #
############

def addJuncEvent2groups(event, groups):

    ASclass, groupID, eventID, gene, chrom, strand, excl_j, incl_j, excl_e, incl_e, data = \
      event[1], event[-2], event[-1], event[2], event[3], event[4], event[5], event[6], event[7], event[8], event[11:-2]

    if ASclass not in groups:
        groups[ASclass] = {}
    if groupID not in groups[ASclass]:
        groups[ASclass][groupID] = {}
    if eventID not in groups[ASclass]:
        groups[ASclass][groupID][eventID] = {}
    else:
        print >>sys.stderr, "%s is repeated in %s" % (eventID, groupID)

    groups[ASclass][groupID][eventID]["gene"] = gene
    groups[ASclass][groupID][eventID]["chrom"] = chrom
    groups[ASclass][groupID][eventID]["strand"] = strand
    groups[ASclass][groupID][eventID]["excl_j"] = excl_j
    groups[ASclass][groupID][eventID]["incl_j"] = incl_j
    groups[ASclass][groupID][eventID]["excl_e"] = excl_e
    groups[ASclass][groupID][eventID]["incl_e"] = incl_e
    groups[ASclass][groupID][eventID]["data"] = data

    return groups


def getData(data, typ, typ2):
    n_data = []
    for d in data:
        if typ == "counts":
            e, i = d.split(";")
            if typ2 == "excl":
                n_data.append(e)
            elif typ2 == "incl":
                n_data.append(i)
        elif typ == "psi":
            if typ2 == "excl":
                n_data.append(str(1 - float(d))) 
            elif typ2 == "incl":
                n_data.append(d)
    return n_data
    
def getJuncSets(groups, juncs_KN):

    complexLines = []
    juncSets = {}
    
    for cl in groups:
        print >>sys.stderr, "working on %s" % cl

        for g in groups[cl]:
            events = groups[cl][g]
                                           
            j_sets, gene, chrom, strand, eventIDs, complx = getJuncSetsOfGroup(events, cl, juncs_KN)
            
            if complx:
                complexLines.extend(eventIDs)
            else:
                if gene not in juncSets: juncSets[gene] = {}
                if cl not in juncSets[gene]: juncSets[gene][cl] = {}
                juncSets[gene][cl][g] = {}
                for j in j_sets:
                    juncSetID, KN, exclusionJuncs, data = \
                      j_sets[j]["juncSetID"], j_sets[j]["KN"], j_sets[j]["exclusionJuncs"], j_sets[j]["data"]
                      
                    juncSets[gene][cl][g][juncSetID] = {}
                    juncSets[gene][cl][g][juncSetID]["chrom"] = chrom
                    juncSets[gene][cl][g][juncSetID]["strand"] = strand
                    juncSets[gene][cl][g][juncSetID]["KN"] = KN
                    juncSets[gene][cl][g][juncSetID]["junctions"] = j
                    juncSets[gene][cl][g][juncSetID]["exclusionJuncs"] = exclusionJuncs
                    juncSets[gene][cl][g][juncSetID]["data"] = data
                                                
    return juncSets, complexLines

def getJuncSetsOfGroup(events, ASclass, juncs_KN):

    i = 1
    sg = 1
    j_sets = {}
    eventIDs = []
    complx = False

    exons = {}
    exons2 = {}
    incl_exons = {}
    grp_non_alt_ss = ''
    
    for e in events:
        gene, chrom, strand = events[e]["gene"], events[e]["chrom"], events[e]["strand"]
        excl_j, incl_j = events[e]["excl_j"], events[e]["incl_j"]
        eventIDs.append(e)

        if ASclass == "cassette" or ASclass == "coord_cassette":
            if events[e]["excl_j"] in j_sets:
                complx = True
            elif events[e]["incl_e"] in exons:
                complx = True
                
            else:
                exons[events[e]["incl_e"]] = 1

                j_sets[excl_j] = {}
                j_sets[excl_j]["KN"] = "K"
                if isNovel(excl_j, juncs_KN):
                    j_sets[excl_j]["KN"] = "N"
                j_sets[excl_j]["juncSetID"] = "sgrp_%s-set_%s" % (sg, i)
                i += 1
                j_sets[excl_j]["exclusionJuncs"] = incl_j
                j_sets[excl_j]["data"] = getData(events[e]["data"], "counts", "excl")
                    
                j_sets[incl_j] = {}
                j_sets[incl_j]["KN"] = "K"
                if isNovel(incl_j, juncs_KN):
                    j_sets[incl_j]["KN"] = "N"
                j_sets[incl_j]["juncSetID"] = "sgrp_%s-set_%s" % (sg, i)
                i += 1
                j_sets[incl_j]["exclusionJuncs"] = excl_j
                j_sets[incl_j]["data"] = getData(events[e]["data"], "counts", "incl")
                sg += 1

        elif ASclass == 'mutually_exclusive' or ASclass == 'alternative_first_exon' or ASclass == 'alternative_last_exon':

            excl_j_list = excl_j.split(";")
            excl_e_count = len(events[e]["excl_e"].split(";"))
            non_exon_excl_juncs = len(excl_j_list) - excl_e_count

            if events[e]["excl_e"] in exons and events[e]["incl_e"] not in exons[events[e]["excl_e"]]:
                complx = True
                break
            else:
                exons[events[e]["excl_e"]] = [events[e]["incl_e"]]
                      
            if events[e]["incl_e"] in exons and events[e]["excl_e"] not in exons[events[e]["incl_e"]]:
                complx = True
                break
            else:
                exons[events[e]["incl_e"]] = events[e]["excl_e"].split(";")
                
            if events[e]["incl_e"] in incl_exons:
                if incl_exons[events[e]["incl_e"]] != events[e]["incl_j"]:
                    complx = True
                    break
            else:
                incl_exons[events[e]["incl_e"]] = events[e]["incl_j"]

            all_juncs = excl_j_list[:]
            all_juncs.append(incl_j)
            for jj in all_juncs:
                for j in jj.split(','):
                    non_alt_ss = ''                     # should be same for all junctions in group
                    chrom, se  = j.split(":")
                    start, end = se.split("-")

                    if ASclass == 'alternative_last_exon':
                        if strand == "+":
                            non_alt_ss = start
                        elif strand == "-":
                            non_alt_ss = end
                    elif ASclass == 'alternative_first_exon':                                
                        if strand == "+":
                            non_alt_ss = end
                        elif strand == "-":
                            non_alt_ss = start

                    if "non_alt_ss" not in exons2:
                        exons2["non_alt_ss"] = non_alt_ss
                    elif exons2["non_alt_ss"] != non_alt_ss:
                        complx = True
                        break

            if excl_e_count == 1:

                if excl_j not in j_sets and incl_j not in j_sets:
                    j_sets[excl_j] = {}
                    j_sets[excl_j]["KN"] = "K"
                    if isNovel(excl_j, juncs_KN):
                        j_sets[excl_j]["KN"] = "N"
                    j_sets[excl_j]["juncSetID"] = "sgrp_1-set_%s" % i
                    i += 1
                    j_sets[excl_j]["exclusionJuncs"] = incl_j
                    j_sets[excl_j]["non_exon_exclJuncs"] = 0
                    j_sets[excl_j]["data"] = getData(events[e]["data"], "counts", "excl")

                    j_sets[incl_j] = {}
                    j_sets[incl_j]["KN"] = "K"
                    if isNovel(incl_j, juncs_KN):
                        j_sets[incl_j]["KN"] = "N"
                    j_sets[incl_j]["juncSetID"] = "sgrp_1-set_%s" % i
                    i += 1
                    j_sets[incl_j]["exclusionJuncs"] = excl_j
                    j_sets[incl_j]["non_exon_exclJuncs"] = non_exon_excl_juncs
                    j_sets[incl_j]["data"] = getData(events[e]["data"], "counts", "incl")

                elif excl_j not in j_sets:
                    incl_j_exclusionJuncs = j_sets[incl_j]["exclusionJuncs"].split(";")
                    if (len(incl_j_exclusionJuncs) - j_sets[incl_j]["non_exon_exclJuncs"]) > 1 and excl_j_list[-1] in incl_j_exclusionJuncs:
                        continue
                    else:
                        complx = True

                else:
                    complx = True

            elif excl_e_count > 1:

                if incl_j not in j_sets:
                    j_sets[incl_j] = {}
                    j_sets[incl_j]["KN"] = "K"
                    if isNovel(incl_j, juncs_KN):
                        j_sets[incl_j]["KN"] = "N"
                    j_sets[incl_j]["juncSetID"] = "sgrp_1-set_%s" % i
                    i += 1
                    j_sets[incl_j]["exclusionJuncs"] = excl_j
                    j_sets[incl_j]["non_exon_exclJuncs"] = non_exon_excl_juncs
                    j_sets[incl_j]["data"] = getData(events[e]["data"], "counts", "incl")

                else:
                    incl_j_exclusionJuncs = j_sets[incl_j]["exclusionJuncs"].split(";")

                    if (len(incl_j_exclusionJuncs) - j_sets[incl_j]["non_exon_exclJuncs"]) < (len(excl_j_list) - non_exon_excl_juncs):
                        fl = True
                        for i_j_e_j in incl_j_exclusionJuncs[j_sets[incl_j]["non_exon_exclJuncs"]:]:
                            if i_j_e_j not in excl_j_list:
                                fl = False
                        if fl:
                            j_sets[incl_j]["KN"] = "K"
                            if isNovel(incl_j, juncs_KN):
                                j_sets[incl_j]["KN"] = "N"
                            j_sets[incl_j]["juncSetID"] = j_sets[incl_j]["juncSetID"]
                            j_sets[incl_j]["exclusionJuncs"] = excl_j
                            j_sets[incl_j]["non_exon_exclJuncs"] = non_exon_excl_juncs
                            j_sets[incl_j]["data"] = getData(events[e]["data"], "counts", "incl")

                        else:
                            complx = True

                    elif (len(incl_j_exclusionJuncs) - j_sets[incl_j]["non_exon_exclJuncs"]) > (len(excl_j_list) - non_exon_excl_juncs):
                        fl = True
                        for e_j in excl_j_list[non_exon_excl_juncs:]:
                            if e_j not in incl_j_exclusionJuncs:
                                fl = False
                        if fl:
                            continue
                        else:
                            complx = True

                    else:
                        complx = True
    

        elif ASclass == 'alternative_acceptor' or ASclass == 'alternative_donor':

            excl_j_list = excl_j.split(";")
            all_juncs = excl_j_list[:]
            all_juncs.append(incl_j)
            
            for j in all_juncs:
                non_alt_ss = ''                     # should be same for all junctions in group
                chrom, se  = j.split(":")
                start, end = se.split("-")

                if ASclass == 'alternative_acceptor':
                    if strand == "+":
                        non_alt_ss = start
                    elif strand == "-":
                        non_alt_ss = end
                else:                                
                    if strand == "+":
                        non_alt_ss = end
                    elif strand == "-":
                        non_alt_ss = start

                if "non_alt_ss" not in exons:
                    exons["non_alt_ss"] = non_alt_ss
                elif exons["non_alt_ss"] != non_alt_ss:
                    complx = True
                    break
            
            if len(excl_j_list) == 1:

                if excl_j not in j_sets and incl_j not in j_sets:
                    j_sets[excl_j] = {}
                    j_sets[excl_j]["KN"] = "K"
                    if isNovel(excl_j, juncs_KN):
                        j_sets[excl_j]["KN"] = "N"
                    j_sets[excl_j]["juncSetID"] = "sgrp_1-set_%s" % i
                    i += 1
                    j_sets[excl_j]["exclusionJuncs"] = incl_j
                    j_sets[excl_j]["data"] = getData(events[e]["data"], "counts", "excl")

                    j_sets[incl_j] = {}
                    j_sets[incl_j]["KN"] = "K"
                    if isNovel(incl_j, juncs_KN):
                        j_sets[incl_j]["KN"] = "N"
                    j_sets[incl_j]["juncSetID"] = "sgrp_1-set_%s" % i
                    i += 1
                    j_sets[incl_j]["exclusionJuncs"] = excl_j
                    j_sets[incl_j]["data"] = getData(events[e]["data"], "counts", "incl")

                elif excl_j not in j_sets:
                    incl_j_exclusionJuncs = j_sets[incl_j]["exclusionJuncs"].split(";")
                    if len(incl_j_exclusionJuncs) > 1 and excl_j_list[0] in incl_j_exclusionJuncs:
                        continue
                    else:
                        complx = True

                else:
                    complx = True

            elif len(excl_j_list) > 1:

                if incl_j not in j_sets:
                    j_sets[incl_j] = {}
                    j_sets[incl_j]["KN"] = "K"
                    if isNovel(incl_j, juncs_KN):
                        j_sets[incl_j]["KN"] = "N"
                    j_sets[incl_j]["juncSetID"] = "sgrp_1-set_%s" % i
                    i += 1
                    j_sets[incl_j]["exclusionJuncs"] = excl_j
                    j_sets[incl_j]["data"] = getData(events[e]["data"], "counts", "incl")

                else:
                    incl_j_exclusionJuncs = j_sets[incl_j]["exclusionJuncs"].split(";")

                    if len(incl_j_exclusionJuncs) < len(excl_j_list):
                        fl = True
                        for i_j_e_j in incl_j_exclusionJuncs:
                            if i_j_e_j not in excl_j_list:
                                fl = False
                        if fl:
                            j_sets[incl_j]["KN"] = "K"
                            if isNovel(incl_j, juncs_KN):
                                j_sets[incl_j]["KN"] = "N"
                            j_sets[incl_j]["juncSetID"] = j_sets[incl_j]["juncSetID"]
                            j_sets[incl_j]["exclusionJuncs"] = excl_j
                            j_sets[incl_j]["data"] = getData(events[e]["data"], "counts", "incl")

                        else:
                            complx = True

                    elif len(incl_j_exclusionJuncs) > len(excl_j_list):
                        fl = True
                        for e_j in excl_j_list:
                            if e_j not in incl_j_exclusionJuncs:
                                fl = False
                        if fl:
                            continue
                        else:
                            complx = True

                    else:
                        complx = True
            
                        
    return j_sets, gene, chrom, strand, eventIDs, complx

    
def isNovel(juncSet, juncs_KN):
    junctions = juncSet.split(";")
    for j in junctions:
        for jj in j.split(","):
            if jj not in juncs_KN:
                continue
            #    print >>sys.stderr, "%s does not have K/N annotation" % jj
            else:
                return False
    return True

##[chrom/strand/KN/junctions/exclusionJuncs/readCounts(list)]
def makeJuncSetBedLines(juncSet, name, asclass, color_ind, jcn_len):

    color = COLORS_SMALL[color_ind]

    chrom = juncSet["chrom"]
    strand = juncSet["strand"]

    line = ''
    
    for jncs in juncSet["junctions"].split(","):
        for jnc in jncs.split(";"):
            chrom, se = jnc.split(":")
            start, end = se.split("-")
            start = int(start) - jcn_len - 1
            end = int(end) + jcn_len
            line += "%s\t%s\t%s\t%s\t0\t%s\t0\t0\t%s\t2\t%s,%s\t0,%s\n" % (chrom, start, end, name, strand,
                                                                   color, jcn_len, jcn_len, end - start - jcn_len)

    
    return line
                             
#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
