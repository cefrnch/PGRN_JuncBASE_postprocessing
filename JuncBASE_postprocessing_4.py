#!usr/bin/env python

#
# this script will pull the pvals values for "junction sets" from collection of pairwise _all_psi_output.txt files
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
   
    opt_parser.add_option("-a",
                          dest="all_psi_out",
                          type="string",
                          help="dir of inputs",
                          default=None)
    opt_parser.add_option("-j",
                          dest="junctionSets",
                          type="string",
                          help="junctionSets",
                          default=None)
    
    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("-a")
    opt_parser.check_required("-j")

    in_dir = options.all_psi_out
    # First read in junctionSet information we already have

    juncSets = {}                     # juncSets[asclass][group][junctions][gene/set/chr/str/KN/exclJuncs/comp]

    j = open(options.junctionSets)
    jh = j.readline()
    jlines = j.readlines()
    for jline in jlines:
        jline = formatLine(jline).split("\t")
        old_group = "_".join(jline[2].split("_")[0:2])
        if jline[1] not in juncSets:
            juncSets[jline[1]] = {}
        if old_group not in juncSets[jline[1]]:
            juncSets[jline[1]][old_group] = {}
        if jline[7] not in juncSets[jline[1]][old_group]:
            juncSets[jline[1]][old_group][jline[7]] = {}
        juncSets[jline[1]][old_group][jline[7]]["gene"] = jline[0]
        juncSets[jline[1]][old_group][jline[7]]["group"] = jline[2]
        juncSets[jline[1]][old_group][jline[7]]["set"] = jline[3]
        juncSets[jline[1]][old_group][jline[7]]["chr"] = jline[4]
        juncSets[jline[1]][old_group][jline[7]]["str"] = jline[5]
        juncSets[jline[1]][old_group][jline[7]]["KN"] = jline[6]
        juncSets[jline[1]][old_group][jline[7]]["exclJuncs"] = jline[8]
    j.close()
        

    # Next read in PSI values for OK junction sets
    comps = []
    for fl in os.listdir(in_dir):
        comp = fl.split("_")[0]
        comps.append(comp)
        f = open(in_dir + fl)
        header = f.readline()
        lines = f.readlines()
        for line in lines:
            line = formatLine(line).split("\t")
            excl_j, incl_j = line[5], line[6]
            old_grp = line[-7]
            if old_grp in juncSets[line[1]]:
                if incl_j in juncSets[line[1]][old_grp] and juncSets[line[1]][old_grp][incl_j]["exclJuncs"] == excl_j:
                    if comp not in juncSets[line[1]][old_grp][incl_j]:
                        juncSets[line[1]][old_grp][incl_j][comp] = line[-3:]
                    elif not isMatch(juncSets[line[1]][old_grp][incl_j][comp], line[-3:]):
                        juncSets[line[1]][old_grp][incl_j][comp] = ["NA","NA","NA"]
                if excl_j in juncSets[line[1]][old_grp] and juncSets[line[1]][old_grp][excl_j]["exclJuncs"] == incl_j:
                    if comp not in juncSets[line[1]][old_grp][excl_j]:
                        juncSets[line[1]][old_grp][excl_j][comp] = flipDeltaPSI(line[-3:])
                    elif not isMatch(juncSets[line[1]][old_grp][excl_j][comp], flipDeltaPSI(line[-3:])):
                        juncSets[line[1]][old_grp][excl_j][comp] = ["NA","NA","NA"]     
                        
    f.close()
    
    # Now print output 
    
    print "gene\tASclass\tgroupID\tjuncSetID\tchr\tstrand\tknown_or_novel\tjunctions\texclusion_junctions",
    for c in comps:
        print "\t%s_deltaPSI\t%s_raw_pval\t%s_corrected_pval" % (c, c, c),
    print
    
    for asclass in juncSets:
        for g in juncSets[asclass]:
            for j in juncSets[asclass][g]:
                print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (juncSets[asclass][g][j]["gene"], asclass, juncSets[asclass][g][j]["group"], 
                                                              juncSets[asclass][g][j]["set"], juncSets[asclass][g][j]["chr"], juncSets[asclass][g][j]["str"],
                                                              juncSets[asclass][g][j]["KN"], j, juncSets[asclass][g][j]["exclJuncs"]),
                for comp in comps:
                    if comp not in juncSets[asclass][g][j]:
                        juncSets[asclass][g][j][comp] = ["NA","NA","NA"]
                    print "\t%s" % ("\t".join(juncSets[asclass][g][j][comp])),
                print
    
            
                    
  
############
# END_MAIN #
############

def flipDeltaPSI(psis):
    psis[0] = str(float(psis[0]) * -1)
    return psis

def isMatch(list1, list2):
    for i in range(len(list1)):
        if list1[i] != list2[i]:
            return False
    return True
        
def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
