#!usr/bin/env python

#
# this script will take a JuncBASE output and add an event_ID and group related events with a group_ID
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

# file to use to make event and group IDs
FIRST_FILE = "AS_exclusion_inclusion_counts.txt"

# use info from FIRST_FILE to edit these files
OTHER_FILES = ["AS_exclusion_inclusion_counts_lenNorm.txt", "left_intron_counts.txt", "right_intron_counts.txt", "left_intron_counts_lenNorm.txt", "right_intron_counts_lenNorm.txt"]

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
   
    opt_parser.add_option("--directory",
                          dest="directory",
                          type="string",
                          help="directory of files to fix",
                          default=None)
    opt_parser.add_option("--prefix",
                          dest="prefix",
                          type="string",
                          help="prefix of files to fix",
                          default=None)

    (options, args) = opt_parser.parse_args()
	
    # validate the command line arguments
    opt_parser.check_required("--directory")
    opt_parser.check_required("--prefix")
    
    directory = options.directory
    prefix = options.prefix
    
    # First read in all events, merging as we go
    f = open(directory + "/" + prefix + FIRST_FILE)
    header = f.readline()
    header = formatLine(header)
    lines = f.readlines()

    as_type2chr2strand2redundantGroup2event = {}
    
    for event in lines:    
        event=formatLine(event)

        buildDictionaries(event, as_type2chr2strand2redundantGroup2event)
        
    f.close()

    # Then, merge redundant events
    changeOccurred = True
    while changeOccurred:
        changeOccurred = mergeRedundantEvents(as_type2chr2strand2redundantGroup2event)


    # next determine group IDs and events IDs, print for FIRST_FILE and store info for OTHER_FILES

    group_ids = {}                           # event_ids[cols 4-9 joined with |] = event_ID  (group ID)
    event_ids = {}                               # e_ids[cols 4-9 joined with |] = event_#   (event #)

    group_id_num = 0
    event_id_num = 0    

    outfile = open(directory + '/' + prefix + 'grouped_' + FIRST_FILE, 'w')
    print >>outfile, header + "\tgroup_ID\tevent_ID"
    for as_type in as_type2chr2strand2redundantGroup2event:
        for chr in as_type2chr2strand2redundantGroup2event[as_type]:
            for strand in as_type2chr2strand2redundantGroup2event[as_type][chr]:
                for group in as_type2chr2strand2redundantGroup2event[as_type][chr][strand]:
                    group_id_num += 1
                    group_id = "group_%d" % group_id_num
                    for event in as_type2chr2strand2redundantGroup2event[as_type][chr][strand][group]:
                        event_id_num += 1
                        event_id = "event_%d" % event_id_num
                        print >>outfile, event + "\t%s\t%s" % (group_id, event_id)
                        ev = event.split("\t")
                        n = "|".join(ev[3:9])
                        group_ids[n] = group_id
                        event_ids[n] = event_id
    outfile.close()
                            
                    
    # then do for other files
    
    for nm in OTHER_FILES:
        f = open(directory + "/" + prefix + nm)
        outfile = open(directory + '/' + prefix + 'grouped_' + nm, 'w')

        header = f.readline()
        header = formatLine(header)
        print >>outfile, header + "\tgroup_ID\tevent_ID"
        
        lines = f.readlines()
        
        for line in lines:
            line = line.replace("\"", "")
            line = formatLine(line)
            ln = line.split("\t")
            n = "|".join(ln[3:9])
            print >>outfile, line + "\t%s\t%s" % (group_ids[n], event_ids[n])

        f.close()
        outfile.close()


############
# END_MAIN #
############

def checkFileType(path_file):
    names = path_file.split(".")
    return names[-1]

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def getFileName(path_file, no_ext='no'):
    filename = path_file.split("/")[-1]
    if no_ext == 'yes':
        names = filename.split(".")
        return ".".join(names[0:-1])
    return filename

def getPath(path_file):
    filename = path_file.split("/")[0:-1] + "/"
    return filename



## FROM ANGELA BROOKS: makeNonRedundantAS_forCourtney_DONTUSE.py

# edited to get rid of pval related stuff
def buildDictionaries(event, as_type2chr2strand2redundantGroup2event):

    line_list = event.split("\t")

    as_type = line_list[1]

    chr = line_list[3]
    strand = line_list[4]

    if as_type == "mutually_exclusive":
        # Largest region will be in exclusion junctions
        group_start, group_end = findLargestRegion(",".join([line_list[5],
                                                             line_list[6]]))    
    else:
        group_start, group_end = findLargestRegion(line_list[5])    

    redundantRegion = (chr, strand, group_start, group_end)

    updateRedundantDictionary(as_type2chr2strand2redundantGroup2event, as_type, redundantRegion,
                              event)


def findLargestRegion(coords_string):

    first_list = coords_string.split(";")
    full_list = []
    for elem1 in first_list:
        for elem2 in elem1.split(","):
            full_list.append(elem2)

    leftmost = INFINITY
    rightmost = -1 

    for coord in full_list:
        chr, start, end = convertCoordStr(coord)

        if start < leftmost:
            leftmost = start

        if end > rightmost:
            rightmost = end

    return leftmost, rightmost 

def mergeRedundantEvents(as_type2chr2strand2redundantGroup2event):
    """
    Because of the construction of the dictionary, some of the redundant
    groups are overlapping.
    """
    for as_type in as_type2chr2strand2redundantGroup2event:
        for chr in as_type2chr2strand2redundantGroup2event[as_type]:
            for strand in as_type2chr2strand2redundantGroup2event[as_type][chr]:
                for first_group in as_type2chr2strand2redundantGroup2event[as_type][chr][strand]:
                    for second_group in as_type2chr2strand2redundantGroup2event[as_type][chr][strand]:
                        if first_group == second_group:
                            continue
                        if coordsOverlap(first_group[2], first_group[3],
                                         second_group[2], second_group[3]):

                            # Add sets together
                            new_set = as_type2chr2strand2redundantGroup2event[as_type][chr][strand][first_group].union(as_type2chr2strand2redundantGroup2event[as_type][chr][strand][second_group])
                            new_start = min(first_group[2], second_group[2])
                            new_end = max(first_group[3], second_group[3])

                            new_region = (first_group[0], first_group[1],
                                          new_start, new_end)

                            if new_region in as_type2chr2strand2redundantGroup2event[as_type][chr][strand]:
                                as_type2chr2strand2redundantGroup2event[as_type][chr][strand][new_region].update(new_set)
                            else:
                                as_type2chr2strand2redundantGroup2event[as_type][chr][strand][new_region] = new_set

                            # Remove old sets and return with flag
                            if new_region != first_group:
                                del as_type2chr2strand2redundantGroup2event[as_type][chr][strand][first_group]
                            if new_region != second_group:
                                del as_type2chr2strand2redundantGroup2event[as_type][chr][strand][second_group]

                            return True

    # If not returned from within the loop, then no change occurred
    return False

def updateRedundantDictionary(as_type2chr2strand2redundantGroup2event, as_type, redundantRegion,
                              event):

    this_chr = redundantRegion[0]
    this_strand = redundantRegion[1]    

    if as_type in as_type2chr2strand2redundantGroup2event:
        if not this_chr in as_type2chr2strand2redundantGroup2event[as_type]:
            as_type2chr2strand2redundantGroup2event[as_type][this_chr] = {this_strand:{redundantRegion: set([event])}}
            return
        if not this_strand in as_type2chr2strand2redundantGroup2event[as_type][this_chr]:
            as_type2chr2strand2redundantGroup2event[as_type][this_chr][this_strand] = {redundantRegion: set([event])}
            return

        foundOverlap = False
        for redundantGroup in as_type2chr2strand2redundantGroup2event[as_type][this_chr][this_strand]:
            if coordsOverlap(redundantGroup[2], redundantGroup[3],
                             redundantRegion[2], redundantRegion[3]):
                foundOverlap = True
                # Pick longest region as the key
                region_start = min(redundantGroup[2], redundantRegion[2])
                region_end = max(redundantGroup[3], redundantRegion[3])

                # Copy the existing set of exon_names
                cur_set = set(as_type2chr2strand2redundantGroup2event[as_type][this_chr][this_strand][redundantGroup]) 
                # Add the current event
                cur_set.add(event)

                del as_type2chr2strand2redundantGroup2event[as_type][this_chr][this_strand][redundantGroup]

                # Check for existing new group
                new_group = (redundantGroup[0], redundantGroup[1],
                             region_start, region_end)
                if new_group in as_type2chr2strand2redundantGroup2event[as_type][this_chr][this_strand]:
                    as_type2chr2strand2redundantGroup2event[as_type][this_chr][this_strand][new_group].update(cur_set)
                else:
                    as_type2chr2strand2redundantGroup2event[as_type][this_chr][this_strand][new_group] = cur_set
                break
        if not foundOverlap:
            # Add this non overlapping region on its own
            as_type2chr2strand2redundantGroup2event[as_type][this_chr][this_strand][redundantRegion] = set([event]) 
    else:
        as_type2chr2strand2redundantGroup2event[as_type] = {this_chr:{this_strand:{redundantRegion: set([event])}}}

                             
#################
# END FUNCTIONS #	
#################	
if __name__ == "__main__": main()
