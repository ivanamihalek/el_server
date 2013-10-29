#!/usr/local/bin/python2.7

import sys, commands
from  el_utils.mysql     import *
from  el_utils.ensembl   import *
from  el_utils.exon      import Exon
from  el_utils.bitstring import Bits
from  el_utils.tree      import species_sort
from  el_utils.map       import Map, get_maps
from  el_utils.utils     import output_fasta

#########################################
def find_human_cognate(cursor, ensembl_db_name, exon_id, exon_known, species_db_id):

    switch_to_db (cursor, ensembl_db_name['homo_sapiens']);
    
    qry  = "select exon_id, exon_known from exon_map where cognate_exon_id = %d " % exon_id
    qry += " and cognate_exon_known = %d "   % exon_known
    qry += " and cognate_genome_db_id = %d " % species_db_id
    rows = search_db (cursor, qry)
     
    if not rows or 'ERROR' in rows[0]:
        return [0,0]

    return rows[0]
    
#########################################
def make_exon_alignment(cursor, ensembl_db_name, human_exon_id, human_exon_known):

    sequence = {}

    # find all other exons that map to the human exon
    maps = get_maps(cursor, ensembl_db_name, human_exon_id, human_exon_known)

    for map in maps:
        # get the raw (unaligned) sequence for the exon that maps onto human
        exon_seqs = get_exon_seqs(cursor, map.exon_id_2, map.exon_known_2, ensembl_db_name[map.species_2])
        if (not exon_seqs ):
            print map
            exit (1)
        exon_seqs = exon_seqs[1:] # the first entry is database id
        [pep_seq, pepseq_transl_start, 
         pepseq_transl_end, left_flank, right_flank, dna_seq] = exon_seqs

        # inflate the compressed sequence
        unaligned_sequence = pep_seq
        if not unaligned_sequence:continue

        bs = Bits(bytes=map.bitmap)
        # check bitmap has correct number of 1s
        if ( not bs.count(1) == len(unaligned_sequence)):
	    print bs.bin
	    print unaligned_sequence
            print "bitmap check fails (?)"
            continue

        # rebuild aligned sequence
        usi = iter(unaligned_sequence)
        reconstructed_sequence = "".join(('-' if c=='0' else next(usi) for c in bs.bin))

        # come up with a unique name for this sequence
        species       = map.species_2                
        sequence_name = species + "_" + str(map.exon_id_2)+"_"+str(map.exon_known_2)

        sequence[sequence_name] = reconstructed_sequence
    
    return sequence


#########################################
def sort_names (sorted_species, alignment):

    sorted_names = []
    for species in sorted_species:
        for seq_name, sequence in alignment.iteritems():
            if (species in seq_name):
                sorted_names.append(seq_name)
    return sorted_names
            
    
#########################################
def main():

    # exon_id comes from the command line
    if len(sys.argv) < 5:
        print "Usage: %s <exon_id>  <exon_known>  <species>  <output name>" % sys.argv[0]
        exit (1)
        
    exon_id    = long(sys.argv[1])
    exon_known = int(sys.argv[2])
    species    = sys.argv[3]
    afa_name   = sys.argv[4]

    db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)
   
    if (not is_coding_exon(cursor, exon_id, exon_known, ensembl_db_name[species] ) ):
        # make an empty file
        print "NOT CODING"
        cmd = "touch " + afa_name
        ret = commands.getoutput(cmd)
        cursor.close()
        db.close()
        return

    ######################################
    if (species == 'homo_sapiens'):
        [human_exon_id, human_exon_known] = [exon_id, exon_known]
        ok = True
    else:
        # find the human exon this guy maps to
        species_db_id = species2genome_db_id(cursor, species)
        if (species_db_id):
            [human_exon_id, human_exon_known] = find_human_cognate(cursor, ensembl_db_name, exon_id, exon_known, species_db_id)
        ok = species_db_id > 0 and human_exon_id>0

    ######################################
    if (ok):
        alignment    = make_exon_alignment(cursor, ensembl_db_name, human_exon_id, human_exon_known)   
        
    if (ok and alignment):
        # sort the remaining species  taxonomically
        sorted_species = species_sort(cursor, all_species, species)
        sorted_names = sort_names (sorted_species, alignment)
        output_fasta (afa_name, sorted_names, alignment)
    else: 
        # make file consisting of the original sequence only
        pepseq = get_exon_pepseq (cursor, exon_id, exon_known, ensembl_db_name[species])
        if pepseq:
            alignment = {}
            seq_name  = species + "_" + str(exon_id) + "_" + str(is_known)
            alignment[seq_name] = pepseq
            output_fasta (afa_name, [seq_name], alignment)
        else:
            # if not even the original sequence can be found, its definitely somebody else's fault;
            # make an empty file
            print "NO SEQUENCE"
            cmd = "touch " + afa_name
            ret = commands.getoutput(cmd)

    cursor.close()
    db.close()
    
    return

 
#########################################
if __name__ == '__main__':
    main()
