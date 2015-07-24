#!/usr/local/bin/python2.7

import sys, commands
from  el_utils.mysql     import *
from  el_utils.ensembl   import *
from  el_utils.exon      import Exon
from  el_utils.bitstring import Bits
from  el_utils.tree      import species_sort
from  el_utils.map       import Map, get_maps
from  el_utils.utils     import output_fasta
# BioPython
from Bio.Seq      import Seq
from Bio.Alphabet import generic_dna

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
def align_nucseq_by_pepseq(aligned_pepseq, nucseq):
    if (not len(aligned_pepseq.replace('-',''))*3 == len(nucseq)):
        print "length mismatch: ", len(aligned_pepseq.replace('-',''))*3, len(nucseq)
        return ""
    codons = iter(map(''.join, zip(*[iter(nucseq)]*3)))
    aligned_nucseq = ''.join(('---' if c=='-' else next(codons) for c in aligned_pepseq))
    return aligned_nucseq

#########################################
def expand_pepseq (reconstructed_pep_sequence, exon_seqs):
    
    dna_aln_expanded = ""
    
    flank_length     = 10 # where should this come from?
    
    [pepseq, pepseq_transl_start, 
     pepseq_transl_end, left_flank, right_flank, dna_seq] = exon_seqs
    # coding dna sequence:
    cds = Seq(dna_seq[pepseq_transl_start:pepseq_transl_end], generic_dna).tostring()
    if not cds: return ""
    aligned_nucseq  = align_nucseq_by_pepseq(reconstructed_pep_sequence, cds)
    if not aligned_nucseq: return ""

    effective_left_flank  = ""
    effective_right_flank = "" 

    #######
    effective_left_flank  = left_flank
    if pepseq_transl_start>0:
        effective_left_flank += dna_seq[:pepseq_transl_start]
    if len(effective_left_flank) > flank_length: 
        effective_left_flank = effective_left_flank[-flank_length:]
    effective_left_flank = effective_left_flank.lower()
    #######
    effective_right_flank = right_flank
    delta = len(dna_seq)-pepseq_transl_end
    if delta>0:
        effective_right_flank = dna_seq[-delta:]+effective_right_flank
    if len(effective_right_flank) > flank_length: 
        effective_right_flank = effective_right_flank[:flank_length]
    effective_right_flank = effective_right_flank.lower()
    
    #######
    # pad the flanking seqs to the needed length
    effective_left_flank  = effective_left_flank.rjust (flank_length, '-')
    effective_right_flank = effective_right_flank.ljust(flank_length, '-')

    #print effective_left_flank, " ** ",  effective_right_flank

    
    dna_aln_expanded = effective_left_flank + aligned_nucseq + effective_right_flank
    return dna_aln_expanded

#########################################
def  make_seq_name (cursor, ensembl_db_name, species, exon_id, exon_known, exon_seqs):    

    sequence_name = ""

    if not exon_seqs:
        # get the raw (unaligned) sequence for the exon that maps onto human
        exon_seqs = get_exon_seqs(cursor, map.exon_id_2, map.exon_known_2, ensembl_db_name[map.species_2])
        if exon_seqs and  len(exon_seqs)>=7:
            [exon_seq_id, pepseq, pepseq_transl_start, pepseq_transl_end, 
             left_flank, right_flank, dna_seq] = exon_seqs

    if exon_known == 1:
        exon_stable_id = exon2stable (cursor, exon_id, ensembl_db_name[species])
        if exon_stable_id:
            sequence_name = exon_stable_id + "_" + species
        else:
            sequence_name = str(exon_id) + "_" + species
    if not sequence_name: 
        gene_id = exon_id2gene_id (cursor, ensembl_db_name[species],  exon_id, exon_known)
        if gene_id:
            gene_stable_id = gene2stable (cursor, gene_id, ensembl_db_name[species])
            if gene_stable_id:
                sequence_name = gene_stable_id  + "_" + species
                if not pepseq_transl_start is None and not  pepseq_transl_end is None:
                    sequence_name += "_" +  pepseq_transl_start+ "_" + pepseq_transl_end
                else:
                    sequence_name += "_unknown_region"

    return sequence_name

#########################################
def make_exon_alignment(cursor, ensembl_db_name, human_exon_id, human_exon_known, nt=False):
    sequence = {}
    shortest_l = -1 # Uninitialized leading padding length
    shortest_r = -1 # Uninitialized trailing padding length

    # find all other exons that map to the human exon
    maps = get_maps(cursor, ensembl_db_name, human_exon_id, human_exon_known)

        
    for map in maps:
        # get the raw (unaligned) sequence for the exon that maps onto human
        exon_seqs = get_exon_seqs(cursor, map.exon_id_2, map.exon_known_2, ensembl_db_name[map.species_2])
        if not exon_seqs or len(exon_seqs)<7:
            #print map
            continue

        [exon_seq_id, pepseq, pepseq_transl_start, pepseq_transl_end, left_flank, right_flank, dna_seq] = exon_seqs

        bs = Bits(bytes=map.bitmap)
        if (not bs.count(1) == len(pepseq)): continue # check bitmap has correct number of 1s
        usi = iter(pepseq)
        reconst_pepseq = "".join(('-' if c=='0' else next(usi) for c in bs.bin))
  
        #########################################################
        # come up with a unique name for this sequence
        species       = map.species_2
        sequence_name = make_seq_name (cursor, ensembl_db_name, species,  map.exon_id_2, map.exon_known_2, exon_seqs)

        if not sequence_name: # for whichever reason we still do not have the name here
            sequence_name = "anon_" + species 
        #########################################################
           
        if nt:
            reconst_ntseq = expand_pepseq (reconst_pepseq, exon_seqs[1:])
            if reconst_ntseq: 
                sequence[sequence_name] = reconst_ntseq
                aln_length = len(reconst_ntseq)
        else:
            if reconst_pepseq: 
                sequence[sequence_name] = reconst_pepseq
                aln_length = len(reconst_pepseq)
                
    # strip common gaps
    all_gaps = {}  
    for pos in range(aln_length):
        all_gaps[pos] = True
        for name, seq in sequence.iteritems():
            if (not seq[pos]=='-'):
                all_gaps[pos] = False
                break

    sequence_stripped = {}
    for name, seq in sequence.iteritems():
        sequence_stripped[name] = ""
        for pos in range(aln_length):
            if all_gaps[pos]: continue
            sequence_stripped[name] += seq[pos]

    return sequence_stripped

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
        print "Usage: %s <exon_id>  <exon_known>  <species>  <output name> [nt]" % sys.argv[0]
        exit (1)
        
    exon_id    = long(sys.argv[1])
    exon_known = int(sys.argv[2])
    species    = sys.argv[3]
    afa_name   = sys.argv[4]

    nt =  len(sys.argv)>5 and sys.argv[5]=='nt'
    
    ######################################
    db     = connect_to_mysql(user="root", passwd="sqljupitersql", host="jupiter.private.bii", port=3307)
    cursor = db.cursor()
    [all_species, ensembl_db_name] = get_species (cursor)
   
    if not is_coding_exon(cursor, exon_id, exon_known, ensembl_db_name[species]) and not nt:
        # make an empty file
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
            [human_exon_id, human_exon_known] = find_human_cognate(cursor, ensembl_db_name, exon_id,
                                                                   exon_known, species_db_id)
        ok = species_db_id > 0 and human_exon_id>0

    ######################################
    if (ok):
        alignment = make_exon_alignment(cursor, ensembl_db_name, human_exon_id, human_exon_known, nt)   
    if (ok and alignment):
        # sort the remaining species  taxonomically
        sorted_species = species_sort(cursor, all_species, species)
        sorted_names = sort_names (sorted_species, alignment)
        output_fasta (afa_name, sorted_names, alignment)
    else: 
        # make file consisting of the original sequence only
        if nt:
            seq = get_exon_seqs (cursor, exon_id, exon_known, ensembl_db_name[species])[-1];
        else:
            seq = get_exon_pepseq (cursor, exon_id, exon_known, ensembl_db_name[species])
        if seq:
            alignment = {}
            sequence_name  = make_seq_name (cursor, ensembl_db_name, species, exon_id, exon_known, [])
            alignment[sequence_name] = seq;
            output_fasta (afa_name, [seq_name], alignment)
        else:
            # if not even the original sequence can be found, its definitely somebody else's fault;
            # make an empty file
            cmd = "touch " + afa_name
            ret = commands.getoutput(cmd)

    cursor.close()
    db.close()
    
    return

 
#########################################
if __name__ == '__main__':
    main()
