package models;
import java.util.*;
import javax.persistence.*;
import play.db.jpa.*;

@MappedSuperclass
//@Inheritance(strategy = InheritanceType.TABLE_PER_CLASS)
public abstract class Exon extends Model implements Comparable<Exon> {
    public String exon_key;
    public String ensembl_gene_id;
    public String ensembl_exon_id;
    public long   start_in_gene;
    public long   end_in_gene;
    public int    strand;
    public int    is_known;	 
    public int    is_coding;	 
    public int    is_canonical;	 
    public int    is_constitutive;	 
    public String species;
    public String source;
    public String protein_seq;
    public String left_flank;
    public String right_flank;
    public String dna_seq;
    //public String maps_to_human_exon_id;

    public int compareTo(Exon other) {
	return (int)(this.start_in_gene - other.start_in_gene);
    }
}
