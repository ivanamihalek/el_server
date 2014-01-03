package models;

import java.util.*;

import javax.persistence.*;

import play.db.jpa.*;

@Entity
@Table(name="ortholog")
public class Ortholog extends Model implements Comparable<Ortholog> {
    public String ensembl_gene_id;
    public String cognate_gene_id;
    public String species;
    public String common_name;
    public Ortholog() {
	ensembl_gene_id = "";
	cognate_gene_id = "";
	species = "";
	common_name = "";
    }

    public int compareTo(Ortholog other) {
	return this.species.compareTo(other.species);
    }
}
