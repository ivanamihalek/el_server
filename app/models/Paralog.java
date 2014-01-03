package models;

import java.util.*;
import javax.persistence.*;
import play.db.jpa.*;

@Entity
@Table(name="paralog")
public class Paralog extends Model {
    public String gene_id1;
    public String gene_id2;
    public Paralog() {
	gene_id1 = "";
	gene_id2 = "";
    }
}
