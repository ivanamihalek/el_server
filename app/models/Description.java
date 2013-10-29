package models;

import java.util.*;

import javax.persistence.*;

import play.db.jpa.*;

@Entity
@Table(name="description")
public class Description extends Model{
    public String ensembl_gene_id;
    public String descr_txt;

    public String getDesc () {
        String desc = descr_txt.split("\\[")[0]; // Remove trailing junk
        if (desc.equals("None")) desc = "";
        return desc;
    }
}
