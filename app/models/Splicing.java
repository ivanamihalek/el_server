package models;

import java.util.*;

import javax.persistence.*;

import play.db.jpa.*;

@Entity
@Table(name="splicing")

public class Splicing extends Model {
    public String gene_id;
    public String transcript_id;
    public String hit_name;

@ElementCollection
@CollectionTable(
  name="splicing_transcript",
  joinColumns=@JoinColumn(referencedColumnName="transcript_id", name="transcript_id")
)
@Column(name="exon_id")
    public List<String> exon_ids;
}
