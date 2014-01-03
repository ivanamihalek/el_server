package models;

import java.util.*;

import javax.persistence.*;

import play.db.jpa.*;

@Entity
@Table(name="name_resolution")
public class Resolution extends Model{
    public String synonym;
    public String stable_id;

    public String getDesc () {
        return synonym;
    }
}
