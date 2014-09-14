package util;

import java.util.*;
import java.util.regex.Pattern;
import java.lang.reflect.Method;
import play.db.jpa.GenericModel;
import models.*;
import play.Logger;

public class GeneQuery {
    private String ensembl_id;
    private String species;

    private static final Map<String, String> speciesMap;
    static {
        Map<String,String> tempMap = new HashMap<String,String>();
	tempMap.put("AME", "ailuropoda_melanoleuca");
	tempMap.put("APL", "anas_platyrhynchos");
	tempMap.put("ACA", "anolis_carolinensis");
	tempMap.put("AMX", "astyanax_mexicanus");
	tempMap.put("BTA", "bos_taurus");
	tempMap.put("CJA", "callithrix_jacchus");
	tempMap.put("CAF", "canis_familiaris");
	tempMap.put("CPO", "cavia_porcellus");
	tempMap.put("CHO", "choloepus_hoffmanni");
	tempMap.put("DAR", "danio_rerio");
	tempMap.put("DNO", "dasypus_novemcinctus");
	tempMap.put("DOR", "dipodomys_ordii");
	tempMap.put("ETE", "echinops_telfairi");
	tempMap.put("ECA", "equus_caballus");
	tempMap.put("EEU", "erinaceus_europaeus");
	tempMap.put("FCA", "felis_catus");
	tempMap.put("FAL", "ficedula_albicollis");
	tempMap.put("GMO", "gadus_morhua");
	tempMap.put("GAL", "gallus_gallus");
	tempMap.put("GAC", "gasterosteus_aculeatus");
	tempMap.put("GGO", "gorilla_gorilla");
	tempMap.put("STO", "ictidomys_tridecemlineatus");
	tempMap.put("LAC", "latimeria_chalumnae");
	tempMap.put("LOC", "lepisosteus_oculatus");
	tempMap.put("LAF", "loxodonta_africana");
	tempMap.put("MMU", "macaca_mulatta");
	tempMap.put("MEU", "macropus_eugenii");
	tempMap.put("MGA", "meleagris_gallopavo");
	tempMap.put("MIC", "microcebus_murinus");
	tempMap.put("MOD", "monodelphis_domestica");
	tempMap.put("MUS", "mus_musculus");
	tempMap.put("MPU", "mustela_putorius_furo");
	tempMap.put("MLU", "myotis_lucifugus");
	tempMap.put("NLE", "nomascus_leucogenys");
	tempMap.put("OPR", "ochotona_princeps");
	tempMap.put("ONI", "oreochromis_niloticus");
	tempMap.put("OAN", "ornithorhynchus_anatinus");
	tempMap.put("OCU", "oryctolagus_cuniculus");
	tempMap.put("ORL", "oryzias_latipes");
	tempMap.put("OAR", "ovis_aries");
	tempMap.put("OGA", "otolemur_garnettii");
	tempMap.put("PTR", "pan_troglodytes");
	tempMap.put("PAN", "papio_anubis");
	tempMap.put("PFO", "poecilia_formosa");
	tempMap.put("PSI", "pelodiscus_sinensis");
	tempMap.put("PMA", "petromyzon_marinus");
	tempMap.put("PPY", "pongo_abelii");
	tempMap.put("PCA", "procavia_capensis");
	tempMap.put("PVA", "pteropus_vampyrus");
	tempMap.put("RNO", "rattus_norvegicus");
	tempMap.put("SHA", "sarcophilus_harrisii");
	tempMap.put("SAR", "sorex_araneus");
	tempMap.put("SSC", "sus_scrofa");
	tempMap.put("TGU", "taeniopygia_guttata");
	tempMap.put("TRU", "takifugu_rubripes");
	tempMap.put("TSY", "tarsius_syrichta");
	tempMap.put("TNI", "tetraodon_nigroviridis");
	tempMap.put("TBE", "tupaia_belangeri");
	tempMap.put("TTR", "tursiops_truncatus");
	tempMap.put("VPA", "vicugna_pacos");
	tempMap.put("XET", "xenopus_tropicalis");
	tempMap.put("XMA", "xiphophorus_maculatus");
        speciesMap = Collections.unmodifiableMap(tempMap);
    }

    public static boolean isEnsembl (String id) {
	return Pattern.matches("^ENS(?:[A-Z]{3})?[GPTE][0-9]{11}$", id);
    }

    public static String findSpecies (String id) {
        String sp = speciesMap.get(id.substring(3,6));
        if ((sp == null) && Pattern.matches("ENS[GPTE]\\d{11}", id))
            sp = "homo_sapiens";
        return sp;
    }

    public GeneQuery (String id) {
        if (id == null || id == "")
            throw new IllegalArgumentException("No id given");
	ensembl_id = id.replaceAll("\\s","").toUpperCase();
	if (!isEnsembl(ensembl_id))
            throw new IllegalArgumentException(id + "is not a valid Ensembl id");
        species = findSpecies(ensembl_id);
	if (species == null)
            throw new IllegalArgumentException(id + "does not belong to a recognizable species");
    }

    public String getId ()      { return ensembl_id; }
    public String getSpecies () { return species; }

    public String getSpeciesCode () {
        String species_code;
        if (species.equals("homo_sapiens"))
            species_code = "HSA";
        else
            species_code = ensembl_id.substring(3,6);
        return species_code;
    }

    public String getDescription () {
        List<Description> descs = Description.find("ensembl_gene_id", ensembl_id).fetch(1);
        String desc = "";
        if (!descs.isEmpty()) desc = descs.get(0).getDesc();
        return desc;
    }

    public List<Exon> getExons () {
	//Logger.info("A log message");
        List<Exon> exons = new ArrayList<Exon>();
        try {
            Class c = Class.forName("models.Exon_" + species);
            Method m = c.getMethod("find", String.class, Object[].class);
	    exons = ((GenericModel.JPAQuery) m.invoke(
                null, "ensembl_gene_id", new Object[] {ensembl_id})).fetch();
            Collections.sort(exons);
        } catch (Exception e) {
            System.out.println(e);
            // Just return empty list if there is an exception
        }
        return exons;
    }

    public List<Ortholog> getOrthologs () {
        List<Ortholog> orthologs = new ArrayList<Ortholog>();
        String human_id;
        if (!species.equals("homo_sapiens")) {
            // Find human ortholog gene id if not human
            orthologs = Ortholog.find("cognate_gene_id", ensembl_id).fetch(1);
            if (orthologs.isEmpty())
                return orthologs; // empty list on fail
            human_id = orthologs.get(0).ensembl_gene_id;
        } else {
            human_id = ensembl_id;
        }
        // Find orthologs of human gene
        orthologs = Ortholog.find("ensembl_gene_id", human_id).fetch();
        // Add human to head of list of orthologs
        if (!orthologs.isEmpty()) {
            Ortholog homo = new Ortholog();
            homo.species = "homo_sapiens";
            homo.common_name = "human";
            homo.ensembl_gene_id = human_id;
            homo.cognate_gene_id = human_id;
            orthologs.add(0,homo);
            Collections.sort(orthologs);
        }
        return orthologs;
    }

    public List<String> getParalogIds () {
        List<String> paralog_ids = new ArrayList<String>();
        List<Paralog> paralogs = new ArrayList<Paralog>();
        String primary_id = ensembl_id;
	paralogs = Paralog.find("gene_id1 = ? or gene_id2 = ?", primary_id, primary_id).fetch();
        if (!paralogs.isEmpty()) {
            if (!paralogs.get(0).gene_id1.equals(primary_id)) {
                primary_id = paralogs.get(0).gene_id1;
                paralogs = Paralog.find("gene_id1", primary_id).fetch();
            }
            paralog_ids.add(primary_id); // First element is primary id
            for (Paralog p: paralogs)
                paralog_ids.add(p.gene_id2);
        }
        return paralog_ids; // empty list on fail
    }

    public List<Splicing> getSplicings () {
        return Splicing.find("gene_id", ensembl_id).fetch();
    }
}
