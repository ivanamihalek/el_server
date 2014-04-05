// FIXME: Exceptions - Handle or throw exceptions properly instead of logging and returning success values
// FIXME: OO - Make query functions (get*) methods in the relevant classes
// FIXME: Script calling - Encapsulate script calling in functions or script object

package controllers;

import play.*;
import play.mvc.*;
import play.data.validation.Required;
import play.libs.Codec;
import play.libs.WS;
import play.db.jpa.GenericModel;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.lang.reflect.Method;

import models.*;

import util.GeneQuery;

public class Application extends Controller {

    @Before
    static void addDefaults() {
	renderArgs.put("applicationTitle", Play.configuration
		       .getProperty("application.title"));
	renderArgs.put("applicationBaseline", Play.configuration
		       .getProperty("application.baseline"));
    }

    /////////////////// Application entry points//////////////////////

    public static void index() {
	// Home page
	render();
    }

    public static void help() {
	// Help page
	render();
    }

    public static void contactPage(String error) {
	// Contact page (with potential error message)
	render(error);
    }

    public static void confirmation(String first_name, String last_name, String telephone,
	@Required String email, @Required String comments) {
	// Confirmation page (called from Contact page)
	// If the required forms are filled in contact page, comment is sent to
	// database and user is directed to the confirmation page
	// for now, comments are collected in /Resource/UserComments.txt
	// FIXME: later they should be changed to e-mails
	
	try {
	    BufferedWriter out = new BufferedWriter(new FileWriter(
		Play.applicationPath + "/Resource/UserComments.txt", true));
	    out.write("Name:     " + first_name);
	    out.newLine();
	    out.write("Surname:  " + last_name);
	    out.newLine();
	    out.write("E-mail:   " + email);
	    out.newLine();
	    out.write("Phone:    " + telephone);
	    out.newLine();
	    out.write("Comments: " + comments);
	    out.newLine();
	    out.newLine();
	    out.close();
	} catch (IOException e) {
	}

	render(first_name, last_name, email, telephone, comments);
    }

    public static void Sequence(String protein_seq) {
	
        String ens_id = "", Final = "";
	Set<Description> hits = new HashSet<Description>();
	
	if (protein_seq == null || !InDatabaseRange(protein_seq)) {
            Final = "Protein sequence is too long or short.";
        } else if (!fasta.matcher(protein_seq).matches()) {
            Final = "Protein sequence should be in fasta format.";
        } else {
            String result[] = blastall(protein_seq,
                Play.configuration.getProperty("all_species")).split("\t");
            if (result.length > 3) {
                ens_id = result[1];
                Final = "The best match is " + ens_id + " with " + result[2] +
                        "% identity over the aligned length of " + result[3] +
                        " positions. ";
		List<Description> descs = Description.find("ensembl_gene_id", ens_id).fetch(1);
		hits.addAll(descs);
           } else {
                Final = "BLAST lookup failed.";
            }
        }
	render(ens_id, protein_seq, Final, hits);
    }

    public static void Search(
	@Required(message = "Ensembl gene id is required") String query) {
	query = query.replaceAll("\\s+","").toUpperCase();
	// Search by ensembl ID, description or synonym
        try {
            GeneQuery q = new GeneQuery(query);
	    Result(q.getId()); // Valid ensembl ID, so give the result page
        } catch (IllegalArgumentException e) {
	    // Not a valid ensembl id, so resolve it as a description or synonym
	    List<Description> hits_desc = 
	    Description.find("descr_txt like ? and ensembl_gene_id like 'ENSG0%'","%"+query+"%").fetch();	    
	    
            List<Resolution> synonyms;
	    synonyms =  Resolution.find("synonym like ?", "%"+query+"%").fetch();
	    List <Resolution> hits_res = new ArrayList<Resolution>();

	    for (Resolution res: synonyms) {
		if (res.synonym.equals(query) ) {
		    hits_res.add(res);
		    break;
		}
	    }
	    for (Resolution res: synonyms) {
		if (res.synonym.equals(query) )  continue;
		hits_res.add(res);
	    }
	    
	    render(hits_desc, hits_res); // Give the matching IDs for the user to choose, if any
	}
    }

    public static void Result(String id) {
        try {
            GeneQuery q = new GeneQuery(id);
            String ensembl_id = q.getId();
            String species = q.getSpecies();
            String description = q.getDescription();
	    render(ensembl_id, species, description);
        } catch (IllegalArgumentException e) {
            index(); // FIXME: Invalid ensembl id
        }
    }

    public static void tabExons(String id) {
        try {
            GeneQuery q = new GeneQuery(id);
            List<Exon> exons = q.getExons();
	    if (exons.isEmpty())
	        renderText("No exons found for " + id);
            String jv_p = "/Resource/scratch/" + q.getId() + ".exon.pep.afa";
            String jv_n = jv_p.replace("pep","dna");
            generateJalviewFile("/lib/pep_assembler.pl", q.getId(), jv_p);
            generateJalviewFile("/lib/nt_assembler.pl", q.getId(), jv_n);
            render(exons, jv_p, jv_n);		
        } catch (IllegalArgumentException e) {
	    renderText(e);
        }
    }

    public static void tabOrthologs(String id) {
        try {
            GeneQuery q = new GeneQuery(id);
            List<Ortholog> orthologs = q.getOrthologs();
            if (orthologs.isEmpty())
                renderText("No orthologs found for " + id);
            String jv_p = "/Best_MSA/pep/" + orthologs.get(0).ensembl_gene_id + ".afa";
            String jv_n = jv_p.replace("pep","dna");
            String jv_notes = jv_p.replace("pep","notes");
	    jv_notes = jv_notes.replace("afa","txt");
	    render(orthologs, jv_p, jv_n, jv_notes);
        } catch (IllegalArgumentException e) {
	    renderText(e);
        }
    }

    public static void tabParalogs(String id) throws IOException {
        try {
            GeneQuery q = new GeneQuery(id);
            List<String> paralog_ids = q.getParalogIds();
            if (paralog_ids.isEmpty())
                renderText("No paralogs found for " + id);
            Map<String,String> desc_map = new HashMap<String,String>();
            for (String pid: paralog_ids) {
                try {
                    desc_map.put(pid, new GeneQuery(pid).getDescription());
                } catch (IllegalArgumentException e) {
                    System.out.println(e); // Bad entry in paralog table!
                }
            }
            String jv_p = "/Best_MSA/para/" + q.getSpeciesCode() + "/pep/" 
                          + paralog_ids.get(0) + ".afa";
            String jv_n = jv_p.replace("pep","dna");
            render(desc_map, jv_p, jv_n);
        } catch (IllegalArgumentException e) {
	    renderText(e);
        }
    }

    public static void tabSplicings(String id) {
        try {
            GeneQuery q = new GeneQuery(id);
            List<Splicing> splicings = q.getSplicings();
	    if (splicings.isEmpty())
	        renderText("No alternative splicings found for " + id + " in CCDS.");
	    String jv_p = "/Best_MSA/alt/" + q.getSpeciesCode() + "/" + q.getId() + ".afa";
	    render (splicings, jv_p);
        } catch (IllegalArgumentException e) {
	    renderText(e);
        }
    }

    public static void Specialization(@Required List<String> afa_names) throws IOException {
        WS.HttpResponse resultpage;
        String ens1 = afa_names.get(0);
        String ens2 = afa_names.get(1);
        File afa1 = new File(Play.applicationPath + "/Best_MSA/pep/" + ens1 + ".afa");
        File afa2 = new File(Play.applicationPath + "/Best_MSA/pep/" + ens2 + ".afa");
	File temp_afa1 = File.createTempFile("specs_",".afa");
	File temp_afa2 = File.createTempFile("specs_",".afa");

	// Write temp_afa files with no Z and modified sequence names using sed
        runAndSaveStdout(new String[]{"/bin/sed", "/^>/s/$/_"+ens1+"/;t;s/Z/-/g", afa1.getPath()}, temp_afa1);
        runAndSaveStdout(new String[]{"/bin/sed", "/^>/s/$/_"+ens2+"/;t;s/Z/-/g", afa2.getPath()}, temp_afa2);

	// Request page from SPECS
	resultpage =  WS.url(Play.configuration.getProperty("specs_url"))
			.files(new WS.FileParam(temp_afa1, "fnm"), new WS.FileParam(temp_afa2, "fnm"))
			.post();

	temp_afa1.delete();
	temp_afa2.delete();

	// Find job page URL on SPECS result and redirect there
	Pattern p = Pattern.compile(Pattern.quote(Play.configuration.getProperty("specs_url")) + "\\?jobID=\\d+");
	Matcher m = p.matcher(resultpage.getString());
	if (m.find()) {
	    redirect(m.group());
	} else {
	    // No match on webpage returned from SPECS!
	    renderHtml(resultpage.getString());
	}
    }

    public static void SpecsRedirect(String afa_name) throws IOException {
	WS.HttpResponse resultpage;
        File afa = new File(Play.applicationPath + afa_name);
	File temp_msf = File.createTempFile("test123_",".msf");

	// Convert temp_afa file to msf using afa2msf.pl script - why are we doing this?
	// to make the server understand that htese are already aligned
	runAndSaveStdout(new String[]{"lib/afa2msf.pl", afa.getPath()}, temp_msf);

	// Request page from SPECS
	resultpage =  WS.url(Play.configuration.getProperty("specs_url"))
			.files(new WS.FileParam(temp_msf, "fnm"))
			.post();

	temp_msf.delete();

	// Find job page URL on SPECS result and redirect there
	Pattern p = Pattern.compile(Pattern.quote(Play.configuration.getProperty("specs_url")) + "\\?jobID=\\d+");
	Matcher m = p.matcher(resultpage.getString());
	if (m.find()) {
	    redirect(m.group());
	} else {
	    // No match on webpage returned from SPECS!
	    renderText("The Cube server reported an error. Conservation visualization failed.\n" + temp_msf.getName() );
	}
    }

    public static void download_exon(String exon_key, String type) {
	String[] key = exon_key.split("_");
	String gene_id = key[0];
	String exon_id = key[1];
	String is_known = key[2];
	if (!is_known.equals("1")) { is_known = "0"; }
	String species = GeneQuery.findSpecies(gene_id);
	String filename = String.format("%s_%s_%s.fa", gene_id, exon_id, type);

	try {
	    File exonfile = File.createTempFile("exon_",".fa");

	    // Run alignment script to produce aligned exons
	    String cmd = String.format("lib/exon_alignment.py %s %s %s %s",
			     exon_id, is_known, species, exonfile.getPath());
	    if (type.equals("N")) { cmd += " nt"; } // Add nt argument for nuc align
	    Process ps = Runtime.getRuntime().exec(cmd);
	    ps.waitFor();

	    InputStream is = new FileInputStream(exonfile);
	    exonfile.delete();

	    renderBinary(is, filename, "text/plain", true);
	} catch (Exception e) { System.out.println(e); }
    }

    //////functions that are called by application entry points//////

    private static final Pattern fasta = Pattern.compile(
	"(?m)\\A(>.*$)?[a-zA-Z-*\\s]+\\Z");

    private static boolean InDatabaseRange(String pro_id) {
	//checks if the length of sequence is in the correct range
	return (pro_id.length() < 33424 && pro_id.length() > 6);
    }

    private static boolean generateJalviewFile(String script, String id, String filename) {
	// Run script against exolocator_db and id to produce output at filename
	boolean success = true;
	String script_path = Play.applicationPath + script;
	String db_name = "exolocator_db";
	String file_path = Play.applicationPath + filename;
	String[] cmd = {script_path, db_name, id, file_path};

	try {
	    Process ps = Runtime.getRuntime().exec(cmd);
	    ps.waitFor();
	} catch (Exception e) {
	    success = false;
	    System.out.println(e);
	}

	return success;
    }

    private static String blastall(String pro_id, String database) {
	String out = null;
	try {
	    // Execute blastall specifying database and tabbed output
	    Process ps = Runtime.getRuntime().exec(
                Play.configuration.getProperty("blastall_path") +
                " -p blastp -m 8 -d " + database);
            // Send pro_id to blastall stdin
            BufferedWriter in = new BufferedWriter(
                new OutputStreamWriter(ps.getOutputStream()));
            in.write(pro_id);
            in.close();
            // Read blastall stdout as a string and return it
            out = new Scanner(ps.getInputStream(),"UTF-8")
            			.useDelimiter("\\A").next();
	} catch (IOException e) {
	    System.out.println(e);
	}
        return out;
    }

    private static boolean runAndSaveStdout(String[] cmd, File target) {
	// Run a command given by cmd and save stdout to target
	boolean success = true;
	try {
	    Process ps = Runtime.getRuntime().exec(cmd);
	    BufferedReader rd = new BufferedReader(new InputStreamReader(ps.getInputStream()));
	    BufferedWriter wr = new BufferedWriter(new FileWriter(target));
	    String line;
	    while ((line = rd.readLine()) != null) {
        	wr.write(line);
        	wr.newLine();
	    }
	    rd.close();
            wr.close();
	} catch (Exception e) {
	    success = false;
	    System.out.println(e);
	}
	return success;
    }
    
} 
