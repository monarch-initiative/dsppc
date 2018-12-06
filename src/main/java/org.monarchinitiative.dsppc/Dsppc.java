package org.monarchinitiative.dsppc;

import com.google.common.collect.ImmutableSet;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.monarchinitiative.phenol.base.PhenolException;
import org.monarchinitiative.phenol.formats.hpo.HpoDisease;
import org.monarchinitiative.phenol.formats.hpo.HpoOntology;
import org.monarchinitiative.phenol.io.obo.hpo.HpOboParser;
import org.monarchinitiative.phenol.io.obo.hpo.HpoDiseaseAnnotationParser;
import org.monarchinitiative.phenol.ontology.data.TermId;
import org.monarchinitiative.phenol.ontology.data.TermPrefix;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

public class Dsppc {
    private static final String GENE_SETS_FILENAME = "src/main/resources/ENTREZ_gene_sets.tsv";
    private static final String HPO_FILENAME = "src/main/resources/hp.obo";
    private static final String HPOA_FILENAME = "/src/main/resources/phenotype.hpoa";
    private final static String MIM2GENE_MEDGEN_FILENAME = "src/main/resources/mim2gene_medgen";

    private static final Logger logger = LogManager.getLogger();

    /*
    collector <- c()
gene.sets <- list()
for (i in 2:length(lines.read)) {
  if (lines.read[i] == c("")) {
    gene.sets[[1]] <- as.set(unlist(lapply(collector, as.integer)))
    collector <- c()
  }
  else {
    collector <- c(collector, unlist(strsplit(lines.read[i], "\t")))
  }
}
gene.sets[[2]] <- as.set(unlist(lapply(collector, as.integer)))

  public Map<TermId, MpGene> parseMarkers() throws IOException, PhenolException {
    ImmutableMap.Builder<TermId, MpGene> bld = ImmutableMap.builder();
    BufferedReader br = new BufferedReader(new FileReader(mgiMarkerPath));
    // skip over first line of file, which is a header line
    String line = br.readLine();
    while ((line=br.readLine()) != null) {
      String[] fields = line.split("\t");
      // first field is MGI Accession ID, seventh is Marker Symbol, tenth is Marker Type
      //String mgiId = fields[0];
      TermId mgiId = TermId.constructWithPrefix(fields[0]);
        bld.put(mgiId, createMpGene(mgiId, fields[6], fields[9]));
    }
    br.close();
    return bld.build();
  }

     */
    private static void parseGeneSets(Set<TermId> pathwayGenes, Set<TermId> anchoredGenes)
            throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(GENE_SETS_FILENAME));
        ArrayList<TermId> geneIds = new ArrayList<>();
        TermPrefix prefix = new TermPrefix("ENTREZ");

        // skip over first line of file, which is a header line
        String line = br.readLine();
        while ((line = br.readLine()) != null) {
            if (line.equals("")) {
                // blank line marks end of the first set of genes
                pathwayGenes.addAll(geneIds);
                geneIds.clear();
            } else {
                for (String s : line.split("\t")) {
                    geneIds.add(new TermId(prefix, s));
                }
            }
        }
        // end of file marks end of the second set of genes
        anchoredGenes.addAll(geneIds);
    }

    private static void parseHPO() {
        final HpoOntology hpo;
        try {
            hpo = new HpOboParser(new File(HPO_FILENAME)).parse();
        } catch (IOException | PhenolException e) {
            e.printStackTrace();
            System.exit(1);
            return; // javac complains otherwise
        }
        System.out.println("DONE: Loading HPO");
        final Map<TermId, HpoDisease> diseaseMap;
        try {
            HpoDiseaseAnnotationParser parser = new HpoDiseaseAnnotationParser(HPOA_FILENAME, hpo);
            diseaseMap = parser.parse();
        } catch (PhenolException e) {
            e.printStackTrace();
            System.exit(1);
            return; // javac complains otherwise
        }

    }

    public static void main (String[] args) {
        final Map<TermId, HpoDisease> diseaseMap;
        final Set<TermId> gpiAnchoredGenes = new TreeSet<>();
        final Set<TermId> gpiPathwayGenes = new TreeSet<>();
        final HpoOntology hpo;
        //final
        try {
            hpo = new HpOboParser(new File(HPO_FILENAME)).parse();
            logger.info("DONE: Loading HPO");
            HpoDiseaseAnnotationParser parser = new HpoDiseaseAnnotationParser(HPOA_FILENAME, hpo);
            diseaseMap = parser.parse();
            parseGeneSets(gpiPathwayGenes, gpiAnchoredGenes);
            new ComputeSimilarityDemo(hpo).run();
        } catch (IOException | PhenolException e) {
            logger.fatal("Fatal error parsing inputs to similarity function. ", e);
        }
    }
}
