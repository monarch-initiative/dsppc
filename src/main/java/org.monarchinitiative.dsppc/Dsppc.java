package org.monarchinitiative.dsppc;

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
import java.util.*;

public class Dsppc {
    private static final String GENE_SETS_FILENAME = "src/main/resources/ENTREZ_gene_sets.tsv";
    private static final String HPO_FILENAME = "src/main/resources/hp.obo";
    private static final String HPOA_FILENAME = "src/main/resources/phenotype.hpoa";
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
        br.close();
    }

    private static Map<TermId, Set<TermId>> parseMedgen() throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(MIM2GENE_MEDGEN_FILENAME));
        TermPrefix diseasePrefix = new TermPrefix("OMIM");
        Set<TermId> diseases;
        String[] fields;
        TermId diseaseId;
        TermId geneId;
        final Map<TermId, Set<TermId>> geneIdToDiseaseIds = new HashMap<>();
        TermPrefix genePrefix = new TermPrefix("ENTREZ");

        // skip over first line of file, which is a header line
        String line = br.readLine();
        while ((line = br.readLine()) != null) {
            fields = line.split("\t");
            if (fields[2].equals("phenotype") && !(fields[1].equals("-"))) {
                diseaseId = new TermId(diseasePrefix, fields[0]);
                geneId = new TermId(genePrefix, fields[1]);
                diseases = geneIdToDiseaseIds.get(geneId);
                if (diseases == null) {
                    // this is the first disease for this gene
                    diseases = new TreeSet<>();
                    diseases.add(diseaseId);
                    geneIdToDiseaseIds.put(geneId, diseases);
                } else {
                    // already saw one or more diseases for this gene, just add the current one
                    diseases.add(diseaseId);
                }
            }
        }
        br.close();
        return geneIdToDiseaseIds;
    }

    /*
     * Optional command line argument to specific minimum number of diseases to decide whether a
     * phenotype should be considered in the similarity function (defaults to 0, no minimum)
     */
    public static void main (String[] args) {
        final Map<TermId, HpoDisease> diseaseMap;
        final Map<TermId, Set<TermId>> geneToDiseasesMap;
        final Set<TermId> gpiAnchoredGenes = new TreeSet<>();
        final Set<TermId> gpiPathwayGenes = new TreeSet<>();
        final HpoOntology hpo;
        final int minDiseases;

        minDiseases = args.length > 0 ? Integer.parseInt(args[0]) : 0;
        try {
            hpo = new HpOboParser(new File(HPO_FILENAME)).parse();
            logger.info("DONE: Loading HPO");
            HpoDiseaseAnnotationParser parser = new HpoDiseaseAnnotationParser(HPOA_FILENAME, hpo);
            diseaseMap = parser.parse();
            logger.info("DONE: Parsing HPO disease annotations");
            geneToDiseasesMap = parseMedgen();
            logger.info("DONE: Parsing gene to disease annotations");
            parseGeneSets(gpiPathwayGenes, gpiAnchoredGenes);
            logger.info("DONE: Parsing gene sets");
            new ComputeSimilarity(hpo, diseaseMap, geneToDiseasesMap,
                    gpiPathwayGenes, gpiAnchoredGenes).run(minDiseases);
        } catch (IOException | PhenolException e) {
            logger.fatal("Fatal error parsing inputs to similarity function. ", e);
        }
    }
}
