package org.monarchinitiative.dsppc;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.monarchinitiative.phenol.base.PhenolException;
import org.monarchinitiative.phenol.formats.hpo.HpoDisease;
import org.monarchinitiative.phenol.io.OntologyLoader;
import org.monarchinitiative.phenol.io.obo.hpo.HpoDiseaseAnnotationParser;
import org.monarchinitiative.phenol.ontology.data.Ontology;
import org.monarchinitiative.phenol.ontology.data.TermId;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class Dsppc {
    static final String ALL_GENES_FILENAME = "human_protein_coding_genes.tsv";
    static final String GENE_SETS_FILENAME = "ENTREZ_gene_sets.tsv";
    static final String HPO_FILENAME = "hp.obo";
    static final String HPOA_FILENAME = "phenotype.hpoa";
    static final String MIM2GENE_MEDGEN_FILENAME = "mim2gene_medgen";
    static final String REPORT_DIRECTORY = "results";
    static final String REPORT_FILENAME = "counterOutput.txt";

    private static final Logger logger = LogManager.getLogger();

    /**
     * Checks path string and adds a final separator character if not already there.
     * @param path       String containing path as user typed it on command line
     * @return String    path with final separator added if it was not already there
     */
    private static String fixFinalSeparator(String path) {
        return path.endsWith(File.separator) ? path : path + File.separator;
    }

    /**
     * Parses the ALL_GENES_FILENAME file to create a set of genes identified by ENTREZ id.
     * @param inputPath    path to ALL_GENES_FILENAME file
     * @return set of ENTREZ gene TermIds
     * @throws IOException if cannot open/read input file
     */
    static List<TermId> parseAllGenes(String inputPath) throws IOException {
        List<TermId> allGenes = new ArrayList<>();
        BufferedReader br = new BufferedReader(new FileReader(inputPath));
        String[] fields;
        String line;
        String entrezprefix = "ENTREZ";

        while ((line = br.readLine()) != null) {
            fields = line.split("\t");
            allGenes.add(TermId.of(entrezprefix, fields[0]));
        }
        br.close();
        return allGenes;
    }

    /**
     * Parses three optional command line arguments:
     * -- directory in which to find input files (defaults to "src/main/resources/")
     * -- minimum number of diseases to decide whether a phenotype should be considered in the
     * similarity function (defaults to 1, no filtering)
     * -- threshold for IC content of phenotype to be included in thresholded similarity calculation
     */
    private static CommandLine parseCommandLineArgs(String[] args) {
        // create the command line Options
        Option helpOpt = Option.builder("h")
                .longOpt("help")
                .required(false)
                .hasArg(false)
                .build();
        Option dataPathOpt = Option.builder("d")
                .longOpt("dataDir")
                .desc("directory containing input files")
                .hasArg()
                .optionalArg(false)
                .argName("directory")
                .required(false)
                .build();
        Option minDisOpt = Option.builder("m")
                .longOpt("minDiseases")
                .desc("minimum number of diseases for phenotype")
                .hasArg()
                .optionalArg(false)
                .argName("min")
                .required(false)
                .build();
        Option thresholdOpt = Option.builder("t")
                .longOpt("ICthreshold")
                .desc("threshold on information content")
                .hasArg()
                .optionalArg(false)
                .argName("threshold")
                .required(false)
                .build();
        Options options = new Options();
        options.addOption(dataPathOpt);
        options.addOption(helpOpt);
        options.addOption(minDisOpt);
        options.addOption(thresholdOpt);

        // create the command line parser and help formatter
        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        try {
            // parse the command line looking for help option
            CommandLine cmdl = parser.parse(options, args);
            if (cmdl.hasOption("h")) {
                // automatically generate usage information, write to System.out
                formatter.printHelp("java -jar dsppc-1.0-SNAPSHOT.jar", options);
                return null;
            } else {
                return cmdl;
            }
        } catch (ParseException e) {
            formatter.printHelp("java -jar dsppc-1.0-SNAPSHOT.jar", options);
            return null;
        }
    }

    static void parseGeneSets(String inputPath, Set<TermId> pathwayGenes, Set<TermId> anchoredGenes)
            throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(inputPath));
        ArrayList<TermId> geneIds = new ArrayList<>();
        String entrezprefix = "ENTREZ";

        // skip over first line of file, which is a header line
        String line = br.readLine();
        while ((line = br.readLine()) != null) {
            if (line.equals("")) {
                // blank line marks end of the first set of genes
                pathwayGenes.addAll(geneIds);
                geneIds.clear();
            } else {
                for (String s : line.split("\t")) {
                    geneIds.add( TermId.of(entrezprefix, s));
                }
            }
        }
        // end of file marks end of the second set of genes
        anchoredGenes.addAll(geneIds);
        br.close();
    }

    /**
     * Parses the phenotype.hpoa file.
     * @param inputPath      path to phenotype.hpoa file
     * @param hpo            ontology object
     * @return map from OMIM TermId to corresponding HpoDisease object
     * @throws PhenolException if error in parsing phenotype.hpoa file
     */
    static Map<TermId, HpoDisease> parseHPOA(String inputPath, Ontology hpo) throws PhenolException {
        HpoDiseaseAnnotationParser parser = new HpoDiseaseAnnotationParser(inputPath, hpo);
        Map<TermId, HpoDisease> diseaseMap = parser.parse();
        Map<TermId, HpoDisease> omimMap = new HashMap<>();
        diseaseMap.forEach((tid, disease) ->
        {
            if (disease.getDatabase().equals("OMIM")) omimMap.put(tid, disease);
        });
        // There should be roughly 7000 OMIM entries left after filtering.
        return omimMap;
    }

    /**
     * Parses the mim2gene_medgen file to find OMIM diseases for each ENTREZ gene mentioned
     * @param inputPath       path to mim2gene_medgen file
     * @return map from ENTREZ gene TermId to set of OMIM disease TermIds
     * @throws IOException if cannot open/read mim2gene_medgen file
     */
    static Map<TermId, Set<TermId>> parseMedgen(String inputPath) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(inputPath));
        String diseasePrefix = "OMIM";
        Set<TermId> diseases;
        String[] fields;
        TermId diseaseId;
        TermId geneId;
        final Map<TermId, Set<TermId>> geneIdToDiseaseIds = new HashMap<>();
        String genePrefix = "ENTREZ";

        // skip over first line of file, which is a header line
        String line = br.readLine();
        while ((line = br.readLine()) != null) {
            fields = line.split("\t");
            if (fields[2].equals("phenotype") &&
                    !(fields[1].equals("-") || fields[5].equals("nondisease"))) {
                diseaseId =  TermId.of(diseasePrefix, fields[0]);
                geneId = TermId.of(genePrefix, fields[1]);
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

    public static void main(String[] args) {
        final CommandLine cmdl = parseCommandLineArgs(args);
        final List<TermId> allGenes;
        final Map<TermId, HpoDisease> diseaseMap;
        final Map<TermId, Set<TermId>> geneToDiseasesMap;
        final Set<TermId> gpiAnchoredGenes = new TreeSet<>();
        final Set<TermId> gpiPathwayGenes = new TreeSet<>();
        final Ontology hpo;
        final String dataDir;
        final int minDiseases;
        final double threshold;

        if (cmdl != null) {
            dataDir = fixFinalSeparator(cmdl.getOptionValue("d", "src/main/resources/"));
            minDiseases = Integer.parseInt(cmdl.getOptionValue("m", "1"));
            threshold = Double.parseDouble(cmdl.getOptionValue("t", "0.5"));
            hpo= OntologyLoader.loadOntology(new File(dataDir + HPO_FILENAME));
            try {
                logger.info("DONE: Loading HPO");
                diseaseMap = parseHPOA(dataDir + HPOA_FILENAME, hpo);
                logger.info("DONE: Parsing HPO disease annotations");
                geneToDiseasesMap = parseMedgen(dataDir + MIM2GENE_MEDGEN_FILENAME);
                logger.info("DONE: Parsing gene to disease annotations");
                allGenes = parseAllGenes(dataDir + ALL_GENES_FILENAME);
                logger.info("DONE: Parsing all human protein-coding genes");
                parseGeneSets(dataDir + GENE_SETS_FILENAME, gpiPathwayGenes, gpiAnchoredGenes);
                logger.info("DONE: Parsing gene sets");

                new ComputeSimilarity(hpo, diseaseMap, geneToDiseasesMap, allGenes,
                        gpiPathwayGenes, gpiAnchoredGenes).run(minDiseases, threshold);
            } catch (IOException e) {
                logger.fatal("Fatal error in file I/O. ", e);
            } catch (PhenolException e) {
                logger.fatal("Fatal error parsing inputs to dsppc. ", e);
            }
        }
    }
}
