package org.monarchinitiative.dsppc;

/*
 * created 06 Feb 2019
 */

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

import org.monarchinitiative.phenol.formats.hpo.HpoDisease;
import org.monarchinitiative.phenol.ontology.data.Ontology;
import org.monarchinitiative.phenol.ontology.data.Term;
import org.monarchinitiative.phenol.ontology.data.TermId;

import static org.monarchinitiative.dsppc.ComputeSimilarity.NUM_ITER;
import static org.monarchinitiative.dsppc.Dsppc.REPORT_DIRECTORY;
import static org.monarchinitiative.dsppc.Dsppc.REPORT_FILENAME;
import static org.monarchinitiative.dsppc.SimFuns.randomSample;

/**
 * The Counter class provides utility functions to compute averages over groups of randomly chosen gene sets.
 * Finds the average for three aspects:
 *     number of disease genes per set
 *     number of diseases related to those disease genes
 *     number of phenotypes related to those diseases
 * Can output the counts and corresponding TermIds to results file specified in Dsppc.java.
 */
class Counter {
    // all human protein-coding genes (excluding GPI pathway genes)
    private final List<TermId> allGenes;
    // map from ENTREZ id to gene name
    private final Map<TermId, String> allGenesMap;
    // average number of diseases associated with genes in random sample
    private double avgDiseases = -1.0;
    // average number of disease genes per random sample
    private double avgDiseaseGenes = -1.0;
    // average number of phenotypes associated with diseases annotated to genes in random sample
    private double avgPhenotypes = -1.0;
    // map from disease ID to HPO disease object
    private final Map<TermId, HpoDisease> diseaseMap;
    // map from gene ID to set of disease IDs
    private final Map<TermId, Set<TermId>> genesToDiseasesMap;
    // ontology to look up names for TermIds
    private final Ontology hpo;
    // output file
    private final File resultsFile;

    // types of things that show up in a report
    private enum EntityType {
        DISEASE, GENE, PHENOTYPE
    }

    // array indices for the array of counts
    static final int NUM_DISEASE_GENES = 0;
    static final int NUM_DISEASES = 1;
    static final int NUM_PHENOTYPES = 2;

    Counter(Ontology hpo, List<TermId> allGenes, Map<TermId, String> allGenesMap,
            Map<TermId, HpoDisease> diseaseMap, Map<TermId, Set<TermId>> genesToDiseasesMap)
            throws IOException {
        this.allGenes = allGenes;
        this.allGenesMap = allGenesMap;
        this.diseaseMap = diseaseMap;
        this.genesToDiseasesMap = genesToDiseasesMap;
        this.hpo = hpo;

        File resultsDir = new File(REPORT_DIRECTORY);
        if (!resultsDir.exists()) {
            resultsDir.mkdirs();
        }
        resultsFile = new File(resultsDir, REPORT_FILENAME);
        if (resultsFile.exists()) {
            // reset file pointer to overwrite previous contents of file
            try (PrintWriter pw = new PrintWriter(resultsFile)) {}
        }
        else {
            resultsFile.createNewFile();
        }
    }

    /**
     * Finds averages over NUM_ITER sets of genes
     * each set chosen at random from the set of all genes
     * average for three counts:
     *     number of disease genes in the set
     *     number of distinct diseases related to those disease genes
     *     number of distinct phenotypes related to those diseases
     * @param cardinality  how many genes per set
     */
    void countAverages(int cardinality) throws IOException {
        int[] counts;
        Set<TermId> randomGenes = null;
        Random rand = new Random();
        int sumDiseases = 0;
        int sumDiseaseGenes = 0;
        int sumPhenotypes = 0;
        for (int i = 0; i < NUM_ITER; i++) {
            randomGenes = randomSample(rand, cardinality, allGenes);
            counts = countOneSet(randomGenes);
            sumDiseaseGenes += counts[NUM_DISEASE_GENES];
            sumDiseases += counts[NUM_DISEASES];
            sumPhenotypes += counts[NUM_PHENOTYPES];
        }
        countAndReport(randomGenes, "Random");
        double denom = (double) NUM_ITER;
        avgDiseaseGenes = sumDiseaseGenes / denom;
        avgDiseases = sumDiseases / denom;
        avgPhenotypes = sumPhenotypes / denom;
    }

    /**
     * For the input set of genes, counts:
     *     number of disease genes in the set
     *     number of diseases related to those disease genes
     *     number of phenotypes related to those diseases
     * @param genes  Set of gene TermIds
     * @return array of int containing the three counts
     */
    int[] countOneSet(Set<TermId> genes) {
        int[] counts = new int[3];
        Set<TermId> diseaseGenes = new HashSet<>(genes);
        diseaseGenes.retainAll(genesToDiseasesMap.keySet());
        counts[NUM_DISEASE_GENES] = diseaseGenes.size();
        Set<TermId> diseases = diseaseGenes.stream()
                .flatMap(g -> genesToDiseasesMap.get(g).stream())
                .collect(Collectors.toSet());
        counts[NUM_DISEASES] = diseases.size();
        Set<TermId> phenotypes = diseases.stream()
                .map(diseaseMap::get)
                .filter(Objects::nonNull)
                .flatMap(d -> d.getPhenotypicAbnormalityTermIdList().stream())
                .collect(Collectors.toSet());
        counts[NUM_PHENOTYPES] = phenotypes.size();
        return counts;
    }

    /**
     * Same functionality as countOneSet except that this method outputs a report of the
     * disease genes, diseases, and phenotypes for the input set of genes. The name of
     * the report file is specified in Dsppc.java
     * @param genes         Set of gene TermIds
     * @param whichGenes    String that identifies which set of genes
     * @return array of int containing the three counts: disease genes, diseases, phenotypes
     */
    int[] countAndReport(Set<TermId> genes, String whichGenes) throws IOException {
        int[] counts = new int[3];
        Set<TermId> diseaseGenes = new HashSet<>(genes);
        diseaseGenes.retainAll(genesToDiseasesMap.keySet());
        counts[NUM_DISEASE_GENES] = diseaseGenes.size();
        Set<TermId> diseases = diseaseGenes.stream()
                .flatMap(g -> genesToDiseasesMap.get(g).stream())
                .collect(Collectors.toSet());
        counts[NUM_DISEASES] = diseases.size();
        Set<TermId> phenotypes = diseases.stream()
                .map(diseaseMap::get)
                .filter(Objects::nonNull)
                .flatMap(d -> d.getPhenotypicAbnormalityTermIdList().stream())
                .collect(Collectors.toSet());
        counts[NUM_PHENOTYPES] = phenotypes.size();

        reportOneSet(resultsFile, genes, EntityType.GENE, whichGenes + " genes");
        reportOneSet(resultsFile, diseaseGenes, EntityType.GENE, whichGenes + " disease genes");
        reportOneSet(resultsFile, diseases, EntityType.DISEASE, whichGenes + " diseases");
        reportOneSet(resultsFile, phenotypes, EntityType.PHENOTYPE, whichGenes + " phenotypes");
        return counts;
    }

    private SortedMap<String, TermId> findNames(Set<TermId> tids, EntityType typ) {
        TreeMap<String, TermId> names = new TreeMap<>();

        switch (typ) {
            case DISEASE:
                for (TermId tid : tids) {
                    if (diseaseMap.containsKey(tid)) {
                        names.put(diseaseMap.get(tid).getName(), tid);
                    } else {
                        names.put(tid.getId(), tid);
                    }
                }
                break;
            case GENE:
                for (TermId tid : tids) {
                    if (allGenesMap.containsKey(tid)) {
                        names.put(allGenesMap.get(tid), tid);
                    } else {
                        names.put(tid.getId(), tid);
                    }
                }
                break;
            case PHENOTYPE:
                Map<TermId, Term> termMap = hpo.getTermMap();
                for (TermId tid : tids) {
                    names.put(termMap.get(tid).getName(), tid);
                }
        }
        return names;
    }

    double getAvgDiseases() {
        return avgDiseases;
    }

    double getAvgDiseaseGenes() {
        return avgDiseaseGenes;
    }

    double getAvgPhenotypes() {
        return avgPhenotypes;
    }

    /**
     * Outputs the term ids in the set to report file specified in Dsppc.java, appending to previous
     * contents of that file.
     * @param file    File for output
     * @param tids    Set of TermIds
     * @param typ     EntityType that corresponds to the tids
     * @param header  String to say what type of TermIds these are (gene, disease, phenotype)
     *                and where they come from (GPI pathway, GPI anchored, ...)
     */
    private void reportOneSet(File file, Set<TermId> tids, EntityType typ, String header)
            throws IOException {
        Map<String, TermId> names = findNames(tids, typ);
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(file, true))) {
            bw.write(header + " : " + tids.size());
            bw.newLine();
            for (Map.Entry<String, TermId> e : names.entrySet()) {
                bw.write(String.format("%10s   %s%n", e.getValue().getId(), e.getKey()));
            }
            bw.newLine();
        }
    }
}
