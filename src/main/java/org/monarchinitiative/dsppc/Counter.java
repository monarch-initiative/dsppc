package org.monarchinitiative.dsppc;

/*
 * created 06 Feb 2019
 */

import java.io.*;
import java.util.*;
import java.util.stream.Collectors;

import org.monarchinitiative.phenol.formats.hpo.HpoDisease;
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
    // output file
    private final File resultsFile;

    // array indices for the array of counts
    static final int NUM_DISEASE_GENES = 0;
    static final int NUM_DISEASES = 1;
    static final int NUM_PHENOTYPES = 2;

    Counter(List<TermId> allGenes, Map<TermId, HpoDisease> diseaseMap,
            Map<TermId, Set<TermId>> genesToDiseasesMap) throws IOException {
        this.allGenes = allGenes;
        this.diseaseMap = diseaseMap;
        this.genesToDiseasesMap = genesToDiseasesMap;

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

        reportOneSet(resultsFile, genes, whichGenes + " genes");
        reportOneSet(resultsFile, diseaseGenes, whichGenes + " disease genes");
        reportOneSet(resultsFile, diseases, whichGenes + " diseases");
        reportOneSet(resultsFile, phenotypes, whichGenes + " phenotypes");
        return counts;
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
     * @param tids    Set of TermIds
     * @param header  String to say what type of TermIds these are (gene, disease, phenotype, ...)
     */
    private void reportOneSet(File file, Set<TermId> tids, String header) throws IOException {
        TermId[] tidsArray = tids.toArray(new TermId[0]);
        Arrays.sort(tidsArray);
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(file, true))) {
            bw.write(header + " : " + tids.size());
            bw.newLine();
            for (TermId tid : tidsArray) {
                bw.write(tid.getId() + "\t");
            }
            bw.newLine();
            bw.newLine();
        }
    }
}
