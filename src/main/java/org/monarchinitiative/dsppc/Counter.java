package org.monarchinitiative.dsppc;

/*
 * created 06 Feb 2019
 */

import java.util.*;
import java.util.stream.Collectors;

import org.monarchinitiative.phenol.formats.hpo.HpoDisease;
import org.monarchinitiative.phenol.ontology.data.TermId;

import static org.monarchinitiative.dsppc.ComputeSimilarity.NUM_ITER;
import static org.monarchinitiative.dsppc.SimFuns.randomSample;

/**
 * The Counter class provides utility functions to compute averages over groups of randomly chosen gene sets.
 * Finds the average for three aspects:
 *     number of disease genes per set
 *     number of diseases related to those disease genes
 *     number of phenotypes related to those diseases
 */
class Counter {
    // all human protein-coding genes (excluding GPI pathway genes)
    private final List<TermId> allGenes;
    // average number of diseases associated with genes in random sample
    private double avgDiseases;
    // average number of disease genes per random sample
    private double avgDiseaseGenes;
    // average number of phenotypes associated with diseases annotated to genes in random sample
    private double avgPhenotypes;
    // map from disease ID to HPO disease object
    private final Map<TermId, HpoDisease> diseaseMap;
    // map from gene ID to set of disease IDs
    private final Map<TermId, Set<TermId>> genesToDiseasesMap;
    // set of genes for proteins that normally are anchored by the GPI pathway anchor

    // array indices for the array of counts
    static final int NUM_DISEASE_GENES = 0;
    static final int NUM_DISEASES = 1;
    static final int NUM_PHENOTYPES = 2;

    Counter(List<TermId> allGenes, Map<TermId, HpoDisease> diseaseMap,
            Map<TermId, Set<TermId>> genesToDiseasesMap, int cardinality) {
        this.allGenes = allGenes;
        this.diseaseMap = diseaseMap;
        this.genesToDiseasesMap = genesToDiseasesMap;
        countAverages(cardinality);
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
    private void countAverages(int cardinality) {
        int[] counts;
        Set<TermId> randomGenes;
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

    double getAvgDiseases() {
        return avgDiseases;
    }

    double getAvgDiseaseGenes() {
        return avgDiseaseGenes;
    }

    double getAvgPhenotypes() {
        return avgPhenotypes;
    }
}
