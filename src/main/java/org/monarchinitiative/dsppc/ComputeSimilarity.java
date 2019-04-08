package org.monarchinitiative.dsppc;

import java.io.IOException;
import java.util.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.monarchinitiative.phenol.formats.hpo.HpoDisease;
import org.monarchinitiative.phenol.ontology.algo.InformationContentComputation;
import org.monarchinitiative.phenol.ontology.data.Ontology;
import org.monarchinitiative.phenol.ontology.data.TermId;
import org.monarchinitiative.phenol.ontology.data.TermIds;
import org.monarchinitiative.phenol.ontology.similarity.PrecomputingPairwiseResnikSimilarity;
import org.monarchinitiative.phenol.ontology.similarity.ResnikSimilarity;

import static org.monarchinitiative.dsppc.Counter.*;
import static org.monarchinitiative.dsppc.SimFuns.*;

// adapted from ComputeSimilarity.java in phenol library

/**
 * App for computing similarity scores between gene (by entrez ID) and HPO term list.
 *
 * @author <a href="mailto:manuel.holtgrewe@bihealth.de">Manuel Holtgrewe</a>
 * @author <a href="mailto:peter.robinson@jax,org">Peter Robinson</a>
 */
class ComputeSimilarity {

    private final List<TermId> allDiseaseGenes;
    private final List<TermId> allGenes;
    private final Map<TermId, HpoDisease> diseaseMap;
    private final Map<TermId, Set<TermId>> genesToDiseasesMap;
    private final Set<TermId> gpiAnchoredGenes;
    private final Set<TermId> gpiPathwayGenes;
    private final Ontology hpo;
    private final Map<TermId, Double> icMap;

    private static final Logger logger = LogManager.getLogger();
    static final int NUM_ITER = 1000;
    private static final int NUM_THREADS = 4;

    ComputeSimilarity(Ontology hpo, Map<TermId, HpoDisease> diseaseMap,
                      Map<TermId, Set<TermId>> geneToDiseasesMap, List<TermId> allGenes,
                      Set<TermId> gpiPathway, Set<TermId> gpiAnchored) throws IOException {
        this.hpo = hpo;
        this.diseaseMap = diseaseMap;
        this.genesToDiseasesMap = geneToDiseasesMap;
        this.allGenes = allGenes;
        gpiAnchoredGenes = new HashSet<>(gpiAnchored);
        gpiPathwayGenes = new HashSet<>(gpiPathway);
        // delete the GPI pathway genes from the set over which we will randomly choose genes
        // to compare to the GPI pathway genes
        this.allGenes.removeAll(gpiPathwayGenes);
        compareCounts();
        // for randomization, consider only disease genes and eliminate genes which are not associated with
        // any disease
        allDiseaseGenes = new ArrayList<>(this.allGenes);
        System.err.println("All genes: " + allGenes.size() + " ; all GPI pathway genes: " +
                gpiPathwayGenes.size() + " ; all GPI anchored genes: " + gpiAnchoredGenes.size());
        allDiseaseGenes.retainAll(genesToDiseasesMap.keySet());
        gpiAnchoredGenes.retainAll(genesToDiseasesMap.keySet());
        gpiPathwayGenes.retainAll(genesToDiseasesMap.keySet());
        System.err.println("All disease genes: " + allDiseaseGenes.size() +
                " ; all GPI pathway disease genes: " + gpiPathwayGenes.size() +
                " ; all GPI anchored disease genes: " + gpiAnchoredGenes.size());
        icMap = computeICmap();
    }

    /**
     * Outputs the number of disease genes, number of diseases, and number of phenotypes for the set of
     * GPI anchored genes and averaged over NUM_ITER sets of the same size chosen at random from allGenes.
     */
    private void compareCounts() throws IOException {
        Counter counter = new Counter(hpo, allGenes, diseaseMap, genesToDiseasesMap);
        int[] pathwayCounts = counter.countAndReport(gpiPathwayGenes, "GPI pathway");
        int[] anchoredCounts = counter.countAndReport(gpiAnchoredGenes, "GPI anchored");
        counter.countAverages(gpiAnchoredGenes.size());
        System.err.println("                   # Disease Genes    # Diseases    # Phenotypes");
        System.err.println(String.format("GPI pathway genes\t\t%6.2f\t\t%6.2f\t\t%6.2f",
                (float) pathwayCounts[NUM_DISEASE_GENES], (float) pathwayCounts[NUM_DISEASES],
                (float) pathwayCounts[NUM_PHENOTYPES]));
        System.err.println(String.format("GPI anchored genes\t\t%6.2f\t\t%6.2f\t\t%6.2f",
                (float) anchoredCounts[NUM_DISEASE_GENES], (float) anchoredCounts[NUM_DISEASES],
                (float) anchoredCounts[NUM_PHENOTYPES]));
        System.err.println(String.format("Average of random\t\t%6.2f\t\t%6.2f\t\t%6.2f",
                counter.getAvgDiseaseGenes(), counter.getAvgDiseases(), counter.getAvgPhenotypes()));
    }

    /**
     * Compute information content of each phenotype term associated with one or more genes.
     * @return map from phenotype TermId to its information content.
     */
    private Map<TermId, Double> computeICmap() {
        /*
        // Compute mapping from OMIM ID to phenotype TermIds and from phenotype to OMIM termId.
        logger.info("Mapping from OMIM ID to HPO phenotype terms and the reverse");
        final Map<TermId, Collection<TermId>> diseaseIdToTermIds = new HashMap<>();
        final Map<TermId, Collection<TermId>> termIdToDiseaseIds = new HashMap<>();

        for (TermId tid0 : diseaseMap.keySet()) {
            HpoDisease disease = diseaseMap.get(tid0);
            List<TermId> hpoTerms = disease.getPhenotypicAbnormalityTermIdList();
            diseaseIdToTermIds.putIfAbsent(tid0, new HashSet<>());
            // add term anscestors
            final Set<TermId> inclAncestorTermIds =
                    TermIds.augmentWithAncestors(hpo, Sets.newHashSet(hpoTerms), true);

            for (TermId tid1 : inclAncestorTermIds) {
                termIdToDiseaseIds.putIfAbsent(tid1, new HashSet<>());
                termIdToDiseaseIds.get(tid1).add(tid0);
                diseaseIdToTermIds.get(tid0).add(tid1);
            }
        }
        */
        final Map<TermId, Collection<TermId>> termIdToGeneIds = new HashMap<>();
        Set<TermId> phenoSet;
        for (TermId tid0 : genesToDiseasesMap.keySet()) {
            phenoSet = diseasesToPhenotypes(1, genesToDiseasesMap.get(tid0));
            // add term anscestors
            final Set<TermId> inclAncestorTermIds =
                    TermIds.augmentWithAncestors(hpo, phenoSet, true);

            for (TermId tid1 : inclAncestorTermIds) {
                termIdToGeneIds.putIfAbsent(tid1, new HashSet<>());
                termIdToGeneIds.get(tid1).add(tid0);
            }
        }

        // Compute information content of HPO terms, given the term-to-gene annotation.
        logger.info("Performing IC precomputation...");
        final Map<TermId, Double> icMap =
                new InformationContentComputation(hpo)
                        .computeInformationContent(termIdToGeneIds);
        logger.info("DONE: Performing IC precomputation");
        return icMap;
    }

    /**
     * Resnik similarity precomputation
     * @return ResnikSimilarity object (asymmetric score)
     */
    private ResnikSimilarity createResnik() {
        logger.info("Performing Resnik precomputation...");
        final PrecomputingPairwiseResnikSimilarity pairwiseResnikSimilarity =
                new PrecomputingPairwiseResnikSimilarity(hpo, icMap, NUM_THREADS);
        logger.info("DONE: Performing Resnik precomputation");

        final ResnikSimilarity resnikSimilarity =
                new ResnikSimilarity(pairwiseResnikSimilarity, false);
        logger.info(String.format("name: %s  params %s",
                resnikSimilarity.getName(),
                resnikSimilarity.getParameters()));
        return resnikSimilarity;
    }

    /**
     * Finds set of phenotype ids that are associated with (minDiseases or more) of the input diseases
     * @param minDiseases  integer to filter phenotypes: only keep those that occur in at least minDiseases
     * @param diseaseIds   set of disease TermIds
     * @return set of TermIds for phenotypes associated with at least minDiseases of the diseases
     */
    private Set<TermId> diseasesToPhenotypes( int minDiseases, Collection<TermId> diseaseIds ) {
        Set<TermId> phenotypes = new HashSet<>();
        // minDiseases = 1 (or any integer < 2) means no filter on phenotypes
        if (minDiseases < 2) {
            for (TermId tid : diseaseIds) {
                HpoDisease disease = diseaseMap.get(tid);
//                if (disease == null) {
//                    logger.warn("Disease listed in mim2gene_medgen but not in phenotype.hpoa: " + tid);
//                } else {
                if (disease != null) {
                    phenotypes.addAll(disease.getPhenotypicAbnormalityTermIdList());
                }
            }
        } else { // filter phenotypes by minDiseases
            Map<TermId, Integer> counts = new HashMap<>();
            for (TermId tid0 : diseaseIds) {
                HpoDisease disease = diseaseMap.get(tid0);
//                if (disease == null) {
//                    logger.warn("Disease listed in mim2gene_medgen but not in phenotype.hpoa: " + tid0);
//                } else {
                if (disease != null) {
                    for (TermId tid1 : disease.getPhenotypicAbnormalityTermIdList()) {
                        counts.putIfAbsent(tid1, 0);
                        counts.put(tid1, counts.get(tid1) + 1);
                    }
                }
            }
            counts.forEach((pheno, count) -> { if (count >= minDiseases) phenotypes.add(pheno); });
        }
        return phenotypes;
    }

    /**
     * Finds all diseases associated with the input collection of gene ids. Ignores any genes
     * for which no diseases are known.
     * @param geneIds Collection of TermIds for genes
     * @return set of TermIds for all diseases associated with the genes
     */
    private Set<TermId> genesToDiseases( Collection<TermId> geneIds ) {
        Set<TermId> retSet = new HashSet<>();
        Collection<TermId> diseasesForGene;

        for (TermId tid : geneIds) {
            // skip over any gene for which there are no diseases in genesToDiseasesMap
            diseasesForGene = genesToDiseasesMap.get(tid);
            if (diseasesForGene != null) {
                retSet.addAll(diseasesForGene);
            }
        }
        return retSet;
    }

    int countTotalHpoTerms(Map<TermId, Set<TermId>> mysetofgenes) {
        Set<TermId> allterms=new HashSet<>();
        for (Set<TermId> hposet : mysetofgenes.values()) {
            allterms.addAll(hposet);
        }
        return allterms.size();
    }

    /**
     * Run randomization.
     * @param minDiseases gpi pathway phenotype must appear in at least this number of diseases to be considered
     * minDiseases = 0 means no filter on phenotypes
     * @param threshold   minimum info content of phenotype to be counted in simfunAboveThreshold calculation
     * threshold = -1 means skip simfunAboveThreshold
     */
    void run( int minDiseases, double threshold ) {
        ResnikSimilarity resnikSimilarity = createResnik();

        // create a set of diseases, a set of the associated phenotypes for genes in the GPI pathway
        // if minDiseases > 1, filter out the phenotypes that occur in fewer diseases
        Set<TermId> gpiPathwayDiseases = genesToDiseases(gpiPathwayGenes);
        Set<TermId> gpiPathwayPhenotypes = diseasesToPhenotypes(minDiseases, gpiPathwayDiseases);
        double sim1;
        double sim2;
        double sim3 = 0.0;

        Map<TermId, Set<TermId>> gpiAnchoredPhenotypes = targetPhenotypes(gpiAnchoredGenes);

        logger.info(String.format("minDiseases = %d", minDiseases));
        logger.info(String.format("threshold = %.2f", threshold));

        sim1 = simfunAllPhenotypes(gpiPathwayPhenotypes, gpiAnchoredPhenotypes.values(), resnikSimilarity);
        logger.info(String.format(
                "Similarity function allPhenotypes applied to GPI pathway genes, GPI anchored genes: %.2f",
                sim1));

        sim2 = simfunBestMatch(gpiPathwayPhenotypes, gpiAnchoredPhenotypes.values(), resnikSimilarity);
        logger.info(String.format(
                "Similarity function bestMatch applied to GPI pathway genes, GPI anchored genes: %.2f",
                sim2));

        if (threshold > -1.0) {
            sim3 = simfunAboveThreshold(gpiPathwayPhenotypes, gpiAnchoredPhenotypes.values(), threshold,
                    resnikSimilarity);
            logger.info(String.format(
                    "Similarity function aboveThreshold applied to GPI pathway genes, GPI anchored genes: %.2f",
                    sim3));
        }

//        int NtermsRealSet=countTotalHpoTerms(gpiAnchoredPhenotypes);

        // initialize counters for the three similarity functions
        int m1 = 0;
        int m2 = 0;
        int m3 = 0;
        Random rand = new Random();
        Set<TermId> sample;
        Map<TermId, Set<TermId>> samplePhenotypes;
        int sampleSize = gpiAnchoredGenes.size();

        /* debugging
        System.err.println("Size of GPI pathway genes: " + gpiPathwayGenes.size());
        System.err.println("Sample size for random sets of genes: " + sampleSize);
        */

        for (int i = 0; i < NUM_ITER; i++) {
            sample = randomSample(rand, sampleSize, allDiseaseGenes);
            samplePhenotypes = targetPhenotypes(sample);
            double simAllPhen=simfunAllPhenotypes(gpiPathwayPhenotypes, samplePhenotypes.values(),
                    resnikSimilarity);
            if (simAllPhen>= sim1)
                m1++;
            if (simfunBestMatch(gpiPathwayPhenotypes, samplePhenotypes.values(), resnikSimilarity) >= sim2)
                m2++;
            if (threshold > -1.0 && simfunAboveThreshold(gpiPathwayPhenotypes, samplePhenotypes.values(),
                    threshold, resnikSimilarity) >= sim3)
                m3++;
//            System.err.println("All Phen: "+simAllPhen);
//            int M = countTotalHpoTerms(samplePhenotypes);
//            System.err.println("Number of HPO terms for sample "+M +" (Real="+NtermsRealSet);
//            System.err.println("Number of HPO terms for sample "+M +" (Real="+NtermsRealSet);
        }

        logger.info(String.format("Number of iterations: %d", NUM_ITER));
        logger.info(String.format("p value for allPhenotypes similarity: %.4f", (double) m1/NUM_ITER));
        logger.info(String.format("p value for bestMatch similarity: %.4f", (double) m2/NUM_ITER));
        if (threshold > -1.0)
            logger.info(String.format("p value for aboveThreshold similarity: %.4f", (double) m3/NUM_ITER));

        /* example of computing score between the sets of HPO terms that annotate two
        // diseases (get the diseases at random)
        logger.info("About to calculate phenotype similarity from two random diseases from a map of size " + diseaseMap.size());
        List<HpoDisease> valuesList = new ArrayList<>(diseaseMap.values());
        int randomIndex1 = new Random().nextInt(valuesList.size());
        HpoDisease randomDisease1 = valuesList.get(randomIndex1);
        int randomIndex2 = new Random().nextInt(valuesList.size());
        HpoDisease randomDisease2 = valuesList.get(randomIndex2);

        List<TermId> phenoAbnormalities1 = randomDisease1.getPhenotypicAbnormalityTermIdList();
        List<TermId> phenoAbnormalities2 = randomDisease2.getPhenotypicAbnormalityTermIdList();

        double similarity = resnikSimilarity.computeScore(phenoAbnormalities1,phenoAbnormalities2);

        logger.info(String.format("Similarity score between the query %s and the target disease %s was %.4f",
                randomDisease1.getName(),randomDisease2.getName(),similarity));

        //public final double computeScore(Collection<TermId> query, Collection<TermId> target) {

        // Temporary storage of term count to score distributions.
        final Map<Integer, ScoreDistribution> scoreDists = new HashMap<>();
        */
    }

    /**
     * For each gene in input set, finds phenotypes of diseases related to that gene
     * @param targetGenes genes for which phenotypes are sought
     * @return map from gene TermId to set of phenotype TermIds for diseases related to that gene
     */
    private Map<TermId, Set<TermId>> targetPhenotypes(Set<TermId> targetGenes) {
        Map<TermId, Set<TermId>> phenotypesByGene = new HashMap<>();
//        targetGenes.forEach(tid -> phenotypesByGene.put(tid, diseasesToPhenotypes(1,
//                genesToDiseasesMap.getOrDefault(tid, new HashSet<>()))));
        for (TermId tid : targetGenes) {
            Set<TermId> diseases = genesToDiseasesMap.getOrDefault(tid, new HashSet<>());
//            if (diseases.size() > 0) System.out.println(tid.getIdWithPrefix() + " " + diseases.size());
            phenotypesByGene.put(tid, diseasesToPhenotypes(1, diseases));
        }
        return phenotypesByGene;
    }
}
