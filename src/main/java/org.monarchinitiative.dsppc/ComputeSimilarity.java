package org.monarchinitiative.dsppc;

import com.google.common.collect.ImmutableSet;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.monarchinitiative.phenol.formats.hpo.HpoDisease;
import org.monarchinitiative.phenol.formats.hpo.HpoOntology;
import org.monarchinitiative.phenol.ontology.algo.InformationContentComputation;
import org.monarchinitiative.phenol.ontology.data.TermId;
import org.monarchinitiative.phenol.ontology.data.TermIds;
import org.monarchinitiative.phenol.ontology.similarity.PrecomputingPairwiseResnikSimilarity;
import org.monarchinitiative.phenol.ontology.similarity.ResnikSimilarity;

import java.util.*;

// adapted from ComputeSimilarity.java in phenol library

/**
 * App for computing similarity scores between gene (by entrez ID) and HPO term list.
 *
 * @author <a href="mailto:manuel.holtgrewe@bihealth.de">Manuel Holtgrewe</a>
 * @author <a href="mailto:peter.robinson@jax,org">Peter Robinson</a>
 */
class ComputeSimilarity {

    private final Map<TermId, HpoDisease> diseaseMap;
    private final HpoOntology hpo;
    private final Map<TermId, Set<TermId>> genesToDiseasesMap;
    private final Set<TermId> gpiAnchoredGenes;
    private final Set<TermId> gpiPathwayGenes;
    private final ResnikSimilarity resnikSimilarity;

    private static final Logger logger = LogManager.getLogger();
    private static final int numThreads = 4;

    ComputeSimilarity(HpoOntology hpo, Map<TermId, HpoDisease> diseaseMap,
                      Map<TermId, Set<TermId>> geneToDiseasesMap,
                      Set<TermId> gpiPathway, Set<TermId> gpiAnchored) {
        this.hpo = hpo;
        this.diseaseMap = diseaseMap;
        this.genesToDiseasesMap = geneToDiseasesMap;
        gpiPathwayGenes = gpiPathway;
        gpiAnchoredGenes = gpiAnchored;
        resnikSimilarity = createResnik();
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
        Map<TermId, Double> icMap = computeICmap();
        logger.info("Performing Resnik precomputation...");
        final PrecomputingPairwiseResnikSimilarity pairwiseResnikSimilarity =
                new PrecomputingPairwiseResnikSimilarity(hpo, icMap, numThreads);
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
                if (disease == null) {
                    logger.warn("Disease listed in mim2gene_medgen but not in phenotype.hpoa: " + tid);
                } else {
                    phenotypes.addAll(disease.getPhenotypicAbnormalityTermIdList());
                }
            }
        } else { // filter phenotypes by minDiseases
            Map<TermId, Integer> counts = new HashMap<>();
            for (TermId tid0 : diseaseIds) {
                HpoDisease disease = diseaseMap.get(tid0);
                if (disease == null) {
                    logger.warn("Disease listed in mim2gene_medgen but not in phenotype.hpoa: " + tid0);
                } else {
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

    /**
     * This function finds the sum of those phenotype of EACH target gene to the set of phenotypes
     * in the source diseases (e.g. GPI) and only counts matches that are above threshold.
     * @param source    set of phenotype ids
     * @param targets   collection of sets of phenotype ids
     * @param threshold minimum IC matching score for a phenotype match between an HPO of targets and the terms in source
     * @return          similarity score
     */
    private double simfunAboveThreshold( Set<TermId> source, Collection<Set<TermId>> targets, double threshold) {
        double score = 0.0;
        for (Set<TermId> targetPhenoSet : targets) {
            for (TermId feature : targetPhenoSet) {
                double tmp = resnikSimilarity.computeScore(source, ImmutableSet.of(feature));
                if (tmp>threshold)  score += tmp;
            }
        }
        return score;
    }

    /**
     * This function compares ALL phenotypes of EACH target gene to the set of phenotypes
     * in the source diseases (e.g. GPI).
     * @param source  set of phenotype ids
     * @param targets collection of sets of phenotype ids
     * @return        similarity score
     */
    private double simfunAllPhenotypes(Set<TermId> source, Collection<Set<TermId>> targets) {
        double score = 0.0;
        for (Set<TermId> targetPheno : targets) {
            score += resnikSimilarity.computeScore(source, targetPheno);
        }
        return score;
    }

    /**
     * This function finds the sum of the best matching phenotype of EACH target gene to
     * the set of phenotypes in the source diseases (e.g. GPI).
     * @param source  set of phenotype ids
     * @param targets collection of sets of phenotype ids
     * @return        similarity score
     */
    private double simfunBestMatch( Set<TermId> source, Collection<Set<TermId>> targets ) {
        double score = 0.0;
        for (Set<TermId> targetPhenoSet : targets) {
            double max=0d;
            for (TermId feature : targetPhenoSet) {
                double tmp = resnikSimilarity.computeScore(source, ImmutableSet.of(feature));
                if (tmp>max) max=tmp;
            }
            score += max;
        }
        return score;
    }

    /**
     * Run randomization.
     * @param minDiseases gpi pathway phenotype must appear in at least this number of diseases to be considered
     * minDiseases = 0 means no filter on phenotypes
     * @param threshold   minimum info content of phenotype to be counted in simfunAboveThreshold calculation
     * threshold = -1 means skip simfunAboveThreshold
     */
    void run( int minDiseases, double threshold ) {
        // create a set of diseases, a set of the associated phenotypes for genes in the GPI pathway
        // if minDiseases > 1, filter out the phenotypes that occur in fewer diseases
        Set<TermId> gpiPathwayDiseases = genesToDiseases(gpiPathwayGenes);
        Set<TermId> gpiPathwayPhenotypes = diseasesToPhenotypes(minDiseases, gpiPathwayDiseases);

        Map<TermId, Set<TermId>> gpiAnchoredPhenotypes = new HashMap<>();
//        gpiAnchoredGenes.forEach(tid -> gpiAnchoredPhenotypes.put(geneId, diseasesToPhenotypes(1,
//                genesToDiseasesMap.getOrDefault(tid, new HashSet<>()))));

        for (TermId tid : gpiAnchoredGenes) {
            Set<TermId> diseases = genesToDiseasesMap.getOrDefault(tid, new HashSet<>());
//            if (diseases.size() > 0) System.out.println(geneId.getIdWithPrefix() + " " + diseases.size());
            gpiAnchoredPhenotypes.put(tid, diseasesToPhenotypes(1, diseases));
        }

        logger.info(String.format(
                "Similarity function allPhenotypes applied to GPI pathway genes, GPI anchored genes: %.2f",
                simfunAllPhenotypes(gpiPathwayPhenotypes, gpiAnchoredPhenotypes.values())));

        logger.info(String.format(
                "Similarity function bestMatch applied to GPI pathway genes, GPI anchored genes: %.2f",
                simfunBestMatch(gpiPathwayPhenotypes, gpiAnchoredPhenotypes.values())));

        if (threshold > -1.0) {
            logger.info(String.format(
                    "Similarity function aboveThreshold (%.2f) applied to GPI pathway genes, GPI anchored genes: %.2f",
                    threshold, simfunAboveThreshold(gpiPathwayPhenotypes, gpiAnchoredPhenotypes.values(), threshold)));
        }
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
}
