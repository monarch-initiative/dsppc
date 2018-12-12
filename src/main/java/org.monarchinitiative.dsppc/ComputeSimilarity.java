package org.monarchinitiative.dsppc;

import com.google.common.collect.Sets;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.monarchinitiative.phenol.formats.hpo.HpoDisease;
import org.monarchinitiative.phenol.formats.hpo.HpoOntology;
import org.monarchinitiative.phenol.ontology.algo.InformationContentComputation;
import org.monarchinitiative.phenol.ontology.data.TermId;
import org.monarchinitiative.phenol.ontology.data.TermIds;
import org.monarchinitiative.phenol.ontology.scoredist.ScoreDistribution;
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
    private final Map<TermId, Collection<TermId>> genesToDiseasesMap;
    private final Set<TermId> gpiAnchoredGenes;
    private final Set<TermId> gpiPathwayGenes;
    private final int numThreads = 4;

    private static final Logger logger = LogManager.getLogger();

    ComputeSimilarity(HpoOntology hpo, Map<TermId, HpoDisease> diseaseMap,
                      Map<TermId, Collection<TermId>> geneToDiseasesMap,
                      Set<TermId> gpiPathway, Set<TermId> gpiAnchored) {
        this.hpo = hpo;
        this.diseaseMap = diseaseMap;
        this.genesToDiseasesMap = geneToDiseasesMap;
        gpiPathwayGenes = gpiPathway;
        gpiAnchoredGenes = gpiAnchored;
    }

    /**
     * Compute information content of each phenotype term associated with one or more diseases.
     * @return map from phenotype TermId to its information content.
     */
    private Map<TermId, Double> computeICmap() {
        // Compute mapping from OMIM ID to phenotype TermIds and from phenotype to OMIM termId.
        logger.info("Mapping from OMIM ID to HPO phenotype terms and the reverse...");
        final Map<TermId, Collection<TermId>> diseaseIdToTermIds = new HashMap<>();
        final Map<TermId, Collection<TermId>> termIdToDiseaseIds = new HashMap<>();

        for (TermId diseaseId : diseaseMap.keySet()) {
            HpoDisease disease = diseaseMap.get(diseaseId);
            List<TermId> hpoTerms = disease.getPhenotypicAbnormalityTermIdList();
            diseaseIdToTermIds.putIfAbsent(diseaseId, new HashSet<>());
            // add term anscestors
            final Set<TermId> inclAncestorTermIds =
                    TermIds.augmentWithAncestors(hpo, Sets.newHashSet(hpoTerms), true);

            for (TermId tid : inclAncestorTermIds) {
                termIdToDiseaseIds.putIfAbsent(tid, new HashSet<>());
                termIdToDiseaseIds.get(tid).add(diseaseId);
                diseaseIdToTermIds.get(diseaseId).add(tid);
            }
        }

        // Compute information content of HPO terms, given the term-to-gene annotation.
        logger.info("Performing IC precomputation...");
        final Map<TermId, Double> icMap =
                new InformationContentComputation(hpo)
                        .computeInformationContent(termIdToDiseaseIds);
        logger.info("DONE: Performing IC precomputation");
        return icMap;
    }

    /**
     * Resnik similarity precomputation
     * @param icMap map from phenotype TermId to information content
     * @return ResnikSimilarity object (asymmetric score)
     */
    private ResnikSimilarity createResnik(Map<TermId, Double> icMap) {
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

    private Set<TermId> diseasesToPhenotypes( int minDiseases, Collection<TermId> diseases ) {
        Set<TermId> phenotypes = new HashSet<>();
        // minDiseases = 0 means no filter on phenotypes
        if (minDiseases < 2) {
            diseases.forEach(d -> phenotypes.addAll(diseaseMap.get(d).getPhenotypicAbnormalityTermIdList()));
        } else {
            Map<TermId, Integer> counts = new HashMap<>();
            for (TermId d : diseases) {
                for (TermId p: diseaseMap.get(d).getPhenotypicAbnormalityTermIdList()) {
                    counts.putIfAbsent(p, 0);
                    counts.put(p, counts.get(p) + 1);
                }
            }
            counts.forEach((pheno, count) -> { if (count >= minDiseases) phenotypes.add(pheno); });
        }
        return phenotypes;
    }

    /**
     * Run application.
     * @param minDiseases gpi pathway phenotype must appear in at least this number of diseases to be considered
     * minDiseases = 0 means no filter on phenotypes
     */
    void run( int minDiseases ) {

        Map<TermId, Double> icMap = computeICmap();
        final ResnikSimilarity resnikSimilarity = createResnik(icMap);
//    for (TermId t:icMap.keySet()) {
//      System.out.println("IC-> "+t.getIdWithPrefix() + ": " + icMap.get(t));
//    }
        // create a set of diseases, a set of the associated phenotypes for genes in the GPI pathway
        // if minDiseases > 0, filter out the phenotypes that occur in fewer diseases
        Set<TermId> gpiPathwayDiseases = new HashSet<>();
        gpiPathwayGenes.forEach(gene -> gpiPathwayDiseases.addAll(genesToDiseasesMap.get(gene)));
        Set<TermId> gpiPathwayPhenotypes = diseasesToPhenotypes(minDiseases, gpiPathwayDiseases);

        Map<TermId, Double> anchoredScores = new HashMap<>();
        gpiAnchoredGenes.forEach(gene ->
            anchoredScores.put(gene, resnikSimilarity.computeScore(gpiPathwayPhenotypes,
                    diseasesToPhenotypes(0, genesToDiseasesMap.get(gene)))));

        // example of computing score between the sets of HPO terms that annotate two
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
    }
}
