package org.monarchinitiative.dsppc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.util.*;
import java.util.stream.Collectors;


import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import org.monarchinitiative.phenol.base.PhenolException;
import org.monarchinitiative.phenol.formats.hpo.HpoGeneAnnotation;
import org.monarchinitiative.phenol.formats.hpo.HpoOntology;
import org.monarchinitiative.phenol.io.base.TermAnnotationParserException;
import org.monarchinitiative.phenol.io.obo.hpo.HpOboParser;
import org.monarchinitiative.phenol.io.obo.hpo.HpoGeneAnnotationParser;
import org.monarchinitiative.phenol.ontology.algo.InformationContentComputation;
import org.monarchinitiative.phenol.ontology.data.*;
import org.monarchinitiative.phenol.ontology.scoredist.ScoreDistribution;
import org.monarchinitiative.phenol.ontology.scoredist.ScoreDistributions;
import org.monarchinitiative.phenol.ontology.scoredist.ScoreSamplingOptions;
import org.monarchinitiative.phenol.ontology.scoredist.SimilarityScoreSampling;
import org.monarchinitiative.phenol.ontology.similarity.PrecomputingPairwiseResnikSimilarity;
import org.monarchinitiative.phenol.ontology.similarity.ResnikSimilarity;
import org.monarchinitiative.phenol.io.obo.hpo.*;
import org.monarchinitiative.phenol.formats.hpo.*;


// adapted from ComputeSimilarityDemo.java in phenol library

/**
 * App for computing similarity scores between gene (by entrez ID) and HPO term list.
 *
 * @author <a href="mailto:manuel.holtgrewe@bihealth.de">Manuel Holtgrewe</a>
 * @author <a href="mailto:peter.robinson@jax,org">Peter Robinson</a>
 */
class ComputeSimilarityDemo {


    /**
     * Number of threads to use.
     */
    private final int numThreads = 4;

    /**
     * Command line arguments.
     */
    private final HpoOntology hpo;

    /**
     * Path to hp.obo file to read.
     */
    private String pathHpObo;

    /**
     * Path to {@code phenotype.hpoa}
     */
    private String pathPhenotypeHpoa;

    /**
     * Path to TSV file to read
     */
    private String pathTsvFile = "resources/omim-example.txt";

    /**
     * Construct with argument list.
     *
     * @param hpo HPO ontology.
     */
    ComputeSimilarityDemo(HpoOntology hpo) {
        this.hpo = hpo;
    }

    /**
     * Run application.
     */
    void run() {
        final Map<TermId, HpoDisease> diseaseMap;
        try {
            HpoDiseaseAnnotationParser parser = new HpoDiseaseAnnotationParser(this.pathPhenotypeHpoa, hpo);
            diseaseMap = parser.parse();
        } catch (PhenolException e) {
            e.printStackTrace();
            System.exit(1);
            return; // javac complains otherwise
        }
        // Compute list of annoations and mapping from OMIM ID to term IDs.
        System.out.println("Loading HPO to disease annotation file...");
        final Map<TermId, Collection<TermId>> diseaseIdToTermIds = new HashMap<>();
        final Map<TermId, Collection<TermId>> termIdToDiseaseIds = new HashMap<>();

        for (TermId diseaseId : diseaseMap.keySet()) {
            HpoDisease disease = diseaseMap.get(diseaseId);
            List<TermId> hpoTerms = disease.getPhenotypicAbnormalityTermIdList();
            diseaseIdToTermIds.putIfAbsent(diseaseId, new HashSet<>());
            // add term anscestors
            final Set<TermId> inclAncestorTermIds = TermIds.augmentWithAncestors(hpo, Sets.newHashSet(hpoTerms), true);

            for (TermId tid : inclAncestorTermIds) {
                termIdToDiseaseIds.putIfAbsent(tid, new HashSet<>());
                termIdToDiseaseIds.get(tid).add(diseaseId);
                diseaseIdToTermIds.get(diseaseId).add(tid);
            }
        }


        // Compute information content of HPO terms, given the term-to-gene annotation.
        System.out.println("Performing IC precomputation...");
        final Map<TermId, Double> icMap =
                new InformationContentComputation(hpo)
                        .computeInformationContent(termIdToDiseaseIds);
        System.out.println("DONE: Performing IC precomputation");
//    int i=0;
//    for (TermId t:icMap.keySet()) {
//      System.out.println("IC-> "+t.getIdWithPrefix() + ": " + icMap.get(t));
//    }


        // Initialize Resnik similarity precomputation
        System.out.println("Performing Resnik precomputation...");
        final PrecomputingPairwiseResnikSimilarity pairwiseResnikSimilarity =
                new PrecomputingPairwiseResnikSimilarity(hpo, icMap, numThreads);
        System.out.println("DONE: Performing Resnik precomputation");
        final ResnikSimilarity resnikSimilarity =
                new ResnikSimilarity(pairwiseResnikSimilarity, false);
        System.out.println(String.format("name: %s  params %s",
                resnikSimilarity.getName(),
                resnikSimilarity.getParameters()));
        // example of computing score between the sets of HPO terms that annotate two
        // diseases (get the diseases at random)
        System.out.println("About to calculate phenotype similarity from two random diseases from a map of size " + diseaseMap.size());
        List<HpoDisease> valuesList = new ArrayList<>(diseaseMap.values());
        int randomIndex1 = new Random().nextInt(valuesList.size());
        HpoDisease randomDisease1 = valuesList.get(randomIndex1);
        int randomIndex2 = new Random().nextInt(valuesList.size());
        HpoDisease randomDisease2 = valuesList.get(randomIndex2);

        List<TermId> phenoAbnormalities1 = randomDisease1.getPhenotypicAbnormalityTermIdList();
        List<TermId> phenoAbnormalities2 = randomDisease2.getPhenotypicAbnormalityTermIdList();


        double similarity = resnikSimilarity.computeScore(phenoAbnormalities1,phenoAbnormalities2);

        System.out.println(String.format("Similarity score between the query %s and the target disease %s was %.4f",
                randomDisease1.getName(),randomDisease2.getName(),similarity));

        //public final double computeScore(Collection<TermId> query, Collection<TermId> target) {


        // Temporary storage of term count to score distributions.
        final Map<Integer, ScoreDistribution> scoreDists = new HashMap<>();

        // Read file line-by line and process.


    }
}
