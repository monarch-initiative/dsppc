package org.monarchinitiative.dsppc;

/*
 * created 04 Feb 2019
 */

import com.google.common.collect.ImmutableSet;
import org.monarchinitiative.phenol.ontology.data.TermId;
import org.monarchinitiative.phenol.ontology.similarity.ResnikSimilarity;

import java.util.*;

class SimFuns {

    /**
     * Selects a random sample (with replacement) of the specified size from termIds
     * @param rand         random number generator
     * @param desiredSize  how big the random sample should be
     * @param termIds      list from which to select the sample
     * @return Set of TermIds selected at random from input termIds
     */
    static Set<TermId> randomSample(Random rand, int desiredSize, List<TermId> termIds) {
        Set<TermId> sample = new HashSet<>();
        int upperLimit = termIds.size();
        while (sample.size() < desiredSize) {
            sample.add(termIds.get(rand.nextInt(upperLimit)));
        }
        return sample;
    }

    /**
     * This function finds the sum of those phenotype of EACH target gene to the set of phenotypes
     * in the source diseases (e.g. GPI) and only counts matches that are above threshold.
     * @param source    set of phenotype ids
     * @param targets   collection of sets of phenotype ids
     * @param threshold minimum IC matching score for a phenotype match between an HPO of targets and the terms in source
     * @param resnikSimilarity ResnikSimilarity object to compute similarity of two sets of TermIds
     * @return          similarity score
     */
    static double simfunAboveThreshold(Set<TermId> source, Collection<Set<TermId>> targets,
                                        double threshold, ResnikSimilarity resnikSimilarity) {
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
     * @param resnikSimilarity ResnikSimilarity object to compute similarity of two sets of TermIds
     * @return        similarity score
     */
    static double simfunAllPhenotypes(
            Set<TermId> source, Collection<Set<TermId>> targets, ResnikSimilarity resnikSimilarity) {
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
     * @param resnikSimilarity ResnikSimilarity object to compute similarity of two sets of TermIds
     * @return        similarity score
     */
    static double simfunBestMatch(
            Set<TermId> source, Collection<Set<TermId>> targets, ResnikSimilarity resnikSimilarity ) {
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
}
