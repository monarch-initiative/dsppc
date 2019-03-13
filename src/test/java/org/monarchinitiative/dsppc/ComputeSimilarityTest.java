package org.monarchinitiative.dsppc;

import org.junit.BeforeClass;
import org.junit.Test;
import org.monarchinitiative.phenol.formats.hpo.HpoDisease;

import org.monarchinitiative.phenol.io.OntologyLoader;
import org.monarchinitiative.phenol.ontology.data.Ontology;
import org.monarchinitiative.phenol.ontology.data.TermId;
//import org.monarchinitiative.phenol.ontology.similarity.AbstractCommonAncestorSimilarity;
import org.monarchinitiative.phenol.ontology.similarity.PairwiseResnikSimilarity;
import org.monarchinitiative.phenol.ontology.similarity.PrecomputingPairwiseResnikSimilarity;
import org.monarchinitiative.phenol.ontology.similarity.ResnikSimilarity;

import java.lang.reflect.*;
import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import static org.junit.Assert.*;
import static org.monarchinitiative.dsppc.Dsppc.*;

/*
 * created 14 Jan 2019
 */
public class ComputeSimilarityTest {
    private static Map<TermId, Double> icMap;
//    private int minDiseases;
    private static PairwiseResnikSimilarity pwResSim;

    // abnormality of the liver
    private static TermId al = TermId.of( "HP:0001392");
    // abnormality of the phalanges of the 5th toe
    private static TermId ap5t = TermId.of( "HP:0010342");
    // phenotypic abnormality
    private static TermId pa = TermId.of("HP:0000118");
    // polycystic liver disease
    private static TermId pld = TermId.of( "HP:0006557");
    // splenic cyst
    private static TermId sc = TermId.of( "HP:0030423");

    private double threshold;

    @BeforeClass
    public static void setUp() throws Exception {
        String dataDir = "src/main/resources/";
        Set<TermId> gpiAnchoredGenes = new TreeSet<>();
        Set<TermId> gpiPathwayGenes = new TreeSet<>();

        List<TermId> allGenes = parseAllGenes(dataDir + ALL_GENES_FILENAME);
        Ontology hpo = OntologyLoader.loadOntology(new File(dataDir + HPO_FILENAME));
        Map<TermId, HpoDisease> diseaseMap = parseHPOA(dataDir + HPOA_FILENAME, hpo);
        Map<TermId, Set<TermId>> geneToDiseasesMap = parseMedgen(dataDir + MIM2GENE_MEDGEN_FILENAME);
        parseGeneSets(dataDir + GENE_SETS_FILENAME, gpiPathwayGenes, gpiAnchoredGenes);

        ComputeSimilarity comSim = new ComputeSimilarity(hpo, diseaseMap, geneToDiseasesMap, allGenes,
                gpiPathwayGenes, gpiAnchoredGenes);
        Class csc = comSim.getClass();
        Field icm = csc.getDeclaredField("icMap");
        icm.setAccessible(true);
        icMap = ( Map<TermId, Double> ) icm.get(comSim);
        pwResSim = new PairwiseResnikSimilarity(hpo, icMap);
    }

    @Test
    public void testIC() {
        System.out.println(String.format(
                "Information content%n%6.4f %s%n%6.4f %s%n%6.4f %s%n%6.4f %s%n%6.4f %s%n",
                icMap.get(pa), "phentypic abnormality",
                icMap.get(al), "abnormality of the liver ", icMap.get(pld), "polycystic liver disease",
                icMap.get(sc), "splenic cyst", icMap.get(ap5t), "abnormality of the phalanges of the 5th toe"));
        // more general concept should have lower information content
        assertTrue("IC of al is >= IC of pld", icMap.get(al) < icMap.get(pld));
    }

    @Test
    public void testRS() {
        double score_al_pld, score_ap5t_pld, score_pld_sc, score_al_sc, score_ap5t_ap5t, score_al_al;

        score_al_pld = pwResSim.computeScore(al, pld);
        score_ap5t_pld = pwResSim.computeScore(ap5t, pld);
        score_pld_sc = pwResSim.computeScore(pld, sc);
        score_al_sc = pwResSim.computeScore(al, sc);
        score_ap5t_ap5t = pwResSim.computeScore(ap5t, ap5t);
        score_al_al = pwResSim.computeScore(al, al);
        // two liver phenotypes should be more similar than a liver phenotype and a bone abnormality
        assertTrue("ap5t is more similar to pld than al is to pld",
                 score_al_pld > score_ap5t_pld);
        // two abdominal diseases involving cysts should be more similar than generic abnormality of liver is
        // to splenic cyst, but it appears to be equal
        assertTrue("al is more similar to sc than pld is to sc", score_pld_sc >= score_al_sc);
        assertEquals("similarity(pld, sc) is not equal to similarity(al, sc)", score_pld_sc, score_al_sc,
                0.00001);

        System.out.println(String.format(
                "Pairwise similarity%n%6.4f %s%n%6.4f %s%n%6.4f %s%n%6.4f %s%n%6.4f %s%n%6.4f %s%n",
                score_al_pld, "sim(abnormality of the liver, polycystic liver disease)",
                score_al_sc, "sim(abnormality of the liver, splenic cyst)",
                score_pld_sc, "sim(polycystic liver disease, splenic cyst)",
                score_ap5t_pld, "sim(abnormality of the phalanges of the 5th toe, polycystic liver disease)",
                score_ap5t_ap5t, "sim(abnormality of the phalanges of the 5th toe, " +
                        "abnormality of the phalanges of the 5th toe)",
                score_al_al, "sim(abnormality of the liver, abnormality of the liver)"));
    }
}