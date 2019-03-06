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
    private Map<TermId, Double> icMap;
//    private int minDiseases;
    private PairwiseResnikSimilarity pwResSim;

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
    public void setUp() throws Exception {
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
        Field rsf = csc.getDeclaredField("resnikSimilarity");
        rsf.setAccessible(true);
        ResnikSimilarity resSim = (ResnikSimilarity) rsf.get(comSim);
        Class<? extends ResnikSimilarity> rsc = resSim.getClass();
//        Class<? extends ResnikSimilarity> rsc = resSim.getClass();
//        Class cl = Class.forName("AbstractCommonAncestorSimilarity");
        Class scl = rsc.getSuperclass();
        Field pwsf = scl.getDeclaredField("pairwiseSimilarity");
//        Field pwsf = rsc.getDeclaredField("pairwiseSimilarity");
        pwsf.setAccessible(true);
//        pwResSim = (PrecomputingPairwiseResnikSimilarity) pwsf.get(resSim);
//        icMap = pwResSim.getTermToIc();
    }

    @Test
    public void testIC() {
        System.out.println(String.format(
                "Information content%n%6.4f %s%n%6.4f %s%n%6.4f %s%n%6.4f %s%n%6.4f %s%n",
                icMap.get(pa), "phentypic abnormality",
                icMap.get(al), "abnormality of the liver ", icMap.get(pld), "polycystic liver disease",
                icMap.get(sc), "splenic cyst", icMap.get(ap5t), "abnormality of the phalanges of the 5th toe"));
        // more general concept should have lower information content
        assertTrue(icMap.get(al) < icMap.get(pld));
    }

    @Test
    public void testRS() {
        double score0, score1;

        // two liver phenotypes should be more similar than a liver phenotype and a bone abnormality
        assertTrue("ap5t is more similar to pld than al is to pld",
                pwResSim.computeScore(al, pld) > pwResSim.computeScore(ap5t, pld));
        // two abdominal diseases involving cysts should be more similar than generic abnormality of liver is
        // to splenic cyst, but it appears to be equal
        score0 = pwResSim.computeScore(pld, sc);
        score1 = pwResSim.computeScore(al, sc);
        assertTrue("al is more similar to sc than pld is to sc", score0 >= score1);
        assertEquals("similarity(pld, sc) is not equal to similarity(al, sc)", score0, score1, 0.00001);
    }
}