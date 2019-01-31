package org.monarchinitiative.dsppc;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.monarchinitiative.phenol.formats.hpo.HpoDisease;

import org.monarchinitiative.phenol.io.OntologyLoader;
import org.monarchinitiative.phenol.ontology.data.Ontology;
import org.monarchinitiative.phenol.ontology.data.TermId;
import org.monarchinitiative.phenol.ontology.similarity.PairwiseResnikSimilarity;

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
    private ComputeSimilarity comSim;
    private String dataDir = "src/main/resources/";
    private List<TermId> allGenes;
    private Map<TermId, HpoDisease> diseaseMap;
    private Map<TermId, Set<TermId>> geneToDiseasesMap;
    private Set<TermId> gpiAnchoredGenes = new TreeSet<>();
    private Set<TermId> gpiPathwayGenes = new TreeSet<>();
    private Ontology hpo;
    private Map<TermId, Double> icMap;
    private int minDiseases;

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

    @Before
    public void setUp() throws Exception {
        allGenes = parseAllGenes(dataDir + ALL_GENES_FILENAME);
        hpo = OntologyLoader.loadOntology(new File(dataDir + HPO_FILENAME));
        diseaseMap = parseHPOA(dataDir + HPOA_FILENAME, hpo);
        geneToDiseasesMap = parseMedgen(dataDir + MIM2GENE_MEDGEN_FILENAME);
        parseGeneSets(dataDir + GENE_SETS_FILENAME, gpiPathwayGenes, gpiAnchoredGenes);

        comSim = new ComputeSimilarity(hpo, diseaseMap, geneToDiseasesMap, allGenes,
                gpiPathwayGenes, gpiAnchoredGenes);
        icMap = comSim.computeICmap();
    }

    @After
    public void tearDown() throws Exception {
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
        PairwiseResnikSimilarity pwResSim = new PairwiseResnikSimilarity(hpo, icMap);
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