package org.monarchinitiative.dsppc;

import org.junit.Before;
import org.junit.Test;

import org.monarchinitiative.phenol.formats.hpo.HpoDisease;
import org.monarchinitiative.phenol.io.OntologyLoader;
import org.monarchinitiative.phenol.ontology.data.Ontology;
import org.monarchinitiative.phenol.ontology.data.TermId;

import java.io.File;
import java.util.*;

import static org.junit.Assert.*;
import static org.monarchinitiative.dsppc.Counter.*;
import static org.monarchinitiative.dsppc.Dsppc.*;

/*
 * created 08 Feb 2019
 */
public class CounterTest {
    private Counter ctr;
    private String dataDir = "src/main/resources/";
    private List<TermId> allGenes;
    private Map<TermId, HpoDisease> diseaseMap;
    private Map<TermId, Set<TermId>> geneToDiseasesMap;
    private Set<TermId> gpiAnchoredGenes = new TreeSet<>();
    private Set<TermId> gpiPathwayGenes = new TreeSet<>();
    private Ontology hpo;

    private static final int CARDINALITY = 5;

    @Before
    public void setUp() throws Exception {
        allGenes = parseAllGenes(dataDir + ALL_GENES_FILENAME);
        hpo = OntologyLoader.loadOntology(new File(dataDir + HPO_FILENAME));
        diseaseMap = parseHPOA(dataDir + HPOA_FILENAME, hpo);
        geneToDiseasesMap = parseMedgen(dataDir + MIM2GENE_MEDGEN_FILENAME);
        parseGeneSets(dataDir + GENE_SETS_FILENAME, gpiPathwayGenes, gpiAnchoredGenes);
        allGenes.removeAll(gpiPathwayGenes);
        ctr = new Counter(allGenes, diseaseMap, geneToDiseasesMap, CARDINALITY);
    }

    /*
    @Test
    public void countOneSetTest() {
        Set<TermId> genes = new HashSet<>(CARDINALITY);
        genes.add(TermId.of("ENTREZ","293"));
        genes.add(TermId.of("ENTREZ","720"));
        genes.add(TermId.of("ENTREZ","4589"));
        genes.add(TermId.of("ENTREZ", "200008"));
        genes.add(TermId.of("ENTREZ", "107986818"));
        int[] counts = ctr.countOneSet(genes);
        assertEquals("NUM_DISEASE_GENES is " + counts[NUM_DISEASE_GENES] + ", should be 2.",
                counts[NUM_DISEASE_GENES], 2);
        assertEquals("NUM_DISEASES is " + counts[NUM_DISEASES] + ", should be 2.",
                counts[NUM_DISEASES], 2);
        assertEquals("NUM_PHENOTYES is " + counts[NUM_PHENOTYPES] + ", should be 10.",
                counts[NUM_PHENOTYPES], 10);
    }
    */
}