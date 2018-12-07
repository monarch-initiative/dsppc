# dsppc: Disease set phenotypic profile comparison

### Goal
Compare two sets of genes to see whether their phenotypic profiles are more similar
than one would expect by chance. The two sets of interest are the GPI pathway genes, and the genes
that code for proteins anchored by the GPI pathway.

### Definitions
* GPI genes -- the list of genes associated with the GPI enzymes
* target set -- the genes that are modified with a GPI anchor from [list](https://docs.google.com/spreadsheets/d/1opZrVNZD5eNSLr0AyxmVSGXBLERF4Bg8NhiKrmRzAKs/edit#gid=0)

### Hypothesis
The phenotypes that characterize the GPI-gene related diseases, such as [Mabry syndrome](https://omim.org/entry/239300), which is caused by mutations in PIGV, are related to defects caused by misfunction of the target genes.

### Algorithm
#### Make a set of all HPOs associated with the GPI diseases
Optionally apply a filter to demand that a phenotype is seen in at least _k_ of the diseases
* For the GPI genes -> merge all HPOs seen in at least _k_ diseases into one set.
* For the target set -> Do not merge the HPOs into one big set. Instead, if a certain target gene is associated with multiple OMIM diseases, merge all of those HPOs into a set. Keep one such set for each target gene associated with one or more disease(s). Ignore target genes associated with no disease.

#### Calculate similarity between the GPI set and the target set
Try two similarity functions:
1. for each gene/disease associated with the target set, choose the best matching HPO term (highest IC score), and form the sum
2. for each gene/disease associated with the target set, take the sum of all matching HPO terms whose IC score exceeds a threshold _t_.
