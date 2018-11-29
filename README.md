# dsppc: Disease set phenotypic profile comparison

### Goal
Compare two sets of genes to see whether their phenotypic profiles are more similar
than one would expect by chance. The two sets of interest are the GPI pathway genes, and the genes
that code for proteins anchored by the GPI pathway. But any two disjoint sets of genes
can be given as input.

### Algorithm
1. For each set of genes (G1, G2), create a set containing all diseases related to those genes (D1, D2).
2. For each of D1, D2 create a set containing all phenotypes related to those diseases (P1, P2).
3. Compute the information content for every HPO phenotype term.
4. Compute the Resnik similarity of P1 and P2. Resnik similarity of two sets defined as the average of all
the pairwise Resnik similarity values for pairs whose first element is drawn from P1 and second element from P2.

To estimate the p value by simulation, hold D1 fixed and choose the second set of diseases at random from all HPO diseases.

5. Iterate _n_ times:
+ Select at random from all HPO diseases a set D3 equal in size to D2.
+ Create the corresponding set of phenotypes P3.
+ Compute the Resnik similarity of P1 and P3.
+ Maintain running count _m_ of how many times ResnikSim(P1, P3) >= ResnikSim(P1, P2).

6. Calculate p value: _m_ / _n_
