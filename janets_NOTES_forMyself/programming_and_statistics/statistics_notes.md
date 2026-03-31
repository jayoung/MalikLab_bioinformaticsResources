# Random older notes

QQ plot:   for normal distribution, should be a straight line.

Robert Gentleman’s comments on modelling. Linear modelling is one option. Lasso is another option – drops non-significant terms

(I think from Charles Cooperberg, not sure): if averaging, better to average after taking log rather than vice versa (trims off skew better)

From Li Hsu:
- K-S test is on cumulative distributions
- the 95% confidence band will be wider than the 95% confidence interval for a specific point on the curve, because it's giving confidence in the whole modelled curve, not just that portion of it.

UpSet plots - an alternative to Venn diagrams. Quite nice. 
http://caleydo.org/tools/upset/



# Seminars 

## Feb 26, 2026  Timothy Barry (candidate for something in PHS)

Postdoctoral Daniel Bauer. PhD in Statistics (Kathryn Roeder and Eugene Katsevich)

Key publications related to the talk:
- T Barry, Z Niu, E Katsevich, X Lin. “The permuted score test for robust differential expression analysis.” ArXiv preprint, 2025. arxiv.org/abs/2501.03530.
- T Barry, K Mason, K Roeder, E Katsevich. “Robust differential expression testing for single-cell CRISPR screens at low multiplicity of infection.” Genome Biology, 2024. link.springer.com/article/10.1186/s13059-024-03254-2.

Also see the work of [Jingyi Jessica Li](https://jsb-lab.org/) (Hutch PI)

He thinks a lot about the statistics of CRISPR screens / single cell RNA-seq, especially in the context of single cells (CRISPR+RNA-seq in single cells = Perturb-seq)

He basically adapts the negative binomial model (commonly used in DESeq etc) to be more robust to violations of models. I think this is done using perturbations.

DESeq etc make heavy use of the Negative Binomial (NB) distribution:
- has a dispersion parameter, phi. 
- when phi=0, Negative Binomial is the same as the Poisson distribution


DESeq etc use Negative Binomial GLMs ("Generalized linear model") for regression, to test whether treatments (etc) have an effect on counts.  Seurat (single cells) uses Mann-Whitney tests and does not adjust for covariates.

Rather than doing simple NB GLMs, Tim proposes PERMUTING the data abd then doing some test. 
- you perform the NB GLM for each gene on the REAL data
- you permute the treatment labels a bunch of times and repeat the NB GLM test, to see how the real score compares to the permuted distribution
- you can make the permutations more computationally EFFICIENT by making them ADAPTIVE (often you don't need to do many permutations to see that a gene is NOT significant, in which case you can stop, but if it seems interesting, then you do more permutations)
- apparently this is robust to batch effects on the dispersion parameter
- permuations give you a synthetic negative control

Research output:
- 'robustified' versions of DESeq2 and MASS (add permutations). Improces performance in various ways, especially the FDR rate.
