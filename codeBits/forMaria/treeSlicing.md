test_tree_thinning
================
Janet Young
2023-10-27

Load a test tree and show some basic information on it. It’s an object
of class `phylo`.

Plot the tree (plot function comes from the ape package, I think - see
?plot.phylo).

``` r
plot(mammal.tree, cex=0.5)
add.scale.bar()
```

![](treeSlicing_files/figure-gfm/simple%20tree%20plot-1.png)<!-- -->

We’ll use the `cutree` function to slice the tree into some number of
clusters (5 in this example): it returns a vector of cluster IDs for
each taxon, but I turn that into a tibble.

``` r
cladeIDs <- dendextend::cutree(mammal.tree, k=5) %>% 
    as_tibble(rownames="taxon") %>% 
    dplyr::rename(clusterID=value) %>% 
    mutate(clusterID=factor(paste("cluster_",clusterID,sep="")))

cladeIDs %>% 
    head()
```

    ## # A tibble: 6 × 2
    ##   taxon         clusterID
    ##   <chr>         <fct>    
    ## 1 U._maritimus  cluster_1
    ## 2 U._arctos     cluster_1
    ## 3 U._americanus cluster_1
    ## 4 N._narica     cluster_1
    ## 5 P._lotor      cluster_1
    ## 6 M._mephitis   cluster_1

Tabulate cluster membership

``` r
cladeIDs %>% 
    count(clusterID)
```

    ## # A tibble: 5 × 2
    ##   clusterID     n
    ##   <fct>     <int>
    ## 1 cluster_1    13
    ## 2 cluster_2     6
    ## 3 cluster_3     3
    ## 4 cluster_4     3
    ## 5 cluster_5    24

`ggtree` is a nice package for tree plotting. The [ggtree ‘book’
website](https://yulab-smu.top/treedata-book/index.html) is helpful.

Here’s that simple plot again using ggtree:

``` r
ggtree(mammal.tree) +
    geom_tiplab(size=2.5) +
    hexpand(0.5)
```

![](treeSlicing_files/figure-gfm/siple%20ggtree%20plot-1.png)<!-- -->

But ggtree has lots of cool ways to display more info about the tree.

Let’s say we have info on the taxa in a tree, and maybe we’ve stored it
in a tibble. To use ggtree’s plotting functions, we want to make a
single object that contains the tree and the tibble together - that will
be an object of class `treedata` (defined by the `tidytree` package)

Maybe there’s an easier way to do that, but I wrote a function called
`addInfoToTree` to merge tree and taxon info that includes a couple of
sanity checks. `left_join` is the function that actually merges the tree
and the information.

``` r
###### addInfoToTree - make a treedata object by combining a phylo object with a tibble that has info on the tips
addInfoToTree <- function(tree, info, colnameForTaxonLabels="taxon") {
    ##### get info in same order as taxa in the tree:
    if (! colnameForTaxonLabels %in% colnames(info)) {
        stop("\n\nERROR - there should be a ",colnameForTaxonLabels," column in the info table\n\n")
    }
    ## check all taxa are in the info table
    tipLabelsInInfoTable <- info %>% select(all_of(colnameForTaxonLabels)) %>% deframe()
    if(length(setdiff(tree$tip.label, tipLabelsInInfoTable))>0) {
        stop("\n\nERROR - there are taxon labels in the tree that are not in the info table\n\n")
    }
    # now get info
    desiredRows <- match(tree$tip.label, tipLabelsInInfoTable)
    info_treeorder <- info[desiredRows,] 
    # add info to tree
    tree_withInfo <-  left_join(
        tree, 
        info_treeorder ,
        by=c("label"=colnameForTaxonLabels))
    return(tree_withInfo)
}
```

Combine the mammal.tree with the cladeID tibble we made using cutree

``` r
mammal.tree.withInfo <- addInfoToTree(mammal.tree, cladeIDs)
```

Plot the tree again and color the clades

``` r
ggtree(mammal.tree.withInfo,aes(color=clusterID) ) +
    geom_tiplab(aes(color=clusterID), size=2.5) + 
    hexpand(0.5) +
    theme(legend.position = "top")
```

![](treeSlicing_files/figure-gfm/plot%20tree%20with%20colored%20clades-1.png)<!-- -->

We could randomly choose one taxon per clade like this:

``` r
cladeMembersList <- split( cladeIDs$taxon, cladeIDs$clusterID)
oneMemberPerClade <- sapply(cladeMembersList, sample, size=1)
```

And we can update the info table associated with the tree to show which
taxa we chose, and label those tips on the tree plot:

``` r
mammal.tree.withInfo <- mammal.tree.withInfo %>% 
    mutate( chosen= label %in% oneMemberPerClade) 

ggtree(mammal.tree.withInfo,aes(color=clusterID) ) +
    geom_tiplab(aes(color=clusterID), size=2.5)+ 
    geom_tippoint(aes(subset= chosen), size=3) +
    hexpand(0.5) +
    theme(legend.position = "top")
```

![](treeSlicing_files/figure-gfm/label%20chosen%20taxa-1.png)<!-- -->
