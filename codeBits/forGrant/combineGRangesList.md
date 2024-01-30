how to combine GrangesList into single set of regions
================
Janet Young

2024-01-30

Task : take a list of regions (e.g. each list item has regions
identified in a single sample) and make a single set of regions found in
any sample

I found a much simpler solution to the thing we struggled with on Friday
while reading the help page `?GRangesList`.

We used that weird `reduce(Reduce())` syntax, but `unlist` is a simpler
way to do it.

The help page says:

    unlist(x, recursive = TRUE, use.names = TRUE):
        Concatenates the elements of x into a single GRanges object.

Make example objects:

``` r
# make single GRanges object
gr1 <- GRanges(seqnames=c("chr1","chr1","chr2","chr2"),
               ranges=IRanges(start=c(1,100,1,100), width=20),
               score1=30)
# make ordinary list
gr_list <- list(first=gr1, second=gr1, third=gr1)
# make GRangesList
gr_GRlist <- GRangesList(gr_list)
```

`unlist`’s default is that it adds names to each region (shown on the
left) to show which of the original objects they came from:

``` r
unlist(gr_GRlist)
```

    ## GRanges object with 12 ranges and 1 metadata column:
    ##          seqnames    ranges strand |    score1
    ##             <Rle> <IRanges>  <Rle> | <numeric>
    ##    first     chr1      1-20      * |        30
    ##    first     chr1   100-119      * |        30
    ##    first     chr2      1-20      * |        30
    ##    first     chr2   100-119      * |        30
    ##   second     chr1      1-20      * |        30
    ##      ...      ...       ...    ... .       ...
    ##   second     chr2   100-119      * |        30
    ##    third     chr1      1-20      * |        30
    ##    third     chr1   100-119      * |        30
    ##    third     chr2      1-20      * |        30
    ##    third     chr2   100-119      * |        30
    ##   -------
    ##   seqinfo: 2 sequences from an unspecified genome; no seqlengths

Notice that this doesn’t work on a regular list object:

``` r
unlist(gr_list)
```

    ## $first
    ## GRanges object with 4 ranges and 1 metadata column:
    ##       seqnames    ranges strand |    score1
    ##          <Rle> <IRanges>  <Rle> | <numeric>
    ##   [1]     chr1      1-20      * |        30
    ##   [2]     chr1   100-119      * |        30
    ##   [3]     chr2      1-20      * |        30
    ##   [4]     chr2   100-119      * |        30
    ##   -------
    ##   seqinfo: 2 sequences from an unspecified genome; no seqlengths
    ## 
    ## $second
    ## GRanges object with 4 ranges and 1 metadata column:
    ##       seqnames    ranges strand |    score1
    ##          <Rle> <IRanges>  <Rle> | <numeric>
    ##   [1]     chr1      1-20      * |        30
    ##   [2]     chr1   100-119      * |        30
    ##   [3]     chr2      1-20      * |        30
    ##   [4]     chr2   100-119      * |        30
    ##   -------
    ##   seqinfo: 2 sequences from an unspecified genome; no seqlengths
    ## 
    ## $third
    ## GRanges object with 4 ranges and 1 metadata column:
    ##       seqnames    ranges strand |    score1
    ##          <Rle> <IRanges>  <Rle> | <numeric>
    ##   [1]     chr1      1-20      * |        30
    ##   [2]     chr1   100-119      * |        30
    ##   [3]     chr2      1-20      * |        30
    ##   [4]     chr2   100-119      * |        30
    ##   -------
    ##   seqinfo: 2 sequences from an unspecified genome; no seqlengths

You probably still want to use the `reduce()` function to merge
overlapping intervals that are present in more than one sample:

``` r
reduce(unlist(gr_GRlist))
```

    ## GRanges object with 4 ranges and 0 metadata columns:
    ##       seqnames    ranges strand
    ##          <Rle> <IRanges>  <Rle>
    ##   [1]     chr1      1-20      *
    ##   [2]     chr1   100-119      *
    ##   [3]     chr2      1-20      *
    ##   [4]     chr2   100-119      *
    ##   -------
    ##   seqinfo: 2 sequences from an unspecified genome; no seqlengths

Exploring `unlist()` a bit more (notice I’m not using `reduce()` now):

If I know I don’t need to keep track of the original objects, I use the
option `use.names = FALSE`. Having duplicate names for regions can
sometimes cause trouble down the road

``` r
unlist(gr_GRlist, use.names = FALSE)
```

    ## GRanges object with 12 ranges and 1 metadata column:
    ##        seqnames    ranges strand |    score1
    ##           <Rle> <IRanges>  <Rle> | <numeric>
    ##    [1]     chr1      1-20      * |        30
    ##    [2]     chr1   100-119      * |        30
    ##    [3]     chr2      1-20      * |        30
    ##    [4]     chr2   100-119      * |        30
    ##    [5]     chr1      1-20      * |        30
    ##    ...      ...       ...    ... .       ...
    ##    [8]     chr2   100-119      * |        30
    ##    [9]     chr1      1-20      * |        30
    ##   [10]     chr1   100-119      * |        30
    ##   [11]     chr2      1-20      * |        30
    ##   [12]     chr2   100-119      * |        30
    ##   -------
    ##   seqinfo: 2 sequences from an unspecified genome; no seqlengths

Sometimes I do want to keep track of the original objects, and I find it
works better to do make a metadata column to contain that information:

``` r
combined_gr_GRlist <- unlist(gr_GRlist)
# put the names in a metadata column
mcols(combined_gr_GRlist)[,"whichOriginalGR"] <- names(combined_gr_GRlist)
# remove the names
names(combined_gr_GRlist) <- NULL
combined_gr_GRlist
```

    ## GRanges object with 12 ranges and 2 metadata columns:
    ##        seqnames    ranges strand |    score1 whichOriginalGR
    ##           <Rle> <IRanges>  <Rle> | <numeric>     <character>
    ##    [1]     chr1      1-20      * |        30           first
    ##    [2]     chr1   100-119      * |        30           first
    ##    [3]     chr2      1-20      * |        30           first
    ##    [4]     chr2   100-119      * |        30           first
    ##    [5]     chr1      1-20      * |        30          second
    ##    ...      ...       ...    ... .       ...             ...
    ##    [8]     chr2   100-119      * |        30          second
    ##    [9]     chr1      1-20      * |        30           third
    ##   [10]     chr1   100-119      * |        30           third
    ##   [11]     chr2      1-20      * |        30           third
    ##   [12]     chr2   100-119      * |        30           third
    ##   -------
    ##   seqinfo: 2 sequences from an unspecified genome; no seqlengths
