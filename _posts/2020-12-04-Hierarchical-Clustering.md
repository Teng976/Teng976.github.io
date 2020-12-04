---
layout: post
title: Hierarchical Clustering of Posterior Means
---

## Hierarchical Clustering

Hierarchical clustering is one of the ubiquitous clustering algorithms,
and is widely used in genomics. We apply hierarchical clustering to
posterior means of DayCent ecosystem model parameters from various
regions around the globe.

## Toy sample illustration

First we generate a toy sample for illustration of the method. `toy` is
a 15 by 2 matrix, simulated from Gaussians with various means
and variances. Every five vectors share the same mean and variance which
make the data consist of three clusters.

``` r
set.seed(2020)
mu1 <- c(1, 3, 5)
sigma1 <- c(1, 0.5, 1)
mu2 <- c(1, 3, 5)*10
sigma2 <- c(1, 0.5, 2)*2
toy <- matrix(c(rnorm(5*3, rep(mu1,each=5), rep(sigma1,each=5)), 
                rnorm(5*3, rep(mu2,each=5), rep(sigma2,each=5))), ncol=2)
```

We scale the data so that each dimension has same weight in determining
the clusters. We calculate the distance between data points by `dist`,
and the default metric is the Euclidean distance. `hclust` does
hierarchical clustering. Various methods exist. Check documentation.

``` r
toy_s <- scale(toy)
d <- dist(toy_s)
hc <- hclust(d, method = 'complete')
```

Plot the cluster dendrogram and pick height `h` to determine the
clusters. Put boxes on the figure.

``` r
plot(hc)
rect.hclust(hc, h=1.5, border=4)
```

![](HierarchicalClustering_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

You can also use `cutree` to cluster either by number of clusters `k` or
by height `h`.

``` r
hc_k3 = cutree(hc, k=3)
table(hc_k3)
```

    ## hc_k3
    ## 1 2 3 
    ## 5 5 5

``` r
hc_h2 = cutree(hc, h=1.5)
table(hc_h2)
```

    ## hc_h2
    ## 1 2 3 
    ## 5 5 5

## Application to Posterior Means of Daycent Model Parameters

We read in the real data and apply hierarchical clustering.

