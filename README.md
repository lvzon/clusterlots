# clusterlots

This is a fast perl/PDL implementation of UPGMA data clustering for very large datasets. The problem is that existing clustering packages have difficulty handling datasets with more than a few thousand data points. Especially the distance matrices tend to become a problem. For example, clustering one field of a 300x300 grid-based simulation (90,000 data points) would require a (non-sparse) distance matrix of 8.1 billion entries (90000^2). This would use over 30 Gb of memory when stored as a block of 4-byte floating point values.

I needed to cluster a lot of such datasets, so to make this manageable I implemented a simple but fast UPGMA clustering algorithm in PDL (the Perl Data Language). To conserve memory, it doesn't store a full hierarchical clustering tree, but rather it partitions the data into clusters based on a pre-defined minimum distance between the cluster centroids. Moreover, to handle large data sets it can compute the clusters based on a random sample of the data points, and then assign all other points to the nearest cluster. This is fairly fast and seems to work very well in practice.

The clustering subroutine has the following parameters:
   -  an input-matrix, which is essentially a concatenated set of vectors. The length of these vectors is the number of points to be clustered, and there is one vector for each variable or dimension.
   -  a threshold, indicating the minimum (euclidian) distance that is required between separate cluster centroids.
   
The algorithm works as follows:

 1. A distance matrix is calculated, which contains the euclidian distance between each point and each other point. 
 2. A physical copy is made of the distance matrix (needed for proper UPGMA centroid calculation, but it does double memory use).
 3. If there are any points/clusters left which have a distance smaller than the threshold distance, the two points with the smallest distance are selected from the distance matrix.
 4. If any or both of these points already belong to a cluster, these points and/or clusters are merged. If neither point belongs to a cluster, a new cluster is created, containing both points.
 5. The working-copy of the distance matrix is updated. The distance of other points to points in the new cluster is recalculated as the mean of all original distances (taken from the copy of the original distance matrix).
 6. The steps 3-5 are repeated until no values below the threshold remain in the distance matrix.

Note that, if there are outliers, the final set of clusters may include clusters of only one point.

When finished, the subroutine returns three variables:
 
 - The number of clusters that were created (as a Perl scalar).
 - The average variable values (average locations in each dimension) for each cluster, as a PDL matrix (clusters x variables).
 - An array with an index vector for each cluster, which contains the indices of the cluster members.

