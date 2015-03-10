#!/usr/bin/perl -w

use PDL;  

# =========================================================================

# The clustering algorithm works great, but clustering a 90000 point field
# is an almost intractable problem (the distance matrix alone will be 30 Gb!)
# So we have to think of another way to do it.
# Perhaps do it in 2 steps:
# - Subsample the original field (select 300 indices or so at random)
# - Create clusters, calculate the cluster averages
# - Assign all 90000 points to a cluster
#	- Create a new membership vector
#	- For each point, calculate the distance to all clusters
#	- Assign to the cluster with the smallest distance
#	- Export the new membership vector

sub cluster_species {
	
	my ($properties, $threshold) = @_;
	
	my $numprop = $properties->getdim(1);
	my $veclen = $properties->getdim(0);
	
	$sample_idx = floor(random(300) * $veclen);
	
	print $sample_idx;
	
	$sample_values = $properties->dice($sample_idx);
	
	print $sample_values;
	
	my ($numclust, $sample_clust_mean, @sample_clust) = cluster_matrix($sample_values, $threshold);
	
	print "$numclust clusters made, with mean values: " . $sample_clust_mean . "\n";
	
	# For each point, calculate the distance to all clusters
	# Assign each point to the cluster with the smallest distance
	
	my $cluster = zeroes($veclen);
	
	for (my $ind = 0 ; $ind < $veclen - 1 ; $ind++) {
		
		my $distance = $sample_clust_mean->copy;
		
		$distance -= $properties->slice("$ind,:");
		$distance *= $distance;
		$distance = sqrt(sumover($distance->xchg(0,1)));
		
		$cluster->slice("$ind") .= minimum_ind($distance);
	}
	
	
#	TODO: If smallest distance exceeds threshold, repeat clustering
#	with a bigger sample
	
	
	# Export the new membership vector and mean values
	
	my $mean_properties = zeroes($numclust, $numprop);
	my @outclusters;
	
	for (my $clusteridx = 0 ; $clusteridx < $numclust ; $clusteridx++) {
		my $clustermap = which($cluster == $clusteridx);
		$outclusters[$clusteridx] = $clustermap;
		
		$mean_properties->slice("($clusteridx)") .= average($properties->dice($clustermap));
	}
	
	return ($numclust, $mean_properties, @outclusters);
}


sub cluster_matrix {
	
	my ($properties, $threshold) = @_;
	
	# properties: a matrix consisting of vectors that contain organism properties
	# threshold: minimum distance between clusters
	
	my $numprop = $properties->getdim(1);		# Check this
	my $veclen = $properties->getdim(0);
	
	# Some settings to control the behaviour of the algorithm
	
	my $keep_orginal_distmat = 1;
	
	print "Clustering $numprop properties, $veclen points\n";
	
	# Maximum distance in property space
	#my $maxdist = sqrt($numprop * ($maxval * $maxval));
	my $mindist;
	my $dummy = -1;
	
	my $distmat = ones($veclen, $veclen) * $dummy;	# Distance matrix
	my $cluster = zeroes($veclen);			# Cluster membership vector
	
	print "Allocated distance matrix and cluster membership vector...\n";
	
	my $clusteridx = 1;
	my $clustercount = 0;
	my @clusteralias;
	my $mask = zeroes($veclen, $veclen);
	my $update_cluster = 0;
	
	print "Allocated mask...\nCalculating distances...\n";
	
	for (my $ind = 0 ; $ind < $veclen - 1 ; $ind++) {
		
		# Calculate euclidian distances between this point
		# and all other points
		
		my $first = $ind + 1;
		my $last = $veclen - 1;
		my $distance = $properties->slice("$first:$last,:")->copy;
		#print $distance . "\n" . $properties->slice("$ind,:");
		$distance -= $properties->slice("$ind,:");
		$distance *= $distance;
		$distance = sqrt(sumover($distance->xchg(0,1)));
		
		
		# Write result to distance matrix
		
		$distmat->slice("($ind),$first:$last") .= $distance;
		#$mask->slice("($ind),:") .= ($distmat->slice("($ind),:") != $dummy);
		#$distmat->slice("$first:$last,($ind)") .= $distance;
		
		# (Note that we only calculate one half of the matrix,
		#  but we fill the entire matrix for later use)
	}
	
	print "Finalising distance matrix...\n";
	
	# Make a mask of the distances that matter for clustering
	
	$mask = ($distmat != $dummy);
	
	for (my $ind = 0 ; $ind < $veclen - 1 ; $ind++) {
		
		#  Fill the entire distance matrix for later use
		
		my $first = $ind + 1;
		my $last = $veclen - 1;
		$distmat->slice("$first:$last,($ind)") .= $distmat->slice("($ind),$first:$last");
	}
	
	
	my $original_distmat = $distmat;
	
	if ($keep_original_distmat) {
		# Make a physical copy of the distance matrix 
		# (this is needed for proper UPGMA centroid calculation,
		# but does double memory use)
		$original_distmat = $distmat->copy;
	}
	
	print "Clustering points...\n";
	
	do {
		
		my $flat_mask = $mask->clump(-1);
		my $valid_indices = which($flat_mask);
		my $valid_values = $distmat->clump(-1)->where($flat_mask);
		my $min_dist_idx = $valid_indices->at(minimum_ind($valid_values));
		
		my $pair_row = $min_dist_idx % $veclen;
		my $pair_col = int($min_dist_idx / $veclen);
		$mindist = $distmat->at($pair_row,$pair_col);
		
		#print "Distance matrix: " . $distmat . "\n";
		#print "Mask: " . $mask . "\n";
		
		print "Minimum distance in matrix: $mindist ($pair_col, $pair_row)\n";
		
		if ($mindist < $threshold) {
			
			my $rowclust = $cluster->at($pair_row);
			my $colclust = $cluster->at($pair_col);
			
			if ($colclust) {
				
				if ($rowclust) {
					
					if ($rowclust == $colclust) {
						print "Both points are already in cluster $rowclust\n";
						
					} else {
						print "Both points are in clusters, adding cluster $rowclust to $colclust\n";
						
						#push @{ $clusteralias[$colclust] }, $rowclust;
						$cluster->where($cluster == $rowclust) .= $colclust;
						$clusteralias[$rowclust] = $colclust;
						$update_cluster = $colclust;
						$clustercount--;
					}
					
				} else {
					
					print "Point $pair_col is already in cluster $colclust, adding $pair_row\n";
					
					$cluster->slice("$pair_row") .= $colclust;
					$update_cluster = $colclust;
					
				}
				
				
			} elsif ($rowclust) {
				
				print "Point $pair_row is already in cluster $rowclust, adding $pair_col\n";
				
				$cluster->slice("$pair_col") .= $rowclust;
				$update_cluster = $rowclust;
				
			} else {
				
				print "Clustering $pair_col and $pair_row in cluster $clusteridx\n";
				
				$clustercount++;
				$cluster->slice("$pair_col") .= $clusteridx;
				$cluster->slice("$pair_row") .= $clusteridx;
				$clusteralias[$clusteridx] = 0;
				$update_cluster = $clusteridx;
				#push @{ $clusteralias[$clusteridx] }, $clusteridx;
				$clusteridx++;
			}
			
			#$distmat->slice("$pair_row,$pair_col") .= $maxdist;
			$mask->slice("$pair_row,$pair_col") .= 0;
			
			# In the simplest case (single linkage) we could
			# just use the minimum distance between points
			# as cluster distances.
			# But for better clustering we have to somehow replace 
			# all values in this row and column of the distance 
			# matrix by either the maximum (complete linkage)
			# or mean (UPGMA) distance to our new cluster.
			
			# So:
			# - Use $cluster->which($cluster == $update_cluster)
			#   to extract an index piddle.
			
			if ($update_cluster) {
				
				# Recalclulate distances to this cluster.
				
				# Construct an index vector of points that need
				# to be updated
				
				print "Updating cluster $update_cluster...\n";
				#print $cluster . "\n";
				#print ($cluster == $update_cluster) . "\n";
				
				my $update_points = which($cluster == $update_cluster);
				
				# Use the index piddle as argument to dice,
				# to extract rows/columns from the distance matrix
				
				my $subdistmat = $original_distmat->dice($update_points);
				
				# Average (or maximise) over the rows/columns
				
				my $clust_average = average($subdistmat);
				
				# Write the resulting vector with the cluster 
				# distance to all other points back into the
				# current row and column of the distance matrix
				
				$distmat->dice(X,$update_points) .= $clust_average;
				$distmat->dice($update_points) .= transpose($clust_average);
				
				# Now we need to clear the mask for all points that
				# are in this cluster, relative to each other.
				# Otherwise we can get strange results (negative distances, etc).
				
				$mask->dice($update_points,$update_points) .= 0;
				
				$update_cluster = 0;
			}
		}
		
	} while ($mindist < $threshold);
	
	print "Number of clusters made: $clustercount\n";
	
	my $residuals = which($cluster == 0);
	my $nr_residuals = $residuals->getdim(0);
	print "Residual points (not assigned to clusters): $nr_residuals\n";
	
	# Now we should loop over the clusters (and residual points), 
	# and return an array of index piddles for each cluster.
	# Also return the number of clusters, and a matrix with the
	# average property values for each cluster.
	# We can use that later to assign points or other clusters
	# to these clusters.
	
	my @finalclusters;
	my $last = $clusteridx;
	my $outclusteridx = 0;
	my $mean_properties = zeroes($clustercount + $nr_residuals, $numprop);
	
	for ($clusteridx = 1 ; $clusteridx < $last ; $clusteridx++) {
		if ($clusteralias[$clusteridx] == 0) {
				print "Cluster $clusteridx is a valid cluster, assigning new index $outclusteridx\n";
				my $clustermap = which($cluster == $clusteridx);
				$finalclusters[$outclusteridx] = $clustermap;
				
				$mean_properties->slice("($outclusteridx)") .= average($properties->dice($clustermap));
				
				$outclusteridx++;
		}
	}
	
	for ($clusteridx = 0 ; $clusteridx < $nr_residuals ; $clusteridx++) {
		print "Assigning residual point " . $residuals->at($clusteridx) .
			" as cluster at index $outclusteridx\n";
		
		$finalclusters[$outclusteridx] = $residuals->slice("$clusteridx");
		my $point_index = $residuals->at($clusteridx);
		$mean_properties->slice("($outclusteridx)") .= $properties->slice("($point_index)");
		
		$outclusteridx++;
	}
	
	print "Total clusters: $outclusteridx\n";
	
	print "Mean properties: " . $mean_properties . "\n";
	
	return ($outclusteridx, $mean_properties, @finalclusters);
}
