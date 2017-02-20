## Parameters used by HINGE

All the parameters below can be set using the .ini file read by the HINGE programs.



###[filter]
- length_threshold = 1000; // Minimum read length
- aln_threshold = 2500; // Minimum alignment length between two reads to be considered when building graph
- min_cov = 5; // Minimum coverage depth for a segment on a read to not be considered erroneous/chimeric
- cut_off = 300; // When looking for chimeric segments, we look for coverage gaps on a read, after reducing all matches by cut_off in the beginning and in the end
- theta = 300; // When classifying a match between two reads as a right/left overlap, internal match, etc., overhangs of length up to theta are ignored
- use_qv = true; // Use qv scores provided by DAligner when creating the read masks (i.e., the part of the read that will actually be used for assembly)
- coverage = true; // Use coverage values when creating the read masks. If both use_qv and coverage are set to true, an intersection of the two masks is taken.
- coverage_frac_repeat_annotation = 3; 
- min_repeat_annotation_threshold = 10;
- max_repeat_annotation_threshold = 20; 

// A repeat annotation is placed on the read at position i+reso if 

``` |coverage[i]-coverage[i+reso]| > min( max( coverage[i+reso]/coverage_frac_repeat_annotation, min_repeat_annotation_threshold), max_repeat_annotation_threshold) ```

- repeat_annotation_gap_threshold = 300; // How far two hinges of the same type can be on a read
- no_hinge_region = 500; // Hinges cannot be placed within no_hinge_region of the start and end of the read
- hinge_min_support = 7; // Minimum number of reads that have to start in a `reso` (default 40) length interval to be considered in hinge calling
- hinge_unbridged = 6; // Number of reads that one has to see before a pileup to declare a potential hinge unbridged
- hinge_bin = 100; // Physical length of the bins considered
- hinge_tolerance_length = 100; // Matches starting within hinge_tolerance_length of a hinge are considered to be starting at the hinge

<!-- - quality_threshold = 0.23; // Quality threshold for edges to be considered in the backbone -->
<!--- n_iter = 2; // iterations of filtering, the filtering needs several iterations, because when filter reads, you got rid of some edges; when filter edges, you got rid of some reads (if the last edge is filtered.) Typically 2-3 iterations will be enough.-->
<!-- - theta2 = 0; // When classifying a match between two reads as a right/left overlap, internal match, etc., an overhang must have length at least theta2 for the match to be seen as internal -->


###[running]
- n_proc = 12; // number of CPUs for layout step




###[layout]

- hinge_tolerance = 150; // This is how far an overlap must start from a hinge to be considered an internal overlap.
- hinge_slack = 1000; // This is the amount by which  a forward overlap must be longer than a forward internal overlap to be preferred while building a graph.

- matching_hinge_slack = 200; // We identify two in-hinges (out-hinges) on two different reads as corresponding to the same repeat event, if the reads match in the repeat part, and the two hinges are within matching_hinge_slack of each other
- min_connected_component_size = 8; // In order to actually add a hinge to a read, we require that at least min_connected_component_size reads have a repeat annotation and they are all identified as the beginning (or end) of the same repeat
- kill_hinge_overlap = 300; 
- kill_hinge_internal = 40; 

// When filtering hinges (so that only one in-hinge and one out-hinge are left for each reapeat), we kill an in-hinge (out-hinge) if there is a forward (backward) extension read that starts at least kill_hinge_overlap before (after) the hinge, 
or if there is a forward_internal (backward_internal) extension read that starts at most kill_hinge_internal after (before) the hinge, as illustrated below. 

<img src="misc/param_description1.png" width=600px/>

- num_events_telomere = 7; 
- del_telomeres = 0; // If set to 1, any read with more than num_events_telomere repeat annotations will be classified as a telomere read and will be deleted.
- aggressive_pruning = 0; //If set to 1, the pruning will be more aggressive. We recommend it be set to 1 for large
genome.
- use_two_matches = 1; // Allow the HINGE algorithm to consider the top two matches between a pair of reads (as opposed to just the longest match)


###[draft]
<!--- min_cov = 10; //obsolete-->
<!--- trim = 200; //obsolete-->
<!--- edge_safe = 100; //obsolete-->
- tspace = 900; //space between new "trace points"
- step = 50;



###[consensus]
- min_length = 4000; // Minimal length of reads used for final consensus
- trim_end = 200; // Trim ends for alignments for final consensus
- best_n = 1; // If one read has multiple alignments with the bacbone assembly, choose the longest n segments for consensus.
- quality_threshold = 0.23; // alignment quality threshold


