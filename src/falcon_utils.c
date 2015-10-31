
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include "common.h"



consensus_data * generate_consensus_large( char ** input_seq,
                                     unsigned int n_seq,
                                     unsigned min_cov,
                                     unsigned K,
                                     unsigned long local_match_count_window,
                                     unsigned long local_match_count_threshold,
                                     double min_idt) {
    // local_match_count_window, local_match_count_threshold obsoleted, keep the interface for a while

    unsigned int j;
    unsigned int seq_count;
    unsigned int aligned_seq_count;
    kmer_lookup * lk_ptr;
    seq_array sa_ptr;
    seq_addr_array sda_ptr;
    kmer_match * kmer_match_ptr;
    aln_range * arange;
    alignment * aln;
    align_tags_t ** tags_list;
    //char * consensus;
    consensus_data * consensus;
    double max_diff;
    max_diff = 1.0 - min_idt;

    seq_count = n_seq; // n sequences
    printf("XX n_seq %d\n", n_seq);
    for (j=0; j < seq_count; j++) {
        printf("seq_len: %u %u\n", j, strlen(input_seq[j]));
    };
    fflush(stdout);

    tags_list = calloc( seq_count, sizeof(align_tags_t *) );
    lk_ptr = allocate_kmer_lookup( 1 << (K * 2) );
    sa_ptr = allocate_seq( (seq_coor_t) strlen( input_seq[0]) ); // seed sequence
    sda_ptr = allocate_seq_addr( (seq_coor_t) strlen( input_seq[0]) );
    printf("before add sequence\n");
    add_sequence( 0, K, input_seq[0], strlen(input_seq[0]), sda_ptr, sa_ptr, lk_ptr);
    //mask_k_mer(1 << (K * 2), lk_ptr, 16);
    // sequence 0 is the base sequence
    // others are pile ups
    printf("after add sequence\n");
    char * base_str = input_seq[0];


    aligned_seq_count = 0;
    for (j=1; j < seq_count; j++) {

        printf("seq_len: %ld %u\n", j, strlen(input_seq[j]));

        kmer_match_ptr = find_kmer_pos_for_seq(input_seq[j], strlen(input_seq[j]), K, sda_ptr, lk_ptr);
#define INDEL_ALLOWENCE_0 6

        arange = find_best_aln_range(kmer_match_ptr, K, K * INDEL_ALLOWENCE_0, 5);  // narrow band to avoid aligning through big indels

        //printf("1:%ld %ld %ld %ld\n", arange_->s1, arange_->e1, arange_->s2, arange_->e2);

        //arange = find_best_aln_range2(kmer_match_ptr, K, K * INDEL_ALLOWENCE_0, 5);  // narrow band to avoid aligning through big indels

        //printf("2:%ld %ld %ld %ld\n\n", arange->s1, arange->e1, arange->s2, arange->e2);

#define INDEL_ALLOWENCE_1 0.10
        if (arange->e1 - arange->s1 < 100 || arange->e2 - arange->s2 < 100 ||
                                                                       abs( (arange->e1 - arange->s1 ) - (arange->e2 - arange->s2) ) >
                    (int) (0.5 * INDEL_ALLOWENCE_1 * (arange->e1 - arange->s1 + arange->e2 - arange->s2))) {
            free_kmer_match( kmer_match_ptr);
            free_aln_range(arange);
            continue;
        }
        //printf("%ld %s\n", strlen(input_seq[j]), input_seq[j]);
        //printf("%ld %s\n\n", strlen(input_seq[0]), input_seq[0]);


#define INDEL_ALLOWENCE_2 150
        // WTF are these heuristics
        //printf("reach\n");
        aln = align(input_seq[j]+arange->s1, arange->e1 - arange->s1 ,
                    input_seq[0]+arange->s2, arange->e2 - arange->s2 ,
                    INDEL_ALLOWENCE_2, 1);
        if (aln->aln_str_size > 500 && ((double) aln->dist / (double) aln->aln_str_size) < max_diff) {
            tags_list[aligned_seq_count] = get_align_tags( aln->q_aln_str,
                                                           aln->t_aln_str,
                                                           aln->aln_str_size,
                                                           arange, j,
                                                           0);
            aligned_seq_count ++;


            char * s = aln->q_aln_str;
            int nq = 0;
            while (*s++) {
                if (*s != '-') nq++;
            }
            s = aln->t_aln_str;
            int nt = 0;
            while (*s++) {
                if (*s != '-') nt++;
            }


            //printf("Q:%s\nT:%s\n%d\n%d %d\n", aln->q_aln_str, aln->t_aln_str, aln->aln_str_size, nq, nt);
            //printf("[%d %d] x [%d %d]\n", arange->s1, arange->e1, arange->s2, arange->e2);
            //printf("j:%s\n0:%s\n", input_seq[j] + arange->s1, input_seq[0] + arange->s2);

        } // get align tags

        /*
        int k;
        int j;
        j = aligned_seq_count - 1;
        for (k = 0; k < tags_list[j]->len; k++) {
            printf("%d %d %ld %d %c %c\n",j, k, tags_list[j]->align_tags[k].t_pos,
                                   tags_list[j]->align_tags[k].delta,
                                   //tags_list[j]->align_tags[k].p_q_base,
                   base_str[tags_list[j]->align_tags[k].t_pos],
                   tags_list[j]->align_tags[k].q_base);
        }
        printf("seq:%s\n",input_seq[j]);
        */
        free_aln_range(arange);
        free_alignment(aln);
        free_kmer_match( kmer_match_ptr);
    }

    if (aligned_seq_count > 0) {
        consensus = get_cns_from_align_tags_large( tags_list, aligned_seq_count, strlen(input_seq[0]), min_cov );
    } else {
        // allocate an empty consensus sequence
        consensus = calloc( 1, sizeof(consensus_data) );
        consensus->sequence = calloc( 1, sizeof(char) );
        consensus->eqv = calloc( 1, sizeof(unsigned int) );
    }
    //free(consensus);
    free_seq_addr_array(sda_ptr);
    free_seq_array(sa_ptr);
    free_kmer_lookup(lk_ptr);
    for (j=0; j < aligned_seq_count; j++) {
        free_align_tags(tags_list[j]);
    }
    free(tags_list);
    return consensus;
}

align_tags_t * get_align_tags2( char * aln_q_seq, 
                               char * aln_t_seq, 
                               seq_coor_t aln_seq_len,
                               aln_range * range,
                               unsigned q_id,
                               seq_coor_t t_offset) {
    char p_q_base;
    align_tags_t * tags;
    seq_coor_t i, j, jj, k, p_j, p_jj;

    tags = calloc( 1, sizeof(align_tags_t) );
    tags->len = aln_seq_len; 
    tags->align_tags = calloc( aln_seq_len + 1, sizeof(align_tag_t) );
    i = range->s1 - 1;
    j = range->s2 - 1;
    jj = 0;
    p_j = -1;
    p_jj = 0;
    p_q_base = '.';

    for (k = 0; k < aln_seq_len; k++) {
        if (aln_q_seq[k] != '-') {
            i ++;
            jj ++;
        } 
        if (aln_t_seq[k] != '-') {
            j ++;
            jj = 0;
        }
        //printf("t %d %d %d %c %c\n", q_id, j, jj, aln_t_seq[k], aln_q_seq[k]);
       
       
        if ( j + t_offset >= 0 && jj < UINT8_MAX && p_jj < UINT8_MAX) {
            (tags->align_tags[k]).t_pos = j + t_offset;
            (tags->align_tags[k]).delta = jj;
            (tags->align_tags[k]).p_t_pos = p_j + t_offset;
            (tags->align_tags[k]).p_delta = p_jj;
            (tags->align_tags[k]).p_q_base = p_q_base;
            (tags->align_tags[k]).q_base = aln_q_seq[k];
            (tags->align_tags[k]).q_id = q_id;
            
            p_j = j;
            p_jj = jj;
            p_q_base = aln_q_seq[k];
        }
    }
    // sentinal at the end
    //k = aln_seq_len;
    tags->len = k; 
    (tags->align_tags[k]).t_pos = UINT_MAX;
    (tags->align_tags[k]).delta = UINT8_MAX;
    (tags->align_tags[k]).q_base = '.';
    (tags->align_tags[k]).q_id = UINT_MAX;
    return tags;
}


consensus_data * get_cns_from_align_tags_large( align_tags_t ** tag_seqs, 
                                          unsigned n_tag_seqs, 
                                          unsigned t_len, 
                                          unsigned min_cov ) {

    seq_coor_t i, j, t_pos;
    unsigned int * coverage;
    unsigned int * local_nbase;

    consensus_data * consensus;
    //char * consensus;
    align_tag_t * c_tag;
    static msa_pos_t * msa_array = NULL;

    coverage = calloc( t_len, sizeof(unsigned int) );
    local_nbase = calloc( t_len, sizeof(unsigned int) );

#ifndef STATIC_ALLOCATE

    msa_array = calloc(t_len, sizeof(msa_pos_t *));

    for (i = 0; i < t_len; i++) {
        msa_array[i] = calloc(1, sizeof(msa_delta_group_t));
        msa_array[i]->size = 8;
        allocate_delta_group(msa_array[i]);
    }

#endif    

#ifdef STATIC_ALLOCATE

    if ( msa_array == NULL) {
        msa_array = get_msa_working_sapce( 6000000 );
    } 

    assert(t_len < 6000000);

#endif    

#define DBG

	
    
    // loop through every alignment
    //printf("XX %d\n", n_tag_seqs);
    for (i = 0; i < n_tag_seqs; i++) {

        // for each alignment position, insert the alignment tag to msa_array
        for (j = 0; j < tag_seqs[i]->len; j++) {
            c_tag = tag_seqs[i]->align_tags + j;
            unsigned int delta;
            delta = c_tag->delta;
            if (delta == 0) {
                t_pos = c_tag->t_pos;
                coverage[ t_pos ] ++;
            }
            // Assume t_pos was set on earlier iteration.
            if (delta > msa_array[t_pos]->max_delta) {
                msa_array[t_pos]->max_delta = delta;
                if (msa_array[t_pos]->max_delta + 4 > msa_array[t_pos]->size ) {
                    realloc_delta_group(msa_array[t_pos], msa_array[t_pos]->max_delta + 8);
                }
            }
            
            unsigned int base;
            switch (c_tag->q_base) {
                case 'A': base = 0; break;
                case 'C': base = 1; break;
                case 'G': base = 2; break;
                case 'T': base = 3; break;
                case '-': base = 4; break;
            }
            // Note: On bad input, base may be uninitialized.
            update_col( &(msa_array[t_pos]->delta[delta].base[base]), c_tag->p_t_pos, c_tag->p_delta, c_tag->p_q_base);
            local_nbase[ t_pos ] ++;
        }
    }

    // propogate score throught the alignment links, setup backtracking information
    align_tag_col_t * g_best_aln_col = 0;
    unsigned int g_best_ck = 0;
    seq_coor_t g_best_t_pos = 0;
    {
        int kk; 
        int ck;
        // char base;
        int best_i;
        int best_j;
        int best_b;
        int best_ck = -1;
        double score;
        double best_score;
        double g_best_score;
        // char best_mark;

        align_tag_col_t * aln_col;
        
        g_best_score = -1;

        for (i = 0; i < t_len; i++) {  //loop through every template base
            //printf("max delta: %d %d\n", i, msa_array[i]->max_delta);
            for (j = 0; j <= msa_array[i]->max_delta; j++) { // loop through every delta position
                for (kk = 0; kk < 5; kk++) {  // loop through diff bases of the same delta posiiton
                    /*
                    switch (kk) {
                        case 0: base = 'A'; break;
                        case 1: base = 'C'; break;
                        case 2: base = 'G'; break;
                        case 3: base = 'T'; break;
                        case 4: base = '-'; break;
                    }
                    */
                    aln_col = msa_array[i]->delta[j].base + kk;
                    if (aln_col->count >= 0) {
                        best_score = -1;
                        best_i = -1;
                        best_j = -1;
                        best_b = -1;

                        for (ck = 0; ck < aln_col->n_link; ck++) { // loop through differnt link to previous column
                            int pi;
                            int pj;
                            int pkk;
                            pi = aln_col->p_t_pos[ck];
                            pj = aln_col->p_delta[ck];
                            switch (aln_col->p_q_base[ck]) {
                                case 'A': pkk = 0; break;
                                case 'C': pkk = 1; break;
                                case 'G': pkk = 2; break;
                                case 'T': pkk = 3; break;
                                case '-': pkk = 4; break;
                                default: pkk = 4;
                            }

                            if (aln_col->p_t_pos[ck] == -1) {
                                score =  (double) aln_col->link_count[ck] - (double) coverage[i] * 0.5;
                            } else {
                                score = msa_array[pi]->delta[pj].base[pkk].score + 
                                        (double) aln_col->link_count[ck] - (double) coverage[i] * 0.5;
                            }
                            // best_mark = ' ';
                            if (score > best_score) {
                                best_score = score;
                                aln_col->best_p_t_pos = best_i = pi;
                                aln_col->best_p_delta = best_j = pj;
                                aln_col->best_p_q_base = best_b = pkk;
                                best_ck = ck;
                                // best_mark = '*';
                            } // find max
                            /*
                            printf("X %d %d %d %c %d %d %d %c %d %lf %c\n", coverage[i], i, j, base, aln_col->count,
                                                                  aln_col->p_t_pos[ck], 
                                                                  aln_col->p_delta[ck], 
                                                                  aln_col->p_q_base[ck], 
                                                                  aln_col->link_count[ck],
                                                                  score, best_mark);
                            */
                        }
                        aln_col->score = best_score;
                        if (best_score > g_best_score) {
                            g_best_score = best_score;
                            g_best_aln_col = aln_col;
                            g_best_ck = best_ck;
                            g_best_t_pos = i;
                            //printf("GB %d %d %d %d\n", i, j, ck, g_best_aln_col);
                        }
                    }
                }
            }
        }
    }
    //assert(g_best_score != -1);

    // reconstruct the sequences
    unsigned int index;
    char bb;
    int ck;
    char * cns_str;
    int * eqv;
    double score0;
    
    consensus = calloc( 1, sizeof(consensus_data) );
    consensus->sequence = calloc( t_len * 2 + 1, sizeof(char) );
    consensus->eqv = calloc( t_len * 2 + 1, sizeof(unsigned int) );
    cns_str = consensus->sequence;
    eqv =  consensus->eqv;

    index = 0;
    ck = g_best_ck;
    i = g_best_t_pos;

    while (1) {
        if (coverage[i] > min_cov) {
            switch (ck) {
                case 0: bb = 'A'; break;
                case 1: bb = 'C'; break;
                case 2: bb = 'G'; break;
                case 3: bb = 'T'; break;
                case 4: bb = '-'; break;
            }
        } else {
            switch (ck) {
                case 0: bb = 'a'; break;
                case 1: bb = 'c'; break;
                case 2: bb = 'g'; break;
                case 3: bb = 't'; break;
                case 4: bb = '-'; break;
            }
        }
        // Note: On bad input, bb will keep previous value, possibly unitialized.

        score0 = g_best_aln_col->score;
        i = g_best_aln_col->best_p_t_pos;
        if (i == -1 || index >= t_len * 2) break;
        j = g_best_aln_col->best_p_delta;
        ck = g_best_aln_col->best_p_q_base;
        g_best_aln_col = msa_array[i]->delta[j].base + ck;

        if (bb != '-') {
            cns_str[index] = bb;
            eqv[index] = (int) score0 - (int) g_best_aln_col->score;
            //printf("C %d %d %c %lf %d %d\n", i, index, bb, g_best_aln_col->score, coverage[i], eqv[index] );
            index ++;
        }
    }
    
    // reverse the sequence
    for (i = 0; i < index/2; i++) {
        cns_str[i] = cns_str[i] ^ cns_str[index-i-1];
        cns_str[index-i-1] = cns_str[i] ^ cns_str[index-i-1];
        cns_str[i] = cns_str[i] ^ cns_str[index-i-1];
        eqv[i] = eqv[i] ^ eqv[index-i-1];
        eqv[index-i-1] = eqv[i] ^ eqv[index-i-1];
        eqv[i] = eqv[i] ^ eqv[index-i-1];
    }

    cns_str[index] = 0;
    //printf("%s\n", cns_str);
#ifndef STATIC_ALLOCATE
    for (i = 0; i < t_len; i++) {
        free_delta_group(msa_array[i]);
        free(msa_array[i]);
    }
    
    free(msa_array);
#endif

#ifdef STATIC_ALLOCATE
    clean_msa_working_space(msa_array, t_len+1);
#endif
    
    free(coverage);
    free(local_nbase);
    return consensus;
}

//const unsigned int K = 8;
