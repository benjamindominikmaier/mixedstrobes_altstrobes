//
//  index.cpp
//  cpp_course
//
//  Created by Kristoffer Sahlin on 4/21/21.
//

#include "index.hpp"
#include <iostream>
#include <math.h>       /* pow */
#include <bitset>
#include <vector>
#include <tuple>
#include <bits/stdc++.h>



/**********************************************************
 *
 * hash kmer into uint64
 *
 * *******************************************************/
// copy from http://www.cse.yorku.ca/~oz/hash.html:

uint64_t hash(std::string kmer)
{
    unsigned long hash = 5381;
    int c;
    for (std::string::size_type i=0; i< kmer.length(); i++) {
        c = kmer[i];
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
    }
    return hash;
}

/**********************************************************
 *
 * hash kmer into uint64
 *
 * *******************************************************/
// copy from minimap2:sketch.c :
static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}//hash64


static unsigned char seq_nt4_table[256] = {
        0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
}; //seq_nt4_table

static inline uint64_t kmer_to_uint64(std::string &kmer, uint64_t kmask)
{
    uint64_t bkmer = 0;

    for (char i : kmer) {
        int c = seq_nt4_table[(uint8_t)i];
        bkmer = (bkmer << 2 | c) & kmask;

    }
    return bkmer;
}


mers_vector construct_flat_vector_three_pos(three_pos_index &tmp_index, uint64_t &unique_elements){
    mers_vector flat_vector;
    for (auto &it : tmp_index)  {
        for (auto &t : it.second) // it.second is the vector of k-mers, t is a tuple
        {
            flat_vector.push_back(t);
        }
    }
    //    flat_array sort
    std::sort(flat_vector.begin(), flat_vector.end());

    uint64_t prev_k;
    std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> t = flat_vector[0];
    prev_k = std::get<0>(t);
    uint64_t curr_k;
    unique_elements = 1;
    for ( auto &t : flat_vector ) {
//        std::cout << t << std::endl;
        curr_k = std::get<0>(t);
        if (curr_k != prev_k){
            unique_elements ++;
        }
        prev_k = curr_k;
    }
    std::cout << "Unique Elements: " << unique_elements << std::endl; 
    return flat_vector;
}


void index_vector_three_pos(mers_vector  &flat_vector, kmer_lookup &mers_index){
//    robin_hood::unordered_map< uint64_t, std::tuple<uint64_t, unsigned int >> mers_index;
    uint64_t offset = 0;
    uint64_t prev_offset = 0;
    unsigned int count = 0;

    uint64_t prev_k;
    std::tuple<uint64_t, unsigned int, unsigned int, unsigned short int, unsigned short int> t = flat_vector[0];
    prev_k = std::get<0>(t);
    uint64_t curr_k;

    for ( auto &t : flat_vector ) {
//        std::cout << t << std::endl;
        curr_k = std::get<0>(t);
        if (curr_k == prev_k){
            count ++;
        }
        else {
            std::tuple<uint64_t, unsigned int> s(prev_offset, count);
            mers_index[prev_k] = s;
            count = 1;
            prev_k = curr_k;
            prev_offset = offset;
        }
        offset ++;
    }

    // last k-mer
    std::tuple<uint64_t, unsigned int> s(prev_offset, count);
    mers_index[curr_k] = s;

//    return mers_index;
}




// initialize queue and current minimum and position
static inline void initialize_window(std::vector<uint64_t> &string_hashes, std::deque <uint64_t> &q, uint64_t &q_min_val, int &q_min_pos, int w_min, int w_max, int k){
    robin_hood::hash<std::string> robin_hash;
    for (int i = w_min; i < w_max; i++) {
//        auto strobe = seq.substr(i, k);
//        uint64_t bstrobe = kmer_to_uint64(strobe, kmask);
//        q.push_back(bstrobe);
        uint64_t hv = string_hashes[i];
        q.push_back(hv);
        if (hv < q_min_val) {
            q_min_val = hv;
            q_min_pos = i;
        }
    }
}

// update queue and current minimum and position
static inline void update_window(std::deque <uint64_t> &q, uint64_t &q_min_val, int &q_min_pos, uint64_t new_strobe_hashval, int w_min, int w_max, int i, int seq_length ){
    uint64_t popped_val;
    popped_val = q.front();
    q.pop_front();
    q.push_back(new_strobe_hashval);
//    std::cout << "LOL " << popped_val << " " << q_min_val  << std::endl;
    if (popped_val == q_min_val){ // we popped the minimum value, find new brute force
//        std::cout << "OK "  << std::endl;
        q_min_val = UINT64_MAX;
//        q_min_pos = -1;
//        q_min_pos = i+w_min;
        q_min_pos = ((i+w_min + seq_length) - abs(i+w_min - seq_length))/2;
//        std::cout << "OK " << q_min_val  << std::endl;
//        std::cout << "OK " << q_min_pos  << std::endl;
        for (int j = 0; j <= q.size()-1; j++) {
//            std::cout << "(" << q[j] << " " << j << " " << i + w_min << "), ";
            if (q[j] < q_min_val) {
                q_min_val = q[j];
                q_min_pos = i + w_min + j + 1;
            }
        }
//        std::cout << "Changed: " << q_min_pos << " " << q_min_val << std::endl;
    }
    else if ( new_strobe_hashval < q_min_val ) { // the new value added to queue is the new minimum
        q_min_val = new_strobe_hashval;
        q_min_pos = i + w_max;
    }
}



static inline void make_string_to_hashvalues2(std::string &seq, std::vector<uint64_t> &string_hashes, std::vector<unsigned int> &pos_to_seq_choord, int k, uint64_t kmask) {
    robin_hood::hash<uint64_t> robin_hash;
//    std::vector<std::tuple<uint64_t, unsigned int, unsigned int> > kmers;
    unsigned int hash_count = 0;
    int l;
    int i;
    uint64_t x = 0;
    for (int i = l = 0; i < seq.length(); i++) {
        int c = seq_nt4_table[(uint8_t) seq[i]];
        if (c < 4) { // not an "N" base
            x = (x << 2 | c) & kmask;                  // forward strand
            if (++l >= k) { // we find a k-mer
//                uint64_t hash_k = x;
//                uint64_t hash_k = robin_hash(x);
                uint64_t hash_k = hash64(x, kmask);
//                std::cout << hash_k << " " << x << std::endl;
                string_hashes.push_back(hash_k);
                pos_to_seq_choord.push_back( i - k + 1);
//                pos_to_seq_choord[hash_count] = i - k + 1;
                hash_count ++;
            }
        } else {
            l = 0, x = 0; // if there is an "N", restart
        }
    }
}


uint64_t make_string_to_hashvalue(const std::string &seq, int k, int start) {
    uint64_t kmask=(1ULL<<2*k) - 1;
    uint64_t x = 0;
    for (int i = start; i < start+k; i++) {
        int c = seq_nt4_table[(uint8_t) seq[i]];
        if (c < 4) { // not an "N" base
            x = (x << 2 | c) & kmask;
            }
        else {
            return 0;
        }
    }
    uint64_t hash_k = hash64(x, kmask);
    return hash_k;
}



static inline void get_next_strobe(std::vector<uint64_t> &string_hashes, uint64_t strobe_hashval, unsigned int &strobe_pos_next, uint64_t &strobe_hashval_next,  unsigned int w_start, unsigned int w_end, uint64_t q){
    uint64_t min_val = UINT64_MAX;
//    unsigned int min_pos;
//    min_pos = -1;
    for (auto i = w_start; i <= w_end; i++) {
//        uint64_t res = (strobe_hashval + string_hashes[i]) & q ;
        uint64_t res = (strobe_hashval ^ string_hashes[i]) ;
        if (res < min_val){
            min_val = res;
            strobe_pos_next = i;
            strobe_hashval_next = string_hashes[i];
        }
    }
}

mers_vector seq_to_kmers(int k, std::string &seq, unsigned int ref_index)
{
    mers_vector kmers;
    int l;
    int i;
    uint64_t mask=(1ULL<<2*k) - 1;
    uint64_t x = 0;
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    int cnt = 0;
    for (int i = l = 0; i <= seq.length()-1; i++) {
        int c = seq_nt4_table[(uint8_t)seq[i]];
        if (c < 4) { // not an "N" base
            x = (x << 2 | c) & mask;                  // forward strand
            if (++l >= k) { // we find a k-mer
                uint64_t hash_k = x;
                std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_k, ref_index, i-k+1, i-k+1, i-k+1);
                kmers.push_back(s);
                cnt ++;
//                if ((cnt % 10000000) == 0 ){
//                    std::cout << cnt << " kmers created." << std::endl;
//                }
            }
        }
        else {
            l = 0, x = 0; // if there is an "N", restart
        }

    }
    return  kmers;
}

mers_vector seq_to_randstrobes2(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index)
{
    mers_vector randstrobes2;

    if (seq.length() < w_max) {
        return randstrobes2;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    uint64_t kmask=(1ULL<<2*k) - 1;
    uint64_t q = pow (2, 16) - 1;
//    std::bitset<64> x(q);
//    std::cout << x << '\n';
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    std::vector<uint64_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
    make_string_to_hashvalues2(seq, string_hashes, pos_to_seq_choord, k, kmask);
    unsigned int seq_length = string_hashes.size();

    if (string_hashes.size() == 0) {
        return randstrobes2;
    }
//    std::cout << seq << std::endl;

    // create the randstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {

//        if ((i % 1000000) == 0 ){
//            std::cout << i << " strobemers created." << std::endl;
//        }
        unsigned int strobe_pos_next;
        uint64_t strobe_hashval_next;

        if (i + w_max <= seq_length - 1){
            unsigned int w_start = i+w_min;
            unsigned int w_end = i+w_max;
            uint64_t strobe_hash;
            strobe_hash = string_hashes[i];
            get_next_strobe(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
        }
        else if ((i + w_min + 1 < seq_length) && (seq_length <= i + w_max) ){
            unsigned int w_start = i+w_min;
            unsigned int w_end = seq_length -1;
            uint64_t strobe_hash;
            strobe_hash = string_hashes[i];
            get_next_strobe(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
        }
        else{
            return randstrobes2;
        }

//        uint64_t hash_randstrobe2 = (string_hashes[i] << k) ^ strobe_hashval_next;
        uint64_t hash_randstrobe2 = (string_hashes[i]/2) + (strobe_hashval_next/3);

        unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
        unsigned int seq_pos_strobe2 = pos_to_seq_choord[strobe_pos_next];
        std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_randstrobe2, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe2);
        randstrobes2.push_back(s);


//        auto strobe1 = seq.substr(i, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe_pos_next - (i+k), ' ') << std::string(k, 'X') << std::endl;

    }
    return randstrobes2;
}


mers_vector seq_to_randstrobes3(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index)
{
    mers_vector randstrobes3;

    if (seq.length() < 2*w_max) {
        return randstrobes3;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    uint64_t kmask=(1ULL<<2*k) - 1;
    uint64_t q = pow (2, 16) - 1;
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    std::vector<uint64_t> string_hashes;
//    std::vector<uint64_t> pos_to_seq_choord;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
    make_string_to_hashvalues2(seq, string_hashes, pos_to_seq_choord, k, kmask);
    unsigned int seq_length = string_hashes.size();

    if (string_hashes.size() == 0) {
        return randstrobes3;
    }
//    std::cout << seq << std::endl;

    // create the randstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {

//        if ((i % 1000000) == 0 ){
//            std::cout << i << " randstrobes created." << std::endl;
//        }
        uint64_t strobe_hash;
        strobe_hash = string_hashes[i];

        unsigned int strobe_pos_next1;
        uint64_t strobe_hashval_next1;
        unsigned int strobe_pos_next2;
        uint64_t strobe_hashval_next2;

        if (i + 2*w_max <= seq_length - 1){
            unsigned int w1_start = i+w_min;
            unsigned int w1_end = i+w_max;
            get_next_strobe(string_hashes, strobe_hash, strobe_pos_next1, strobe_hashval_next1, w1_start, w1_end, q);

            unsigned int w2_start = i+w_max + w_min;
            unsigned int w2_end = i+2*w_max;
            uint64_t conditional_next = (strobe_hash ^ strobe_hashval_next1);
            get_next_strobe(string_hashes, conditional_next, strobe_pos_next2, strobe_hashval_next2, w2_start, w2_end, q);
        }
        else if ((i + 2*w_min + 1 < seq_length) && (seq_length <= i + 2*w_max) ){
//            unsigned int w_start = i+w_min;
//            unsigned int w_end = seq_length -1;
//            uint64_t strobe_hash;
//            strobe_hash = string_hashes[i];
//            get_next_strobe(string_hashes, strobe_hash, strobe_pos_next1, strobe_hashval_next1, w_start, w_end, q);

            int overshot;
            overshot = i + 2*w_max - seq_length;
            unsigned int w1_start = i+w_min;
            unsigned int w1_end = i+w_max - overshot/2;
            get_next_strobe(string_hashes, strobe_hash, strobe_pos_next1, strobe_hashval_next1, w1_start, w1_end, q);

            unsigned int w2_start = i+w_max - overshot/2 + w_min;
            unsigned int w2_end = i+2*w_max - overshot;
            uint64_t conditional_next = (strobe_hash ^ strobe_hashval_next1);
            get_next_strobe(string_hashes, conditional_next, strobe_pos_next2, strobe_hashval_next2, w2_start, w2_end, q);
        }
        else{
            return randstrobes3;
        }

//        uint64_t hash_randstrobe3 = (((strobe_hash << k) ^ strobe_hashval_next1) << k) ^ strobe_hashval_next2;
        uint64_t hash_randstrobe3 = (strobe_hash/3) + (strobe_hashval_next1/4) + (strobe_hashval_next2/5);

        unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
        unsigned int seq_pos_strobe2 =  pos_to_seq_choord[strobe_pos_next1]; //seq_pos_strobe1 + (strobe_pos_next1 - i); //
        unsigned int seq_pos_strobe3 =  pos_to_seq_choord[strobe_pos_next2]; //seq_pos_strobe1 + (strobe_pos_next2 - i); //
//        std::cout << i << " " << strobe_pos_next1 << " " << strobe_pos_next2 << " " << seq_pos_strobe1 << " " << seq_pos_strobe2 << " " << seq_pos_strobe3 << " " << pos_to_seq_choord.size() << std::endl;

        // TODO: Take care of corner case (tmep if statement below. Some values in end of string produce a cororidnate of 0 for the last strobe. Probably an off-by-one error in the calculation of the strobe coord in the last strobe window
        if (strobe_pos_next2 ==  seq_length){
//            std::cout << "OMGGGGGGG " << i << " " << seq_pos_strobe1 << " " << seq_pos_strobe2 << " " << seq_pos_strobe3 << std::endl;
            seq_pos_strobe3 = seq_length-1;
        }
        std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_randstrobe3, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe3);
        randstrobes3.push_back(s);
        // std::cout << "randstrobes " <<  hash_randstrobe3 << " " << ref_index << " " << seq_pos_strobe1 << " " << seq_pos_strobe2 << " " << seq_pos_strobe3 << " " << std::endl;


//        auto strobe1 = seq.substr(i, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe_pos_next1 - (i+k), ' ') << std::string(k, 'X') << std::string(strobe_pos_next2 - strobe_pos_next1 - k, ' ') << std::string(k, 'X') << std::endl;
//        std::cout << i << " " << strobe_pos_next1 << " " << strobe_pos_next2 << " " << seq_length << std::endl;
    }
    return randstrobes3;
}


mers_vector seq_to_randstrobes(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index)
{
    mers_vector randstrobes;

    if (seq.length() < (n-1)*w_max) {
        return randstrobes;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    uint64_t kmask=(1ULL<<2*k) - 1;
    uint64_t q = pow (2, 16) - 1;
//    std::bitset<64> x(q);
//    std::cout << x << '\n';
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    std::vector<uint64_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
    make_string_to_hashvalues2(seq, string_hashes, pos_to_seq_choord, k, kmask);
    unsigned int seq_length = string_hashes.size();

    if (string_hashes.size() == 0) {
        return randstrobes;
    }
//    std::cout << seq << std::endl;

    // create the randstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {

//        if ((i % 1000000) == 0 ){
//            std::cout << i << " strobemers created." << std::endl;
//        }
        uint64_t strobe_hash;
        strobe_hash = string_hashes[i];

        unsigned int strobe_pos_next;
        uint64_t strobe_hashval_next;
        std::vector<unsigned int> strobe_positions(n, i);
        std::vector<uint64_t> strobe_hashvalues(n, strobe_hash);

        if (i + (n-1)*w_max <= seq_length - 1){
            for (unsigned int order = 1; order < n; order++){
              unsigned int w_start = i+w_min+(order-1)*w_max;
              unsigned int w_end = i+order*w_max;
              get_next_strobe(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
              strobe_positions[order] = strobe_pos_next;
              strobe_hashvalues[order] = strobe_hashval_next;
              strobe_hash = (strobe_hash ^ strobe_hashval_next);
            }
        }
        else if ((i + (n-1)*w_min + 1 < seq_length) && (seq_length <= i + (n-1)*w_max) ){
            int overshot;
            overshot = i + (n-1)*w_max - seq_length;

            for (unsigned int order = 1; order < n; order++){
              unsigned int w_start = i+w_min+(order-1)*w_max-(order-1)*overshot/2;
              unsigned int w_end = i+order*w_max-order*overshot/2;
              get_next_strobe(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
              strobe_positions[order] = strobe_pos_next;
              strobe_hashvalues[order] = strobe_hashval_next;
              strobe_hash = (strobe_hash ^ strobe_hashval_next);
            }
        }
        else{
            return randstrobes;
        }

//        uint64_t hash_randstrobe2 = (string_hashes[i] << k) ^ strobe_hashval_next;
        uint64_t hash_randstrobe = 0;
        for (unsigned int j = 0; j < n; j++) {
          hash_randstrobe += strobe_hashvalues[j]/(n+j);
        }

        // TODO: Take care of corner case (tmep if statement below. Some values in end of string produce a coordinate of 0 for the last strobe. Probably an off-by-one error in the calculation of the strobe coord in the last strobe window
        if (strobe_positions[n-1] ==  seq_length){
            strobe_positions[n-1] = seq_length-1;
        }

        std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_randstrobe, ref_index, strobe_positions[0], strobe_positions[1], strobe_positions[n-1]);
        randstrobes.push_back(s);
        // std::cout << "randstrobes " <<  hash_randstrobe << " " << ref_index << " " << strobe_positions[0] << " " << strobe_positions[1] << " " << strobe_positions[n-1] << " " << std::endl;


//        auto strobe1 = seq.substr(i, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe_pos_next - (i+k), ' ') << std::string(k, 'X') << std::endl;

    }
    return randstrobes;
}


mers_vector seq_to_mixedrandstrobes2(int n, int k, int w_min, int w_max, int strobe_fraction, std::string &seq, unsigned int ref_index)
{
    mers_vector mixedrandstrobes2;

    if (seq.length() < w_max) {
        return mixedrandstrobes2;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    uint64_t kmask=(1ULL<<2*k) - 1;
    uint64_t q = pow (2, 16) - 1;
//    std::bitset<64> x(q);
//    std::cout << x << '\n';
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    std::vector<uint64_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
    make_string_to_hashvalues2(seq, string_hashes, pos_to_seq_choord, k, kmask);
    unsigned int seq_length = string_hashes.size();

    if (string_hashes.size() == 0) {
        return mixedrandstrobes2;
    }
//    std::cout << seq << std::endl;

    // create the randstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {

//        if ((i % 1000000) == 0 ){
//            std::cout << i << " strobemers created." << std::endl;
//        }
        uint64_t strobe_hash;
        strobe_hash = string_hashes[i];

        if (strobe_hash % 100 < strobe_fraction){

            unsigned int strobe_pos_next;
            uint64_t strobe_hashval_next;

            if (i + w_max <= seq_length - 1){
                unsigned int w_start = i+w_min;
                unsigned int w_end = i+w_max;

                get_next_strobe(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
            }
            else if ((i + w_min + 1 < seq_length) && (seq_length <= i + w_max) ){
                unsigned int w_start = i+w_min;
                unsigned int w_end = seq_length -1;

                get_next_strobe(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
            }
            else{
                return mixedrandstrobes2;
            }

    //        uint64_t hash_randstrobe2 = (string_hashes[i] << k) ^ strobe_hashval_next;
            uint64_t hash_randstrobe2 = (string_hashes[i]/2) + (strobe_hashval_next/3);

            unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
            unsigned int seq_pos_strobe2 = pos_to_seq_choord[strobe_pos_next];
            std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_randstrobe2, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe2);
            mixedrandstrobes2.push_back(s);


    //        auto strobe1 = seq.substr(i, k);
    //        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe_pos_next - (i+k), ' ') << std::string(k, 'X') << std::endl;

        }
        else {
            // strobe_hash = ((strobe_hash << k) ^ string_hashes[i+k]);
            strobe_hash = ((strobe_hash/2) + (string_hashes[i+k]/3));

            std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (strobe_hash, ref_index, i, i+k, i+k);
            mixedrandstrobes2.push_back(s);
        }
    }
    return mixedrandstrobes2;
}


mers_vector seq_to_mixedrandstrobes3(int n, int k, int w_min, int w_max, int strobe_fraction, std::string &seq, unsigned int ref_index)
{
    mers_vector mixedrandstrobes3;

    if (seq.length() < 2*w_max) {
        return mixedrandstrobes3;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    uint64_t kmask=(1ULL<<2*k) - 1;
    uint64_t q = pow (2, 16) - 1;
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    std::vector<uint64_t> string_hashes;
//    std::vector<uint64_t> pos_to_seq_choord;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
    make_string_to_hashvalues2(seq, string_hashes, pos_to_seq_choord, k, kmask);
    unsigned int seq_length = string_hashes.size();

    if (string_hashes.size() == 0) {
        return mixedrandstrobes3;
    }
//    std::cout << seq << std::endl;

    // create the randstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {

//        if ((i % 1000000) == 0 ){
//            std::cout << i << " randstrobes created." << std::endl;
//        }
        uint64_t strobe_hash;
        strobe_hash = string_hashes[i];

        if (strobe_hash % 100 < strobe_fraction){

            unsigned int strobe_pos_next1;
            uint64_t strobe_hashval_next1;
            unsigned int strobe_pos_next2;
            uint64_t strobe_hashval_next2;

            if (i + 2*w_max <= seq_length - 1){
                unsigned int w1_start = i+w_min;
                unsigned int w1_end = i+w_max;
                get_next_strobe(string_hashes, strobe_hash, strobe_pos_next1, strobe_hashval_next1, w1_start, w1_end, q);

                unsigned int w2_start = i+w_max + w_min;
                unsigned int w2_end = i+2*w_max;
                uint64_t conditional_next = (strobe_hash ^ strobe_hashval_next1);
                get_next_strobe(string_hashes, conditional_next, strobe_pos_next2, strobe_hashval_next2, w2_start, w2_end, q);
            }
            else if ((i + 2*w_min + 1 < seq_length) && (seq_length <= i + 2*w_max) ){
    //            unsigned int w_start = i+w_min;
    //            unsigned int w_end = seq_length -1;
    //            uint64_t strobe_hash;
    //            strobe_hash = string_hashes[i];
    //            get_next_strobe(string_hashes, strobe_hash, strobe_pos_next1, strobe_hashval_next1, w_start, w_end, q);

                int overshot;
                overshot = i + 2*w_max - seq_length;
                unsigned int w1_start = i+w_min;
                unsigned int w1_end = i+w_max - overshot/2;
                get_next_strobe(string_hashes, strobe_hash, strobe_pos_next1, strobe_hashval_next1, w1_start, w1_end, q);

                unsigned int w2_start = i+w_max - overshot/2 + w_min;
                unsigned int w2_end = i+2*w_max - overshot;
                uint64_t conditional_next = (strobe_hash ^ strobe_hashval_next1);
                get_next_strobe(string_hashes, conditional_next, strobe_pos_next2, strobe_hashval_next2, w2_start, w2_end, q);
            }
            else{
                return mixedrandstrobes3;
            }

    //        uint64_t hash_randstrobe3 = (((strobe_hash << k) ^ strobe_hashval_next1) << k) ^ strobe_hashval_next2;
            uint64_t hash_randstrobe3 = (strobe_hash/3) + (strobe_hashval_next1/4) + (strobe_hashval_next2/5);

            unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
            unsigned int seq_pos_strobe2 =  pos_to_seq_choord[strobe_pos_next1]; //seq_pos_strobe1 + (strobe_pos_next1 - i); //
            unsigned int seq_pos_strobe3 =  pos_to_seq_choord[strobe_pos_next2]; //seq_pos_strobe1 + (strobe_pos_next2 - i); //
    //        std::cout << i << " " << strobe_pos_next1 << " " << strobe_pos_next2 << " " << seq_pos_strobe1 << " " << seq_pos_strobe2 << " " << seq_pos_strobe3 << " " << pos_to_seq_choord.size() << std::endl;

            // TODO: Take care of corner case (tmep if statement below. Some values in end of string produce a cororidnate of 0 for the last strobe. Probably an off-by-one error in the calculation of the strobe coord in the last strobe window
            if (strobe_pos_next2 ==  seq_length){
    //            std::cout << "OMGGGGGGG " << i << " " << seq_pos_strobe1 << " " << seq_pos_strobe2 << " " << seq_pos_strobe3 << std::endl;
                seq_pos_strobe3 = seq_length-1;
            }
            std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_randstrobe3, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe3);
            mixedrandstrobes3.push_back(s);
            // std::cout << "randstrobes " <<  hash_randstrobe3 << " " << ref_index << " " << seq_pos_strobe1 << " " << seq_pos_strobe2 << " " << seq_pos_strobe3 << " " << std::endl;


    //        auto strobe1 = seq.substr(i, k);
    //        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe_pos_next1 - (i+k), ' ') << std::string(k, 'X') << std::string(strobe_pos_next2 - strobe_pos_next1 - k, ' ') << std::string(k, 'X') << std::endl;
    //        std::cout << i << " " << strobe_pos_next1 << " " << strobe_pos_next2 << " " << seq_length << std::endl;
        }
        else {
            // strobe_hash = ((strobe_hash << k) ^ string_hashes[i+k]);
            // strobe_hash = ((strobe_hash << k) ^ string_hashes[i+2*k]);

            strobe_hash = ((strobe_hash/3) + (string_hashes[i+k]/4) + (string_hashes[i+2*k]/5));

            std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (strobe_hash, ref_index, i, i+k, i+2*k);
            mixedrandstrobes3.push_back(s);
        }
    }
    return mixedrandstrobes3;
}


mers_vector seq_to_mixedrandstrobes(int n, int k, int w_min, int w_max, int strobe_fraction, std::string &seq, unsigned int ref_index)
{
    mers_vector mixedrandstrobes;

    if (seq.length() < (n-1)*w_max) {
        return mixedrandstrobes;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    uint64_t kmask=(1ULL<<2*k) - 1;
    uint64_t q = pow (2, 16) - 1;
//    std::bitset<64> x(q);
//    std::cout << x << '\n';
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    std::vector<uint64_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
    make_string_to_hashvalues2(seq, string_hashes, pos_to_seq_choord, k, kmask);
    unsigned int seq_length = string_hashes.size();

    if (string_hashes.size() == 0) {
        return mixedrandstrobes;
    }
//    std::cout << seq << std::endl;

    // create the randstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {

//        if ((i % 1000000) == 0 ){
//            std::cout << i << " strobemers created." << std::endl;
//        }
        uint64_t strobe_hash;
        strobe_hash = string_hashes[i];

        if (strobe_hash % 100 < strobe_fraction){
            unsigned int strobe_pos_next;
            uint64_t strobe_hashval_next;
            std::vector<unsigned int> strobe_positions(n, i);
            std::vector<uint64_t> strobe_hashvalues(n, strobe_hash);

            if (i + (n-1)*w_max <= seq_length - 1){
                for (unsigned int order = 1; order < n; order++){
                unsigned int w_start = i+w_min+(order-1)*w_max;
                unsigned int w_end = i+order*w_max;
                get_next_strobe(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
                strobe_positions[order] = strobe_pos_next;
                strobe_hashvalues[order] = strobe_hashval_next;
                strobe_hash = (strobe_hash ^ strobe_hashval_next);
                }
            }
            else if ((i + (n-1)*w_min + 1 < seq_length) && (seq_length <= i + (n-1)*w_max) ){
                int overshot;
                overshot = i + (n-1)*w_max - seq_length;

                for (unsigned int order = 1; order < n; order++){
                unsigned int w_start = i+w_min+(order-1)*w_max-(order-1)*overshot/2;
                unsigned int w_end = i+order*w_max-order*overshot/2;
                get_next_strobe(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
                strobe_positions[order] = strobe_pos_next;
                strobe_hashvalues[order] = strobe_hashval_next;
                strobe_hash = (strobe_hash ^ strobe_hashval_next);
                }
            }
            else{
                return mixedrandstrobes;
            }

    //        uint64_t hash_randstrobe2 = (string_hashes[i] << k) ^ strobe_hashval_next;
            uint64_t hash_randstrobe = 0;
            for (unsigned int j = 0; j < n; j++) {
            hash_randstrobe += strobe_hashvalues[j]/(n+j);
            }

            // TODO: Take care of corner case (tmep if statement below. Some values in end of string produce a cororidnate of 0 for the last strobe. Probably an off-by-one error in the calculation of the strobe coord in the last strobe window
            if (strobe_positions[n-1] ==  seq_length){
    //            std::cout << "OMGGGGGGG " << i << " " << seq_pos_strobe1 << " " << seq_pos_strobe2 << " " << seq_pos_strobe3 << std::endl;
                strobe_positions[n-1] = seq_length-1;
            }

            std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_randstrobe, ref_index, strobe_positions[0], strobe_positions[1], strobe_positions[n-1]);
            mixedrandstrobes.push_back(s);
            // std::cout << "randstrobes " <<  hash_randstrobe << " " << ref_index << " " << strobe_positions[0] << " " << strobe_positions[1] << " " << strobe_positions[n-1] << " " << std::endl;


    //        auto strobe1 = seq.substr(i, k);
    //        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe_pos_next - (i+k), ' ') << std::string(k, 'X') << std::endl;
        }
        else {
            for (unsigned int j = 0; j < n; j++) {
                // strobe_hash = ((strobe_hash << k) ^ string_hashes[i+j*k]);
                strobe_hash += string_hashes[i+j*k]/(n+j);
            }
            std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (strobe_hash, ref_index, i, i+k, i+k*(n-1));
            mixedrandstrobes.push_back(s);
        }
    }
    return mixedrandstrobes;
}


mers_vector seq_to_altstrobes2(int n, int k, int k2, int w_min, int w_max, std::string &seq, unsigned int ref_index)
{
    mers_vector altstrobes2;

    if (seq.length() < w_max) {
        return altstrobes2;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    uint64_t kmask=(1ULL<<2*k) - 1;
    uint64_t kmask2=(1ULL<<2*k2) - 1;
    uint64_t q = pow (2, 16) - 1;
//    std::bitset<64> x(q);
//    std::cout << x << '\n';
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    std::vector<uint64_t> string_hashes;
    std::vector<uint64_t> string_hashes2;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
    make_string_to_hashvalues2(seq, string_hashes, pos_to_seq_choord, k, kmask);
    make_string_to_hashvalues2(seq, string_hashes2, pos_to_seq_choord, k2, kmask2);
    unsigned int seq_length = string_hashes.size();

    if (string_hashes.size() == 0) {
        return altstrobes2;
    }
//    std::cout << seq << std::endl;

    int w_min1 = w_min - (k2-k)/2; 
    int w_min2 = w_min + (k2-k)/2;
    int w_max1 = w_max - (k2-k)/2;
    int w_max2 = w_max + (k2-k)/2;

    // create the altstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {

//        if ((i % 1000000) == 0 ){
//            std::cout << i << " strobemers created." << std::endl;
//        }
        unsigned int strobe_pos_next;
        uint64_t strobe_hashval_next;

        uint64_t strobe_hash;
        strobe_hash = string_hashes[i];

        if (strobe_hash % 2 == 0){ // k - k2

            if (i + w_max <= seq_length - 1){
                unsigned int w_start = i+w_min2;
                unsigned int w_end = i+w_max2;
                get_next_strobe(string_hashes2, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
            }
            else if ((i + w_min2 + 1 < seq_length) && (seq_length <= i + w_max2) ){
                unsigned int w_start = i+w_min2;
                unsigned int w_end = seq_length -1;
                get_next_strobe(string_hashes2, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
            }
            else{
                return altstrobes2;
            }

//          uint64_t hash_altstrobe2 = (string_hashes[i] << k) ^ strobe_hashval_next;
            uint64_t hash_altstrobe2 = (string_hashes[i]/2) + (strobe_hashval_next/3);

            unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
            unsigned int seq_pos_strobe2 = pos_to_seq_choord[strobe_pos_next];
            unsigned int seq_pos_strobe3 = pos_to_seq_choord[strobe_pos_next + k];
            std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_altstrobe2, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe3);
            altstrobes2.push_back(s);
        }

        else { // k2 - k

            if (i + w_max <= seq_length - 1){
                unsigned int w_start = i+w_min1;
                unsigned int w_end = i+w_max1;
                uint64_t strobe_hash;
                strobe_hash = string_hashes2[i];
                get_next_strobe(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
            }
            else if ((i + w_min + 1 < seq_length) && (seq_length <= i + w_max) ){
                unsigned int w_start = i+w_min1;
                unsigned int w_end = seq_length -1;
                uint64_t strobe_hash;
                strobe_hash = string_hashes2[i];
                get_next_strobe(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
            }
            else{
                return altstrobes2;
            }

//          uint64_t hash_altstrobe2 = (string_hashes[i] << k) ^ strobe_hashval_next;
            uint64_t hash_altstrobe2 = (string_hashes2[i]/2) + (strobe_hashval_next/3);

            unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
            unsigned int seq_pos_strobe2 = pos_to_seq_choord[i+k];
            unsigned int seq_pos_strobe3 = pos_to_seq_choord[strobe_pos_next];
            std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_altstrobe2, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe3);
            altstrobes2.push_back(s);
        }


//        auto strobe1 = seq.substr(i, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe_pos_next - (i+k), ' ') << std::string(k, 'X') << std::endl;

    }
    return altstrobes2;
}


mers_vector seq_to_altstrobes(int n, int k, int k2, int w_min, int w_max, std::string &seq, unsigned int ref_index)
{
    mers_vector altstrobes;

    if (seq.length() < (n-1) * w_max) {
        return altstrobes;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    uint64_t kmask=(1ULL<<2*k) - 1;
    uint64_t kmask2=(1ULL<<2*k2) - 1;
    uint64_t q = pow (2, 16) - 1;
//    std::bitset<64> x(q);
//    std::cout << x << '\n';
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    std::vector<uint64_t> string_hashes;
    std::vector<uint64_t> string_hashes2;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
    make_string_to_hashvalues2(seq, string_hashes, pos_to_seq_choord, k, kmask);
    make_string_to_hashvalues2(seq, string_hashes2, pos_to_seq_choord, k2, kmask2);
    unsigned int seq_length = string_hashes.size();

    if (string_hashes.size() == 0) {
        return altstrobes;
    }
//    std::cout << seq << std::endl;

    int w_min1 = w_min - (k2-k)/2; 
    int w_min2 = w_min + (k2-k)/2;
    int w_max1 = w_max - (k2-k)/2;
    int w_max2 = w_max + (k2-k)/2;

    bool seed_k1 = true;
    unsigned int required_length = 0;
    unsigned int required_length_1 = n*w_max-w_max1;
    unsigned int required_length_2 = n*w_max-w_max2;

    // create the altstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {

//        if ((i % 1000000) == 0 ){
//            std::cout << i << " strobemers created." << std::endl;
//        }
        uint64_t strobe_hash;
        strobe_hash = string_hashes[i];

        unsigned int strobe_pos_next;
        uint64_t strobe_hashval_next;

        std::vector<unsigned int> strobe_positions(n, i);
        std::vector<uint64_t> strobe_hashvalues(n, strobe_hash);

        if (strobe_hash % 2 == 0){ // k - k2
            bool seed_k1 = false;
            unsigned int required_length = required_length_1;
        }
        else{ // k2 - k
            bool seed_k1 = true;
            strobe_hash = string_hashes2[i];
            strobe_hashvalues[0] = strobe_hash;
            unsigned int required_length = required_length_2;
        }

        if (i + required_length <= seq_length - 1){
            unsigned int current_offset = i;
            for (unsigned int order = 1; order < n; order++){
                if (seed_k1 == true){ // seed strobe of length k
                    unsigned int w_start = w_min2+current_offset;
                    current_offset += w_max1;
                    unsigned int w_end = std::min(current_offset, seq_length);
                    get_next_strobe(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
                }
                else {// seed strobe of length k2
                    unsigned int w_start = w_min1+current_offset;
                    current_offset += w_max2;
                    unsigned int w_end = std::min(current_offset, seq_length);
                    get_next_strobe(string_hashes2, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
                }

                strobe_positions[order] = strobe_pos_next;
                strobe_hashvalues[order] = strobe_hashval_next;
                strobe_hash = (strobe_hash ^ strobe_hashval_next);
                seed_k1 = !seed_k1;
            }
            
        }
        else if ((i + (n-1)*w_min + 1 < seq_length) && (seq_length <= i + (n-1)*w_max) ){
            int overshot = i + (n-1)*w_max - seq_length;

            for (unsigned int order = 1; order < n; order++){
                if (seed_k1 == true){ // seed strobe of length k
                    unsigned int w_start = i+w_min1+(order-1)*w_max-(order-1)*overshot/2;
                    unsigned int w_end = i+order*w_max-order*overshot/2;
                    get_next_strobe(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
                }
                else { // seed strobe of length k2
                    unsigned int w_start = i+w_min2+(order-1)*w_max-(order-1)*overshot/2;
                    unsigned int w_end = i+order*w_max-order*overshot/2;
                    get_next_strobe(string_hashes2, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);   
                }
                strobe_positions[order] = strobe_pos_next;
                strobe_hashvalues[order] = strobe_hashval_next;
                strobe_hash = (strobe_hash ^ strobe_hashval_next);
                seed_k1 = !seed_k1;
            }
        }
        
        else{
            return altstrobes;
        }

        uint64_t hash_altstrobe = 0;
        for (unsigned int j = 0; j < n; j++) {
          hash_altstrobe += strobe_hashvalues[j]/(n+j);
        }

        // TODO: Take care of corner case (tmep if statement below. Some values in end of string produce a coordinate of 0 for the last strobe. Probably an off-by-one error in the calculation of the strobe coord in the last strobe window
        if (strobe_positions[n-1] ==  seq_length){
            strobe_positions[n-1] = seq_length-1;
        }

        std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_altstrobe, ref_index, strobe_positions[0], strobe_positions[1], strobe_positions[n-1]);
        altstrobes.push_back(s);
    }

    return altstrobes;
}


mers_vector seq_to_minstrobes2(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index)
{
    mers_vector minstrobes2;
    if (seq.length() < w_max) {
        return minstrobes2;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    uint64_t kmask=(1ULL<<2*k) - 1;
    std::vector<uint64_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
    make_string_to_hashvalues2(seq, string_hashes, pos_to_seq_choord, k, kmask);
    unsigned int seq_length = string_hashes.size();
//    unsigned int seq_length = seq.length();

    // initialize the deque
    std::deque <uint64_t> q;
    uint64_t q_min_val = UINT64_MAX;
    int q_min_pos = -1;
    initialize_window(string_hashes, q, q_min_val, q_min_pos, w_min, w_max, k);

//    std::cout << seq << std::endl;

    // create the minstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {
//        auto strobe1 = seq.substr(i, k);
//        uint64_t bstrobe = kmer_to_uint64(strobe1, kmask);

        uint64_t bstrobe = string_hashes[i];
//        uint64_t hash_minstrobe2 = (bstrobe << k) ^ q_min_val;
        uint64_t hash_minstrobe2 = (bstrobe/2) + (q_min_val/3);

        unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
        unsigned int seq_pos_strobe2 = pos_to_seq_choord[q_min_pos];
        std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_minstrobe2, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe2);
        minstrobes2.push_back(s);

        // update queue and current minimum and position
        if (i + w_max <= seq_length){
//            auto new_strobe = seq.substr(i + w_max, k);
//            uint64_t new_bstrobe = kmer_to_uint64(new_strobe, kmask);
//            uint64_t new_strobe_hashval = hash64(new_bstrobe, kmask);
            uint64_t new_strobe = string_hashes[i + w_max];

            update_window(q, q_min_val, q_min_pos, new_strobe, w_min, w_max, i, seq_length );
        }
        else if ((i + w_min + 1 < seq_length) && (seq_length < i + w_max) ){
//            uint64_t new_strobe_hashval =  UINT64_MAX;
            uint64_t new_strobe = UINT64_MAX;
            update_window(q, q_min_val, q_min_pos, new_strobe, w_min, w_max, i, seq_length );
        }
        else{
            return minstrobes2;
        }

//        auto strobe1 = seq.substr(i, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(q_min_pos - (i+k), ' ') << std::string(k, 'X') << std::endl;

    }
    return minstrobes2;
}



mers_vector seq_to_mixedminstrobes2(int n, int k, int w_min, int w_max, int strobe_fraction, std::string &seq, unsigned int ref_index)
{
    mers_vector mixedminstrobes2;
    if (seq.length() < w_max) {
        return mixedminstrobes2;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    uint64_t kmask=(1ULL<<2*k) - 1;
    std::vector<uint64_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
    make_string_to_hashvalues2(seq, string_hashes, pos_to_seq_choord, k, kmask);
    unsigned int seq_length = string_hashes.size();
//    unsigned int seq_length = seq.length();

    // initialize the deque
    std::deque <uint64_t> q;
    uint64_t q_min_val = UINT64_MAX;
    int q_min_pos = -1;
    initialize_window(string_hashes, q, q_min_val, q_min_pos, w_min, w_max, k);

//    std::cout << seq << std::endl;

    // create the minstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {
//        auto strobe1 = seq.substr(i, k);
//        uint64_t bstrobe = kmer_to_uint64(strobe1, kmask);

        uint64_t bstrobe = string_hashes[i];

        if (bstrobe % 100 < strobe_fraction){
    //        uint64_t hash_minstrobe2 = (bstrobe << k) ^ q_min_val;
            uint64_t hash_minstrobe2 = (bstrobe/2) + (q_min_val/3);

            unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
            unsigned int seq_pos_strobe2 = pos_to_seq_choord[q_min_pos];
            std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_minstrobe2, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe2);
            mixedminstrobes2.push_back(s);

            // update queue and current minimum and position
            if (i + w_max <= seq_length){
    //            auto new_strobe = seq.substr(i + w_max, k);
    //            uint64_t new_bstrobe = kmer_to_uint64(new_strobe, kmask);
    //            uint64_t new_strobe_hashval = hash64(new_bstrobe, kmask);
                uint64_t new_strobe = string_hashes[i + w_max];

                update_window(q, q_min_val, q_min_pos, new_strobe, w_min, w_max, i, seq_length );
            }
            else if ((i + w_min + 1 < seq_length) && (seq_length < i + w_max) ){
    //            uint64_t new_strobe_hashval =  UINT64_MAX;
                uint64_t new_strobe = UINT64_MAX;
                update_window(q, q_min_val, q_min_pos, new_strobe, w_min, w_max, i, seq_length );
            }
            else{
                return mixedminstrobes2;
            }
        }
        else {
                // bstrobe = ((bstrobe << k) ^ string_hashes[i+k]);
                bstrobe = ((bstrobe/2) + (string_hashes[i+k]/3));
                std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (bstrobe, ref_index, i, i+k, i+k);
                mixedminstrobes2.push_back(s);
            }
    }

//        auto strobe1 = seq.substr(i, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(q_min_pos - (i+k), ' ') << std::string(k, 'X') << std::endl;

    return mixedminstrobes2;
}


mers_vector seq_to_hybridstrobes2(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index)
{
    mers_vector hybridstrobes2;

    // TODO: This if-statement leads to downstream bug:
    //  Process finished with exit code 139 (interrupted by signal 11: SIGSEGV). Fix this
    if (seq.length() < w_max) {
        return hybridstrobes2;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    uint64_t kmask=(1ULL<<2*k) - 1;
    std::vector<uint64_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
    make_string_to_hashvalues2(seq, string_hashes, pos_to_seq_choord, k, kmask);
//    std::cout << seq.length() << " " << string_hashes.size() << " " << k << std::endl;
    unsigned int seq_length = string_hashes.size();

    int x = (w_max - w_min) /3;

    // initialize deque1
    std::deque <uint64_t> q1;
    uint64_t q1_min_val = UINT64_MAX;
    int q1_min_pos = -1;
    initialize_window(string_hashes, q1, q1_min_val, q1_min_pos, w_min, w_min+ x, k);

    // initialize deque2
    std::deque <uint64_t> q2;
    uint64_t q2_min_val = UINT64_MAX;
    int q2_min_pos = -1;
    initialize_window(string_hashes, q2, q2_min_val, q2_min_pos, w_min+x, w_min+2*x, k);

    // initialize deque3
    std::deque <uint64_t> q3;
    uint64_t q3_min_val = UINT64_MAX;
    int q3_min_pos = -1;
    initialize_window(string_hashes, q3, q3_min_val, q3_min_pos, w_min+2*x, w_max, k);

//    std::cout << seq << std::endl;

    // create the hybridstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {
//        auto strobe1 = seq.substr(i, k);
//        uint64_t bstrobe = kmer_to_uint64(strobe1, kmask);
//        uint64_t strobe_hashval = hash64(bstrobe, kmask);

        uint64_t bstrobe = string_hashes[i];

        int strobe2_pos;
        uint64_t strobe2_val;
        unsigned int r =  bstrobe % 3;
        if (r == 0){
            strobe2_val = q1_min_val;
            strobe2_pos = q1_min_pos;
        }
        else if (r == 1){
            strobe2_val = q2_min_val;
            strobe2_pos = q2_min_pos;
        }
        else{
            strobe2_val = q3_min_val;
            strobe2_pos = q3_min_pos;
        }

        uint64_t hash_hybridstrobe2;
//        hash_hybridstrobe2 = (bstrobe << k) ^ strobe2_val;
        hash_hybridstrobe2 = (bstrobe/2) + (strobe2_val/3);

        unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
        unsigned int seq_pos_strobe2 = pos_to_seq_choord[strobe2_pos];
        std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_hybridstrobe2, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe2);
        hybridstrobes2.push_back(s);


        // update queue and current minimum and position
        if (i + w_max < seq_length){
//            auto new_strobe1 = seq.substr(i + w_min+x, k);
//            uint64_t new_bstrobe1 = kmer_to_uint64(new_strobe1, kmask);
//            uint64_t new_strobe_hashval1 = hash64(new_bstrobe1, kmask);
            uint64_t new_strobe_hashval1 = string_hashes[i + w_min+x];
            update_window(q1, q1_min_val, q1_min_pos, new_strobe_hashval1, w_min, w_min+x, i, seq_length);

//            auto new_strobe2 = seq.substr(i + w_min+2*x, k);
//            uint64_t new_bstrobe2 = kmer_to_uint64(new_strobe2, kmask);
//            uint64_t new_strobe_hashval2 = hash64(new_bstrobe2, kmask);
            uint64_t new_strobe_hashval2 = string_hashes[i + w_min+2*x];
            update_window(q2, q2_min_val, q2_min_pos, new_strobe_hashval2, w_min+x, w_min+2*x, i, seq_length);

//            auto new_strobe3 = seq.substr(i + w_max, k);
//            uint64_t new_bstrobe3 = kmer_to_uint64(new_strobe3, kmask);
//            uint64_t new_strobe_hashval3 = hash64(new_bstrobe3, kmask);
            uint64_t new_strobe_hashval3 = string_hashes[i + w_max];
            update_window(q3, q3_min_val, q3_min_pos, new_strobe_hashval3, w_min+2*x, w_max, i, seq_length);
        }
        else if ((i + w_min + 1 < seq_length) && (seq_length <= i + w_min+x) ){
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(q1, q1_min_val, q1_min_pos, new_strobe_hashval, w_min, w_min+x, i, seq_length);
            update_window(q2, q2_min_val, q2_min_pos, new_strobe_hashval, w_min+x, w_min+2*x, i, seq_length);
            update_window(q3, q3_min_val, q3_min_pos, new_strobe_hashval, w_min+2*x, w_max, i, seq_length);
        }
        else if ((i + w_min+x < seq_length) && (seq_length <= i + w_min+2*x)  ){
            uint64_t new_strobe_hashval1 = string_hashes[i + w_min+x];
            update_window(q1, q1_min_val, q1_min_pos, new_strobe_hashval1, w_min, w_min+x, i, seq_length);
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(q2, q2_min_val, q2_min_pos, new_strobe_hashval, w_min+x, w_min+2*x, i, seq_length);
            update_window(q3, q3_min_val, q3_min_pos, new_strobe_hashval, w_min+2*x, w_max, i, seq_length);
        }
        else if ((i + w_min+2*x < seq_length) && (seq_length <= i + w_max) ){
            uint64_t new_strobe_hashval1 = string_hashes[i + w_min+x];
            update_window(q1, q1_min_val, q1_min_pos, new_strobe_hashval1, w_min, w_min+x, i, seq_length);
            uint64_t new_strobe_hashval2 = string_hashes[i + w_min+2*x];
            update_window(q2, q2_min_val, q2_min_pos, new_strobe_hashval2, w_min+x, w_min+2*x, i, seq_length);
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(q3, q3_min_val, q3_min_pos, new_strobe_hashval, w_min+2*x, w_max, i, seq_length);
        }
        else{
            return hybridstrobes2;
        }

//        auto strobe1 = seq.substr(i, k);
//        auto strobe2 = seq.substr(strobe2_pos, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe2_pos - (i+k)-1, ' ') << std::string(k, 'X') << std::endl;
//        std::cout << i << " " << strobe2_pos << " " << seq_length << std::endl;
//        std::cout << i << ": " << strobe1 << strobe2 << std::endl;


    }
    return hybridstrobes2;
}




mers_vector seq_to_hybridstrobes3(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index)
{
    mers_vector hybridstrobes3;

    // TODO: This if-statement leads to downstream bug:
    //  Process finished with exit code 139 (interrupted by signal 11: SIGSEGV). Fix this
    if (seq.length() < 2*w_max) {
        return hybridstrobes3;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    uint64_t kmask=(1ULL<<2*k) - 1;
    std::vector<uint64_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
    make_string_to_hashvalues2(seq, string_hashes, pos_to_seq_choord, k, kmask);
//    std::cout << seq.length() << " " << string_hashes.size() << " " << k << std::endl;
    unsigned int seq_length = string_hashes.size();

    int x = (w_max - w_min) /3;

    // initialize deque1 window1
    std::deque <uint64_t> q1;
    uint64_t q1_min_val = UINT64_MAX;
    int q1_min_pos = -1;
    initialize_window(string_hashes, q1, q1_min_val, q1_min_pos, w_min, w_min+ x, k);

    // initialize deque2 window1
    std::deque <uint64_t> q2;
    uint64_t q2_min_val = UINT64_MAX;
    int q2_min_pos = -1;
    initialize_window(string_hashes, q2, q2_min_val, q2_min_pos, w_min+x, w_min+2*x, k);

    // initialize deque3 window1
    std::deque <uint64_t> q3;
    uint64_t q3_min_val = UINT64_MAX;
    int q3_min_pos = -1;
    initialize_window(string_hashes, q3, q3_min_val, q3_min_pos, w_min+2*x, w_max, k);

    // initialize deque1 window2
    std::deque <uint64_t> q4;
    uint64_t q4_min_val = UINT64_MAX;
    int q4_min_pos = -1;
    initialize_window(string_hashes, q4, q4_min_val, q4_min_pos, w_max+ w_min, w_max + w_min+ x, k);

    // initialize deque2 window2
    std::deque <uint64_t> q5;
    uint64_t q5_min_val = UINT64_MAX;
    int q5_min_pos = -1;
    initialize_window(string_hashes, q5, q5_min_val, q5_min_pos, w_max+ w_min+x, w_max+ w_min+2*x, k);

    // initialize deque3 window2
    std::deque <uint64_t> q6;
    uint64_t q6_min_val = UINT64_MAX;
    int q6_min_pos = -1;
    initialize_window(string_hashes, q6, q6_min_val, q6_min_pos, w_max+ w_min+2*x, 2*w_max, k);

//    std::cout << seq << std::endl;

    // create the hybridstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {
//        auto strobe1 = seq.substr(i, k);
//        uint64_t bstrobe = kmer_to_uint64(strobe1, kmask);
//        uint64_t strobe_hashval = hash64(bstrobe, kmask);

        uint64_t bstrobe = string_hashes[i];

        int strobe2_pos;
        uint64_t strobe2_val;
        unsigned int r =  bstrobe % 3;
        if (r == 0){
            strobe2_val = q1_min_val;
            strobe2_pos = q1_min_pos;
        }
        else if (r == 1){
            strobe2_val = q2_min_val;
            strobe2_pos = q2_min_pos;
        }
        else{
            strobe2_val = q3_min_val;
            strobe2_pos = q3_min_pos;
        }

        uint64_t hash_hybridstrobe2;
//        hash_hybridstrobe2 = (bstrobe << k) ^ strobe2_val;
        hash_hybridstrobe2 = (bstrobe/3) + (strobe2_val/4);

        int strobe3_pos;
        uint64_t strobe3_val;
        unsigned int r2 =  hash_hybridstrobe2 % 3;
        if (r2 == 0){
            strobe3_val = q4_min_val;
            strobe3_pos = q4_min_pos;
        }
        else if (r2 == 1){
            strobe3_val = q5_min_val;
            strobe3_pos = q5_min_pos;
        }
        else{
            strobe3_val = q6_min_val;
            strobe3_pos = q6_min_pos;
        }

        uint64_t hash_hybridstrobe3;
//        hash_hybridstrobe2 = (bstrobe << k) ^ strobe2_val;
        hash_hybridstrobe3 = hash_hybridstrobe2 + (strobe3_val/5);

        unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
        unsigned int seq_pos_strobe2 = pos_to_seq_choord[strobe2_pos];
        unsigned int seq_pos_strobe3 = pos_to_seq_choord[strobe3_pos];
        std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_hybridstrobe3, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe3);
        hybridstrobes3.push_back(s);


        // update queue and current minimum and position
        uint64_t new_strobe_hashval1 = string_hashes[i + w_min+x];
        update_window(q1, q1_min_val, q1_min_pos, new_strobe_hashval1, w_min, w_min+x, i, seq_length);

        uint64_t new_strobe_hashval2 = string_hashes[i + w_min+2*x];
        update_window(q2, q2_min_val, q2_min_pos, new_strobe_hashval2, w_min+x, w_min+2*x, i, seq_length);

        uint64_t new_strobe_hashval3 = string_hashes[i + w_max];
        update_window(q3, q3_min_val, q3_min_pos, new_strobe_hashval3, w_min+2*x, w_max, i, seq_length);

        if (i + 2*w_max < seq_length){
            uint64_t new_strobe_hashval4 = string_hashes[i + w_max + w_min+x];
            update_window(q4, q4_min_val, q4_min_pos, new_strobe_hashval4, w_max + w_min, w_max + w_min+x, i, seq_length);

            uint64_t new_strobe_hashval5 = string_hashes[i + w_max + w_min+2*x];
            update_window(q5, q5_min_val, q5_min_pos, new_strobe_hashval5, w_max + w_min+x, w_max + w_min+2*x, i, seq_length);

            uint64_t new_strobe_hashval6 = string_hashes[i + 2*w_max];
            update_window(q6, q6_min_val, q6_min_pos, new_strobe_hashval6, w_max + w_min+2*x, 2*w_max, i, seq_length);
        }
        else if ((i + w_max + w_min + 1 < seq_length) && (seq_length <= i + w_max + w_min+x) ){
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(q4, q4_min_val, q4_min_pos, new_strobe_hashval, w_max + w_min, w_max + w_min+x, i, seq_length);
            update_window(q5, q5_min_val, q5_min_pos, new_strobe_hashval, w_max + w_min+x, w_max + w_min+2*x, i, seq_length);
            update_window(q6, q6_min_val, q6_min_pos, new_strobe_hashval, w_max + w_min+2*x, 2*w_max, i, seq_length);
        }
        else if ((i + w_max + w_min+x < seq_length) && (seq_length <= i + w_max + w_min+2*x)  ){
            uint64_t new_strobe_hashval4 = string_hashes[i + w_max + w_min+x];
            update_window(q4, q4_min_val, q4_min_pos, new_strobe_hashval4, w_max + w_min, w_max + w_min+x, i, seq_length);
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(q5, q5_min_val, q5_min_pos, new_strobe_hashval, w_max + w_min+x, w_max + w_min+2*x, i, seq_length);
            update_window(q6, q6_min_val, q6_min_pos, new_strobe_hashval, w_max + w_min+2*x, 2*w_max, i, seq_length);
        }
        else if ((i + w_max + w_min+2*x < seq_length) && (seq_length <= i + 2*w_max) ){
            uint64_t new_strobe_hashval4 = string_hashes[i + w_max + w_min+x];
            update_window(q4, q4_min_val, q4_min_pos, new_strobe_hashval4, w_max + w_min, w_max + w_min+x, i, seq_length);
            uint64_t new_strobe_hashval5 = string_hashes[i + w_max + w_min+2*x];
            update_window(q5, q5_min_val, q5_min_pos, new_strobe_hashval5, w_max + w_min+x, w_max + w_min+2*x, i, seq_length);
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(q6, q6_min_val, q6_min_pos, new_strobe_hashval, w_max + w_min+2*x, 2*w_max, i, seq_length);
        }
        else{
            return hybridstrobes3;
        }

//        auto strobe1 = seq.substr(i, k);
//        auto strobe2 = seq.substr(strobe2_pos, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe2_pos - (i+k) - 1, ' ') << std::string(k, 'X') << std::string(strobe3_pos - strobe2_pos - k, ' ') << std::string(k, 'X') << std::endl;

//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe2_pos - (i+k)-1, ' ') << std::string(k, 'X') << std::endl;
//        std::cout << i << " " << strobe2_pos << " " << seq_length << std::endl;
//        std::cout << i << ": " << strobe1 << strobe2 << std::endl;


    }
    return hybridstrobes3;
}




mers_vector seq_to_mixedhybridstrobes2(int n, int k, int w_min, int w_max, int strobe_fraction, std::string &seq, unsigned int ref_index)
{
    mers_vector mixedhybridstrobes2;

    // TODO: This if-statement leads to downstream bug:
    //  Process finished with exit code 139 (interrupted by signal 11: SIGSEGV). Fix this
    if (seq.length() < w_max) {
        return mixedhybridstrobes2;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    uint64_t kmask=(1ULL<<2*k) - 1;
    std::vector<uint64_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
    make_string_to_hashvalues2(seq, string_hashes, pos_to_seq_choord, k, kmask);
//    std::cout << seq.length() << " " << string_hashes.size() << " " << k << std::endl;
    unsigned int seq_length = string_hashes.size();

    int x = (w_max - w_min) /3;

    // initialize deque1
    std::deque <uint64_t> q1;
    uint64_t q1_min_val = UINT64_MAX;
    int q1_min_pos = -1;
    initialize_window(string_hashes, q1, q1_min_val, q1_min_pos, w_min, w_min+ x, k);

    // initialize deque2
    std::deque <uint64_t> q2;
    uint64_t q2_min_val = UINT64_MAX;
    int q2_min_pos = -1;
    initialize_window(string_hashes, q2, q2_min_val, q2_min_pos, w_min+x, w_min+2*x, k);

    // initialize deque3
    std::deque <uint64_t> q3;
    uint64_t q3_min_val = UINT64_MAX;
    int q3_min_pos = -1;
    initialize_window(string_hashes, q3, q3_min_val, q3_min_pos, w_min+2*x, w_max, k);

//    std::cout << seq << std::endl;

    // create the hybridstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {
//        auto strobe1 = seq.substr(i, k);
//        uint64_t bstrobe = kmer_to_uint64(strobe1, kmask);
//        uint64_t strobe_hashval = hash64(bstrobe, kmask);

        uint64_t bstrobe = string_hashes[i];

        if (bstrobe % 100 < strobe_fraction){

            int strobe2_pos;
            uint64_t strobe2_val;
            unsigned int r =  bstrobe % 3;
            if (r == 0){
                strobe2_val = q1_min_val;
                strobe2_pos = q1_min_pos;
            }
            else if (r == 1){
                strobe2_val = q2_min_val;
                strobe2_pos = q2_min_pos;
            }
            else{
                strobe2_val = q3_min_val;
                strobe2_pos = q3_min_pos;
            }

            uint64_t hash_hybridstrobe2;
    //        hash_hybridstrobe2 = (bstrobe << k) ^ strobe2_val;
            hash_hybridstrobe2 = (bstrobe/2) + (strobe2_val/3);

            unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
            unsigned int seq_pos_strobe2 = pos_to_seq_choord[strobe2_pos];
            std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_hybridstrobe2, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe2);
            mixedhybridstrobes2.push_back(s);
        }

        else {
            // bstrobe = ((bstrobe << k) ^ string_hashes[i+k]);
            bstrobe = ((bstrobe/2) + (string_hashes[i+k]/3));
            std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (bstrobe, ref_index, i, i+k, i+k);
            mixedhybridstrobes2.push_back(s);
        }


        // update queue and current minimum and position
        if (i + w_max < seq_length){
//            auto new_strobe1 = seq.substr(i + w_min+x, k);
//            uint64_t new_bstrobe1 = kmer_to_uint64(new_strobe1, kmask);
//            uint64_t new_strobe_hashval1 = hash64(new_bstrobe1, kmask);
            uint64_t new_strobe_hashval1 = string_hashes[i + w_min+x];
            update_window(q1, q1_min_val, q1_min_pos, new_strobe_hashval1, w_min, w_min+x, i, seq_length);

//            auto new_strobe2 = seq.substr(i + w_min+2*x, k);
//            uint64_t new_bstrobe2 = kmer_to_uint64(new_strobe2, kmask);
//            uint64_t new_strobe_hashval2 = hash64(new_bstrobe2, kmask);
            uint64_t new_strobe_hashval2 = string_hashes[i + w_min+2*x];
            update_window(q2, q2_min_val, q2_min_pos, new_strobe_hashval2, w_min+x, w_min+2*x, i, seq_length);

//            auto new_strobe3 = seq.substr(i + w_max, k);
//            uint64_t new_bstrobe3 = kmer_to_uint64(new_strobe3, kmask);
//            uint64_t new_strobe_hashval3 = hash64(new_bstrobe3, kmask);
            uint64_t new_strobe_hashval3 = string_hashes[i + w_max];
            update_window(q3, q3_min_val, q3_min_pos, new_strobe_hashval3, w_min+2*x, w_max, i, seq_length);
        }
        else if ((i + w_min + 1 < seq_length) && (seq_length <= i + w_min+x) ){
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(q1, q1_min_val, q1_min_pos, new_strobe_hashval, w_min, w_min+x, i, seq_length);
            update_window(q2, q2_min_val, q2_min_pos, new_strobe_hashval, w_min+x, w_min+2*x, i, seq_length);
            update_window(q3, q3_min_val, q3_min_pos, new_strobe_hashval, w_min+2*x, w_max, i, seq_length);
        }
        else if ((i + w_min+x < seq_length) && (seq_length <= i + w_min+2*x)  ){
            uint64_t new_strobe_hashval1 = string_hashes[i + w_min+x];
            update_window(q1, q1_min_val, q1_min_pos, new_strobe_hashval1, w_min, w_min+x, i, seq_length);
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(q2, q2_min_val, q2_min_pos, new_strobe_hashval, w_min+x, w_min+2*x, i, seq_length);
            update_window(q3, q3_min_val, q3_min_pos, new_strobe_hashval, w_min+2*x, w_max, i, seq_length);
        }
        else if ((i + w_min+2*x < seq_length) && (seq_length <= i + w_max) ){
            uint64_t new_strobe_hashval1 = string_hashes[i + w_min+x];
            update_window(q1, q1_min_val, q1_min_pos, new_strobe_hashval1, w_min, w_min+x, i, seq_length);
            uint64_t new_strobe_hashval2 = string_hashes[i + w_min+2*x];
            update_window(q2, q2_min_val, q2_min_pos, new_strobe_hashval2, w_min+x, w_min+2*x, i, seq_length);
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(q3, q3_min_val, q3_min_pos, new_strobe_hashval, w_min+2*x, w_max, i, seq_length);
        }
        else{
            return mixedhybridstrobes2;
        }

//        auto strobe1 = seq.substr(i, k);
//        auto strobe2 = seq.substr(strobe2_pos, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe2_pos - (i+k)-1, ' ') << std::string(k, 'X') << std::endl;
//        std::cout << i << " " << strobe2_pos << " " << seq_length << std::endl;
//        std::cout << i << ": " << strobe1 << strobe2 << std::endl;


    }
    return mixedhybridstrobes2;
}




mers_vector seq_to_mixedhybridstrobes3(int n, int k, int w_min, int w_max, int strobe_fraction, std::string &seq, unsigned int ref_index)
{
    mers_vector mixedhybridstrobes3;

    // TODO: This if-statement leads to downstream bug:
    //  Process finished with exit code 139 (interrupted by signal 11: SIGSEGV). Fix this
    if (seq.length() < 2*w_max) {
        return mixedhybridstrobes3;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    uint64_t kmask=(1ULL<<2*k) - 1;
    std::vector<uint64_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
    make_string_to_hashvalues2(seq, string_hashes, pos_to_seq_choord, k, kmask);
//    std::cout << seq.length() << " " << string_hashes.size() << " " << k << std::endl;
    unsigned int seq_length = string_hashes.size();

    int x = (w_max - w_min) /3;

    // initialize deque1 window1
    std::deque <uint64_t> q1;
    uint64_t q1_min_val = UINT64_MAX;
    int q1_min_pos = -1;
    initialize_window(string_hashes, q1, q1_min_val, q1_min_pos, w_min, w_min+ x, k);

    // initialize deque2 window1
    std::deque <uint64_t> q2;
    uint64_t q2_min_val = UINT64_MAX;
    int q2_min_pos = -1;
    initialize_window(string_hashes, q2, q2_min_val, q2_min_pos, w_min+x, w_min+2*x, k);

    // initialize deque3 window1
    std::deque <uint64_t> q3;
    uint64_t q3_min_val = UINT64_MAX;
    int q3_min_pos = -1;
    initialize_window(string_hashes, q3, q3_min_val, q3_min_pos, w_min+2*x, w_max, k);

    // initialize deque1 window2
    std::deque <uint64_t> q4;
    uint64_t q4_min_val = UINT64_MAX;
    int q4_min_pos = -1;
    initialize_window(string_hashes, q4, q4_min_val, q4_min_pos, w_max+ w_min, w_max + w_min+ x, k);

    // initialize deque2 window2
    std::deque <uint64_t> q5;
    uint64_t q5_min_val = UINT64_MAX;
    int q5_min_pos = -1;
    initialize_window(string_hashes, q5, q5_min_val, q5_min_pos, w_max+ w_min+x, w_max+ w_min+2*x, k);

    // initialize deque3 window2
    std::deque <uint64_t> q6;
    uint64_t q6_min_val = UINT64_MAX;
    int q6_min_pos = -1;
    initialize_window(string_hashes, q6, q6_min_val, q6_min_pos, w_max+ w_min+2*x, 2*w_max, k);

//    std::cout << seq << std::endl;

    // create the hybridstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {
//        auto strobe1 = seq.substr(i, k);
//        uint64_t bstrobe = kmer_to_uint64(strobe1, kmask);
//        uint64_t strobe_hashval = hash64(bstrobe, kmask);

        uint64_t bstrobe = string_hashes[i];

        if (bstrobe % 100 < strobe_fraction){

            int strobe2_pos;
            uint64_t strobe2_val;
            unsigned int r =  bstrobe % 3;
            if (r == 0){
                strobe2_val = q1_min_val;
                strobe2_pos = q1_min_pos;
            }
            else if (r == 1){
                strobe2_val = q2_min_val;
                strobe2_pos = q2_min_pos;
            }
            else{
                strobe2_val = q3_min_val;
                strobe2_pos = q3_min_pos;
            }

            uint64_t hash_hybridstrobe2;
    //        hash_hybridstrobe2 = (bstrobe << k) ^ strobe2_val;
            hash_hybridstrobe2 = (bstrobe/3) + (strobe2_val/4);

            int strobe3_pos;
            uint64_t strobe3_val;
            unsigned int r2 =  hash_hybridstrobe2 % 3;
            if (r2 == 0){
                strobe3_val = q4_min_val;
                strobe3_pos = q4_min_pos;
            }
            else if (r2 == 1){
                strobe3_val = q5_min_val;
                strobe3_pos = q5_min_pos;
            }
            else{
                strobe3_val = q6_min_val;
                strobe3_pos = q6_min_pos;
            }

            uint64_t hash_hybridstrobe3;
    //        hash_hybridstrobe2 = (bstrobe << k) ^ strobe2_val;
            hash_hybridstrobe3 = hash_hybridstrobe2 + (strobe3_val/5);

            unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
            unsigned int seq_pos_strobe2 = pos_to_seq_choord[strobe2_pos];
            unsigned int seq_pos_strobe3 = pos_to_seq_choord[strobe3_pos];
            std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_hybridstrobe3, ref_index, seq_pos_strobe1, seq_pos_strobe2, seq_pos_strobe3);
            mixedhybridstrobes3.push_back(s);
        }
        else {
            // bstrobe = ((bstrobe << k) ^ string_hashes[i+k]);
            // bstrobe = ((bstrobe << k) ^ string_hashes[i+2*k]);
            bstrobe = ((bstrobe/3) + (string_hashes[i+k]/4) + (string_hashes[i+2*k]/5));
            std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (bstrobe, ref_index, i, i+k, i+2*k);
            mixedhybridstrobes3.push_back(s);
        }


        // update queue and current minimum and position
        uint64_t new_strobe_hashval1 = string_hashes[i + w_min+x];
        update_window(q1, q1_min_val, q1_min_pos, new_strobe_hashval1, w_min, w_min+x, i, seq_length);

        uint64_t new_strobe_hashval2 = string_hashes[i + w_min+2*x];
        update_window(q2, q2_min_val, q2_min_pos, new_strobe_hashval2, w_min+x, w_min+2*x, i, seq_length);

        uint64_t new_strobe_hashval3 = string_hashes[i + w_max];
        update_window(q3, q3_min_val, q3_min_pos, new_strobe_hashval3, w_min+2*x, w_max, i, seq_length);

        if (i + 2*w_max < seq_length){
            uint64_t new_strobe_hashval4 = string_hashes[i + w_max + w_min+x];
            update_window(q4, q4_min_val, q4_min_pos, new_strobe_hashval4, w_max + w_min, w_max + w_min+x, i, seq_length);

            uint64_t new_strobe_hashval5 = string_hashes[i + w_max + w_min+2*x];
            update_window(q5, q5_min_val, q5_min_pos, new_strobe_hashval5, w_max + w_min+x, w_max + w_min+2*x, i, seq_length);

            uint64_t new_strobe_hashval6 = string_hashes[i + 2*w_max];
            update_window(q6, q6_min_val, q6_min_pos, new_strobe_hashval6, w_max + w_min+2*x, 2*w_max, i, seq_length);
        }
        else if ((i + w_max + w_min + 1 < seq_length) && (seq_length <= i + w_max + w_min+x) ){
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(q4, q4_min_val, q4_min_pos, new_strobe_hashval, w_max + w_min, w_max + w_min+x, i, seq_length);
            update_window(q5, q5_min_val, q5_min_pos, new_strobe_hashval, w_max + w_min+x, w_max + w_min+2*x, i, seq_length);
            update_window(q6, q6_min_val, q6_min_pos, new_strobe_hashval, w_max + w_min+2*x, 2*w_max, i, seq_length);
        }
        else if ((i + w_max + w_min+x < seq_length) && (seq_length <= i + w_max + w_min+2*x)  ){
            uint64_t new_strobe_hashval4 = string_hashes[i + w_max + w_min+x];
            update_window(q4, q4_min_val, q4_min_pos, new_strobe_hashval4, w_max + w_min, w_max + w_min+x, i, seq_length);
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(q5, q5_min_val, q5_min_pos, new_strobe_hashval, w_max + w_min+x, w_max + w_min+2*x, i, seq_length);
            update_window(q6, q6_min_val, q6_min_pos, new_strobe_hashval, w_max + w_min+2*x, 2*w_max, i, seq_length);
        }
        else if ((i + w_max + w_min+2*x < seq_length) && (seq_length <= i + 2*w_max) ){
            uint64_t new_strobe_hashval4 = string_hashes[i + w_max + w_min+x];
            update_window(q4, q4_min_val, q4_min_pos, new_strobe_hashval4, w_max + w_min, w_max + w_min+x, i, seq_length);
            uint64_t new_strobe_hashval5 = string_hashes[i + w_max + w_min+2*x];
            update_window(q5, q5_min_val, q5_min_pos, new_strobe_hashval5, w_max + w_min+x, w_max + w_min+2*x, i, seq_length);
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(q6, q6_min_val, q6_min_pos, new_strobe_hashval, w_max + w_min+2*x, 2*w_max, i, seq_length);
        }
        else{
            return mixedhybridstrobes3;
        }

//        auto strobe1 = seq.substr(i, k);
//        auto strobe2 = seq.substr(strobe2_pos, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe2_pos - (i+k) - 1, ' ') << std::string(k, 'X') << std::string(strobe3_pos - strobe2_pos - k, ' ') << std::string(k, 'X') << std::endl;

//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe2_pos - (i+k)-1, ' ') << std::string(k, 'X') << std::endl;
//        std::cout << i << " " << strobe2_pos << " " << seq_length << std::endl;
//        std::cout << i << ": " << strobe1 << strobe2 << std::endl;


    }
    return mixedhybridstrobes3;
}






mers_vector seq_to_hybridstrobes(int n, int k, int w_min, int w_max, std::string &seq, unsigned int ref_index)
{
    mers_vector hybridstrobes;

    // TODO: This if-statement leads to downstream bug:
    //  Process finished with exit code 139 (interrupted by signal 11: SIGSEGV). Fix this
    if (seq.length() < (n-1)*w_max) {
        return hybridstrobes;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    uint64_t kmask=(1ULL<<2*k) - 1;
    std::vector<uint64_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
    make_string_to_hashvalues2(seq, string_hashes, pos_to_seq_choord, k, kmask);
//    std::cout << seq.length() << " " << string_hashes.size() << " " << k << std::endl;
    unsigned int seq_length = string_hashes.size();

    int x = (w_max - w_min) /3;
    std::vector<std::tuple <std::deque<uint64_t>, uint64_t, int> > windows;

    for (unsigned int order = 0; order < (n-1); order++){
      // initialize deque1
      std::deque <uint64_t> q1;
      uint64_t q1_min_val = UINT64_MAX;
      int q1_min_pos = -1;
      initialize_window(string_hashes, q1, q1_min_val, q1_min_pos, w_max*order+ w_min, w_max*order+ w_min+ x, k);
      windows.emplace_back(q1, q1_min_val, q1_min_pos);

      // initialize deque2
      std::deque <uint64_t> q2;
      uint64_t q2_min_val = UINT64_MAX;
      int q2_min_pos = -1;
      initialize_window(string_hashes, q2, q2_min_val, q2_min_pos, w_max*order+ w_min+x, w_max*order+ w_min+2*x, k);
      windows.emplace_back(q2, q2_min_val, q2_min_pos);

      // initialize deque3
      std::deque <uint64_t> q3;
      uint64_t q3_min_val = UINT64_MAX;
      int q3_min_pos = -1;
      initialize_window(string_hashes, q3, q3_min_val, q3_min_pos, w_max*order+ w_min+2*x, w_max*(order+1), k);
      windows.emplace_back(q3, q3_min_val, q3_min_pos);
    }

//    std::cout << seq << std::endl;

    // create the hybridstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {
//        auto strobe1 = seq.substr(i, k);
//        uint64_t bstrobe = kmer_to_uint64(strobe1, kmask);
//        uint64_t strobe_hashval = hash64(bstrobe, kmask);

        uint64_t bstrobe = string_hashes[i]/n;
        std::vector<unsigned int> strobe_positions(n, pos_to_seq_choord[i]);

        for (unsigned int strobe_num = 1; strobe_num < n; strobe_num++){
          uint64_t strobe_next_val;
          unsigned int r = bstrobe % 3;
          if (r == 0){
              bstrobe += std::get<1>(windows[strobe_num*3-3])/(n+strobe_num);
              strobe_positions[strobe_num] = pos_to_seq_choord[std::get<2>(windows[strobe_num*3-3])];
          }
          else if (r == 1){
            bstrobe += std::get<1>(windows[strobe_num*3-2])/(n+strobe_num);
            strobe_positions[strobe_num] = pos_to_seq_choord[std::get<2>(windows[strobe_num*3-2])];
          }
          else{
            bstrobe += std::get<1>(windows[strobe_num*3-1])/(n+strobe_num);
            strobe_positions[strobe_num] = pos_to_seq_choord[std::get<2>(windows[strobe_num*3-1])];
          }
        }

        std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (bstrobe, ref_index, strobe_positions[0], strobe_positions[1], strobe_positions[n-1]);
        hybridstrobes.push_back(s);


        // update queue and current minimum and position
        for (unsigned int order = 0; order < (n-2); order++){
            uint64_t new_strobe_hashval1 = string_hashes[i + order*w_max + w_min+x];
            update_window(std::get<0>(windows[order*3]), std::get<1>(windows[order*3]), std::get<2>(windows[order*3]), new_strobe_hashval1, order*w_max+w_min, order*w_max+w_min+x, i, seq_length);

            uint64_t new_strobe_hashval2 = string_hashes[i + order*w_max + w_min+2*x];
            update_window(std::get<0>(windows[order*3+1]), std::get<1>(windows[order*3+1]), std::get<2>(windows[order*3+1]), new_strobe_hashval2, order*w_max+w_min+x, order*w_max+w_min+2*x, i, seq_length);

            uint64_t new_strobe_hashval3 = string_hashes[i + order*w_max + w_max];
            update_window(std::get<0>(windows[order*3+2]), std::get<1>(windows[order*3+2]), std::get<2>(windows[order*3+2]), new_strobe_hashval3, order*w_max+w_min+2*x, order*w_max+w_max, i, seq_length);
        }

        if (i + 2*w_max < seq_length){
            uint64_t new_strobe_hashval4 = string_hashes[i + w_max + w_min+x];
            update_window(std::get<0>(windows[(n-2)*3]), std::get<1>(windows[(n-2)*3]), std::get<2>(windows[(n-2)*3]), new_strobe_hashval4, (n-2)*w_max + w_min, (n-2)*w_max + w_min+x, i, seq_length);

            uint64_t new_strobe_hashval5 = string_hashes[i + w_max + w_min+2*x];
            update_window(std::get<0>(windows[(n-2)*3+1]), std::get<1>(windows[(n-2)*3+1]), std::get<2>(windows[(n-2)*3+1]), new_strobe_hashval5, (n-2)*w_max + w_min+x, (n-2)*w_max + w_min+2*x, i, seq_length);

            uint64_t new_strobe_hashval6 = string_hashes[i + 2*w_max];
            update_window(std::get<0>(windows[(n-2)*3+2]), std::get<1>(windows[(n-2)*3+2]), std::get<2>(windows[(n-2)*3+2]), new_strobe_hashval6, (n-2)*w_max + w_min+2*x, (n-1)*w_max, i, seq_length);
        }
        else if ((i + w_max + w_min + 1 < seq_length) && (seq_length <= i + w_max + w_min+x) ){
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(std::get<0>(windows[(n-2)*3]), std::get<1>(windows[(n-2)*3]), std::get<2>(windows[(n-2)*3]), new_strobe_hashval, (n-2)*w_max + w_min, (n-2)*w_max + w_min+x, i, seq_length);
            update_window(std::get<0>(windows[(n-2)*3+1]), std::get<1>(windows[(n-2)*3+1]), std::get<2>(windows[(n-2)*3+1]), new_strobe_hashval, (n-2)*w_max + w_min+x, (n-2)*w_max + w_min+2*x, i, seq_length);
            update_window(std::get<0>(windows[(n-2)*3+2]), std::get<1>(windows[(n-2)*3+2]), std::get<2>(windows[(n-2)*3+2]), new_strobe_hashval, (n-2)*w_max + w_min+2*x, (n-1)*w_max, i, seq_length);
        }
        else if ((i + w_max + w_min+x < seq_length) && (seq_length <= i + w_max + w_min+2*x)  ){
            uint64_t new_strobe_hashval4 = string_hashes[i + w_max + w_min+x];
            update_window(std::get<0>(windows[(n-2)*3]), std::get<1>(windows[(n-2)*3]), std::get<2>(windows[(n-2)*3]), new_strobe_hashval4, (n-2)*w_max + w_min, (n-2)*w_max + w_min+x, i, seq_length);
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(std::get<0>(windows[(n-2)*3+1]), std::get<1>(windows[(n-2)*3+1]), std::get<2>(windows[(n-2)*3+1]), new_strobe_hashval, (n-2)*w_max + w_min+x, (n-2)*w_max + w_min+2*x, i, seq_length);
            update_window(std::get<0>(windows[(n-2)*3+2]), std::get<1>(windows[(n-2)*3+2]), std::get<2>(windows[(n-2)*3+2]), new_strobe_hashval, (n-2)*w_max + w_min+2*x, (n-1)*w_max, i, seq_length);
        }
        else if ((i + w_max + w_min+2*x < seq_length) && (seq_length <= i + 2*w_max) ){
            uint64_t new_strobe_hashval4 = string_hashes[i + w_max + w_min+x];
            update_window(std::get<0>(windows[(n-2)*3]), std::get<1>(windows[(n-2)*3]), std::get<2>(windows[(n-2)*3]), new_strobe_hashval4, (n-2)*w_max + w_min, (n-2)*w_max + w_min+x, i, seq_length);
            uint64_t new_strobe_hashval5 = string_hashes[i + w_max + w_min+2*x];
            update_window(std::get<0>(windows[(n-2)*3+1]), std::get<1>(windows[(n-2)*3+1]), std::get<2>(windows[(n-2)*3+1]), new_strobe_hashval5, (n-2)*w_max + w_min+x, (n-2)*w_max + w_min+2*x, i, seq_length);
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(std::get<0>(windows[(n-2)*3+2]), std::get<1>(windows[(n-2)*3+2]), std::get<2>(windows[(n-2)*3+2]), new_strobe_hashval, (n-2)*w_max + w_min+2*x, (n-1)*w_max, i, seq_length);
        }
        else{
            return hybridstrobes;
        }

//        auto strobe1 = seq.substr(i, k);
//        auto strobe2 = seq.substr(strobe2_pos, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe2_pos - (i+k) - 1, ' ') << std::string(k, 'X') << std::string(strobe3_pos - strobe2_pos - k, ' ') << std::string(k, 'X') << std::endl;

//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe2_pos - (i+k)-1, ' ') << std::string(k, 'X') << std::endl;
//        std::cout << i << " " << strobe2_pos << " " << seq_length << std::endl;
//        std::cout << i << ": " << strobe1 << strobe2 << std::endl;


    }
    return hybridstrobes;
}





mers_vector seq_to_mixedhybridstrobes(int n, int k, int w_min, int w_max, int strobe_fraction, std::string &seq, unsigned int ref_index)
{
    mers_vector mixedhybridstrobes;

    // TODO: This if-statement leads to downstream bug:
    //  Process finished with exit code 139 (interrupted by signal 11: SIGSEGV). Fix this
    if (seq.length() < (n-1)*w_max) {
        return mixedhybridstrobes;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    uint64_t kmask=(1ULL<<2*k) - 1;
    std::vector<uint64_t> string_hashes;
    std::vector<unsigned int> pos_to_seq_choord;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
    make_string_to_hashvalues2(seq, string_hashes, pos_to_seq_choord, k, kmask);
//    std::cout << seq.length() << " " << string_hashes.size() << " " << k << std::endl;
    unsigned int seq_length = string_hashes.size();

    int x = (w_max - w_min) /3;
    std::vector<std::tuple <std::deque<uint64_t>, uint64_t, int> > windows;

    for (unsigned int order = 0; order < (n-1); order++){
      // initialize deque1
      std::deque <uint64_t> q1;
      uint64_t q1_min_val = UINT64_MAX;
      int q1_min_pos = -1;
      initialize_window(string_hashes, q1, q1_min_val, q1_min_pos, w_max*order+ w_min, w_max*order+ w_min+ x, k);
      windows.emplace_back(q1, q1_min_val, q1_min_pos);

      // initialize deque2
      std::deque <uint64_t> q2;
      uint64_t q2_min_val = UINT64_MAX;
      int q2_min_pos = -1;
      initialize_window(string_hashes, q2, q2_min_val, q2_min_pos, w_max*order+ w_min+x, w_max*order+ w_min+2*x, k);
      windows.emplace_back(q2, q2_min_val, q2_min_pos);

      // initialize deque3
      std::deque <uint64_t> q3;
      uint64_t q3_min_val = UINT64_MAX;
      int q3_min_pos = -1;
      initialize_window(string_hashes, q3, q3_min_val, q3_min_pos, w_max*order+ w_min+2*x, w_max*(order+1), k);
      windows.emplace_back(q3, q3_min_val, q3_min_pos);
    }

//    std::cout << seq << std::endl;

    // create the hybridstrobes
    for (unsigned int i = 0; i <= seq_length; i++) {
//        auto strobe1 = seq.substr(i, k);
//        uint64_t bstrobe = kmer_to_uint64(strobe1, kmask);
//        uint64_t strobe_hashval = hash64(bstrobe, kmask);

        uint64_t bstrobe = string_hashes[i]/n;

        if (bstrobe % 100 < strobe_fraction){
            std::vector<unsigned int> strobe_positions(n, pos_to_seq_choord[i]);

            for (unsigned int strobe_num = 1; strobe_num < n; strobe_num++){
            uint64_t strobe_next_val;
            unsigned int r = bstrobe % 3;
            if (r == 0){
                bstrobe += std::get<1>(windows[strobe_num*3-3])/(n+strobe_num);
                strobe_positions[strobe_num] = pos_to_seq_choord[std::get<2>(windows[strobe_num*3-3])];
            }
            else if (r == 1){
                bstrobe += std::get<1>(windows[strobe_num*3-2])/(n+strobe_num);
                strobe_positions[strobe_num] = pos_to_seq_choord[std::get<2>(windows[strobe_num*3-2])];
            }
            else{
                bstrobe += std::get<1>(windows[strobe_num*3-1])/(n+strobe_num);
                strobe_positions[strobe_num] = pos_to_seq_choord[std::get<2>(windows[strobe_num*3-1])];
            }
            }

            std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (bstrobe, ref_index, strobe_positions[0], strobe_positions[1], strobe_positions[n-1]);
            mixedhybridstrobes.push_back(s);
        }
        else {
            for (unsigned int j = 0; j < n; j++) {
                // bstrobe = ((bstrobe << k) ^ string_hashes[i+j*k]);
                bstrobe += string_hashes[i+j*k]/(n+j);
            }
            std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (bstrobe, ref_index, i, i+k, i+k*(n-1));
            mixedhybridstrobes.push_back(s);
        }


        // update queue and current minimum and position
        for (unsigned int order = 0; order < (n-2); order++){
            uint64_t new_strobe_hashval1 = string_hashes[i + order*w_max + w_min+x];
            update_window(std::get<0>(windows[order*3]), std::get<1>(windows[order*3]), std::get<2>(windows[order*3]), new_strobe_hashval1, order*w_max+w_min, order*w_max+w_min+x, i, seq_length);

            uint64_t new_strobe_hashval2 = string_hashes[i + order*w_max + w_min+2*x];
            update_window(std::get<0>(windows[order*3+1]), std::get<1>(windows[order*3+1]), std::get<2>(windows[order*3+1]), new_strobe_hashval2, order*w_max+w_min+x, order*w_max+w_min+2*x, i, seq_length);

            uint64_t new_strobe_hashval3 = string_hashes[i + order*w_max + w_max];
            update_window(std::get<0>(windows[order*3+2]), std::get<1>(windows[order*3+2]), std::get<2>(windows[order*3+2]), new_strobe_hashval3, order*w_max+w_min+2*x, order*w_max+w_max, i, seq_length);
        }

        if (i + 2*w_max < seq_length){
            uint64_t new_strobe_hashval4 = string_hashes[i + w_max + w_min+x];
            update_window(std::get<0>(windows[(n-2)*3]), std::get<1>(windows[(n-2)*3]), std::get<2>(windows[(n-2)*3]), new_strobe_hashval4, (n-2)*w_max + w_min, (n-2)*w_max + w_min+x, i, seq_length);

            uint64_t new_strobe_hashval5 = string_hashes[i + w_max + w_min+2*x];
            update_window(std::get<0>(windows[(n-2)*3+1]), std::get<1>(windows[(n-2)*3+1]), std::get<2>(windows[(n-2)*3+1]), new_strobe_hashval5, (n-2)*w_max + w_min+x, (n-2)*w_max + w_min+2*x, i, seq_length);

            uint64_t new_strobe_hashval6 = string_hashes[i + 2*w_max];
            update_window(std::get<0>(windows[(n-2)*3+2]), std::get<1>(windows[(n-2)*3+2]), std::get<2>(windows[(n-2)*3+2]), new_strobe_hashval6, (n-2)*w_max + w_min+2*x, (n-1)*w_max, i, seq_length);
        }
        else if ((i + w_max + w_min + 1 < seq_length) && (seq_length <= i + w_max + w_min+x) ){
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(std::get<0>(windows[(n-2)*3]), std::get<1>(windows[(n-2)*3]), std::get<2>(windows[(n-2)*3]), new_strobe_hashval, (n-2)*w_max + w_min, (n-2)*w_max + w_min+x, i, seq_length);
            update_window(std::get<0>(windows[(n-2)*3+1]), std::get<1>(windows[(n-2)*3+1]), std::get<2>(windows[(n-2)*3+1]), new_strobe_hashval, (n-2)*w_max + w_min+x, (n-2)*w_max + w_min+2*x, i, seq_length);
            update_window(std::get<0>(windows[(n-2)*3+2]), std::get<1>(windows[(n-2)*3+2]), std::get<2>(windows[(n-2)*3+2]), new_strobe_hashval, (n-2)*w_max + w_min+2*x, (n-1)*w_max, i, seq_length);
        }
        else if ((i + w_max + w_min+x < seq_length) && (seq_length <= i + w_max + w_min+2*x)  ){
            uint64_t new_strobe_hashval4 = string_hashes[i + w_max + w_min+x];
            update_window(std::get<0>(windows[(n-2)*3]), std::get<1>(windows[(n-2)*3]), std::get<2>(windows[(n-2)*3]), new_strobe_hashval4, (n-2)*w_max + w_min, (n-2)*w_max + w_min+x, i, seq_length);
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(std::get<0>(windows[(n-2)*3+1]), std::get<1>(windows[(n-2)*3+1]), std::get<2>(windows[(n-2)*3+1]), new_strobe_hashval, (n-2)*w_max + w_min+x, (n-2)*w_max + w_min+2*x, i, seq_length);
            update_window(std::get<0>(windows[(n-2)*3+2]), std::get<1>(windows[(n-2)*3+2]), std::get<2>(windows[(n-2)*3+2]), new_strobe_hashval, (n-2)*w_max + w_min+2*x, (n-1)*w_max, i, seq_length);
        }
        else if ((i + w_max + w_min+2*x < seq_length) && (seq_length <= i + 2*w_max) ){
            uint64_t new_strobe_hashval4 = string_hashes[i + w_max + w_min+x];
            update_window(std::get<0>(windows[(n-2)*3]), std::get<1>(windows[(n-2)*3]), std::get<2>(windows[(n-2)*3]), new_strobe_hashval4, (n-2)*w_max + w_min, (n-2)*w_max + w_min+x, i, seq_length);
            uint64_t new_strobe_hashval5 = string_hashes[i + w_max + w_min+2*x];
            update_window(std::get<0>(windows[(n-2)*3+1]), std::get<1>(windows[(n-2)*3+1]), std::get<2>(windows[(n-2)*3+1]), new_strobe_hashval5, (n-2)*w_max + w_min+x, (n-2)*w_max + w_min+2*x, i, seq_length);
            uint64_t new_strobe_hashval =  UINT64_MAX;
            update_window(std::get<0>(windows[(n-2)*3+2]), std::get<1>(windows[(n-2)*3+2]), std::get<2>(windows[(n-2)*3+2]), new_strobe_hashval, (n-2)*w_max + w_min+2*x, (n-1)*w_max, i, seq_length);
        }
        else{
            return mixedhybridstrobes;
        }

//        auto strobe1 = seq.substr(i, k);
//        auto strobe2 = seq.substr(strobe2_pos, k);
//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe2_pos - (i+k) - 1, ' ') << std::string(k, 'X') << std::string(strobe3_pos - strobe2_pos - k, ' ') << std::string(k, 'X') << std::endl;

//        std::cout << std::string(i, ' ') << strobe1 << std::string(strobe2_pos - (i+k)-1, ' ') << std::string(k, 'X') << std::endl;
//        std::cout << i << " " << strobe2_pos << " " << seq_length << std::endl;
//        std::cout << i << ": " << strobe1 << strobe2 << std::endl;


    }
    return mixedhybridstrobes;
}


mers_vector seq_to_multistrobes2(int n, int k, int k_b, int w_min, int w_max, std::string &seq, unsigned int ref_index)
{
    mers_vector multistrobes2;

    if (seq.length() < w_max) {
        return multistrobes2;
    }

    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    uint64_t kmask=(1ULL<<2*k_b) - 1;
    uint64_t q = pow (2, 16) - 1;
//    std::bitset<64> x(q);
//    std::cout << x << '\n';
    // make string of strobes into hashvalues all at once to avoid repetitive k-mer to hash value computations
    std::vector<uint64_t> string_hashes;
    std::vector<uint64_t> string_hashes2;
    std::vector<unsigned int> pos_to_seq_choord;
    uint64_t hash_multistrobe2;
//    robin_hood::unordered_map< unsigned int, unsigned int>  pos_to_seq_choord;
    make_string_to_hashvalues2(seq, string_hashes, pos_to_seq_choord, k_b, kmask);
    unsigned int seq_length = string_hashes.size();

    if (string_hashes.size() == 0) {
        return multistrobes2;
    }
//    std::cout << seq << std::endl;

    // compute windows for all strobe length combinations
    std::vector<uint64_t> lengths;
    std::vector<uint64_t> w_mins;
    std::vector<uint64_t> w_maxs;

    for (unsigned int i = k_b; i < 2*k-k_b; i++){
        lengths.push_back(i);
        w_mins.push_back(w_min + i - k);
        w_maxs.push_back(w_max + i - k);
    }

    unsigned int options = lengths.size();

    // create the multistrobes
    for (unsigned int i = 0; i <= seq_length; i++) {

//        if ((i % 1000000) == 0 ){
//            std::cout << i << " strobemers created." << std::endl;
//        }
        unsigned int strobe_pos_next;
        uint64_t strobe_hashval_next;

        uint64_t strobe_hash = string_hashes[i];
        unsigned int option = strobe_hash % options;
        unsigned int k1 = lengths[option]; // determine length of first strobe
        unsigned int k2 = 2*k-k1; // determine length of second strobe
        unsigned int w_start = i + w_mins[option];
        strobe_hash = make_string_to_hashvalue(seq, k1, i);

        if (i + w_max <= seq_length - 1){
            unsigned int w_end = i + w_maxs[option];
            get_next_strobe(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
        }
        else if ((w_start + 1 < seq_length) && (seq_length <= i + w_maxs[option]) ){
            unsigned int w_end = seq_length -1;
            get_next_strobe(string_hashes, strobe_hash, strobe_pos_next, strobe_hashval_next, w_start, w_end, q);
        }
        else{
            return multistrobes2;
        }
        
        unsigned int seq_pos_strobe1 = pos_to_seq_choord[i];
        unsigned int seq_pos_strobe2 = pos_to_seq_choord[strobe_pos_next];
        unsigned int seq_pos_strobe3 = pos_to_seq_choord[strobe_pos_next + k2 - k_b];
        strobe_hashval_next = make_string_to_hashvalue(seq, k2, seq_pos_strobe2);
        hash_multistrobe2 = (strobe_hash/2) + (strobe_hashval_next/3);

        /* if (k1 > k2){
            uint64_t hash1 = strobe_hash & ((1 << 2*k) - 1);
            uint64_t hash2 = (strobe_hashval_next << 2*k2) ^ (strobe_hash >> 2*k);
            hash_multistrobe2 = (hash1 + hash2);
        } 
        else if (k2 > k1){
            uint64_t hash1 = (strobe_hash << 2*k1) ^ (strobe_hashval_next >> 2*k);
            uint64_t hash2 = strobe_hashval_next & ((1 << 2*k) - 1);
            hash_multistrobe2 = (hash1 + hash2);
        }
        else{
            // hash_multistrobe2 = (strobe_hash + strobe_hashval_next);
            hash_multistrobe2 = ((strobe_hash/2) + (strobe_hashval_next/3);
        } */

        std::tuple<uint64_t, unsigned int, unsigned int, unsigned int, unsigned int> s (hash_multistrobe2, ref_index, seq_pos_strobe1, seq_pos_strobe3, seq_pos_strobe3);
        multistrobes2.push_back(s);
    }
    return multistrobes2;
}


mers_vector seq_to_multistrobes(int n, int k, int k_b, int w_min, int w_max, std::string &seq, unsigned int ref_index)
{
    std::cout << "Multistrobes are just implemeted for order 2." << std::endl;
    exit(0);
}