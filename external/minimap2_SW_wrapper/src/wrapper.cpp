#include "ksw2.h"
#include "kalloc.h"
#include <string>
#include "wrapper.hpp"
#include <assert.h>
#include <cstring>
#include <iostream>

SWWrapper::Aligner::Aligner(int8_t match_score, int8_t mismatch_score, int8_t gap_opening, int8_t gap_extension, int8_t gap_opening2, int8_t gap_extension2, int w, int zdrop) {
    // mm should be negatives, gaps should be positives
    
    mismatch_score = (mismatch_score > 0) ? -1 * mismatch_score : mismatch_score;
    gap_opening = (gap_opening > 0) ? gap_opening : -1 * gap_opening;
    gap_opening2 = (gap_opening2 > 0) ? gap_opening2 : -1 * gap_opening2;
    gap_extension = (gap_extension > 0) ? gap_extension : -1 * gap_extension;
    gap_extension2 = (gap_extension2 > 0) ? gap_extension2 : -1 * gap_extension2;
    
    num_char = 5; // A, C, G, T, N
    
    unsigned char seq_nt4_table[] = {
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
    };
    
    characters.resize(256);
    for(int i = 0; i <= 255; i++) characters[i] = seq_nt4_table[i];
    
    score_matrix = new int8_t[num_char*num_char];
    for(int i = 0; i < num_char; i++) for(int j = 0; j < num_char; j++) {
        int get_idx = i * num_char + j;
        if(i == j) score_matrix[get_idx] = match_score;
        else score_matrix[get_idx] = mismatch_score;
    }
    gapo1 = gap_opening;
    gape1 = gap_extension;
    gapo2 = gap_opening2;
    gape2 = gap_extension2;
    
    this-> w = w;
    this->zdrop = zdrop;
    
    km = km_init();
    
    memset(&result, 0, sizeof(ksw_extz_t));
}

SWWrapper::Aligner::Aligner(const std::string & chars, int8_t match_score, int8_t mismatch_score, int8_t gap_opening, int8_t gap_extension, int8_t gap_opening2, int8_t gap_extension2, int w, int zdrop) {
    // mm should be negatives, gaps should be positives
    mismatch_score = (mismatch_score > 0) ? -1 * mismatch_score : mismatch_score;
    gap_opening = (gap_opening > 0) ? gap_opening : -1 * gap_opening;
    gap_opening2 = (gap_opening2 > 0) ? gap_opening2 : -1 * gap_opening2;
    gap_extension = (gap_extension > 0) ? gap_extension : -1 * gap_extension;
    gap_extension2 = (gap_extension2 > 0) ? gap_extension2 : -1 * gap_extension2;
    
    num_char = chars.size(); // A, C, G, T, N
    
    characters.resize(256, num_char);
    
    int current_indexing = 0;
    for(int i = 0; i < num_char; i++) characters[chars[i]] = current_indexing++;
    
    score_matrix = new int8_t[num_char*num_char];
    for(int i = 0; i < num_char; i++) for(int j = 0; j < num_char; j++) {
        int get_idx = i * num_char + j;
        if(i == j) score_matrix[get_idx] = match_score;
        else score_matrix[get_idx] = mismatch_score;
    }
    gapo1 = gap_opening;
    gape1 = gap_extension;
    gapo2 = gap_opening2;
    gape2 = gap_extension2;
    
    this-> w = w;
    this->zdrop = zdrop;
    
    km = km_init();
    
    memset(&result, 0, sizeof(ksw_extz_t));
}

SWWrapper::Aligner::~Aligner() {
    delete [] score_matrix;
}

bool SWWrapper::Aligner::align(const std::string & ref, const std::string & query, SWAlignment * align) {
    ksw_reset_extz(&result);
    
    int tlen = ref.size();
    uint8_t * target = new uint8_t[tlen];
    for(int i = 0; i < tlen; i++) {
        target[i] = characters[ref[i]];
        assert(target[i] < num_char);
    }
    
    int qlen = query.size();
    uint8_t * qr = new uint8_t[qlen];
    for(int i = 0; i < qlen; i++) {
        qr[i] = characters[query[i]];
        assert(qr[i] < num_char);
    }
    
    // we want indels
    int flag = KSW_EZ_GENERIC_SC;
    ksw_extd2_sse(km, qlen, qr, tlen, target, num_char, score_matrix, gapo1, gape1, gapo2, gape2, w, zdrop, 0, flag, &result);
    
    align->score = result.score;
    align->cigar.clear();
    //for(uint32_t i = 0; i < result.m_cigar; i++) align->cigar.push_back(result.cigar[i]);
    for(uint32_t i = 0; i < result.n_cigar; i++) align->cigar.push_back(result.cigar[i]);
    
    delete [] target;
    delete [] qr;
    
    return true;
}

bool SWWrapper::Aligner::align_prefix(const std::string & ref, const std::string & query, SWAlignment * align) {
    ksw_reset_extz(&result);
    
    int tlen = ref.size();
    uint8_t * target = new uint8_t[tlen];
    for(int i = 0; i < tlen; i++) {
        target[i] = characters[ref[i]];
        assert(target[i] < num_char);
    }
    
    int qlen = query.size();
    uint8_t * qr = new uint8_t[qlen];
    for(int i = 0; i < qlen; i++) {
        qr[i] = characters[query[i]];
        assert(qr[i] < num_char);
    }
    
    // we want indels + extension mode, enabling free ends and trimming cigars
    int flag = KSW_EZ_GENERIC_SC ^ KSW_EZ_EXTZ_ONLY;
    // reads should be fully aligned, so we want to assign such end bonus
    ksw_extd2_sse(km, qlen, qr, tlen, target, num_char, score_matrix, gapo1, gape1, gapo2, gape2, w, zdrop, 1000, flag, &result);
    
    align->score = result.score;
    align->cigar.clear();
    for(uint32_t i = 0; i < result.n_cigar; i++) align->cigar.push_back(result.cigar[i]);
    
    delete [] target;
    delete [] qr;
    
    return true;
}

bool SWWrapper::Aligner::align_suffix(const std::string & ref, const std::string & query, SWAlignment * align) {
    ksw_reset_extz(&result);
    
    // minimap2 only partial aligns from left to right. Thus, suffix alignment should be reversed
    int tlen = ref.size();
    uint8_t * target = new uint8_t[tlen];
    for(int i = 0; i < tlen; i++) {
        target[i] = characters[ref[tlen - i - 1]];
        assert(target[i] < num_char);
    }
    
    int qlen = query.size();
    uint8_t * qr = new uint8_t[qlen];
    for(int i = 0; i < qlen; i++) {
        qr[i] = characters[query[qlen - i - 1]];
        assert(qr[i] < num_char);
    }
    
    // we want indels + make minimap2 right (effectively, left) aligns gaps + extension mode, enabling free ends and trimming cigars, also reverse the cigar since the seqs are also reversed
    int flag = KSW_EZ_GENERIC_SC ^ KSW_EZ_RIGHT ^ KSW_EZ_REV_CIGAR ^ KSW_EZ_EXTZ_ONLY;
    // reads should be fully aligned, so we want to assign such end bonus
    ksw_extd2_sse(km, qlen, qr, tlen, target, num_char, score_matrix, gapo1, gape1, gapo2, gape2, w, zdrop, 1000, flag, &result);
    
    align->score = result.score;
    align->cigar.clear();
    //for(uint32_t i = 0; i < result.m_cigar; i++) align->cigar.push_back(result.cigar[i]);
    for(uint32_t i = 0; i < result.n_cigar; i++) align->cigar.push_back(result.cigar[i]);
    
    delete [] target;
    delete [] qr;
    
    return true;
}

bool SWWrapper::Aligner::align_suffix_rev_cigar(const std::string & ref, const std::string & query, SWAlignment * align) {
    ksw_reset_extz(&result);
    
    // minimap2 only partial aligns from left to right. Thus, suffix alignment should be reversed
    int tlen = ref.size();
    uint8_t * target = new uint8_t[tlen];
    for(int i = 0; i < tlen; i++) {
        target[i] = characters[ref[tlen - i - 1]];
        assert(target[i] < num_char);
    }
    
    int qlen = query.size();
    uint8_t * qr = new uint8_t[qlen];
    for(int i = 0; i < qlen; i++) {
        qr[i] = characters[query[qlen - i - 1]];
        assert(qr[i] < num_char);
    }
    
    // we want indels + make minimap2 right (effectively, left) aligns gaps + extension mode, enabling free ends and trimming cigars, don't reverse the cigar: start from end
    int flag = KSW_EZ_GENERIC_SC ^ KSW_EZ_RIGHT ^ KSW_EZ_EXTZ_ONLY;
    // reads should be fully aligned, so we want to assign such end bonus
    ksw_extd2_sse(km, qlen, qr, tlen, target, num_char, score_matrix, gapo1, gape1, gapo2, gape2, w, zdrop, 1000, flag, &result);
    
    align->score = result.score;
    align->cigar.clear();
    //for(uint32_t i = 0; i < result.m_cigar; i++) align->cigar.push_back(result.cigar[i]);
    for(uint32_t i = 0; i < result.n_cigar; i++) align->cigar.push_back(result.cigar[i]);
    
    delete [] target;
    delete [] qr;
    
    return true;
}
