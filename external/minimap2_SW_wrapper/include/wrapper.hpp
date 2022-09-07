#include "ksw2.h"
#include "kalloc.h"
#include <string>
#include <vector>

#ifndef SW_WRAPPER_HEADER_INCLUDED
#define SW_WRAPPER_HEADER_INCLUDED

namespace SWWrapper {
    class Aligner;
    struct SWAlignment;
}

struct SWWrapper::SWAlignment {
    int score; // alignment score
    std::vector<uint32_t> cigar; // cigar in standard htslib format
};

class SWWrapper::Aligner {
    public:
        Aligner(int8_t match_score, int8_t mismatch_score, int8_t gap_opening, int8_t gap_extension, int8_t gap_opening2, int8_t gap_extension2, int w, int zdrop); //ACGTN default
        Aligner(int8_t match_score, int8_t mismatch_score, int8_t gap_opening, int8_t gap_extension) :
            Aligner(match_score, mismatch_score, gap_opening, gap_extension, gap_opening, gap_extension, -1, -1) {} // ACGTN default
        Aligner(const std::string & chars, int8_t match_score, int8_t mismatch_score, int8_t gap_opening, int8_t gap_extension, int8_t gap_opening2, int8_t gap_extension2, int w, int zdrop); //custom chars
        Aligner(const std::string & chars, int8_t match_score, int8_t mismatch_score, int8_t gap_opening, int8_t gap_extension):
            Aligner(chars, match_score, mismatch_score, gap_opening, gap_extension, gap_opening, gap_extension, -1, -1) {} //custom chars
        bool align(const std::string & ref, const std::string & query, SWAlignment * align);
        bool align_prefix(const std::string & ref, const std::string & query, SWAlignment * align);
        bool align_suffix(const std::string & ref, const std::string & query, SWAlignment * align);
        bool align_suffix_rev_cigar(const std::string & ref, const std::string & query, SWAlignment * align);
        
        ~Aligner();
        
    private:
        uint8_t num_char;
    
        int8_t * score_matrix;
        
        int8_t gapo1;
        int8_t gape1;
        int8_t gapo2;
        int8_t gape2;
        
        std::vector<uint8_t> characters; //mapping from all char to its indices
        
        int w;
        int zdrop;
        
        void *km; // memory space
        
        ksw_extz_t result;
};

#endif
