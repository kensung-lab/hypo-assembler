/*
 * 
 * Copyright (c) 2019, Ritu Kundu and Joshua Casey
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/** Class Alignment.
 * It represents an alignment and contains functionality for kmer -support update; minimizer-support update, and arms finding
 */
#pragma once
#include <sdsl/bit_vectors.hpp>
#include "globalDefs.hpp"
#include "PackedSeq.hpp"
#include "wrapper.hpp"

//#include "Contig.hpp"

namespace hypo
{
class Contig;
const UINT32 INDEL_CLUSTER_THRESHOLD = 30;
const UINT MAP_QUAL_TH = 0;
const float CLIP_FRAC =0.1;
enum class ArmType: UINT8 {
    INTERNAL, 
    PREFIX, 
    SUFFIX,
    EMPTY
};

using Arm = struct SArm{
    const UINT32 windex;
    const PackedSeq<2> arm;
    const ArmType armtype;
    std::string basequal; // Used for long arms
    bool is_high_qual; // Used for short arms
    std::string cigstr;
    UINT32 terminal; // Used for Pref/Suff typeof short arm
    // Used for short arms
    SArm (const UINT32 ind, const PackedSeq<2>& ps, UINT32 left, UINT32 right, const ArmType at, bool is_high, std::string& cigs): 
    windex(ind), arm(std::move(PackedSeq<2>(ps,left,right))), armtype(at),basequal(""), is_high_qual(is_high), cigstr(std::move(cigs))
    {}
    // Used for long arms
    SArm (const UINT32 ind, const PackedSeq<2>& ps, UINT32 left, UINT32 right, const ArmType at,const std::string&bq, bool is_high=true): 
    windex(ind), arm(std::move(PackedSeq<2>(ps,left,right))), armtype(at), basequal(bq,left,right-left), is_high_qual(is_high){}
    SArm (const UINT32 ind, const PackedSeq<2>& ps, UINT32 left, UINT32 right, const ArmType at,const std::string&bq, std::string cigs): 
    windex(ind), arm(std::move(PackedSeq<2>(ps,left,right))), armtype(at), basequal(bq,left,right-left), cigstr(std::move(cigs)), is_high_qual(true){}
    SArm (const UINT32 ind): windex(ind), arm(std::move(PackedSeq<2>())), armtype(ArmType::EMPTY){}
  };

class Alignment {
public:
    Alignment(Contig& contig, bam1_t *hts_align, UINT8 norm_edit_th=100);
    // USed for merging gapped aln
    Alignment(std::unique_ptr<Alignment> const& aln1, std::unique_ptr<Alignment>const&  aln2, UINT32 missing_len, UINT32 gap_len);

    // USed for merging overlapping aln
    Alignment(std::unique_ptr<Alignment> const& aln1, std::unique_ptr<Alignment>const&  aln2, UINT32 missing_len);

    Alignment(const Alignment &) = delete;
    Alignment &operator=(const Alignment &) = delete;
	Alignment(Alignment&&) = delete;
    Alignment& operator=(Alignment&&) = delete;
     ~Alignment()= default;

    // Used to pop_back invalid alignment (a long read alignment may be invalid) from the alignment_store by hypo
    bool is_valid; // Can be false only for a long read

    void update_solidkmers_support (const UINT k, Contig& pc);
    void update_minimisers_support (Contig& contig);

    // For adding to short windows
    void find_short_arms(const UINT k, const Mode mode, Contig& contig);
    /** This will add arms into short windows;
     * */
    void add_short_arms(const Mode mode, const Contig& contig);

    // For adding to long windows and adding long-arms in LO mode
    void find_long_arms(Contig& contig);    
    void add_long_arms(const Contig& contig);

    inline void invalidate() {if (is_valid) {_cigar.clear();_basequal.clear();_apseq.clear();reset_head_tail();} is_valid=false;}
    inline void reset_head_tail() {_head.clear();_tail.clear();_head_basequal.clear(); _tail_basequal.clear();}
    inline UINT32 get_rb() const {return _rb;}
    inline UINT32 get_re() const {return _re;}
    inline std::string get_rname() const {return _qname;}
    inline UINT32 get_qlen() const {return _qae+(_left_clip+_right_clip);}
    
    // SUPPLEMENTARY JOINING: some utilities
    inline UINT32 get_real_length() const {return _qae-_qab+(_left_clip+_right_clip);}
    inline bool is_complete() const {return get_real_length() == _head.get_seq_size() + _apseq.get_seq_size() + _tail.get_seq_size();}
    inline std::string get_complete_seq() const {
        return _head.unpack() + _apseq.unpack() + _tail.unpack();
    }
    inline std::string get_complete_qual() const {
        return _head_basequal + _basequal + _tail_basequal;
    }
    inline bool is_softclipped() const {
        return _head.get_seq_size() + _tail.get_seq_size() > 0;
    }
    inline bool is_unclipped() const {
        return _left_clip + _right_clip == 0;
    }
    static void initiate_sw_engines(const UINT32 num_threads);
    bool join_supplementary(Alignment & pair, std::string full_read, std::string full_qual, Contig & contig, int _engine_idx);
    
    
    inline bool is_rev() const {return _is_rev;}
    // Returns absolute pos in query
    inline void get_qcoords(UINT32& b, UINT32& e) const {
        // Here _qae represents the aligned segment len (same as _apseq.get_seq_size())
         UINT32 q_abs_beg = 0;
         UINT32 qlen = _qae+(_left_clip+_right_clip);
         // Left
         for(UINT32 i = 0; i < _num_cigar_op; i++) {        
            UINT32 current_cigar = _cigar[i];
            UINT32 cigar_op = bam_cigar_op(current_cigar);
            UINT32 oplen = bam_cigar_oplen(current_cigar);
            if(cigar_op == BAM_CSOFT_CLIP || cigar_op == BAM_CHARD_CLIP) {
                q_abs_beg+=oplen;
            }
            else {break;}
         }
         // Right
         UINT32 q_abs_end = q_abs_beg + _qae;
         if (_is_rev){
            b = qlen-q_abs_end;
            e = qlen-q_abs_beg;
         }
         else {
             b = q_abs_beg;
             e = q_abs_end;
         }
    }
    

private:

    UINT32 _rb; // ref beginning/start pos
    UINT32 _re; // ref ending pos (points to 1 past the index of last aligned base)
    std::string _qname; 
    UINT32 _qab; // aligned query beg pos (0 based) in the sequence stored
    UINT32 _qae; // aligned query endpos (0 based) in the sequence stored (points to 1 past the index of last aligned base)
    PackedSeq<2> _apseq; // aligned seq of query in Packed_Seq<2> format
    std::vector<UINT32> _cigar; // cigar string in htslib cigar format
    UINT32 _num_cigar_op; // number of cigar operations
    std::string _basequal; // base quality of sequence; Phred scale with +33 offset
    bool _is_rev;// is mapped to reverse strand
    UINT8 _mq; //mapping quality
    bool _is_properly_paired;

    // Head/tail are soft-clips; reset once long term patterns are done (to save space)
    PackedSeq<4> _head;
    PackedSeq<4> _tail;
    std::string _head_basequal;
    std::string _tail_basequal;

    /* Remembers if hard/soft clipped ends; false initially */
    UINT32 _left_clip;
    UINT32 _right_clip;
     
    std::vector<Arm> _arms;

    // Used for remembering minimsers pos used as anchors (succeeding of pvs window will be the same as preceding of next window)
    // Initially: -2
    // If not found: -1
    // Else: starting position of the anchor in the read
    INT64 _anchor_pos;

    /* Indices of windows falling (completely) in this read (valid only when _is_wind_crossed is set)*/
    // Used only for short-reads.
    bool _is_wind_crossed; // Does this read cross a window fully? (Initially false)
    UINT32 _beg_windex;
    UINT32 _end_windex;

    
    
    void initialise_pos(const bam1_t *hts_align);
    void copy_data(const bam1_t *hts_align);
    void copy_head_tail(const bam1_t *hts_align);

    std::vector<UINT32> find_bp(const sdsl::bit_vector::select_1_type& Sreg_pos, const std::vector<RegionType>& reg_type, const UINT32 beg_ind, const UINT32 end_ind);
    std::vector<UINT32> find_bp_with_cigar_long(const sdsl::bit_vector::select_1_type& Sreg_pos, const std::vector<RegionType>& reg_type, const UINT32 beg_ind, const UINT32 end_ind, std::vector<int> & cigar_cuts);
    std::vector<UINT32> find_bp_with_cigar(const sdsl::bit_vector::select_1_type& Sreg_pos, const std::vector<RegionType>& reg_type, const UINT32 beg_ind, const UINT32 end_ind, const Contig & contig, std::vector<std::string> & cigar_cuts);
    
    inline UINT32 to_most_significant(const UINT32 digit) {
        UINT32 temp = digit;
        UINT32 div = 1;
        while(temp >= 10) {
            temp /= 10;
            div *= 10;
        }
        return digit - (digit % div);
    }
    
    // Returns whether armis too short to be added or not
    bool prepare_short_arm(const UINT k, const UINT32 windex, const UINT32 qb, const UINT32 qe, const ArmType origarmtype, Contig& contig, std::string& cigstr);

    // Assumes ref_pos falls in query alignment; For del/ins in query, returns next base
    inline UINT32 get_query_pos(const UINT32 refpos) {  
        UINT32 current_reference_index = _rb;
        UINT32 current_query_index = 0;    
        UINT32 current_cigar;
        UINT32 get_op;
        UINT32 get_oplen;
        UINT32 qe=0;
        for(UINT32 i = 0; i < _num_cigar_op; i++) {        
            current_cigar = _cigar[i];
            get_op = bam_cigar_op(current_cigar);
            get_oplen = bam_cigar_oplen(current_cigar);
            if(get_op == BAM_CSOFT_CLIP || get_op == BAM_CHARD_CLIP) { //softclips are already discarded
                continue;
            }
            bool should_break = false;
            if((bam_cigar_type(get_op) & 3) == 3) { // consume query and reference
                while(get_oplen > 0) {
                    if (current_reference_index==refpos) {qe=current_query_index;should_break=true;break;}
                    get_oplen--;
                    current_reference_index++;
                    current_query_index++;
                }
                if (should_break) {break;}
            } else if(bam_cigar_type(get_op) & 2) { // consume reference only: deletions
                while(get_oplen > 0) {
                    if (current_reference_index==refpos) {qe=current_query_index;should_break=true;break;}
                    get_oplen--;
                    current_reference_index++;
                }
                if (should_break) {break;}
            } else if(bam_cigar_type(get_op) & 1) { // consume query only: insertions
                // get substring of the inserted bases
                while(get_oplen > 0) {
                    ++current_query_index;
                    get_oplen--;
                }
            }
        }
        return qe;    
    }
    
    
    // SUPPLEMENTARY JOINING: similar to in window, a global pool of alignment engines
    static std::vector<std::shared_ptr<SWWrapper::Aligner>> _minimap_sw_engines;
    static std::vector<std::shared_ptr<SWWrapper::SWAlignment>> _minimap_sw_results;
    
    
    
}; // Alignment
} // namespace hypo

