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
/** Class Contig.
 * It represents a contig and contains the functionality for different processing that will be carried out on a contig.
 */
#pragma once
#include <memory>
#include <mutex>
#include <unordered_map>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp>
#include "suk/SolidKmers.hpp"

#include "globalDefs.hpp"
#include "PackedSeq.hpp"
#include "Alignment.hpp"
#include "Window.hpp"

namespace hypo
{
#define MAX_U16_LIMIT 0xffffu
#define MIN_SR_LEN_LSA 50
using KmerInfo = struct Skmerinfo {
    Skmerinfo(UINT64 kmer):kid(kmer),coverage(0),support(0){}
    Skmerinfo(UINT64 kmer, UINT16 c, UINT16 s):kid(kmer),coverage(c),support(s){}
    UINT64 kid;
    UINT16 coverage;
    UINT16 support;
    std::mutex k_mutex;
};

using MWMinimiserInfo = struct Sminimiserinfo {
    std::vector<UINT32> minimisers; // minimser (2 bit based); The MSB represents validity (0 for valid, 1 for invalid)
    std::vector<UINT32> rel_pos; // Relative position of each minimiser wrt pvs one' first is wrt the start of the window
    std::vector<UINT16> support;
    std::vector<UINT16> coverage;
    std::mutex m_mutex;
};

using DivisionPoint = struct Sdivisionpoint {
    UINT32 start_pos; // position of the first base of the window
    UINT32 minimiser;
};


using SCovInfo = struct SscovInfo {
    SscovInfo():coverage(0), pcoverage(0), scoverage(0) {}
    UINT16 coverage;
    UINT16 pcoverage;
    UINT16 scoverage;
    std::mutex c_mutex;
};

using LCovInfo = struct SlcovInfo {
    SlcovInfo():coverage(0){}
    UINT16 coverage;
    std::mutex c_mutex;
};



class Contig {
public:
    Contig(const UINT32 id, const std::string& name, const std::string seq);
    Contig(const UINT32 id, const kseq_t * ks);
    ~Contig() = default;
    // Fills _solid_pos and _kmerinfo
    void find_solid_pos(const std::unique_ptr<suk::SolidKmers>& pSK);

    /** It will mark the SR and Megawindows.
     * It will also prepare the contig for minimiser based window-divsion in SO mode; called after kmer-support has been updated.
     * */
    void prepare_for_division(const UINT k, const std::string wdir);

    // This will divide the contig into regions (windows); in SO mode, called after minimser-support has been updated .
    void divide_into_regions();
    
    /** This will curate the short windows; resetting those with insufficient short-arms; 
     * Called after preparing  and filling arms.
     * */
    void prune_short_windows();
    
    /** This will first fill arms into short windows;
     * In LSA mode, arms will not be filled in High Cov windows; 
     * In SO/SECOND mode, all arms will be filled.
     * Alignments will be destroyed (except in LSA mode because LSA requres to write reads for re-mapping) 
     * */
    inline void fill_short_windows(std::vector<std::unique_ptr<Alignment>>& alignments, UINT64 initial_aid, UINT64 final_aid) {
        /* Fill Short Windows */
        for (UINT i=initial_aid; i < final_aid;++i) {
            alignments[i]->add_short_arms(_mode, *this);
            alignments[i].reset();
        }
        prune_short_windows();
    }

    /** This will fill long arms into pseudo/long windows
     * Then destroy the ds associated with pseudo-windows (in LSA mode); Called after pseudo-arm/long-arm preparing.
     * Alignments will be destroyed
     * Also calls handling Long Del by long windows
     * */
    inline void fill_long_windows(std::vector<std::unique_ptr<Alignment>>& alignments, UINT64 initial_aid, UINT64 final_aid) {
        /* Fill Long Windows */
        for (UINT i=initial_aid; i < final_aid;++i) {
            if (alignments[i]==nullptr) {continue;}
            alignments[i]->add_long_arms(*this);
            alignments[i].reset();
        }
        handle_old();
        _lcovinfo.clear();
        _lcovinfo.shrink_to_fit();
    }

    
    inline void handle_old() {
        for (auto d:_ovl_ld) {
            auto pos = std::get<0>(d);
            auto del_len = std::get<1>(d);
            // Find if part of a long window; Assumes only one long del in a window
            if (_reg_pos[pos]==0) {
                UINT64 ind = _Rreg_pos(pos)-1; // rank returns number of set bits in [0,pos);so subtract 1 for index
                if (_reg_type[ind]==RegionType::LONG) {
                    auto pindex = _reg_info[ind];
                    _plwindows[pindex]->add_longdel(del_len);
                }
            }
        }
    }
    void handle_ovl_li();
    inline void add_ovl(UINT32 b, UINT32 l) {_ovl_li.emplace_back(std::make_pair(b,l));}
    inline void add_gap(UINT32 b, UINT32 e) {_gap_li.emplace_back(std::make_pair(b,e));}
    inline void add_old(UINT32 b, UINT32 l) {_ovl_ld.emplace_back(std::make_pair(b,l));}
    inline void add_gld(UINT32 b, UINT32 e) {_gap_ld.emplace_back(std::make_pair(b,e));}

    
    inline UINT64 get_num_regions() const {return _reg_type.size()-1;}
    inline UINT64 get_num_sr() const {return _numSR;}
    inline UINT64 get_len_sr() const {return _lenSR;}
    inline UINT32 get_contig_len() const {return _len;}

    /* These are used for poa in windows */
    // Returns the number of short windows initially created
    inline UINT32 get_num_swindows() const {return _pswindows.size();}
    inline UINT32 get_num_lwindows() const {return _plwindows.size();}
    inline bool is_valid_swindow(UINT32 ind) {assert(ind<_pswindows.size()); return (_pswindows[ind]!=nullptr);}
    inline bool is_valid_lwindow(UINT32 ind) {assert(ind<_plwindows.size()); return (_plwindows[ind]!=nullptr);}
    // Assumes ind will be valid
    inline void generate_short_consensus(const UINT64 ind, const UINT32 th) {if(_pswindows[ind]!=nullptr){_pswindows[ind]->generate_consensus(th);}} 
    inline void generate_long_consensus(const UINT64 ind, const UINT32 th) {if(_plwindows[ind]!=nullptr){_plwindows[ind]->generate_consensus(th);}}

    void generate_inspect_file(std::string wdir, std::ofstream& bedfile);

    static void set_mode(Mode mode) {_mode=mode;} 

    friend std::ostream &operator<<(std::ostream &, const Contig &);
    friend class Alignment;

    inline std::string get_name() {return _name;};

private:
    const UINT32 _id;
    const std::string _name;
    const UINT32 _len;
    PackedSeq<4> _pseq;
    // Total Number and the length of SR
    UINT64 _numSR;
    UINT64 _lenSR;
    UINT32 _num_wind; // number of windows initially created; (only long windows in LSA mode)

    static Mode _mode;

    // Created by find_solid_pos; Will be used for finding SR; Will be discarded after that
    sdsl::bit_vector _solid_pos;
    sdsl::bit_vector::rank_1_type _Rsolid_pos;
    sdsl::bit_vector::select_1_type _Ssolid_pos;
    std::vector<std::unique_ptr<KmerInfo>> _kmerinfo; 
    

    /** Stores the first and the last kmer-id for each SR (starting from index 1). Index 0 is dummy.
     * Every i^th SR (1-based) will have its first and the last KID stored at (2i-1) and (2i).
     * Each window of type SWS, SW, WS, SWM, or MWS will make use of kmer in this.
     * Will be discarded after filling Short arms.
     * 
     * */
    std::vector<UINT64> _anchor_kmers;

    /** DS used for dividing MegaWindows (Weak Region) (Window between SRs) into windows
     * _reg_pos: bit-vector marking the beginning of each region (SR or MegaWindow); has a dummy at the end (= len of contig)
     *      - Later when MW will be divided into windows, more bits will be set in it (corresponding to each window and Minimiser-SR)
     * rank-select: supporting DS
     * 
     * _is_win_even: Marks whether MW are even numbered(0,2,4,..) or not (1,3,5,..)
     *   - We do not need to save type of each region as it can be deduced.
     *     * if (_is_win_even): i%2==0 => MW; otherwise => SR
     *     * if (!_is_win_even): i%2==0 => SR; otherwise => MW
     *   - Index of Hashtable for each MW can also be deduced
     *     * if (_is_win_even): i/2 is the index
     *     * Else: (i-1)/2 is the index
     * 
     * _minimserinfo: Vector containing, for each MW, its minimser info (hastable for support and corresponding mutex)
     * 
     * Except _reg_pos, rest will be destroyed after division into short-windows
     * */
    sdsl::bit_vector _reg_pos;
    sdsl::bit_vector::rank_1_type _RMreg_pos;
    sdsl::bit_vector::select_1_type _SMreg_pos;
    bool _is_win_even;
    std::vector<std::unique_ptr<MWMinimiserInfo>> _minimserinfo;

    /** Together, the following 3 vectors represent region-division map of the contig
     * A region is either an SR (/MSR) or a window
     * 
     * -reg_type: Stores the type of region
     * 
     * _reg_info: Stores 32-bit unsigned integer showing different entities for different region types;
     *      - minimser value for MSR 
     *      - rank (1-based) for SR (i^th (1-based) in the contig SR will have i stored as info): 
     *          + To find its corresponding kmers in _anchor_kmers
     *          + This is done so that SW, WS , SWS, SWM, MWS know their preceding/succeeding  SR and hence the anchor-KID (through _anchor_kmers)
     *      - Other Region types (i.e. windows) index to their window-pointer in corresponding array;
     *          + Invalid/NoPol windows will have 0; Therefore, window-pointer arrays will have null pointers at 0th index.
     *          + Short windows (all types) as well as long windows (and Long Ins i.e. LI window) will be indexed to _pswindows/_plwindows/_pliwindows resp.
     *          + Short/Long windows also have it as the index of corresponding _cov_info (therefore, _cov_info)
     *
     * 
     * _pswindows/_plwindows: stores pointers to Short/Long Windows
     *      - _pwindows will be a nullptr at 0th index and for invalid windows 
     *      - (short windows will be invalid when merged as long windows)
     * */
    
    sdsl::bit_vector::rank_1_type _Rreg_pos;
    sdsl::bit_vector::select_1_type _Sreg_pos;
    std::vector<RegionType> _reg_type;
    std::vector<UINT32> _reg_info;
    std::vector<std::unique_ptr<Window>> _plwindows;
    std::vector<std::unique_ptr<Window>> _pswindows;
    std::vector<std::unique_ptr<Window>> _pliwindows;
    std::vector<std::unique_ptr<SCovInfo>> _scovinfo; // Needed for SECOND/SHORT
    std::vector<std::unique_ptr<LCovInfo>> _lcovinfo; // Needed for LO and LSA

    // Handling Lon INS
    std::vector<std::pair<UINT32,UINT32>> _ovl_li; // beg and repeat-len of LI
    std::vector<std::pair<UINT32,UINT32>> _gap_li; // beg and end of gap of LI
    std::vector<std::pair<UINT32,UINT32>> _ovl_ld; // beg and del-len of LD
    std::vector<std::pair<UINT32,UINT32>> _gap_ld; // beg and end of gap of LD

    inline void increment_kmer_support(const UINT32 ind) { std::lock_guard<std::mutex> guard(_kmerinfo[ind]->k_mutex); ++(_kmerinfo[ind]->support);}
    inline void increment_kmer_coverage(const UINT32 ind) { std::lock_guard<std::mutex> guard(_kmerinfo[ind]->k_mutex); ++(_kmerinfo[ind]->coverage);}
    inline void increment_minimser_support(const UINT32 mininfo_idx, const UINT32 mindex) { 
        {std::lock_guard<std::mutex> guard(_minimserinfo[mininfo_idx]->m_mutex); ++(_minimserinfo[mininfo_idx]->support[mindex]);}
    }
    inline void increment_minimser_coverage(const UINT32 mininfo_idx, const UINT32 mindex) { 
        {std::lock_guard<std::mutex> guard(_minimserinfo[mininfo_idx]->m_mutex); ++(_minimserinfo[mininfo_idx]->coverage[mindex]);}
    }
    inline void increment_lwindow_coverage(const UINT32 ind) { std::lock_guard<std::mutex> guard(_lcovinfo[ind]->c_mutex); ++(_lcovinfo[ind]->coverage);}
    inline void increment_swindow_coverage(const UINT32 ind) { std::lock_guard<std::mutex> guard(_scovinfo[ind]->c_mutex); ++(_scovinfo[ind]->coverage);}
    inline void increment_swindow_pcoverage(const UINT32 ind) { std::lock_guard<std::mutex> guard(_scovinfo[ind]->c_mutex); ++(_scovinfo[ind]->pcoverage);}
    inline void increment_swindow_scoverage(const UINT32 ind) { std::lock_guard<std::mutex> guard(_scovinfo[ind]->c_mutex); ++(_scovinfo[ind]->scoverage);}
    void initialise_minimserinfo(const std::string& draft_seq, UINT32 minfoind);
    void divide(const UINT32 reg_index, const UINT32 beg, const UINT32 end, char pvs, char nxt);
    void force_divide(const UINT32 beg, const UINT32 end, char pvs, char nxt);
    void fixed_divide(const UINT32 beg, const UINT32 end);
    inline bool is_periodic(const std::string s) {
        auto m = s.size();
        std::vector<UINT16> ft(m,0);
        UINT16 ln = 0; //length of pvs pref-suff
        UINT16 i=1; //
        while (i < m) {
            if (s[i] == s[ln]) {
                ++ln;
                ft[i]=ln;
                ++i;
            }
            else {
                if (ln!=0) {ln = ft[ln-1];}
                else {ft[i]=0;++i;}
            }
        }
        INT32 per = m-ft[m-1];
        bool is_periodic = (per<m) && (m%per==0);
	return is_periodic;
    }
    inline void load_vector(const std::string ifn, std::vector<UINT32>& datavec) {
        std::ifstream infile;        
        infile.open (ifn, std::ios::in | std::ios::binary);
        if (!infile.is_open()) {
            fprintf(stderr, "[Hypo::Contig] Error: File open error: Vector File (%s) could not be opened!\n",ifn.c_str());
            exit(1);
        }

        while (infile) {
            UINT32 val;
            infile.read(reinterpret_cast<char *>(&val), sizeof(UINT32));
            if (infile.bad()) {
                fprintf(stderr, "[Hypo::Contig] Error: File read error: Failed to read from Vector File (%s)!\n",ifn.c_str());
                exit(1);
            }
            if (infile.eof()) break;
            datavec.push_back(val);
        }
        infile.close();
    }
    inline void save_vector(const std::string ofn, const std::vector<UINT32>& datavec) {
        std::ofstream outfile;
        outfile.open (ofn, std::ios::out | std::ios::trunc | std::ios::binary);

        if (!outfile.is_open()) {
            fprintf(stderr, "[Hypo::Contig] Error: File open error: Vector File (%s) could not be opened for writing!\n",ofn.c_str());
            exit(1);
        }

        for (auto val : datavec) {
            outfile.write(reinterpret_cast<const char *>(&val), sizeof(UINT32));
            if (outfile.bad()) {
                fprintf(stderr, "[Hypo::Contig] Error: File write error: Failed to rewritead to Vector File (%s)!\n",ofn.c_str());
                exit(1);
            }
        }
        outfile.close();
    }
}; // Polisher
} // namespace hypo
