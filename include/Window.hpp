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
/** Class Window.
 * It contains the functionality for POA.
 */
#pragma once
#include <cmath>
#include <mutex>
#include <algorithm>
#include <unordered_map>
#include "spoa/spoa.hpp"
#include "globalDefs.hpp"
#include "PackedSeq.hpp"
#include "wrapper.hpp"

namespace hypo
{
const std::string cMarkerLetters = "JO";
const std::string cHead = "JJJJJ";
const std::string cTail = "OOOOO";
const UINT cMarkerLen = cHead.size();
#define HIGHEST_PHRED "~~~~~"
const UINT64 cNumDigits[9] = {0,1,10,100,1000,10000,100000,1000000,10000000};
const UINT64 cRange[9] =     {0,4,9, 90,400,9000,90000,900000,9000000};

const double CONSENSUS_SUPPORT_FRACTION = 0.80;
const double MOST_FREQUENT_FRACTION1 = 0.40;
const double MOST_FREQUENT_FRACTION2 = 0.80;
const int CONSENSUS_CONSTANT_THRESHOLD = 5;
const UINT READS_CLUSTER_THRESHOLD = 20;

const int CONSENSUS_CONSTANT_THRESHOLD2 = 2;
const int RUNLEN_TH = 5;
const double PS_RUN_FRAC = 0.3;

const float LONG_TH = 0.5;

const float LONG_CLUS_LEN_TH = 0.2;
const int VERY_LONG_INDEL = 2000;
const float LONG_CLUS_SECOND_TH = 0.2;
const float MIN_CLUSTER_STRENGTH = 5;

const int RELIABLE_DRAFT_LTH = 25;
const int RELIABLE_DRAFT_STH = 25;

const int DOWNSAMPLE_TH = 30;

const int M = 0;
const int X = 1;
const int I = 2;
const int D = 3;
const int DIST = 4;
const int ALLOWED_ERR = 1;

//-A2 -B8 -O12,32 -E2,1
const int sw_match_score = 2;
const int sw_mismatch_score = 8;
const int sw_gap_open1 = 12;
const int sw_gap_open2 = 32;
const int sw_gap_ext1 = 2;
const int sw_gap_ext2 = 1;
const int lr_match_score = 3;
const int lr_misMatch_score = -5;
const int lr_gap_penalty = -4;

enum class WindowType : UINT8
{
    SHORT,
    LONG
};

enum class AlnType : char
{
    P,
    S,
    I
};

class Window
{
public:
    Window() : _wtype(WindowType::SHORT),
               _num_internal(0), _num_pre(0), _num_suf(0), _num_empty(0), _longest_int_len(0), is_internal_only(false), _empty_maxalen(0), 
               _longest_pre_len(0), _longest_suf_len(0), _num_internal_low_qual(0), _long_del(0) {}

    Window(const PackedSeq<4> &ps, const size_t left_ind, const size_t right_ind, WindowType wt) : _wtype(wt), _longest_int_len(0), _long_del(0), _empty_maxalen(0),
                                                                                                   is_internal_only(false), _num_internal(0), _num_pre(0), _num_suf(0), _num_empty(0), _longest_pre_len(0), _longest_suf_len(0), _draft(ps, left_ind, right_ind), _num_internal_low_qual(0)
    {
    }

    Window(const Window &) = delete;
    Window &operator=(const Window &) = delete;
    Window(Window &&) = delete;
    Window &operator=(Window &&) = delete;
    ~Window() = default;

    static void prepare_for_poa(const UINT32 num_threads);
    static void prepare_for_cluster(const UINT32 num_threads);
    void generate_consensus(const UINT32 engine_idx);
    inline std::string get_consensus() const { return _consensus; }
    inline size_t get_consensus_len() const { return _consensus.size(); }
    inline size_t get_window_len() const { return _draft.get_seq_size(); }

    inline void add_prefix(const PackedSeq<2> &ps)
    {
        {//std::lock_guard<std::mutex> guard(w_mutex);
            UINT arm_len = ps.get_seq_size();
            if (arm_len > _longest_pre_len)
            {
                _longest_pre_len = arm_len;
            }
            ++_num_pre;
            _pre_arms.emplace_back(std::move(ps));
        }
    }

    inline void add_suffix(const PackedSeq<2> &ps)
    {
        {//std::lock_guard<std::mutex> guard(w_mutex);
            UINT arm_len = ps.get_seq_size();        
            if (arm_len > _longest_suf_len)
            {
                _longest_suf_len = arm_len;
            }
            ++_num_suf;
            _suf_arms.emplace_back(std::move(ps));
        }
    }

    inline void add_internal(const PackedSeq<2> &ps)
    {
        {//std::lock_guard<std::mutex> guard(w_mutex);
            UINT arm_len = ps.get_seq_size();
            if (_longest_int_len < arm_len) {
                _longest_int_len = arm_len;
            }
            ++_num_internal;
            _internal_arms.emplace_back(std::move(ps));
        }
    }

    inline void add_internal_with_cigar(const PackedSeq<2> &ps, const std::string &cig)
    {
        {//std::lock_guard<std::mutex> guard(w_mutex);
            UINT arm_len = ps.get_seq_size();
            if (_longest_int_len < arm_len) {
                _longest_int_len = arm_len;
            }
            ++_num_internal;
            _internal_arms.emplace_back(std::move(ps));
            _internal_cigs.emplace_back(std::move(cig));
        }
    }

    inline void add_internal_cig(const std::string &cig)
    {
        {//std::lock_guard<std::mutex> guard(w_mutex);
            _internal_cigs.emplace_back(std::move(cig));
        }
    }

    inline void add_pref_end(const UINT32 pos)
    {
        {//std::lock_guard<std::mutex> guard(w_mutex);
            _pends.emplace_back(pos);
        }
    }

    inline void add_suff_beg(const UINT32 pos)
    {
        {//std::lock_guard<std::mutex> guard(w_mutex);
            _sbegs.emplace_back(pos);
        }
    }

    inline void add_pref_clip(const UINT32 c) {_pclips.emplace_back(c);}
    inline void add_suff_clip(const UINT32 c) {_sclips.emplace_back(c);}

   
    inline void add_pref_wt(const UINT8 w)
    {
            _pre_awt.emplace_back(w);
    }

    inline void add_suff_wt(const UINT8 w)
    {
            _suf_awt.emplace_back(w);
    }

    inline void add_internal_wt(const UINT8 w)
    {
            _internal_awt.emplace_back(w);
    }

    inline void add_max_areflen(const UINT64 l)
    { // Called after adding an internal (long) arm; Update the max len if this has it larger
        _areflen.emplace_back(l);
    }

    inline void add_internal(const PackedSeq<2> &ps, const std::string &basequal)
    {
        {//std::lock_guard<std::mutex> guard(w_mutex);
            UINT arm_len = ps.get_seq_size();
            if (_longest_int_len < arm_len) {
                _longest_int_len = arm_len;
            }
            ++_num_internal;
            _internal_arms.emplace_back(std::move(ps));
            _internal_qual.emplace_back(std::move(basequal));
        }
    }

    inline void add_prefix(const PackedSeq<2> &ps, const std::string &basequal)
    {
        {//std::lock_guard<std::mutex> guard(w_mutex);
            UINT arm_len = ps.get_seq_size();
            if (arm_len > _longest_pre_len)
            {
                _longest_pre_len = arm_len;
            }
            ++_num_pre;
            _pre_arms.emplace_back(std::move(ps));
            _pre_qual.emplace_back(std::move(basequal));
        }
    }

    inline void add_suffix(const PackedSeq<2> &ps, const std::string &basequal)
    {
        {//std::lock_guard<std::mutex> guard(w_mutex);
            UINT arm_len = ps.get_seq_size();        
            if (arm_len > _longest_suf_len)
            {
                _longest_suf_len = arm_len;
            }
            ++_num_suf;
            _suf_arms.emplace_back(std::move(ps));
            _suf_qual.emplace_back(std::move(basequal));  
        }  
    }

    inline void add_empty(UINT64 alen=0) {
        {//std::lock_guard<std::mutex> guard(w_mutex); 
            ++_num_empty; 
        }
        if(alen>_empty_maxalen) {_empty_maxalen=alen;}    
    }
    inline void add_longdel(UINT32 ldlen) {_long_del=ldlen;}

    void reset_consensus(const std::string& s1,const std::string& s2);
       
    inline UINT32 get_num_pre() const { return _num_pre; }
    inline UINT32 get_num_suf() const { return _num_suf; }
    inline UINT32 get_num_internal() const { return _num_internal + _num_empty; }
    inline UINT32 get_num_total() const { return _num_internal + _num_empty + _num_pre + _num_suf; }
    inline UINT32 get_maxlen_pre() const { return _longest_pre_len; }
    inline UINT32 get_maxlen_suf() const { return _longest_suf_len; }
    inline bool is_sh_high() const { return (_num_internal > SWindow_cov_settings.high_th); }
    inline bool is_sh_low() const { return ((_num_internal+std::max(_num_pre,_num_suf)) < SWindow_cov_settings.low_th); }
    inline bool is_superlow_ps() const {UINT MIN_SHORT_NUM = Arms_settings.min_short_num;
     return ((_num_internal < MIN_SHORT_NUM) && (_num_pre < MIN_SHORT_NUM || _num_suf < MIN_SHORT_NUM)); }
    inline void set_low_qual_internal() { 
        {
            //std::lock_guard<std::mutex> guard(w_mutex); 
            ++_num_internal_low_qual; 
        }
    }
    inline void set_internal_only() { is_internal_only=true; }
    inline bool has_low_qual_internal() const { return _num_internal_low_qual>std::ceil(Arms_settings.lq_frac_th*_num_internal); }
    inline void clear_pre_suf()
    {
        _num_pre = 0;
        _num_suf = 0;
        _pre_arms.clear();
        _suf_arms.clear();
        _pre_arms.shrink_to_fit();
        _suf_arms.shrink_to_fit();
        _pre_qual.clear();
        _suf_qual.clear();
        _pre_awt.clear();
        _suf_awt.clear();
        _pends.clear();
        _sbegs.clear();
        _pclips.clear();
        _sclips.clear();
    }

    static void set_mode(Mode mode) {_mode=mode;} 
    
    friend std::ostream &operator<<(std::ostream &, const Window &);

protected:
    //std::mutex w_mutex;
    static Mode _mode;
    UINT32 _engine_idx;
    std::string info;
    WindowType _wtype;
    UINT32 _num_internal;
    UINT32 _num_pre;
    UINT32 _num_suf;
    UINT32 _num_empty;
    PackedSeq<4> _draft;
    std::vector<PackedSeq<2>> _internal_arms;
    std::vector<PackedSeq<2>> _pre_arms;
    std::vector<PackedSeq<2>> _suf_arms;
    std::string _consensus;
    /** For storing multiple consensus
     * Used only in spcl cases (Internal only window with 'good' clusters);
     * Has marked consensuses (Actual string is from (begin()+cMarkerLen) to (end()-cMarkerLen) )
     * Main consensus will be chosen from it.
     * 
     */
    std::vector<std::string> _multi_consensus;
    // Base quality
    std::vector<std::string> _internal_qual;
    std::vector<std::string> _pre_qual;
    std::vector<std::string> _suf_qual;

    // IF an arm is reliable (i.e. properly paired read); A reliable arms will have twice wt (2) than that of unreliable (1)
    std::vector<UINT8> _internal_awt;
    std::vector<UINT8> _pre_awt;
    std::vector<UINT8> _suf_awt;
    // Used only in short window
    UINT32 _num_internal_low_qual; // Number of internal arms having 2 or more bases of low quality;
    std::vector<std::string>_internal_cigs;
    bool is_internal_only; // should window be tried to get polished using internalreads only (Used only for SHORT windows)

    // Keep track of longest arms of each type (len)
    UINT32 _longest_pre_len;
    UINT32 _longest_suf_len;
    UINT32 _longest_int_len;

    // Ending/Starting pos of each prefix/suffix; used in PS windows to find if it needs patching
    std::vector<UINT32> _pends;
    std::vector<UINT32> _sbegs;
    // Clipped lengths of each prefix/suffix; Needed where there is no spanning P/S to determining patching units
    std::vector<UINT32> _pclips;
    std::vector<UINT32> _sclips;

    // Len of any long del in this window, initialisedto 0
    UINT32 _long_del;

    // Aligned reference length; Considers the whole alignment of the read not just this window;
    // Used to select the cluster in long windows. 
    std::vector<UINT64> _areflen;
    UINT64 _empty_maxalen;

    static std::vector<std::shared_ptr<spoa::AlignmentEngine>> _alignment_engines_long;
    static std::vector<std::shared_ptr<spoa::AlignmentEngine>> _alignment_engine_pslong;
    static std::vector<std::shared_ptr<SWWrapper::Aligner>> _minimap_sw_engines;
    static std::vector<std::shared_ptr<SWWrapper::SWAlignment>> _minimap_sw_results;

    // Setting consensus
    inline void set_consensus(std::string const &con) { _consensus = std::move(con); }
    inline void set_marked_consensus(std::string const &con) { 
        if (con.size() <2*cMarkerLen) {_consensus.assign(_draft.unpack());}
        else {
            auto boffset = (con[0]==cHead[0])?(cMarkerLen):(0);
            auto eoffset = (con[con.size()-1]==cTail[0])?(cMarkerLen):(0);
            _consensus.assign(con.begin() + boffset, con.end() - eoffset);
        } 
    }
    inline void set_marked3_consensus(std::string const &con) { 
        if (con.size() <2*cMarkerLen || con[0]!=cHead[0] || con[con.size()-1]!=cTail[0]) {_consensus.assign(_draft.unpack());}
        else {_consensus.assign(con.begin() + cMarkerLen, con.end() - cMarkerLen);} 
    }
    // For N's as head and tail (coming from Short Pileup approach used in PSR/PSA)
    inline void set_marked2_consensus(std::string const &con) { 
        if (con.size() <2*cMarkerLen || con[0]!='N' || con[con.size()-1]!='N') {_consensus.assign(_draft.unpack());}
        else {_consensus.assign(con.begin() + cMarkerLen, con.end() - cMarkerLen);} 
    } 

    // Clearing
    inline void clear_arms()
    {
        _num_internal = 0;
        _internal_arms.clear();
        _internal_arms.shrink_to_fit();
        _internal_qual.clear();
        _internal_awt.clear();
        _internal_cigs.clear();

        clear_pre_suf();
    }
    
    // Consensus generation in Long windows using POA
    void generate_consensus_long();
    std::string curate_long(const std::string &con, const std::vector<UINT32> &dst, const UINT32 num_arms);
    bool has_long_indel(std::string& cig);

    // Consensus generation in short window; Used for short windows with internal arms only;
    void generate_consensus_short();

    /* Most frequent arm method */
    bool select_most_frequent();
    void select_consensus_from_draft(bool is_marked);
    bool is_reliable(std::string& cig);
    
    /* PScluster method; Used in short window when most frequent internal arm method fails or window has too many pref/suff */
    void generate_pscluster_consensus_short();
    std::string get_SW_pileup_consensus(const std::string& backbone, const std::vector<std::string>& internal, const std::vector<std::vector<UINT16>>& iw, const std::vector<std::string>& pre, const std::vector<std::vector<UINT16>>& pw, const std::vector<std::string>& suf, const std::vector<std::vector<UINT16>>& sw);
    void sw_align_internal(const std::string& bb, const std::string& seq, const std::vector<UINT16>& w, std::vector<UINT16>& cov, std::vector<std::vector<UINT16>>& bases, std::vector<std::unordered_map<UINT16, UINT16>> & num_dels, std::vector<std::unordered_map<std::string,UINT16>>& ins);
    void sw_align_prefix(const std::string& bb, const std::string& seq, const std::vector<UINT16>& w, std::vector<UINT16>& cov, std::vector<std::vector<UINT16>>& bases, std::vector<std::unordered_map<UINT16, UINT16>> & num_dels, std::vector<std::unordered_map<std::string,UINT16>>& ins);
    void sw_align_suffix(const std::string& bb, const std::string& seq, const std::vector<UINT16>& w, std::vector<UINT16>& cov, std::vector<std::vector<UINT16>>& bases, std::vector<std::unordered_map<UINT16, UINT16>> & num_dels, std::vector<std::unordered_map<std::string,UINT16>>& ins);
    // Patching related
    bool identify_patch_region(UINT32& patchb, UINT32& patche, UINT32& patch_unitlen);
    bool get_query_pos_left(const std::string& bb, const std::string& seq, const UINT32 refb, const UINT32 refe, UINT32& qb, UINT32& qe);
    bool get_query_pos_right(const std::string& bb, const std::string& seq, const UINT32 refb, const UINT32 refe, UINT32& qb, UINT32& qe);
    void sw_diff(const std::string& bb, const std::string& seq, std::vector<UINT16>& diff);    
    
    // unit is reverseif given string is not forward (prefix); else forward (suffix)
    inline UINT32 longest_periodic_pref(std::string s, std::string& unit, bool is_fwd=true) {
        if (!is_fwd) {std::reverse(s.begin(),s.end());}
        UINT32 m = s.size();
        std::vector<UINT16> ft(m,0);
        UINT16 ln = 0; //length of pvs pref-suff
        UINT32 i=1; //
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
        INT per = 0;
        UINT32 plen = 0;
        int unit_len = 0;
        for (UINT32 i=1; i <= m; ++i) {per = i - ft[i-1]; if (per<i && (i%per==0)){plen=i; unit_len=per;}}
        unit = s.substr(0,unit_len);
        UINT final = plen;
        for (int i=0; i < unit_len-1; ++i) {if (s[plen+i]==s[i]) {++final;}else{break;}}
	    return final;
    }

    // Whether b supports a
    inline bool is_supported (const AlnType at, const std::string& a, const std::string& b) {
        UINT l1 = std::min(a.size(),b.size());
        UINT l2 = std::max(a.size(),b.size());
        if (at==AlnType::I && l1!=l2) {
            return false;
        }
        bool same = true;
        for (UINT i=0; i< l1; ++i) {
            if (at==AlnType::S) {
                if (a[a.size()-1-i]!=b[b.size()-1-i]) {
                    same = false;
                    break;
                }
            }
            else {
                if (a[i]!=b[i]) { // I or P
                    same = false;
                    break;
                }
            }            
        }
        return same;
    }

    // Returns first base after the run
    inline UINT32 extend_right(const std::string& modified, const std::string& unit, const UINT32 pos) {
        auto unit_len = unit.size();
        UINT32 runlen = 0;
        UINT32 endpos = pos;
        while (endpos < modified.size()) {
            if (modified[endpos]!= unit[runlen]) {break;}
            ++endpos;
            ++runlen;
            if (runlen>=unit_len){runlen=0;}
        }
        return endpos;
    }

    // Returns pos of the first base of the the run; unit is in reverse direction
    inline UINT32 extend_left(const std::string& modified, const std::string& unit, const UINT32 pos) {
        auto unit_len = unit.size();
        UINT32 runlen = 0;
        int endpos = pos;
        while (endpos >= 0) {
            if (modified[endpos]!= unit[runlen]) {
                break;
            }
            if (endpos==0) {return endpos;}
            --endpos;
            ++runlen;
            if (runlen>=unit_len){runlen=0;}
        }
        return UINT32(endpos+1);
    }

}; // Window
} // namespace hypo
