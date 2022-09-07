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
/** Defines the class Alignment.
 * It represents an alignment and contains functionality for kmer -support update; minimizer-support update, and arms finding
 */
#include <cmath>
#include "Alignment.hpp"
#include "Contig.hpp"
#include "MinimizerDeque.hpp"

namespace hypo
{
    
Alignment::Alignment(Contig& contig, bam1_t *hts_align, UINT8 norm_edit_th): is_valid(true), _anchor_pos(-2), _is_wind_crossed(false),_beg_windex(0),_end_windex(0), _left_clip(0), _right_clip(0) {
    initialise_pos(hts_align);
    UINT32 clen = contig._len;
    if (_rb >= clen || _re > clen) {
        _qname = std::string(bam_get_qname(hts_align), hts_align->core.l_qname);
        fprintf(stdout, "[Hypo::Alignment] Error: Alignment File error: Looks like the reference in the alignment file is different from the draft. Contig (%s): Read (%s): rb (%u): re (%u): clen (%u)\n",contig._name.c_str(), _qname.c_str(), _rb, _re, clen); 
        exit(1);
    }

    std::string name = std::string(bam_get_qname(hts_align), hts_align->core.l_qname);
    INT e = name.size();
    while (name[e-1]=='\0') {--e;}
    _qname = name.substr(0,e);

    copy_data(hts_align); // may unset is_valid
    /*
    auto nmp = bam_aux_get(hts_align,"NM");
    if (nmp) {
        INT64 edit_dist = bam_aux2i(nmp);
        UINT32 rlen = _re-_rb;
        auto norm_edit_dist = std::ceil(edit_dist*100/rlen);
        if (norm_edit_dist > norm_edit_th) { // edit dist above tolerance
            is_valid = false;
        }
    }
    */
   if (is_valid && norm_edit_th < 100){
        UINT ln = std::max(_qae,_re-_rb);
        double diff = ln - std::min(_qae,_re-_rb);
        auto err = std::ceil(diff*100/ln);
        if (err > norm_edit_th) {
                is_valid = false;
        }
   }   
}

Alignment::Alignment(std::unique_ptr<Alignment> const& aln1,std::unique_ptr<Alignment> const& aln2, UINT32 missing_len, UINT32 gap_len): is_valid(true), _anchor_pos(-2), _is_wind_crossed(false),_beg_windex(0),_end_windex(0), _left_clip(0), _right_clip(0) {
    //Assumes aln1 is left, aln2 is right and missing is in between
    // Discard the clip related information of the original alns;
    std::string missing_qual="";
    std::string missing_seq = "";
    if (aln1->_tail.get_seq_size() >= missing_len) {
        missing_seq = aln1->_tail.unpack().substr(0,missing_len);
        missing_qual=aln1->_tail_basequal.substr(0,missing_len);
    }
    else if (aln2->_head.get_seq_size() >= missing_len) {
        missing_seq = aln2->_head.unpack().substr(aln2->_head.get_seq_size()-missing_len);
        missing_qual=aln2->_head_basequal.substr(aln2->_head.get_seq_size()-missing_len);
    }
    if (missing_seq.empty()) {is_valid = false;}
    else {    
        _rb=aln1->_rb;
        _re=aln2->_re;
        _qname=aln1->_qname;
        _qab=aln1->_qab;
        _qae = (aln1->_qae+aln2->_qae+missing_len);
        _is_rev  =aln1->_is_rev;
        _mq = aln1->_mq;
        std::string combined = (aln1->_apseq.unpack() + missing_seq + aln2->_apseq.unpack());
        PackedSeq<2> ps (combined);
        if (ps.is_valid()) {
            _apseq = std::move(ps);
        }
        else {
            is_valid = false;
        }
        _basequal = aln1->_basequal + missing_qual + aln2->_basequal;

        // Save cigar (without changing HTS format)
        // Cigar data is encoded 4 bytes per CIGAR operation.
        _num_cigar_op = aln1->_num_cigar_op + aln2->_num_cigar_op + 2; // added del and ins op
        _cigar.insert(std::end(_cigar),std::begin(aln1->_cigar),std::end(aln1->_cigar));
        // Del the gaplen; insertthe missing len
        UINT32 del = ((gap_len <<BAM_CIGAR_SHIFT) | BAM_CDEL);
        UINT32 ins = ((missing_len <<BAM_CIGAR_SHIFT) | BAM_CINS);
        _cigar.emplace_back(del);
        _cigar.emplace_back(ins);
        _cigar.insert(std::end(_cigar),std::begin(aln2->_cigar),std::end(aln2->_cigar));
    }

}

Alignment::Alignment(std::unique_ptr<Alignment> const& aln1,std::unique_ptr<Alignment> const& aln2, UINT32 missing_len): is_valid(true), _anchor_pos(-2), _is_wind_crossed(false),_beg_windex(0),_end_windex(0), _left_clip(0), _right_clip(0) {
    //Assumes aln1 is left, aln2 is right and missing is in between
    // Discard the clip related information of the original alns;
    // Attaches missing seq as an INS in the first; last base of missing seq is treated as M/X in the corresponding cigar.
    std::string missing_qual="";
    std::string missing_seq = "";
    if (aln1->_tail.get_seq_size() >= missing_len) {
        missing_seq = aln1->_tail.unpack().substr(0,missing_len);
        missing_qual=aln1->_tail_basequal.substr(0,missing_len);
    }
    else if (aln2->_head.get_seq_size() >= missing_len) {
        missing_seq = aln2->_head.unpack().substr(aln2->_head.get_seq_size()-missing_len);
        missing_qual=aln2->_head_basequal.substr(aln2->_head.get_seq_size()-missing_len);
    }
    if (missing_seq.empty()) {is_valid = false;}
    else {    
        _rb=aln1->_rb;
        _re=aln1->_re+1;
        _qname=aln1->_qname;
        _qab=aln1->_qab;
        _qae = (aln1->_qae+missing_len);
        _is_rev  =aln1->_is_rev;
        _mq = aln1->_mq;
        std::string combined = (aln1->_apseq.unpack() + missing_seq);
        PackedSeq<2> ps (combined);
        if (ps.is_valid()) {
            _apseq = std::move(ps);
        }
        else {
            is_valid = false;
        }
        _basequal = aln1->_basequal + missing_qual;

        // Save cigar (without changing HTS format)
        // Cigar data is encoded 4 bytes per CIGAR operation.
        _num_cigar_op = aln1->_num_cigar_op + 2; // added ins and M op
        _cigar.insert(std::end(_cigar),std::begin(aln1->_cigar),std::end(aln1->_cigar));
        // Del the gaplen; insertthe missing len
        UINT32 mat = ((1 <<BAM_CIGAR_SHIFT) | BAM_CMATCH);
        UINT32 ins = (((missing_len-1) <<BAM_CIGAR_SHIFT) | BAM_CINS);
        _cigar.emplace_back(ins);
        _cigar.emplace_back(mat);
    }

}


void Alignment::update_solidkmers_support (const UINT k, Contig& contig) {
    /* Find kmer range */
    auto first = contig._Rsolid_pos(_rb);
    auto last = contig._Rsolid_pos(_re);
    // discard those which do not wholly fall in the alignment
    // NOTE SELECT query's index is 1-based
    for (auto i = last; i >first; --i) {
        auto pos = contig._Ssolid_pos(i);
        if (pos + k <= _re ) { // falls wholly within
            last = i;
            break;
        }
    }
    // Here first+1 th and last th pos fall in the read; index falling in the read: first and last-1; lst points to index of kmer outside (or not fully in)
    if  (last>first) { // some kmer found
        std::unordered_multimap<UINT64,UINT32> kmap;
        kmap.reserve(last-first);
        std::vector<UINT64> spos;
        spos.reserve(last-first);
        // update coverage; Add kid to map; add pos of solid kmer to spos
        for (auto i = first; i <last; ++i) {
            contig.increment_kmer_coverage(i);
            kmap.insert({contig._kmerinfo[i]->kid,i-first});
            spos.emplace_back(contig._Ssolid_pos(i+1));
        }
        // Find kmers in the read
        UINT64 kmer=0;
        UINT kmer_len=0;
        const UINT64 kmask = (1ULL<<2*k) - 1;
        size_t num_qbases = _apseq.get_seq_size();
        UINT32 num_cbases = _re-_rb;
        INT64 pvs_supp_kpos = -1; 
        UINT32 pvs_supp_r_bind = 0;
        for (size_t r_ind=0; r_ind < num_qbases; ++r_ind) { //r_ind is index of the end of kmer
            BYTE b = _apseq.enc_base_at(r_ind);
            //Always ACGT for reads
            kmer = (kmer << 2 | b) & kmask;
            if (kmer_len<k) { ++kmer_len;}   
            // Update support if conditions meet           
            if (kmer_len==k) {
                UINT32 r_bind = r_ind+1-k;
                auto kit = kmap.equal_range(kmer);
                for  (auto itr=kit.first; itr!=kit.second;++itr) {
                    UINT32 c_ind = itr->second;
                    INT64 c_dist = spos[c_ind] - _rb;
                    UINT32 srange_left = ((c_dist>k)?(c_dist-k):(0));
                    UINT32 srange_right = std::min(INT64(num_cbases),c_dist+(k));
                    // WARNING: The same read-kmer can provide support to many in the contig
                    // (if it is not unique and repeat close to each other). But such cases  are rare. So ignore.
                    if (r_bind >= srange_left && r_bind <= srange_right) { // supports
                        bool should_update = true;
                        /*
                        // Check for special case of pvs adjacent/overlapping nbr (This is just heuristics; to ignore some of the possible wrong cases)
                        if  (pvs_supp_kpos > -1 && (spos[c_ind]<=k+UINT64(pvs_supp_kpos))) { // spcl case: pvs kmer exists in the range and is overlapping/adjacent
                            if ( (r_bind-pvs_supp_r_bind) != (spos[c_ind]-pvs_supp_kpos)) { // insertion in reads between adjacent; do not support
                                should_update = false;
                            }
                        } 
                        */
                         // Check for special case of pvs adjacent/overlapping nbr (This is just heuristics; to ignore some of the possible wrong cases)
                        if  (pvs_supp_kpos > -1 ) { // spcl case: pvs kmer exists in the range
                            if (spos[c_ind]<=k+UINT64(pvs_supp_kpos)) { // pvs kmer is overlapping/adjacent
                                if ( (r_bind-pvs_supp_r_bind) != (spos[c_ind]-pvs_supp_kpos)) { // insertion in reads between adjacent; do not support
                                    should_update = false;
                                }
                            }
                            else { // a window between pvs and this (supported) kmer
                                if ( r_bind< k + pvs_supp_r_bind) { // insertion in reads between adjacent; do not support
                                    should_update = false;
                                }
                            }
                        } 
                        if (should_update) {
                            pvs_supp_kpos = spos[c_ind];
                            pvs_supp_r_bind = r_bind;
                            contig.increment_kmer_support(first+c_ind);
                        }                         
                    }
                }
            }
        }
    }
}

void Alignment::update_minimisers_support (Contig& contig) {
    UINT MINIMIZER_K = Minimizer_settings.k;
    UINT MINIMIZER_W = Minimizer_settings.w;
    /* Find starting and ending region */
    // Rank(i) returns number of set bits in [0, i). Therefore look for _rb+1 so that index of starting region can be deduced by subtracting 1.
    auto first = contig._RMreg_pos(_rb+1)-1;
    auto last = contig._RMreg_pos(_re);

    auto first_windex = ((contig._is_win_even && first%2==0) || (!contig._is_win_even && first%2==1)) ? (first) : (first+1);
    auto last_windex = ((contig._is_win_even && last%2==0) || (!contig._is_win_even && last%2==1)) ? (last) : (last-1);

    // Update minimser support
    if (last_windex >= first_windex) {
        // Find minimsers
        UINT32 last_found_position = _apseq.get_seq_size() + 1; //a unique identifier for 'first minimizer'
        
        //UINT32 shift = 2 * (MINIMIZER_K - 1);
        UINT32 mask = (1ULL<<2*MINIMIZER_K) - 1;
        UINT32 kmer[2] = {0,0};
        MinimizerDeque<UINT32> minimizer_window(MINIMIZER_W + 1);
        UINT32 count_not_N = 0;
        UINT32 processed_kmer = 0;
        UINT32 current_start_position = 0;
        
        //we have to make a vector of found minimizers and a map to validity, to handle duplicates
        //is there a better solution?
        std::unordered_multimap<UINT32,UINT32> found_minimizers; // maps read minimser to its position in the read
        
        for(size_t i = 0; i < _apseq.get_seq_size(); ++i) {
            BYTE c = _apseq.enc_base_at(i);
            if(c < 4) {
                ++count_not_N;                
                kmer[0] = (kmer[0] << 2ull | c) & mask;           // forward k-mer
                //kmer[1] = (kmer[1] >> 2ull) | (3ULL^c) << shift; // reverse k-mer
                //int z = kmer[0] < kmer[1] ? 0 : 1;
                int z = 0;
                if(count_not_N >= MINIMIZER_K) {
                    while(!minimizer_window.empty() && std::get<0>(minimizer_window.back()) > kmer[z]) minimizer_window.pop_back();
                    minimizer_window.push_back(std::make_tuple(kmer[z], i));
                    while(std::get<1>(minimizer_window.front()) + MINIMIZER_W <= i) minimizer_window.pop_front();
                    ++processed_kmer;
                    if(processed_kmer >= MINIMIZER_W) {
                        current_start_position = std::get<1>(minimizer_window.front()) - MINIMIZER_K + 1;
                        if(current_start_position != last_found_position) { //first minimizer
                            found_minimizers.insert({std::get<0>(minimizer_window.front()),current_start_position});
                        }
                        last_found_position = current_start_position;
                    }                    
                }
            } else {
                count_not_N = 0;
            }
        }

        // Check minimser in each window; incrementing support if it falls within reasonable range (2*mk on either side)
        UINT16 num_cbases = _re-_rb;
        for (UINT32 i=first_windex; i<=last_windex; i=i+2) {
            UINT32 minfoidx = (contig._is_win_even) ? (i/2) : ((i-1)/2);
            auto num_minimsers = contig._minimserinfo[minfoidx]->rel_pos.size();
            // NOTE SELECT query's index is 1-based
            auto minimiser_pos = contig._SMreg_pos(i+1);         
            for (UINT32 mi=0; mi< num_minimsers; ++mi) {
                minimiser_pos+=contig._minimserinfo[minfoidx]->rel_pos[mi];
                UINT32 c_dist = minimiser_pos - _rb;
                UINT32 range_left = ((c_dist>(2*MINIMIZER_K))?(c_dist-(2*MINIMIZER_K)):(0));
                UINT32 range_right = std::min(num_cbases,(UINT16)(c_dist+(3*MINIMIZER_K)));
                if (minimiser_pos>=_rb && minimiser_pos<_re) {
                    // coverage
                    contig.increment_minimser_coverage(minfoidx,mi);
                    // support
                    auto mit = found_minimizers.equal_range(contig._minimserinfo[minfoidx]->minimisers[mi]);
                    for  (auto itr=mit.first; itr!=mit.second;++itr) {
                       
                        // WARNING: The same read-kmer can provide support to many in the contig
                        // (if it is not unique and repeat close to each other). But such cases  are rare. So ignore.
                        if (itr->second >= range_left && itr->second <= range_right) { // supports
                            contig.increment_minimser_support(minfoidx,mi);
                        }
                    }
                }
                if (minimiser_pos>=_re) {
                    break;
                }
            }
        }
    }    
}
// For adding to short windows
void Alignment::find_short_arms(const UINT k, const Mode mode, Contig& contig) {
    /* Find the range */
    auto b_ind = contig._Rreg_pos(_rb);
    if (contig._reg_pos[_rb]==0) { // starts before
        --b_ind;
    }
    // e_ind will point to 1 past the last region falling in the read (i.e. starting pos of the next reg)
    auto e_ind = contig._Rreg_pos(_re); 
    // Ignore if the whole read falls into a single SR or window
    if (e_ind-b_ind > 1) { 
        // Update coverage
        
        if (mode==Mode::SECOND || mode==Mode::SO) {
            for (auto ind=b_ind+1; ind < e_ind-1; ++ind) {
                if (_mq>MAP_QUAL_TH && contig._reg_type[ind] != RegionType::SR && contig._reg_type[ind] != RegionType::MSR) { // a window
                    auto pindex = contig._reg_info[ind];
                    contig.increment_swindow_coverage(pindex);
                }
            }
        }
        
        _is_wind_crossed = true;
        _beg_windex = b_ind;
        _end_windex = e_ind;
        // Find breaking-points of the read; at least 1 will exist
        std::vector<std::string> cig_strings;
        std::vector<UINT32> bp = find_bp_with_cigar(contig._Sreg_pos,contig._reg_type,b_ind,e_ind, contig, cig_strings );
        _arms.reserve(bp.size());
        // first arm may be suffix or internal or in an SR
        ArmType armtype = (contig._reg_pos[_rb]==0) ? ArmType::SUFFIX: ArmType::INTERNAL;
        if (contig._reg_type[b_ind] != RegionType::SR && contig._reg_type[b_ind] != RegionType::MSR) { // a window
            prepare_short_arm(k, b_ind, _qab, bp[0], armtype, contig, cig_strings[0]);
            auto pindex = contig._reg_info[b_ind];
            if (_mq>MAP_QUAL_TH) {
                if (armtype==ArmType::INTERNAL) {contig.increment_swindow_coverage(pindex);} // first window lies wholly in the read
                else {contig.increment_swindow_scoverage(pindex);}
            }                      
        }
        // All others are internal arms (to be added) if in a window; Check for empty arm; If empty just incease the pointer
        UINT bp_ind = 0;
        for (auto ind=b_ind+1; ind < e_ind-1; ++ind,++bp_ind) {
            if (contig._reg_type[ind] != RegionType::SR && contig._reg_type[ind] != RegionType::MSR) { // a window
                if (bp[bp_ind+1] == bp[bp_ind]) {// Empty seq in this window
                    _arms.emplace_back(Arm(ind));
                }
                else {
                    prepare_short_arm(k, ind, bp[bp_ind], bp[bp_ind+1], ArmType::INTERNAL, contig, cig_strings[bp_ind+1]);
                }
            }
        }
        // last window may be suffix or internal or in an SR
        armtype = (contig._reg_pos[_re]==0) ? ArmType::PREFIX: ArmType::INTERNAL;
        if (contig._reg_type[e_ind-1] != RegionType::SR && contig._reg_type[e_ind-1] != RegionType::MSR) { // a window
            prepare_short_arm(k, e_ind-1,  bp[bp_ind], _qae, armtype, contig, cig_strings[bp_ind+1]);
            auto pindex = contig._reg_info[e_ind-1];
            if (_mq>MAP_QUAL_TH) {
                if (armtype==ArmType::INTERNAL) {contig.increment_swindow_coverage(pindex);} // first window lies wholly in the read
                else {contig.increment_swindow_pcoverage(pindex);}
            }            
        }
    }
    _arms.shrink_to_fit(); 
}

// For adding to long windows
void Alignment::find_long_arms(Contig& contig) {
    if (!is_valid) {return;}
    /* Find the range */
    auto b_ind = contig._Rreg_pos(_rb);
    if (contig._reg_pos[_rb]==0) { // starts before
        --b_ind;
    }
    // e_ind will point to 1 past the last region falling in the read (i.e. starting pos of the next reg)
    auto e_ind = contig._Rreg_pos(_re); 

    
    // Update coverage
    if (_mq>0) {
        for (auto ind=b_ind; ind < e_ind; ++ind) {
            if (contig._reg_type[ind] != RegionType::SR && contig._reg_type[ind] != RegionType::INVALID) { // a long window
                auto pindex = contig._reg_info[ind];
                contig.increment_lwindow_coverage(pindex);
            }
        }
    }
    
    
    // Ignore if the whole read falls into a single SR or window
    if (e_ind-b_ind > 1) { 
        // Find breaking-points of the read; at least 1 will exist
        std::vector<int> cigs;
        std::vector<UINT32> bp = find_bp_with_cigar_long(contig._Sreg_pos,contig._reg_type,b_ind,e_ind, cigs );
        for(auto & i : bp) std::cout << i << std::endl;
        _arms.reserve(bp.size());
        // first arm may be suffix or internal or in an SR
        ArmType armtype = (contig._reg_pos[_rb]==0) ? ArmType::SUFFIX: ArmType::INTERNAL;
        if (contig._reg_type[b_ind] != RegionType::SR && contig._reg_type[b_ind] != RegionType::INVALID) { // a valid window
            if ((bp[0]) >= Arms_settings.long_arm_min_len) {
                _arms.emplace_back(Arm(b_ind, _apseq, _qab, bp[0], armtype,_basequal,  std::to_string((cigs[0]+cigs[1]))));
            }            
        }
        // All others are internal arms (to be added) if in a window; Check for empty arm; If empty just incease the pointer
        UINT bp_ind = 0;
        for (auto ind=b_ind+1; ind < e_ind-1; ++ind,++bp_ind) {
            if (contig._reg_type[ind] != RegionType::SR && contig._reg_type[ind] != RegionType::INVALID) { // a valid window
                if (bp[bp_ind+1] == bp[bp_ind]) {// Empty seq in this window
                    _arms.emplace_back(Arm(ind));
                }
                else {
                    _arms.emplace_back(Arm(ind, _apseq,bp[bp_ind], bp[bp_ind+1], ArmType::INTERNAL,_basequal,  std::to_string((cigs[bp_ind]+cigs[bp_ind+1]+cigs[bp_ind+2]))));                   
                }
            }
        }
        // last window may be suffix or internal or in an SR
        armtype = (contig._reg_pos[_re]==0) ? ArmType::PREFIX: ArmType::INTERNAL;
        if (contig._reg_type[e_ind-1] != RegionType::SR && contig._reg_type[e_ind-1] != RegionType::INVALID) { // a valid window
            if (_qae-bp[bp_ind] >= Arms_settings.long_arm_min_len) {
                _arms.emplace_back(Arm(e_ind-1, _apseq, bp[bp_ind], _qae, armtype,_basequal,  std::to_string((cigs[bp_ind]+cigs[bp_ind+1]))));
            }            
        }
    } 
    _arms.shrink_to_fit(); 
}

void Alignment::add_short_arms(const Mode mode, const Contig& contig) {
    //std::cout << _qname <<" " << _arms.size()<<std::endl;
    // In LSA mode, arms are only for SWS windows
    for (const auto& a: _arms) {
        bool should_add = true;
        auto pindex = contig._reg_info[a.windex];
        
        auto pscov = std::min(contig._scovinfo[pindex]->pcoverage,contig._scovinfo[pindex]->scoverage);
        auto intcov = contig._scovinfo[pindex]->coverage;
        auto app_cov = std::max(pscov,intcov);
        // If coverage of window is too low, use mq=o; otherwise not
        if ( _mq  <= MAP_QUAL_TH && app_cov > Arms_settings.min_short_num) {
            should_add = false;
        }
        
        if (should_add) {
            UINT8 w = (_is_properly_paired)?(2):(1);
            if (a.armtype==ArmType::PREFIX) {
                contig._pswindows[pindex]->add_prefix(a.arm);
                contig._pswindows[pindex]->add_pref_end(a.terminal);
                contig._pswindows[pindex]->add_pref_wt(w);
                contig._pswindows[pindex]->add_pref_clip(_right_clip);
            }
            else if (a.armtype==ArmType::SUFFIX) {
                contig._pswindows[pindex]->add_suffix(a.arm);
                contig._pswindows[pindex]->add_suff_beg(a.terminal);
                contig._pswindows[pindex]->add_suff_wt(w);
                contig._pswindows[pindex]->add_suff_clip(_left_clip);
            }
            else if (a.armtype==ArmType::INTERNAL) {
                contig._pswindows[pindex]->add_internal_with_cigar(a.arm,a.cigstr);
                if (!a.is_high_qual) {contig._pswindows[pindex]->set_low_qual_internal();}
                contig._pswindows[pindex]->add_internal_wt(w);
            }
            else {
                contig._pswindows[pindex]->add_empty();
            }
        }
    }
    _arms.clear();
}

void Alignment::add_long_arms(const Contig& contig) {
    if (!is_valid) {return;}
    //std::cout << _qname <<" " << _arms.size()<<std::endl;    
    // Add only if the read spans a high-cov region=> starts and ends in a normal coverage
    bool should_add = true;
    auto hic_thres = LWindow_cov_settings.high_th;
    auto loc_thres = LWindow_cov_settings.low_th;
    auto clip_frac_th = Filtering_settings.clip_ratio;
    auto contig_end_threshold = Filtering_settings.contig_end_leniency;
    if (_arms.size() > 1) {
        auto bindex = contig._reg_info[_arms[0].windex]; // first
        auto eindex = contig._reg_info[_arms[_arms.size()-1].windex]; // last
        auto bcov = contig._lcovinfo[bindex]->coverage;
        auto ecov = contig._lcovinfo[eindex]->coverage;
        
        // only discard unnatural clip values
        double current_clip_fraction = (double)((_left_clip+_right_clip)) / double(_apseq.get_seq_size());
        
        // disregard clips, thus don't discard alignment, in case of reads being aligned to start or end of contig
        if (_rb > contig_end_threshold && _re < (contig._len - contig_end_threshold) && current_clip_fraction > clip_frac_th && bcov > hic_thres && ecov > hic_thres ) {
            should_add = false;
        }

        // If coverage of window is too low, use mq=o; otherwise not
        if ( _mq == 0 && bcov > loc_thres && ecov > loc_thres) {
            should_add = false;
        }

    }
    if (should_add) {
        for (const auto& a: _arms) {
            auto pindex = contig._reg_info[a.windex];
            if (a.armtype==ArmType::PREFIX) {
                contig._plwindows[pindex]->add_prefix(a.arm, a.basequal);
            }
            else if (a.armtype==ArmType::SUFFIX) {
                contig._plwindows[pindex]->add_suffix(a.arm, a.basequal);
            }
            else if (a.armtype==ArmType::INTERNAL) {
                contig._plwindows[pindex]->add_internal(a.arm, a.basequal);
                contig._plwindows[pindex]->add_internal_cig(a.cigstr);
                UINT alen = _re-_rb;
                if (_left_clip > CLIP_FRAC*alen || _right_clip > CLIP_FRAC*alen) {
                    alen=0;
                }
                contig._plwindows[pindex]->add_max_areflen(alen);
            }
            else {
                contig._plwindows[pindex]->add_empty(_re-_rb);
            }
        }
    }
    _arms.clear();
}

// beg_ind points to the starting pos of the window the read starts in; end_ind points to the starting pos of the window following the one in which read ends
std::vector<UINT32> Alignment::find_bp(const sdsl::bit_vector::select_1_type& Sreg_pos, const std::vector<RegionType>& reg_type, const UINT32 beg_ind, const UINT32 end_ind) {    
    std::vector<UINT32> results;
    UINT32 num_bp = end_ind-beg_ind-1; // ending and beginning themselves do not fall within the read; so no bp
    results.reserve(num_bp);
    UINT32 current_reference_pos = _rb;
    UINT32 current_processed_index = beg_ind + 1;
    UINT32 next_ref_pos = Sreg_pos(current_processed_index + 1);
    UINT32 current_query_pos = 0;

    UINT32 current_cigar;
    UINT32 get_op;
    UINT32 get_oplen;
    UINT32 len_diff;
    INT8 is_corner = 0;
    for(UINT32 i = 0; i < _num_cigar_op; i++) {        
        current_cigar = _cigar[i];
        get_op = bam_cigar_op(current_cigar);
        get_oplen = bam_cigar_oplen(current_cigar);
        
        if(get_op == BAM_CSOFT_CLIP || get_op == BAM_CHARD_CLIP) { //softclips are already discarded
            continue;
        }
        
        if((bam_cigar_type(get_op) & 3) == 3) { //bit 1 and bit 2 set, consumes query and reference
            if(is_corner) { //a corner but no insertions, so just insert the previous position
                results.push_back(current_query_pos);
                is_corner = 0;
                ++current_processed_index;
                next_ref_pos = Sreg_pos(current_processed_index + 1);
            }
            while(current_reference_pos + get_oplen >= next_ref_pos && !is_corner) {
                len_diff = next_ref_pos - current_reference_pos;
                current_reference_pos = next_ref_pos;
                current_query_pos += len_diff;
                get_oplen -= len_diff;
                if(get_oplen > 0) {
                    results.push_back(current_query_pos);
                    ++current_processed_index;
                    next_ref_pos = Sreg_pos(current_processed_index + 1);
                }
                else is_corner = 1;
            }
            if(get_oplen > 0) { // exhaust this cigar op
                current_reference_pos += get_oplen;
                current_query_pos += get_oplen;
            }
        } else if(bam_cigar_type(get_op) & 2) { //only bit 2 set, consumes reference
            if(is_corner) { //a corner but no insertions, so just insert the previous position
                results.push_back(current_query_pos);
                is_corner = 0;
                ++current_processed_index;
                next_ref_pos = Sreg_pos(current_processed_index + 1);
            }
            while(current_reference_pos + get_oplen >= next_ref_pos && !is_corner) {
                len_diff = next_ref_pos - current_reference_pos;
                current_reference_pos = next_ref_pos;
                get_oplen -= len_diff;
                if(get_oplen > 0) {
                    results.push_back(current_query_pos);
                    ++current_processed_index;
                    next_ref_pos = Sreg_pos(current_processed_index + 1);
                }
                else is_corner = 1;
            }
            if(get_oplen > 0) { // exhaust this cigar op
                current_reference_pos += get_oplen;
            }
        } else if(bam_cigar_type(get_op) & 1) { //only bit 1 set, consumes query
            if(is_corner) { //this is a corner, now add the position
                if(reg_type[current_processed_index-1] == RegionType::SR || reg_type[current_processed_index-1] == RegionType::MSR) { // if this is SR, we want the current_query_pos to be included in the right window
                    results.push_back(current_query_pos);
                } else { // otherwise we include it here
                    results.push_back(current_query_pos + get_oplen);
                }
                ++current_processed_index;
                next_ref_pos = Sreg_pos(current_processed_index + 1);
                is_corner = 0;
            }
            current_query_pos += get_oplen;
        }
        
        if(current_processed_index == end_ind) break; //this should not be needed if end_ind is correct
    }

    return results;
}

// beg_ind points to the starting pos of the window the read starts in; end_ind points to the starting pos of the window following the one in which read ends
std::vector<UINT32> Alignment::find_bp_with_cigar_long(const sdsl::bit_vector::select_1_type& Sreg_pos, const std::vector<RegionType>& reg_type, const UINT32 beg_ind, const UINT32 end_ind, std::vector<int> & cigar_cuts) {    
    std::vector<UINT32> results;
    UINT32 num_bp = end_ind-beg_ind-1; // ending and beginning themselves do not fall within the read; so no bp
    results.reserve(num_bp);
    UINT32 current_reference_pos = _rb;
    UINT32 current_processed_index = beg_ind + 1;
    UINT32 next_ref_pos = Sreg_pos(current_processed_index + 1);
    UINT32 current_query_pos = 0;

    UINT32 current_cigar;
    UINT32 get_op;
    UINT32 get_oplen;
    UINT32 orig_oplen;
    UINT32 len_diff;
    INT8 is_corner = 0;

    INT64 extra = 0;
    
    
    for(UINT32 i = 0; i < _num_cigar_op; i++) {        
        current_cigar = _cigar[i];
        get_op = bam_cigar_op(current_cigar);
        get_oplen = bam_cigar_oplen(current_cigar);
        
        // IT is MUST to ignore soft/hard clips for the sake of gapped long-del handling
        if(get_op == BAM_CSOFT_CLIP || get_op == BAM_CHARD_CLIP) { //softclips are already discarded
            continue;
        }
        
        if((bam_cigar_type(get_op) & 3) == 3) { //bit 1 and bit 2 set, consumes query and reference
            //original_cigar += std::to_string(get_oplen) + "M";
            
            if(is_corner) { //a corner but no insertions, so just insert the previous position
                results.push_back(current_query_pos);
                cigar_cuts.push_back(extra);
                extra=0;
                is_corner = 0;
                ++current_processed_index;
                next_ref_pos = Sreg_pos(current_processed_index + 1);
            }
            while(current_reference_pos + get_oplen >= next_ref_pos && !is_corner) {
                len_diff = next_ref_pos - current_reference_pos;
                current_reference_pos = next_ref_pos;
                current_query_pos += len_diff;
                get_oplen -= len_diff;
                if(get_oplen > 0) {
                    results.push_back(current_query_pos);
                    cigar_cuts.push_back(extra);
                    extra=0;
                    ++current_processed_index;
                    next_ref_pos = Sreg_pos(current_processed_index + 1);
                }
                else {is_corner = 1;}
            }
            if(get_oplen > 0) { // exhaust this cigar op
                current_reference_pos += get_oplen;
                current_query_pos += get_oplen;
            }
        } else if(bam_cigar_type(get_op) & 2) { //only bit 2 set, consumes reference
            //original_cigar += std::to_string(get_oplen) + "D";
            
            if(is_corner) { //a corner but no insertions, so just insert the previous position
                results.push_back(current_query_pos);
                cigar_cuts.push_back(extra);
                extra=0;
                is_corner = 0;
                ++current_processed_index;
                next_ref_pos = Sreg_pos(current_processed_index + 1);
            }
            if(get_oplen >= INDEL_CLUSTER_THRESHOLD) {
                //add_cigar += std::to_string(to_most_significant(get_oplen)) + "D";
                extra -= (INT64(get_oplen/10));
            }
            
            orig_oplen = get_oplen;
            while(current_reference_pos + get_oplen >= next_ref_pos && !is_corner) {
                len_diff = next_ref_pos - current_reference_pos;
                current_reference_pos = next_ref_pos;
                get_oplen -= len_diff;
                if(get_oplen > 0) {
                    results.push_back(current_query_pos);
                    cigar_cuts.push_back(extra);
                    extra=0;
                    if(orig_oplen >= INDEL_CLUSTER_THRESHOLD) {
                        //add_cigar += std::to_string(to_most_significant(orig_oplen)) + "D";
                        extra -= (INT64(orig_oplen/10));
                    }
                    ++current_processed_index;
                    next_ref_pos = Sreg_pos(current_processed_index + 1);
                }
                else {is_corner = 1;}
            }
            if(get_oplen > 0) { // exhaust this cigar op
                current_reference_pos += get_oplen;
            }
        } else if(bam_cigar_type(get_op) & 1) { //only bit 1 set, consumes query
            //original_cigar += std::to_string(get_oplen) + "I";
            
            if(is_corner) { //this is a corner, now add the position
                if(reg_type[current_processed_index-1] == RegionType::SR || reg_type[current_processed_index-1] == RegionType::MSR) { // if this is SR, we want the current_query_pos to be included in the right window
                    results.push_back(current_query_pos);
                    cigar_cuts.push_back(extra);
                    extra=0;
                    if(get_oplen >= INDEL_CLUSTER_THRESHOLD) {
                        //add_cigar += std::to_string(to_most_significant(get_oplen)) + "I";
                        extra += (INT64(get_oplen/10));
                    }
                } else { // otherwise we include it here
                    results.push_back(current_query_pos + get_oplen);
                    if(get_oplen >= INDEL_CLUSTER_THRESHOLD) {
                        //add_cigar += std::to_string(to_most_significant(get_oplen)) + "I";
                        extra += (INT64(get_oplen/10));
                    }
                    cigar_cuts.push_back(extra);
                    extra=0;
                }
                ++current_processed_index;
                next_ref_pos = Sreg_pos(current_processed_index + 1);
                is_corner = 0;
            } else {
                if(get_oplen >= INDEL_CLUSTER_THRESHOLD) {
                    //add_cigar += std::to_string(to_most_significant(get_oplen)) + "I";
                    extra += (INT64(get_oplen/10));
                }
            }
            current_query_pos += get_oplen;
        }
        
        if(current_processed_index == end_ind) {break;} //this should not be needed if end_ind is correct
    }
    if(is_corner) { //a corner but no insertions, so just insert the previous position
        results.push_back(current_query_pos);
        cigar_cuts.push_back(extra);
        extra=0;
    }
    cigar_cuts.push_back(extra);
    return results;
}


// beg_ind points to the starting pos of the window the read starts in; end_ind points to the starting pos of the window following the one in which read ends
std::vector<UINT32> Alignment::find_bp_with_cigar(const sdsl::bit_vector::select_1_type& Sreg_pos, const std::vector<RegionType>& reg_type, const UINT32 beg_ind, const UINT32 end_ind, const Contig & contig, std::vector<std::string> & cigar_cuts) {
    std::vector<UINT32> results;
    UINT32 num_bp = end_ind-beg_ind-1; // ending and beginning themselves do not fall within the read; so no bp
    results.reserve(num_bp);
    UINT32 current_reference_pos = _rb;
    UINT32 current_processed_index = beg_ind + 1;
    UINT32 next_ref_pos = Sreg_pos(current_processed_index + 1);
    UINT32 current_query_pos = 0;

    UINT32 current_cigar;
    UINT32 get_op;
    UINT32 get_oplen;
    UINT32 len_diff;
    INT8 is_corner = 0;
    
    // cigar op parts of each section
    std::vector<std::vector< std::tuple<UINT32, UINT32> > > cigar_parts(num_bp + 1);
    UINT32 cigar_parts_idx = 0;
    
    INT8 finished_breaking = 0;
    
    for(UINT32 i = 0; i < _num_cigar_op; i++) {        
        current_cigar = _cigar[i];
        get_op = bam_cigar_op(current_cigar);
        get_oplen = bam_cigar_oplen(current_cigar);
        
        if(get_op == BAM_CSOFT_CLIP || get_op == BAM_CHARD_CLIP) { //softclips are already discarded
            continue;
        }
        
        if(finished_breaking) { // once we get to end_ind, just add all the cigar ops
            cigar_parts[cigar_parts_idx].push_back(std::make_tuple(get_op, get_oplen));
            continue;
        }
        
        if((bam_cigar_type(get_op) & 3) == 3) { //bit 1 and bit 2 set, consumes query and reference
            if(is_corner) { //a corner but no insertions, so just insert the previous position
                results.push_back(current_query_pos);
                is_corner = 0;
                ++current_processed_index;
                next_ref_pos = Sreg_pos(current_processed_index + 1);
                ++cigar_parts_idx;
            }
            while(current_reference_pos + get_oplen >= next_ref_pos && !is_corner) {
                len_diff = next_ref_pos - current_reference_pos;
                current_reference_pos = next_ref_pos;
                current_query_pos += len_diff;
                get_oplen -= len_diff;
                cigar_parts[cigar_parts_idx].push_back(std::make_tuple(get_op, len_diff));
                if(get_oplen > 0) {
                    results.push_back(current_query_pos);
                    ++current_processed_index;
                    next_ref_pos = Sreg_pos(current_processed_index + 1);
                    ++cigar_parts_idx;
                }
                else {is_corner = 1;}
            }
            if(get_oplen > 0) { // exhaust this cigar op
                current_reference_pos += get_oplen;
                current_query_pos += get_oplen;
                cigar_parts[cigar_parts_idx].push_back(std::make_tuple(get_op, get_oplen));
            }
        } else if(bam_cigar_type(get_op) & 2) { //only bit 2 set, consumes reference
            if(is_corner) { //a corner but no insertions, so just insert the previous position
                results.push_back(current_query_pos);
                is_corner = 0;
                ++current_processed_index;
                next_ref_pos = Sreg_pos(current_processed_index + 1);
                ++cigar_parts_idx;
            }
            while(current_reference_pos + get_oplen >= next_ref_pos && !is_corner) {
                len_diff = next_ref_pos - current_reference_pos;
                current_reference_pos = next_ref_pos;
                get_oplen -= len_diff;
                cigar_parts[cigar_parts_idx].push_back(std::make_tuple(get_op, len_diff));
                if(get_oplen > 0) {
                    results.push_back(current_query_pos);
                    ++current_processed_index;
                    next_ref_pos = Sreg_pos(current_processed_index + 1);
                    ++cigar_parts_idx;
                }
                else {is_corner = 1;}
            }
            if(get_oplen > 0) { // exhaust this cigar op
                current_reference_pos += get_oplen;
                cigar_parts[cigar_parts_idx].push_back(std::make_tuple(get_op, get_oplen));
            }
        } else if(bam_cigar_type(get_op) & 1) { //only bit 1 set, consumes query
            if(is_corner) { //this is a corner, now add the position
                if(reg_type[current_processed_index-1] == RegionType::SR || reg_type[current_processed_index-1] == RegionType::MSR) { // if this is SR, we want the current_query_pos to be included in the right window
                    results.push_back(current_query_pos);
                    cigar_parts[cigar_parts_idx+1].push_back(std::make_tuple(get_op, get_oplen));
                } else { // otherwise we include it here
                    results.push_back(current_query_pos + get_oplen);
                    cigar_parts[cigar_parts_idx].push_back(std::make_tuple(get_op, get_oplen));
                }
                ++current_processed_index;
                next_ref_pos = Sreg_pos(current_processed_index + 1);
                is_corner = 0;
                ++cigar_parts_idx;
            } else {cigar_parts[cigar_parts_idx].push_back(std::make_tuple(get_op, get_oplen));}
            current_query_pos += get_oplen;
        }
        
        if(current_processed_index >= end_ind) {finished_breaking = 1;} // once we get to end_ind, just add all the cigar ops
    }
    
    // go through results to construct the adjusted cigars
    current_reference_pos = _rb;
    current_query_pos = 0;
    
    for(UINT32 i = 0; i < results.size(); i++) {
        std::string constructed_string = "";
        UINT32 current_match = 0;
        for(UINT32 j = 0; j < cigar_parts[i].size(); j++) {
            get_op = std::get<0>(cigar_parts[i][j]);
            get_oplen = std::get<1>(cigar_parts[i][j]);
            if((bam_cigar_type(get_op) & 3) == 3) { // consume query and reference: check for mismatches
                current_match = 0;
                for(UINT32 k = 0; k < get_oplen; k++) {
                    if(contig._pseq.base_at(current_reference_pos) == _apseq.base_at(current_query_pos)) { // equal, extend the M
                        current_match++;
                    } else {
                        if(current_match > 0) constructed_string += std::to_string(current_match) + "M";
                        constructed_string += _apseq.base_at(current_query_pos);
                        current_match = 0;
                    }
                    current_reference_pos++;
                    current_query_pos++;
                }
                
                // remainding matches
                if(current_match > 0) constructed_string += std::to_string(current_match) + "M";
            } else if(bam_cigar_type(get_op) & 2) { // consume reference only: deletions. Just add the cigar
                constructed_string += std::to_string(get_oplen) + "D";
                current_reference_pos += get_oplen;
            } else if(bam_cigar_type(get_op) & 1) { // consume query only: insertions. Add the query bases
                for(UINT32 k = 0; k < get_oplen; k++) {
                    constructed_string += std::tolower(_apseq.base_at(current_query_pos));
                    current_query_pos++;
                }
                
            }
        }
        assert(current_query_pos == results[i]);
        cigar_cuts.push_back(constructed_string);
    }
    
    // last positions
    std::string constructed_string = "";
    UINT32 current_match = 0;
    UINT32 i = num_bp;
    for(UINT32 j = 0; j < cigar_parts[i].size(); j++) {
        get_op = std::get<0>(cigar_parts[i][j]);
        get_oplen = std::get<1>(cigar_parts[i][j]);
        if((bam_cigar_type(get_op) & 3) == 3) { // consume query and reference: check for mismatches
            current_match = 0;
            for(UINT32 k = 0; k < get_oplen; k++) {
                if(contig._pseq.base_at(current_reference_pos) == _apseq.base_at(current_query_pos)) { // equal, extend the M
                    current_match++;
                } else {
                    if(current_match > 0) constructed_string += std::to_string(current_match) + "M";
                    constructed_string += _apseq.base_at(current_query_pos);
                    current_match = 0;
                }
                current_reference_pos++;
                current_query_pos++;
            }
            
            // remainding matches
            if(current_match > 0) constructed_string += std::to_string(current_match) + "M";
        } else if(bam_cigar_type(get_op) & 2) { // consume reference only: deletions. Just add the cigar
            constructed_string += std::to_string(get_oplen) + "D";
            current_reference_pos += get_oplen;
        } else if(bam_cigar_type(get_op) & 1) { // consume query only: insertions. Add the query bases
            for(UINT32 k = 0; k < get_oplen; k++) {
                constructed_string += std::tolower(_apseq.base_at(current_query_pos));
                current_query_pos++;
            }
            
        }
    }
    cigar_cuts.push_back(constructed_string);
    
    return results;
}


bool Alignment::prepare_short_arm(const UINT k, const UINT32 windex, const UINT32 qb, const UINT32 qe, const ArmType origarmtype, Contig& contig, std::string& cigstr) {
    // Ignore if the arm is too short as compared to window
    auto mk = Minimizer_settings.k;
    auto curr_pos = contig._Sreg_pos(windex+1);
    auto next_pos = contig._Sreg_pos(windex+2);
    UINT32 window_len = next_pos - curr_pos;
    bool part_anc = false;
    ArmType armtype = origarmtype;
    //std:: cout << "wind: qb: qe : " << std::dec<< windex << " " << qb << " " << qe <<" "<<curr_pos<<" "<<next_pos<<std::endl;
    

    UINT32 terminal = 0;
    if (armtype==ArmType::PREFIX) {terminal = _re-curr_pos;}
    else if (armtype==ArmType::SUFFIX) {terminal = _rb-curr_pos;}
    
    // pref/suff should be sufficiently long to be a useful arm 
    if ((armtype!=ArmType::INTERNAL) && (next_pos-curr_pos > (Arms_settings.short_arm_coef*(qe-qb)))) {return true;}    
    
    RegionType wtype = contig._reg_type[windex];
    bool valid = true;
    auto q_beg = qb;
    auto q_end = qe;
    auto hq_th = Arms_settings.base_qual_th;
    
    // Find preceding anchor SR kmer; range =+/- k on either side of expected
    if ((wtype==RegionType::SWS || wtype==RegionType::SW || wtype==RegionType::SWM) && origarmtype!=ArmType::SUFFIX) { // preceded by SR and int or pre type
        // pa->_qab is always 0 (as clipped bases are discarded)
        {
            bool found = false;
            if (q_beg >= k) {
                // Get preceeding kmer ID
                UINT32 prec_SR_rank = contig._reg_info[windex-1];
                UINT64 anchor_kmer = contig._anchor_kmers[(prec_SR_rank<<1)]; // index of last kmer of SR is 2i
                UINT32 expected = q_beg-k;
                found = _apseq.check_kmer(anchor_kmer,k,expected);
            
                // Check if the kmer is at expected position. If not, check if kmer exists in tolerable range; If yes, Find its rightmost position
                if (!found) {
                    UINT32 search_start = (q_beg < (2*k)) ? (0) : (q_beg - (2*k));
                    UINT32 search_end = (q_end < (q_beg+k)) ? (q_end) : (q_beg+k);
                    size_t kmer_ind=0;
                    if (_apseq.find_kmer(anchor_kmer,k, search_start,search_end,expected,false,kmer_ind)) { // found
                        q_beg = kmer_ind + k;
                        found = true;
                    }
                } 
                // Probably a long INS/DEL; Try harder
                if (!found) {
                    expected = get_query_pos(curr_pos-k);
                    if (expected!=q_beg-k && expected <=q_end && expected +k <=_qae ) {
                        found = _apseq.check_kmer(anchor_kmer,k,expected);
                        if (!found) {
                            UINT32 search_start = (expected < (k)) ? (0) : (expected - (k));
                            UINT32 search_end = (q_end < (expected+2*k)) ? (q_end) : (expected+2*k);
                            size_t kmer_ind=0;
                            //std::cout << "A " << expected << " "<<search_start << " " << search_end << " " << _qae<<std::endl;
                            if (_apseq.find_kmer(anchor_kmer,k, search_start,search_end,expected,false,kmer_ind)) { // found
                                q_beg = kmer_ind + k;
                                found = true;
                            }
                        }
                        else {q_beg=expected+k;}
                    }
                } 
            } 
            if (!found) { // preceeding kmer does not exist; 
                if (origarmtype == ArmType::PREFIX || window_len <2 || part_anc) {valid=false;}
                else {
                    part_anc = true;
                    armtype = ArmType::SUFFIX;
                }
            }
        }
    }
    // Find succeeding anchor SR kmer; range =+/- k on either side of expected
    if ((wtype==RegionType::SWS || wtype==RegionType::WS || wtype==RegionType::MWS) && origarmtype!=ArmType::PREFIX) { // succeeded by SR and int or suff type
        // pa->_qab is always 1 (as clipped bases rae discarded)
        {            
            bool found = false;
            if (q_end + k <= _qae) {
                // Get succeeding kmer ID
                UINT32 succ_SR_rank = contig._reg_info[windex+1];
                UINT64 anchor_kmer = contig._anchor_kmers[(succ_SR_rank<<1)-1]; // index of first kmer of SR is 2i-1
                found = _apseq.check_kmer(anchor_kmer,k,q_end);
                UINT32 expected = q_end;
                // Check if the kmer is at expected position. If not, check if kmer exists in tolerable range; If yes, Find its leftmost position
                if (!found) {
                    UINT32 search_start = (q_end < (q_beg+k)) ? (q_beg) : (q_end-k);
                    UINT32 search_end = std::min(_qae,q_end+(2*k));
                    size_t kmer_ind=0;
                    if (_apseq.find_kmer(anchor_kmer, k,search_start,search_end,expected,true,kmer_ind)) { // found
                        q_end = kmer_ind;
                        found = true;
                    }
                }
                // Probably a long INS/DEL; Try harder
                if (!found) {
                    expected = get_query_pos(next_pos);
                    if (expected !=q_end && expected >=q_beg && expected +k <=_qae ) {
                        found = _apseq.check_kmer(anchor_kmer,k,expected);
                        if (!found) {
                            UINT32 search_start = (expected < (q_beg+k)) ? (q_beg) : (expected - (k));
                            UINT32 search_end = std::min(_qae,expected+(2*k));
                            size_t kmer_ind=0;
                            //std::cout << "B " << expected << " "<<search_start << " " << search_end << " " << _qae<<std::endl;
                            if (_apseq.find_kmer(anchor_kmer,k, search_start,search_end,expected,true,kmer_ind)) { // found
                                q_end = kmer_ind;
                                found = true;
                            }
                        }
                        else {q_end=expected;}
                    }
                } 
            }
            if (!found) { // succeeding kmer does not exist
                if (origarmtype == ArmType::SUFFIX|| window_len <2 || part_anc) {valid=false;}
                else {
                    part_anc = true;
                    armtype = ArmType::PREFIX;
                }
            }
        }
    }
    // Find preceding anchor minimiser kmer; range =+/- 2*mk on either side of expected
    if ((wtype==RegionType::MWM || wtype==RegionType::MW || wtype==RegionType::MWS) && origarmtype!=ArmType::SUFFIX) { // preceded by MINI and int or pre type
        // pa->_qab is always 0 (as clipped bases are discarded)
        {
            /* Left minimser should reset anchor pos;It is for the right to set which is to be used by the following window; 
            If the left does not reset, it may be used by any other window down the line. */
            // Use the pvs minimser-kmer as the anchor (We may lose the read-segment starting exactly at the minimser; but such cases should be rare)
            // Check if the kmer is found in previous window. If not, this read-segment is unusable
            if (_anchor_pos==-2) { // not yet set
                bool found = false;
                if (q_beg >= mk) {
                    UINT32 prec_MIN = contig._reg_info[windex-1];
                    // Check if the kmer is at expected position. If not, check if kmer exists in tolerable range; If yes, Find its leftmost position
                    found = _apseq.check_kmer(prec_MIN,mk,q_beg-mk);
                    UINT32 expected = q_beg-mk;
                
                    if (!found) {
                        UINT32 search_start = (q_beg < (3*mk)) ? (0) : (q_beg - (3*mk));
                        UINT32 search_end = (q_end < (q_beg+(2*mk))) ? (q_end) : (q_beg+(2*mk));
                        size_t kmer_ind=0;
                        if (_apseq.find_kmer(prec_MIN,mk, search_start,search_end,expected,false,kmer_ind)) { // found
                            q_beg = kmer_ind + mk;
                            found = true;
                            //_anchor_pos = kmer_ind;
                        }
                    }
                    // Probably a long INS/DEL; Try harder
                    if (!found) {
                        expected = get_query_pos(curr_pos-mk);
                        if (expected!=q_beg-mk && expected <=q_end && expected +mk <=_qae ) {
                            found = _apseq.check_kmer(prec_MIN,mk,expected);
                            if (!found) {
                                UINT32 search_start = (expected < (2*mk)) ? (0) : (expected - (2*mk));
                                UINT32 search_end = (q_end < (expected+3*mk)) ? (q_end) : (expected+3*mk);
                                size_t kmer_ind=0;
                                //std::cout << "C " << expected <<" "<< search_start << " " << search_end << " " << _qae<<std::endl;
                                if (_apseq.find_kmer(prec_MIN,mk, search_start,search_end,expected,false,kmer_ind)) { // found
                                    q_beg = kmer_ind + mk;
                                    found = true;
                                }
                            }
                            else {q_beg=expected+mk;}
                        }
                    } 
                }
                if (!found) { // preceeding kmer does not exist; do base by base match
                    if (origarmtype == ArmType::PREFIX|| window_len <2 || part_anc) {valid=false;}
                    else {
                        part_anc = true;armtype = ArmType::SUFFIX;
                    }
                }
            }
            else if (_anchor_pos==-1) { // set but not found
                if (origarmtype == ArmType::PREFIX|| window_len <2 || part_anc) {valid=false;}
                else {
                    part_anc = true;armtype = ArmType::SUFFIX;
                }
            }
            else {
                q_beg = _anchor_pos + mk;
            }
            _anchor_pos = -2;
        }
    }
    // Find succeeding minimiser kmer; range =+/- 2*mk on either side of expected
    if ((wtype==RegionType::MWM || wtype==RegionType::WM || wtype==RegionType::SWM) && origarmtype!=ArmType::PREFIX) { // succeeded by MINI and int or suff type
        // pa->_qab is always 1 (as clipped bases rae discarded)
        /* Set anchor posonly if the following window is going to use and later resett it; 
          This is to avoid using same anchor by another wnidow down the line
        */
        RegionType nwtype = contig._reg_type[windex+1]; // windex +1 always existsfor a window because of dummy SR in the end
        bool update_anc = (nwtype==RegionType::MWM || nwtype==RegionType::MW || nwtype==RegionType::MWS);
        if (update_anc) {_anchor_pos = q_end;}
        {
            bool found = false;
            if (q_end + mk <= _qae) {
                // Get succeeding kmer ID
                UINT32 succ_MIN = contig._reg_info[windex+1];
                UINT32 expected = q_end;
                // Check if the kmer is at expected position. If not, check if kmer exists in tolerable range; If yes, Find its leftmost position
                found = _apseq.check_kmer(succ_MIN,mk,q_end);               
                if (!found) {
                    UINT32 search_start = (q_end < (q_beg+(2*mk))) ? (q_beg) : (q_end-(2*mk));
                    UINT32 search_end = std::min(_qae,q_end+(3*mk));
                    size_t kmer_ind=0;
                    if (_apseq.find_kmer(succ_MIN, mk,search_start,search_end,expected,true,kmer_ind)) { // found
                        q_end = kmer_ind;
                        found = true;
                    }
                }
                // Probably a long INS/DEL; Try harder
                if (!found) {
                    expected = get_query_pos(next_pos);
                    if (expected!=q_end && expected >=q_beg && expected +mk <=_qae ) {
                        found = _apseq.check_kmer(succ_MIN,mk,expected);
                        if (!found) {
                            UINT32 search_start = (expected < (q_beg+2*mk)) ? (q_beg) : (expected - (2*mk));
                            UINT32 search_end = std::min(_qae,expected+(3*mk));
                            size_t kmer_ind=0;
                            //std::cout << "D " << expected <<" "<< search_start << " " << search_end << " " << _qae<<std::endl;
                            if (_apseq.find_kmer(succ_MIN,mk, search_start,search_end,expected,true,kmer_ind)) { // found
                                q_end = kmer_ind;
                                found = true;
                            }
                        }
                        else {q_end=expected;}
                    }
                } 
            }
            if (!found) { // succeeding kmer does not exist
                if (armtype == ArmType::SUFFIX|| window_len <2 || part_anc) {valid=false;}
                else {
                    part_anc = true; armtype = ArmType::PREFIX;
                }
                if (update_anc){_anchor_pos = -1;}
            }
            else {
                if (update_anc){_anchor_pos = q_end;}
            }
        }
    }
    if (valid) {
        if (q_beg < q_end){
            if (part_anc) {
                if (armtype == ArmType::SUFFIX) {++q_beg; terminal=1;}
                else {
                    --q_end; terminal = window_len-1;
                }
            }
        
            if (q_beg < q_end){
                UINT32 num_lq = 0;
                for (auto q: _basequal.substr(q_beg,q_end-q_beg)) {
                    if (static_cast<UINT>(q) <  hq_th+33) {++num_lq;}
                }
                _arms.emplace_back(Arm(windex, _apseq, q_beg,q_end,armtype,num_lq<=Arms_settings.allowed_num_lq_bases, cigstr));
                if (armtype==ArmType::PREFIX || armtype==ArmType::SUFFIX) {_arms[_arms.size()-1].terminal = terminal;}
            }
        }
        else {
            _arms.emplace_back(Arm(windex));
        }
    }
    return false;
}


void Alignment::initialise_pos(const bam1_t *hts_align) {
    //_qname = std::string(bam_get_qname(hts_align), hts_align->core.l_qname);
    _rb = hts_align->core.pos;
    _qab = 0;
    // Get ending pos and update _qab (if soft-clipped)
    
    UINT32 curr_qp = _qab;
    UINT32 curr_rp = _rb;
    bool clip_before = true;
    UINT32 * pcigar = bam_get_cigar(hts_align);
    _num_cigar_op = hts_align->core.n_cigar;
    UINT32 count_clip_end = 0;
       
    for(UINT32 i = 0; i < _num_cigar_op; ++i) {
        UINT32 current_cigar = pcigar[i];
        UINT32 cigar_op = bam_cigar_op(current_cigar);
        UINT32 cigar_oplen = bam_cigar_oplen(current_cigar);
        
        if(cigar_op == BAM_CSOFT_CLIP || cigar_op == BAM_CHARD_CLIP) {
            if (i==0) {_left_clip=cigar_oplen;}
            else {_right_clip=cigar_oplen;}
        }
        
        if(clip_before) {
            if(cigar_op == BAM_CSOFT_CLIP) {
                _qab += cigar_oplen;
            } else if(cigar_op != BAM_CHARD_CLIP) {
                clip_before = false;
            }
        }        
        if((bam_cigar_type(cigar_op) & 3) == 3) { //bit 1 and bit 2 set, consumes query and reference
            curr_rp += cigar_oplen;
            curr_qp += cigar_oplen;
        } else if(bam_cigar_type(cigar_op) & 2) { //bit 2 set, consumes reference
            curr_rp += cigar_oplen;
        } else if(bam_cigar_type(cigar_op) & 1) { //bit 1 set, consumes query
            if(!clip_before && cigar_op == BAM_CSOFT_CLIP) {count_clip_end += cigar_oplen;}
            curr_qp += cigar_oplen;
        }
    }
    _qae = curr_qp - count_clip_end;
    _re = curr_rp;     
}

void Alignment::copy_data(const bam1_t *hts_align) {
    //_qname = std::string(bam_get_qname(hts_align), hts_align->core.l_qname);
    copy_head_tail(hts_align);
    _mq = hts_align->core.qual;
    UINT32 * pcigar = bam_get_cigar(hts_align);
    _is_rev = bam_is_rev(hts_align);
    _is_properly_paired = hts_align->core.flag & BAM_FPROPER_PAIR;
    // Get seq
    UINT32 qlen = _qae-_qab;
    UINT32 offset = _qab;
    PackedSeq<2> ps (qlen,offset,bam_get_seq(hts_align));
    if (ps.is_valid()) {
        _apseq = std::move(ps);
        // Reset query indices in accordance with only the stored (aligned) part; Clipped part already discarded
        _qae = _qae-_qab;
        _qab = 0;
 
        // Save cigar (without changing HTS format)
        // Cigar data is encoded 4 bytes per CIGAR operation.
        //UINT64 cigar_data_size = _num_cigar_op << 2; // multiply by 4
        UINT64 cigar_data_size = _num_cigar_op ;
        _cigar.reserve(cigar_data_size);
        _cigar.insert(_cigar.begin(),pcigar,pcigar+cigar_data_size); // copy not move; else htslib will destroy it  
        // Quality score; Phred scale with no +33 offset
        UINT8* pbasequal = bam_get_qual(hts_align);
        if (pbasequal!=nullptr)  {
            _basequal.resize(qlen);                
            for ( UINT32 i = 0; i < qlen; ++i ) {
                    _basequal[i] = static_cast<char>(pbasequal[i+offset]+33);
            }
        } 
    }
    else {is_valid=false;}      
}

void Alignment::copy_head_tail(const bam1_t *hts_align) {    
    // Get seq
    UINT32 qlen = hts_align->core.l_qseq;
    if (_qab > 0) {
        PackedSeq<4> hps (_qab,0,bam_get_seq(hts_align));
        _head = std::move(hps);
    }

    if (qlen > _qae) {
        PackedSeq<4> tps (qlen-_qae,_qae,bam_get_seq(hts_align));
        _tail = std::move(tps);
    }
    
    UINT8* pbasequal = bam_get_qual(hts_align);
    if (pbasequal!=nullptr)  {
        _head_basequal.resize(_qab);                
        for ( UINT32 i = 0; i < _qab; ++i ) {
                _head_basequal[i] = static_cast<char>(pbasequal[i]+33);
        }
        _tail_basequal.resize(qlen-_qae);                
        for ( UINT32 i = _qae; i < qlen; ++i ) {
                _tail_basequal[i-_qae] = static_cast<char>(pbasequal[i]+33);
        }
    }     
}

// SUPPLEMENTARY JOINING: try to join this alignment with argument
bool Alignment::join_supplementary(Alignment & pair, std::string full_read, std::string full_qual, Contig & contig, int _engine_idx) {
    
    // sanity checks
    assert(_qname == pair._qname);
    assert(get_real_length() == full_read.size());
    assert(pair.get_real_length() == full_read.size());
    
    // can't join if different orientation
    if(_is_rev != pair._is_rev) return false;
    
    // can't join if distance in contig is too large
    if(std::abs(int(_re) - int(pair._rb)) >= 10000) return false;
    
    // in case second alignment is fully contained in the first
    if(_re >= pair._re) return false; // may want to invalidate
    if (_left_clip + _qae - _qab >= pair._left_clip + pair._qae - pair._qab) return false;
    
    // in case first alignment is actually subset of second
    if (pair._left_clip < _left_clip) return false;
    
    // try joining
    
    // reverse full_read if needed
    if(_is_rev) {
        std::reverse(full_read.begin(), full_read.end());
        std::reverse(full_qual.begin(), full_qual.end());
        for(int j = 0; j < full_read.size(); j++) {
            if(full_read[j] == 'A') full_read[j] = 'T';
            else if(full_read[j] == 'C') full_read[j] = 'G';
            else if(full_read[j] == 'G') full_read[j] = 'C';
            else if(full_read[j] == 'T') full_read[j] = 'A';
        }
    }
    
    // note important positions:
    
    // segment 1 end / segment 2 start, in terms of read pos
    UINT32 r1 = _left_clip + _qae - _qab;
    UINT32 r2 = pair._left_clip;
    
    // segment 1 end / segment 2 start, in terms of contig pos
    UINT32 c1 = _re;
    UINT32 c2 = pair._rb;
    
    UINT32 curr_qp = _left_clip;
    UINT32 curr_rp = _rb;
    
    // read position of c2
    UINT32 d2 = full_read.size();
    // reference position of r2
    UINT32 r2_s1 = _re;
    
    UINT32 last_aligned_query = curr_qp;
    UINT32 last_aligned_contig = curr_rp;
    
    UINT32 c2_adjusted = c2;
    UINT32 r2_adjusted = r2;
    for(UINT32 i = 0; i < _num_cigar_op; ++i) {
        UINT32 current_cigar = _cigar[i];
        UINT32 cigar_op = bam_cigar_op(current_cigar);
        UINT32 cigar_oplen = bam_cigar_oplen(current_cigar);
        
        if(cigar_op == BAM_CSOFT_CLIP || cigar_op == BAM_CHARD_CLIP) {
            continue;
        }
                
        if((bam_cigar_type(cigar_op) & 3) == 3) { //bit 1 and bit 2 set, consumes query and reference
            if(curr_rp <= c2 && curr_rp + cigar_oplen > c2) { // found c2
                d2 = curr_qp; // + (c2 - curr_rp);
                c2_adjusted = curr_rp;
                //std::cout << "A set d2 " << curr_qp << " " << d2 << std::endl;
            }
            if(curr_qp <= r2 && curr_qp + cigar_oplen > r2) { // found r2
                r2_s1 = curr_rp; // + (r2 - curr_qp);
                r2_adjusted = curr_qp;
                //std::cout << "A set r2 " << curr_rp << " " << r2_s1 << std::endl;
            }
            
            last_aligned_query = curr_qp + cigar_oplen - 1;
            last_aligned_contig = curr_rp + cigar_oplen - 1;
            
            curr_rp += cigar_oplen;
            curr_qp += cigar_oplen;
        } else if(bam_cigar_type(cigar_op) & 2) { //bit 2 set, consumes reference
            if(curr_rp <= c2 && curr_rp + cigar_oplen > c2) {
                //std::cout << "B set d2 " << curr_qp << " " << d2 << " " << c2 << " " << cigar_oplen << " " << std::endl;
                // d2 = curr_qp - 1; 
                // c2 = curr_rp - 1;
                // in this case, change c2 to reflect
                d2 = last_aligned_query;
                c2_adjusted = last_aligned_contig;
                //std::cout << "B set d2 " << curr_qp << " " << d2 << " " << c2 << " " << cigar_oplen << " " << std::endl;
            }
            
            curr_rp += cigar_oplen;
        } else if(bam_cigar_type(cigar_op) & 1) { //bit 1 set, consumes query
            if(curr_qp <= r2 && curr_qp + cigar_oplen > r2) { // found r2
                // r2_s1 = curr_rp - 1; // in this case, change r2 to reflect
                // r2 = curr_qp - 1;
                r2_s1 = last_aligned_contig;
                r2_adjusted = last_aligned_query;
                //std::cout << "B set r2 " << curr_rp << " " << r2_s1 << std::endl;
            }
            curr_qp += cigar_oplen;
        }
    }
    
    // find closest read position of c1
    // find reference position of r1
    curr_qp = pair._left_clip;
    curr_rp = pair._rb;
    
    // read posiiton of c1
    UINT32 d1 = 0;
    // reference position of r1
    UINT32 r1_s2 = pair._re;
    
    last_aligned_query = curr_qp;
    last_aligned_contig = curr_rp;
    
    UINT32 c1_adjusted = c1;
    UINT32 r1_adjusted = r1;
    for(UINT32 i = 0; i < pair._num_cigar_op; ++i) {
        UINT32 current_cigar = pair._cigar[i];
        UINT32 cigar_op = bam_cigar_op(current_cigar);
        UINT32 cigar_oplen = bam_cigar_oplen(current_cigar);
        
        if(cigar_op == BAM_CSOFT_CLIP || cigar_op == BAM_CHARD_CLIP) {
            continue;
        }
                
        if((bam_cigar_type(cigar_op) & 3) == 3) { //bit 1 and bit 2 set, consumes query and reference
            if(curr_rp <= c1 && curr_rp + cigar_oplen >= c1) { // found c1
                d1 = curr_qp; // + (c1 - curr_rp);
                c1_adjusted = curr_rp;
                //std::cout << "A set d1 " << curr_qp << " " << d1 << std::endl;
            }
            if(curr_qp <= r1 && curr_qp + cigar_oplen >= r1) { // found r1
                
                r1_s2 = curr_rp; // + (r1 - curr_qp);
                r1_adjusted = curr_qp;
                //std::cout << "A set r1_s2 " << curr_qp << " " << d2 << std::endl;
            }
            
            last_aligned_query = curr_qp + cigar_oplen - 1;
            last_aligned_contig = curr_rp + cigar_oplen - 1;
            
            curr_rp += cigar_oplen;
            curr_qp += cigar_oplen;
            
            
        } else if(bam_cigar_type(cigar_op) & 2) { //bit 2 set, consumes reference
            if(curr_rp <= c1 && curr_rp + cigar_oplen >= c1) { // found c1 on an insertion
                //d1 = curr_qp - 1;
                //c1 = curr_rp - 1;
                d1 = last_aligned_query;
                c1_adjusted = last_aligned_contig;
            }
            curr_rp += cigar_oplen;
        } else if(bam_cigar_type(cigar_op) & 1) { //bit 1 set, consumes query
            if(curr_qp <= r1 && curr_qp + cigar_oplen >= r1) { // found r1 but on a deletion
                // r1_s2 = curr_rp - 1; // put after the deletion
                // r1 = curr_qp - 1;
                
                r1_s2 = last_aligned_contig;
                r1_adjusted = last_aligned_query;
            }
            curr_qp += cigar_oplen;
        }
    }
    
    // c2 in segment 1 and c1 in segment 2
    UINT32 left_end = std::min(std::min(r1, r2_adjusted), d2);
    UINT32 right_end = std::max(std::max(r1_adjusted, r2), d1);
    
    assert(right_end >= left_end); // sanity check
    
    // consistently cut on left if possible
    UINT32 contig_left = contig._pseq.get_seq_size(); 
    if(left_end == r1) contig_left = std::min(c1, contig_left);
    if(left_end == r2_adjusted) contig_left = std::min(r2_s1, contig_left);
    if(left_end == d2) contig_left = std::min(c2_adjusted, contig_left);
    
    UINT32 contig_right = contig._pseq.get_seq_size(); 
    if(right_end == r1_adjusted) contig_right = std::min(r1_s2, contig_right);
    if(right_end == r2) contig_right = std::min(c2, contig_right);
    if(right_end == d1) contig_right = std::min(c1_adjusted, contig_right);
    
    assert(contig_right >= contig_left); // sanity check
    /*
    std::cout << "Possible left cuts: " << r1 << " " << r2 << " " << d2 << std::endl;
    std::cout << "Possible left cuts: " << c1 << " " << r2_s1 << " " << c2 << std::endl;
    std::cout << "Possible right cuts: " << r1 << " " << r2 << " " << d1 << std::endl;
    std::cout << "Possible right cuts: " << r1_s2 << " " << c2 << " " << c1 << std::endl;
    std::cout << "Read cuts: " << left_end << " " << right_end << std::endl;
    std::cout << "Contig cuts: " << contig_left << " " << contig_right << std::endl;
    *///
    
    /*
    // adjust ends until it is aligned
    while(left_end > 0 && left_read_to_ctg.find(left_end) == left_read_to_ctg.end()) left_end--;
    while(right_end <= full_read.size() && right_read_to_ctg.find(right_end) == right_read_to_ctg.end()) right_end++;
    
    if(left_end == 0 || right_end >= full_read.size()) {
        //std::cout << "strange case, unaligned parts" << std::endl;
        return false;
    }
    
    UINT32 contig_left = left_read_to_ctg[left_end];
    UINT32 contig_right = right_read_to_ctg[right_end];
    */
    
    // force alignment of contig contig_left:contig_right with read left_end:right_end
    
    //std::cout << "start aligning" << std::endl;
    
    //SWWrapper::Aligner aligner("ACGTN", 2, 8, 12, 2, 32, 1, -1, -1);
    //SWWrapper::SWAlignment result;
    
    //std::cout << contig_left << " " << contig_right << " " <<contig._pseq.get_seq_size() << std::endl;
    
    std::string contig_part = contig._pseq.unpack(contig_left, contig_right);
    std::string read_part = full_read.substr(left_end, right_end - left_end);
             
    // cases of empty parts?
    if(contig_part.size() == 0) {
        _minimap_sw_results[_engine_idx]->cigar.clear();
        if(read_part.size() > 0) _minimap_sw_results[_engine_idx]->cigar.push_back(bam_cigar_gen(read_part.size(), BAM_CINS));
    } else if(read_part.size() == 0) {
        _minimap_sw_results[_engine_idx]->cigar.clear();
        _minimap_sw_results[_engine_idx]->cigar.push_back(bam_cigar_gen(contig_part.size(), BAM_CDEL));
    } else {
        _minimap_sw_engines[_engine_idx]->align(contig_part, read_part, _minimap_sw_results[_engine_idx].get());
    }
    
    // join the alignments: 
    
    //std::cout << "updating" << std::endl;
    
    // CIGAR: first CIGAR up to left_end
    std::vector<UINT32> new_cigar;
    curr_qp = _left_clip;
    curr_rp = _rb;
    
    UINT32 global_rp = _rb;
    UINT32 global_qp = _left_clip;
    
    for(UINT32 i = 0; i < _num_cigar_op; ++i) {
        UINT32 current_cigar = _cigar[i];
        UINT32 cigar_op = bam_cigar_op(current_cigar);
        UINT32 cigar_oplen = bam_cigar_oplen(current_cigar);
        
        if(cigar_op == BAM_CSOFT_CLIP || cigar_op == BAM_CHARD_CLIP) {
            // if(i == 0) new_cigar.push_back(bam_cigar_gen(cigar_oplen, BAM_CSOFT_CLIP)); // discard clips
            continue;
        }
                
        if((bam_cigar_type(cigar_op) & 3) == 3) { //bit 1 and bit 2 set, consumes query and reference
            if(curr_qp + cigar_oplen <= left_end) new_cigar.push_back(current_cigar);
            else {
                if(left_end - curr_qp > 0) {
                    new_cigar.push_back(bam_cigar_gen(left_end - curr_qp, cigar_op));
                    global_rp += left_end - curr_qp;
                    global_qp += left_end - curr_qp;
                }
                break;
            }
            curr_qp += cigar_oplen;
            global_qp += cigar_oplen;
            
            global_rp += cigar_oplen;
            curr_rp += cigar_oplen;
            
        } else if(bam_cigar_type(cigar_op) & 2) { //bit 2 set, consumes reference
            new_cigar.push_back(current_cigar);
            
            global_rp += cigar_oplen;
            curr_rp += cigar_oplen;
        } else if(bam_cigar_type(cigar_op) & 1) { //bit 1 set, consumes query
            if(curr_qp + cigar_oplen <= left_end) new_cigar.push_back(current_cigar);
            else {
                if(left_end - curr_qp > 0) { 
                    new_cigar.push_back(bam_cigar_gen(left_end - curr_qp, cigar_op));
                    global_qp += left_end - curr_qp;
                }
                break;
            }
            curr_qp += cigar_oplen;
            global_qp += cigar_oplen;
        }
    }
    //std::cout << "L " << global_qp << " " << left_end << " " << curr_qp << std::endl;
    //std::cout << "L " << global_rp << " " << contig_left << std::endl;
    assert(global_qp == left_end);
    assert(global_rp == contig_left);
    
    // CIGAR: alignment results
    
    // any clips on new cigar?
    std::string construct_cigar = "";
    for(int i = 0; i < _minimap_sw_results[_engine_idx]->cigar.size(); i++) {
        UINT32 current_cigar = _minimap_sw_results[_engine_idx]->cigar[i];
        UINT32 cigar_op = bam_cigar_op(current_cigar);
        UINT32 cigar_oplen = bam_cigar_oplen(current_cigar);
        
        assert(cigar_op != BAM_CSOFT_CLIP && cigar_op != BAM_CHARD_CLIP);
        
        if((bam_cigar_type(cigar_op) & 3) == 3) { //bit 1 and bit 2 set, consumes query and reference
            global_rp += cigar_oplen;
            global_qp += cigar_oplen;
            
            construct_cigar += std::to_string(cigar_oplen) + "M";
        } else if(bam_cigar_type(cigar_op) & 2) { //bit 2 set, consumes reference
            global_rp += cigar_oplen;
            
            construct_cigar += std::to_string(cigar_oplen) + "D";
        } else if(bam_cigar_type(cigar_op) & 1) { //bit 1 set, consumes query
            global_qp += cigar_oplen;
            
            construct_cigar += std::to_string(cigar_oplen) + "I";
        }
    }
    
    /*
     std::cout << _rb << " " << _re << " " << _left_clip << " " << _right_clip << std::endl;
    std::cout << pair._rb << " " << pair._re << " " << pair._left_clip << " " << pair._right_clip << std::endl;
    std::cout << construct_cigar << std::endl;
    std::cout << contig_left << " " << contig_right << std::endl;
    std::cout << left_end << " " << right_end << " " << contig_part.size() << " " << read_part.size() << std::endl;
    std::cout << "I " << global_rp << " " << contig_right << " " << global_qp << " " << right_end << std::endl; */
    assert(global_rp == contig_right);
    assert(global_qp == right_end);
    
    
    new_cigar.insert(new_cigar.end(), _minimap_sw_results[_engine_idx]->cigar.begin(), _minimap_sw_results[_engine_idx]->cigar.end());
    
    // CIGAR: starting from right_end
    curr_qp = pair._left_clip;
    curr_rp = pair._rb;
    
    for(UINT32 i = 0; i < pair._num_cigar_op; ++i) {
        UINT32 current_cigar = pair._cigar[i];
        UINT32 cigar_op = bam_cigar_op(current_cigar);
        UINT32 cigar_oplen = bam_cigar_oplen(current_cigar);
        
        if(cigar_op == BAM_CSOFT_CLIP || cigar_op == BAM_CHARD_CLIP) {
            // if(i > 0) new_cigar.push_back(bam_cigar_gen(cigar_oplen, BAM_CSOFT_CLIP)); // discard clips
            continue;
        }
                
        if((bam_cigar_type(cigar_op) & 3) == 3) { //bit 1 and bit 2 set, consumes query and reference
            if(curr_qp <= right_end) {
                if(curr_qp + cigar_oplen > right_end) {
                    if(curr_qp + cigar_oplen - right_end > 0) new_cigar.push_back(bam_cigar_gen(curr_qp + cigar_oplen - right_end, cigar_op));
                    
                    global_rp += curr_qp + cigar_oplen - right_end;
                    global_qp += curr_qp + cigar_oplen - right_end;
                }
                
                
            } else {
                new_cigar.push_back(current_cigar);
                
                global_rp += cigar_oplen;
                global_qp += cigar_oplen;
            }
            curr_qp += cigar_oplen;
            
            
        } else if(bam_cigar_type(cigar_op) & 2) { //bit 2 set, consumes reference
            if(curr_qp > right_end) { new_cigar.push_back(current_cigar); global_rp += cigar_oplen;}
        
            
            
        } else if(bam_cigar_type(cigar_op) & 1) { //bit 1 set, consumes query
            if(curr_qp <= right_end) {
                if(curr_qp + cigar_oplen > right_end) {
                    if(curr_qp + cigar_oplen - right_end > 0) {
                        new_cigar.push_back(bam_cigar_gen(curr_qp + cigar_oplen - right_end, cigar_op));
                        global_qp += curr_qp + cigar_oplen - right_end;
                    }
                    
                }
            } else {
                new_cigar.push_back(current_cigar);
                global_qp += cigar_oplen;
            }
            curr_qp += cigar_oplen;
            
        }
    }
    //std::cout << "R " << global_rp << " " << pair._re << " " << right_end << " " << std::endl;
    //std::cout << "R " << curr_qp << " " << global_qp << " " << pair._qae + pair._left_clip << " " << global_rp << " " << pair._re << std::endl;
    assert(global_rp == pair._re);
    assert(global_qp == pair._qae + pair._left_clip);
    
    // update variables
    _rb = _rb;
    _re = pair._re;
    
    
    _qab = 0; // _left_clip;
    _qae = full_read.size() - _left_clip - pair._right_clip;
    
    _cigar = new_cigar;
    _num_cigar_op = _cigar.size();
    
    _basequal = full_qual.substr(_left_clip, _qae);
    _apseq = PackedSeq<2>(full_read.substr(_left_clip, _qae));
    
    _right_clip = pair._right_clip;
    
    assert(_qae - _qab + _right_clip + _left_clip == full_read.size());
    
    // custom mq, > 60 to show this is custom alignment
    _mq = 111;
    
    // invalidate second alignment
    pair.invalidate();
    
    // CIGAR sanity check
    
    
    curr_qp = 0;
    curr_rp = _rb;
    
    int get_left_clip = 0;
    int get_right_clip = 0;
    
    for(UINT32 i = 0; i < _num_cigar_op; ++i) {
        UINT32 current_cigar = _cigar[i];
        UINT32 cigar_op = bam_cigar_op(current_cigar);
        UINT32 cigar_oplen = bam_cigar_oplen(current_cigar);
        
        if(cigar_op == BAM_CSOFT_CLIP || cigar_op == BAM_CHARD_CLIP) {
            assert(i == 0 || i == _num_cigar_op - 1);
            if(i == 0) get_left_clip = cigar_oplen;
            else get_right_clip = cigar_oplen;
        }
                
        if((bam_cigar_type(cigar_op) & 3) == 3) { //bit 1 and bit 2 set, consumes query and reference
            curr_qp += cigar_oplen;
            curr_rp += cigar_oplen;
        } else if(bam_cigar_type(cigar_op) & 2) { //bit 2 set, consumes reference
            curr_rp += cigar_oplen;
        } else if(bam_cigar_type(cigar_op) & 1) { //bit 1 set, consumes query
            curr_qp += cigar_oplen;
        }
    }
    
    assert(get_left_clip == 0);
    assert(get_right_clip == 0);
    assert(curr_rp == _re);
    assert(curr_qp == _qae);
    
    //std::cout << "End join" << std::endl;
    
    return true;
}

std::vector<std::shared_ptr<SWWrapper::Aligner>> Alignment::_minimap_sw_engines;
std::vector<std::shared_ptr<SWWrapper::SWAlignment>> Alignment::_minimap_sw_results;

void Alignment::initiate_sw_engines(const UINT32 num_threads) {
    if(!_minimap_sw_engines.empty()) return;
    for (UINT32 i = 0; i < num_threads; ++i) {
        _minimap_sw_engines.emplace_back(std::make_unique<SWWrapper::Aligner>("ACGTN", 2, 8, 12, 2, 32, 1, -1, -1));
        _minimap_sw_results.emplace_back(std::make_unique<SWWrapper::SWAlignment>());
    }
}


} // namespace hypo

