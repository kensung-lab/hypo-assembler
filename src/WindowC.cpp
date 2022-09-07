/*
 *
 * Copyright (c) 2019, Ritu Kundu  and Joshua Casey
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
/** Defines the class Window.
 * It contains the functionality for POA.
 */
#include "Window.hpp"
#include <cmath>
#include <set>
#include <numeric>


namespace hypo
{
Mode Window::_mode = Mode::LSA;
std::vector<std::shared_ptr<spoa::AlignmentEngine>>
    Window::_alignment_engines_long;
std::vector<std::shared_ptr<spoa::AlignmentEngine>>
    Window::_alignment_engine_pslong;       
std::vector<std::shared_ptr<SWWrapper::Aligner>> Window::_minimap_sw_engines;
std::vector<std::shared_ptr<SWWrapper::SWAlignment>> Window::_minimap_sw_results;
void Window::prepare_for_poa(const UINT32 num_threads)
{
    if (!_alignment_engines_long.empty()) {return;} // Already initialised
    for (UINT32 i = 0; i < num_threads; ++i)
    {
        _alignment_engines_long.emplace_back(
            spoa::createAlignmentEngine(spoa::AlignmentType::kNW, lr_match_score,
                                        lr_misMatch_score, lr_gap_penalty));
    }
    _alignment_engines_long.shrink_to_fit();
    _alignment_engine_pslong.emplace_back(
            spoa::createAlignmentEngine(spoa::AlignmentType::kOV, lr_match_score,
                                        lr_misMatch_score, lr_gap_penalty));
}

void Window::prepare_for_cluster(const UINT32 num_threads)
{
    if (!_minimap_sw_engines.empty()) {return;} // Already initialised
    for (UINT32 i = 0; i < num_threads; ++i)
    {
        std::string possible_chars = cMarkerLetters + "ACGTN";
        //_minimap_sw_engines.emplace_back(std::make_unique<SWWrapper::Aligner>(possible_chars, sp.sr_match_score, sp.sr_misMatch_score, sp.sr_gap_penalty, sp.sr_gap_penalty));
        _minimap_sw_engines.emplace_back(std::make_unique<SWWrapper::Aligner>(possible_chars, sw_match_score, sw_mismatch_score, sw_gap_open1, sw_gap_ext1, sw_gap_open2, sw_gap_ext2, -1, -1));
        _minimap_sw_results.emplace_back(std::make_unique<SWWrapper::SWAlignment>());
    }
    _minimap_sw_engines.shrink_to_fit();
    _minimap_sw_results.shrink_to_fit();
}

void Window::generate_consensus(const UINT32 engine_idx)
{
    _engine_idx = engine_idx;
    auto wt = _wtype;
    auto num_non_empty_arms = _num_internal + _num_pre + _num_suf;
    //if (_num_empty > num_non_empty_arms)
    if ((wt == WindowType::LONG || (is_internal_only && _draft.get_seq_size() < RELIABLE_DRAFT_STH)) && _num_empty > _num_internal)
    {
        set_consensus(""); // empty sequence is consensus
    }
    else if (num_non_empty_arms >=2)
    { // There should be at least 2 arms to do consensus
        if (wt == WindowType::SHORT)
        {
            generate_consensus_short();
        }
        else
        {
            generate_consensus_long();
        }
    }
    else
    {
        set_consensus(_draft.unpack());
    }
    //clear_arms();
}

void Window::reset_consensus(const std::string& s1,const std::string& s2) {
    assert(!s1.empty() && !s2.empty());
    auto graph = spoa::createGraph();
    auto alignment = _alignment_engine_pslong[0]->align(s1, graph);
    graph->add_alignment(alignment, s1);
    alignment = _alignment_engine_pslong[0]->align(s2, graph);
    graph->add_alignment(alignment, s2);
    //_consensus = graph->generate_consensus();

    std::vector<std::string> msa;
    graph->generate_multiple_sequence_alignment(msa, true);
    std::cout << "================================="<<std::endl;
    for (const auto &it : msa){
        std::cout << it.c_str() << std::endl;
    }
    auto mid = msa[0].size()/2;
    // Left string dominates untilmiddle, then right takes over
    for (UINT i=0; i < mid; ++i) {
        if (msa[0][i]!='-') {_consensus+=msa[0][i];}
    }
    for (UINT i=mid; i < msa[0].size(); ++i) {
        if (msa[1][i]!='-') {_consensus+=msa[1][i];}
    }
    std::cout << _consensus << std::endl;
}

// Only for debugging ()
std::ostream &operator<<(std::ostream &os, const Window &wnd)
{
    // Write numbers
    os << wnd._num_internal << "\t" << wnd._num_pre << "\t" << wnd._num_suf
       << "\t" << wnd._num_empty << "\t" << wnd.info << std::endl;
    // Write draft sequence
    os << "++\t" << wnd._draft.unpack() << std::endl;

    // Write consensus sequence
    os << "++\t" << wnd._consensus << std::endl;
    // Write internal arms

    for (UINT32 i = 0; i < wnd._num_internal; ++i)
    {
        //os << wnd._internal_arms[i].unpack() << std::endl;
    }
    // Write pref arms
    for (UINT32 i = 0; i < wnd._num_pre; ++i)
    {
        //os << wnd._pre_arms[i].unpack() << std::endl;
    }
    // Write suff arms
    for (UINT32 i = 0; i < wnd._num_suf; ++i)
    {
        //os << wnd._suf_arms[i].unpack() << std::endl;
    }

    return os;
}


void Window::select_consensus_from_draft(bool is_marked)
{
    std::string modified_str;
    if (is_marked) {
        modified_str = cHead + _draft.unpack() + cTail;
    }
    else {
        modified_str = _draft.unpack();
    }

    UINT32 min_id = 0;
    UINT32 min_distance = 2*modified_str.size();
    UINT32 min_mis = 2*modified_str.size();
    for (UINT32 j = 0; j < _multi_consensus.size(); j++)
    {
        std::vector<UINT16> diff(DIST+1,0);
        std::cout << "Choosing from multiple\n";
        sw_diff(modified_str, _multi_consensus[j], diff);
        if (diff[X] < min_mis)
        {
            min_id = j;
            min_mis = diff[X];
            min_distance = diff[DIST];
        }
        else if (diff[X] == min_mis && diff[DIST] < min_distance)
        {
            min_id = j;
            min_distance = diff[DIST];
        }
    }
    std::cout << "mismatch: " << min_mis <<"  dis:" << min_distance << " " << _multi_consensus[min_id]<<std::endl;
    if (is_marked) {
        set_marked_consensus(_multi_consensus[min_id]);
    }
    else {
        set_consensus(_multi_consensus[min_id]);
    }
    
}

bool Window::is_reliable(std::string& cig){
    if (is_sh_high() && has_long_indel(cig)) {return false;}
    UINT num_mis = 0;
    UINT num_ins = 0;
    UINT num_del = 0;
    bool is_pvs_ins = false;
    for (auto c: cig) {
        if (c=='A' || c=='C'|| c=='G'|| c=='T') {
            ++num_mis;
            is_pvs_ins=false;
        }
        else if (c=='a' || c=='c'|| c=='g'|| c=='t') {
            if (!is_pvs_ins) { // first ins base
                ++num_ins;
            }
            is_pvs_ins=true;
        }
        else if (c=='D') {
            ++num_del;
            is_pvs_ins=false;
        }
        else {
            is_pvs_ins=false;
        }
    }
    return (num_mis<=(num_ins+num_del));
}

bool Window::has_long_indel(std::string& cig){
    std::cout<<"Checking has long indel: "<<cig<<std::endl;
    UINT len_ins = 0;
    std::string del_len = "0";
    bool is_pvs_ins = false;
    for (auto c: cig) {
        if (c=='A' || c=='C'|| c=='G'|| c=='T') {
            if (len_ins >= RELIABLE_DRAFT_STH ) {std::cout<<c<<" "<<len_ins<<std::endl;return true;}
            len_ins = 0;
            del_len = "";
        }
        else if (c=='a' || c=='c'|| c=='g'|| c=='t') {
            ++len_ins;
            del_len = "0";
        }
        else if (c=='D') {
            if (len_ins >= RELIABLE_DRAFT_STH || std::stoi(del_len)  >= RELIABLE_DRAFT_STH) {std::cout<<"Del:INS "<<del_len<<":"<<len_ins<<std::endl; return true;}
            len_ins = 0;
            del_len = "0";
        }
        else if (c=='M') {
            if (len_ins >= RELIABLE_DRAFT_STH ) {return true;}
            len_ins = 0;
            del_len = "";
        }
        else {
            del_len += c;
        }
    }
    if (len_ins >= RELIABLE_DRAFT_STH ) {std::cout<<"LAst: "<<len_ins<<std::endl;return true;}
    std::cout<<"No Long Indel\n";
    return false;
}

bool Window::select_most_frequent() {
// Will always have at least one arm
    // dictionary of unique cigar string with  vector of id that correspind to it.
    /* Why we need vector of id?
     * If the pvs/nxt anchor has insertion/del (i.e. approx position of anchor is used rather than exact), cigar does not reflect the actual sequence.
     * Therefore, lengths of armswith that cigar should also be considered in identifying the correspondence between sequence and cigar.
     * 
    */
    std::unordered_map<std::string,std::vector<UINT32>> cigseqs;
    std::vector<UINT32> counts;  
    std::vector<UINT32> internal_arm_index; // index of the first arm with that cig
    std::vector<UINT32> internal_arm_len; // length of the first arm with that cig
    UINT32 cig_id = 0;

    for (UINT32 i = 0; i < _internal_arms.size(); ++i)
    {
        // Ignore read seqs with big indels
        if (has_long_indel(_internal_cigs[i])) {std::cout<<"ignore\n"; continue;}
        auto search = cigseqs.find(_internal_cigs[i]);
        bool added_new = true;
        auto seqlen = _internal_arms[i].get_seq_size();
        if (search != cigseqs.end()) {
            for (auto cigid: search->second) {
               if (seqlen == internal_arm_len[cigid]) {
                   counts[cigid]+=_internal_awt[i];
                   added_new = false;
                   break;
               }
            } 
            if (added_new) { // new id needs to be added for this cig
                (search->second).emplace_back(cig_id);
            }           
        }
        else {
            std::vector<UINT32>dummy{cig_id};
            cigseqs.emplace(std::make_pair(_internal_cigs[i],dummy));
        }
        if (added_new) {            
            counts.emplace_back(_internal_awt[i]);
            internal_arm_index.emplace_back(i);
            internal_arm_len.emplace_back(seqlen);
            ++cig_id;
        }        
    }

    // Return if no sequence added
    if (counts.size()==0) {return false;}

    // Check if majority supports
    // Find most frequent
    auto result = std::max_element(counts.begin(), counts.end());
    auto max_count_id = std::distance(counts.begin(), result);
    /// * //////////// For debugging
    std::cout << "*****\t"<<_draft.unpack()<<std::endl;
    for (auto it=cigseqs.begin();it!=cigseqs.end();++it) {
        for (auto ci:it->second) {
            std::cout << it->first << "\t\t\t"<< ci <<"\t"<<counts[ci]<<"\t"<<internal_arm_index[ci]<<std::endl;
            }
    }
    //*/
    if (counts[max_count_id] < 2*CONSENSUS_CONSTANT_THRESHOLD) {
        return false;
    }

    if (counts.size()>1) { // more than 1 sequence present
        // sort in decreasing order of counts
        std::vector<uint32_t> rank;
        rank.reserve(counts.size());
        for (uint32_t i = 0; i < counts.size(); ++i) {
            rank.emplace_back(i);
        }
        std::sort(rank.begin(), rank.end(), [&](uint32_t l, uint32_t r) {
            return counts[l] > counts[r]; });

        INT64 remaining = _num_internal;
        for (auto i =0; i< rank.size();++i) {
            auto cig_ind = rank[i];
            auto seq_ind = internal_arm_index[cig_ind];
            _multi_consensus.emplace_back(_internal_arms[seq_ind].unpack());
            std::cout << "Added " << seq_ind << " " << _internal_arms[seq_ind].unpack() << std::endl;
            // Check if next should be added
            if (i+1<rank.size()) {
                // If heaviest and the second heaviest bundles are close, use evidence from P/S
                if (i < 2 && ((double)counts[cig_ind]*MOST_FREQUENT_FRACTION2 <counts[rank[i+1]])) {
                    std::cout << "Too close. Shifting to PS: " << std::endl;
                     return false;
                }
                // next has frequency < threshold (5)
                bool cond1 = counts[rank[i+1]] < 2*CONSENSUS_CONSTANT_THRESHOLD;
                // this one makes > th (0.8) of the remaining
                bool cond2 = ((double)counts[cig_ind] / (double)remaining) > 2*CONSENSUS_SUPPORT_FRACTION;
                // next one is less than th (40%) of this one
                bool cond3 = (double)counts[cig_ind]*MOST_FREQUENT_FRACTION1 >counts[rank[i+1]];
                // we trust draft only if read sets have more mismatches (draft can not be trusted with read-sets having indels)
                bool cond4 = is_reliable(_internal_cigs[seq_ind]);
                if (cond1 || cond2 || cond3 || cond4) {
                    break;
                }
            }
            remaining -= counts[cig_ind];
        }
        if (_multi_consensus.size()>1) {
            select_consensus_from_draft(false);
        }
        else {
            set_consensus(_multi_consensus[0]);
        }
        
    }
    else { // only one sequence
        auto most_freq_seq_ind = internal_arm_index[max_count_id];
        set_consensus(_internal_arms[most_freq_seq_ind].unpack());
    }
    return true;
}

void Window::generate_consensus_short()
{
    std::cout << "++++++++++++++++++++++++\t++++++++++++++++++++++++++++++++\t" << std::endl;
    std::cout << "++++++++\t" << _draft.unpack() <<std::endl;
    std::cout << _internal_arms.size() <<" " << _pre_arms.size() <<" "<<_suf_arms.size() <<std::endl;
    for (UINT32 i = 0; i < _num_internal; ++i)
    {
        std::cout <<"I" << i << "\t" <<_internal_awt[i]<<"\t"<<_internal_arms[i].unpack() << std::endl;
    }
    // Write pref arms
    for (UINT32 i = 0; i < _num_pre; ++i)
    {
        std::cout <<"P" << i << "\t" <<_pre_awt[i]<<"\t"<<_pre_arms[i].unpack() << std::endl;
    }
    // Write suff arms
    for (UINT32 i = 0; i < _num_suf; ++i)
    {
        std::cout <<"S" << i <<"\t" << _suf_awt[i]<<"\t"<<_suf_arms[i].unpack() << std::endl;
    }

    bool do_cluster = true;
    //if (_num_pre==0 && _num_suf==0 && _num_internal>0) { // internal only window
    if (is_internal_only) {
        do_cluster = !select_most_frequent();
    }
    /* POA only for pref-suf short windows */
    if (do_cluster) {
        generate_pscluster_consensus_short();   
    }
    else {
        //clear_arms();
        clear_pre_suf();
    }
}

std::string Window::get_SW_pileup_consensus(const std::string& backbone, const std::vector<std::string>& internal, const std::vector<std::vector<UINT16>>& iw, const std::vector<std::string>& pre, const std::vector<std::vector<UINT16>>& pw, const std::vector<std::string>& suf, const std::vector<std::vector<UINT16>>& sw) {
    auto conssize = backbone.size();
    std::string consensus;
    consensus.reserve(conssize);
    
    std::vector<UINT16> cov(backbone.size(),0);
    //std::vector<UINT16> num_dels(backbone.size(),0);
    std::vector<std::unordered_map<UINT16, UINT16>> dels(backbone.size(), std::unordered_map<UINT16, UINT16>());
    std::vector<std::unordered_map<std::string,UINT16>> ins(backbone.size(),std::unordered_map<std::string,UINT16>());

    // Whether a window with only very few prefix/suffix 
    bool is_super_low_ps = is_superlow_ps();
    
    std::vector<std::vector<UINT16>> bases_cover(backbone.size(), std::vector<UINT16>(cCode.size(), 0));
    std::cout << "PSA"<<std::endl;
    std::cout << backbone << std::endl;
    // Call SW aln for each
    for (auto i=0; i <internal.size(); ++i) {
        sw_align_internal(backbone, internal[i], iw[i], cov, bases_cover, dels, ins);
    }
    for (auto i=0; i < pre.size(); ++i) {
        //std::string unit;
        //auto trim_len = longest_periodic_pref(pre[i],unit,false);
        //sw_align_prefix(backbone, pre[i], pre[i].size()-trim_len, pw[i], cov, bases_cover, dels, ins);
        sw_align_prefix(backbone, pre[i], pw[i], cov, bases_cover, dels, ins);
    }
    for (auto i=0; i < suf.size(); ++i) {
        //std::string unit;
        //auto trim_len = longest_periodic_pref(suf[i],unit,true);
        //sw_align_suffix(backbone, suf[i], trim_len, sw[i], cov, bases_cover, dels, ins);
        sw_align_suffix(backbone, suf[i], sw[i], cov, bases_cover, dels, ins);
    }   

    // Debug
    /// *
    std::cout << "++++++\t"<<_draft.unpack()<<std::endl;
    std::cout << "Dels: " ;
    for (UINT i = 0; i < backbone.size(); ++i) {
        //for (auto j=0; j<6; ++j) {std::cout << bases_cover[i][j] << " ";}
        //std::cout << std::endl;
        for (auto it=dels[i].begin();it!=dels[i].end();++it) {
            std::cout << i << " " << it->first << ":\t"<< it->second<<" "<<cov[i]<<std::endl;
        }
        for (auto it=ins[i].begin();it!=ins[i].end();++it) {
            std::cout << i << " " << it->first << "\t\t"<< it->second<<" "<<cov[i]<<std::endl;
        }
    }
    // */

    // Rectify backbone 
    for (UINT i = 0; i < backbone.size();){
        if (cov[i]==0) { // not covered by any sequence; just copy the base
            consensus.push_back(backbone[i]);
            ++i;
            continue;
        }
        if (is_super_low_ps && cov[i]==1) { // in case of super low cov p/s, trust the draft
            std::cout << "Superlow\n";
            consensus.push_back(backbone[i]);
            ++i;
            continue;
        }
        // Find most frequent base
        UINT16 max_base_count = 2; // 2 is the max wt (properly paired)
        UINT16 max_base_ind = cNt4Table[backbone[i]];
        for(auto j=0;j<cCode.size();++j) {
            if (bases_cover[i][j]>max_base_count){max_base_count=bases_cover[i][j];max_base_ind=j;}
        }
        if (is_super_low_ps){max_base_ind = cNt4Table[backbone[i]];std::cout<<"Trusting draft\n";} // trust draft base in suoerlow
        // Find blank count
        UINT16 del_len = 0;
        UINT16 num_blank = 0;
        if (!dels[i].empty()) { // at least one del
            auto x = std::max_element(dels[i].begin(), dels[i].end(),
                [](const std::pair<UINT16,UINT16>& p1, const std::pair<UINT16,UINT16>& p2) {
                    return p1.second < p2.second; });
                del_len = x->first;
                num_blank = std::accumulate(dels[i].begin(), dels[i].end(), 0, [](UINT16 v, const std::pair<UINT16,UINT16>& pair){
                    return v + pair.second;


                });
        }
        bool choose_del = false;
        if (del_len > 0 && num_blank > max_base_count) { // this base  and following bases (up until dellen) should be deleted
            if (del_len>= RELIABLE_DRAFT_STH) {
                if (is_super_low_ps) {
                    std::cout << "Replaced with draft " << i << ":" << del_len<<":"<<backbone[i] << "  num-blank:base-cnt" << num_blank <<":"<<max_base_count<< " ";
                    consensus.append(backbone.substr(i,del_len)) ;
                    i+=del_len;
                    choose_del = true;
                }
            }
            else {
                std::cout << "Deleted " << i << ":" << del_len<<":"<<backbone[i] << "  num-blank:base-cnt" << num_blank <<":"<<max_base_count<< " ";
                choose_del = true;
                i+=del_len;                
            }            
        }
        if (!choose_del) {
            char ch = (cCode[max_base_ind]!='N') ? (cCode[max_base_ind]) : ((i<cMarkerLen)?(cHead[0]):(cTail[0]));
            if (max_base_count==0) {ch = backbone[i];} // No base here; take backbone
            consensus.push_back(ch);  
            // Check ins after base
            std::string ins_seq = "";
            UINT16 num_ins = 0;
            if (!ins[i].empty()) {
                auto x = std::max_element(ins[i].begin(), ins[i].end(),
                [](const std::pair<std::string,UINT16>& p1, const std::pair<std::string,UINT16>& p2) {
                    return p1.second < p2.second; });
                ins_seq = x->first;
                //num_ins = x->second;

                num_ins = std::accumulate(ins[i].begin(), ins[i].end(), 0, [](UINT16 v, const std::pair<std::string,UINT16>& pair){
                    return v + pair.second;
                });
            }
            if (num_ins>0 && 2*num_ins > cov[i+1]) { // add insertion (i+1 works because last base is tail and should have no ins)
                if (is_super_low_ps && ins_seq.size()>= RELIABLE_DRAFT_STH ) {;}// Do not trust long indel here
                else {
                consensus.append(ins_seq);
                std::cout << "Inserted " << i << ":" << ins_seq << ":" << num_ins<<":" << cov[i+1]<< "  ";
                }
            }                      
            ++i;
        }
    }
    
    std::cout << "\nC:\t"<<consensus<<std::endl;
    return consensus;
}

// patche-patchb are proper periodic (No leftover); patchb is the first base of the run; patche is the last base
bool Window::identify_patch_region(UINT32& patchb, UINT32& patche, UINT32& patch_unitlen) {
    // 0: pref;1: suff; counts 0/1 has pos of first/last base in run
    std::vector<std::unordered_map<UINT32, UINT16>> counts(2, std::unordered_map<UINT32, UINT16>());
    std::cout <<" IDentifying PAtch region\n";
    // Lambda for finding most supported pos
    auto best_pos = [&](const UINT16 i) -> UINT32 {
        UINT32 pos = 0;
        if (!counts[i].empty()) { // at least one del
            auto x = std::max_element(counts[i].begin(), counts[i].end(),
                [](const std::pair<UINT32,UINT16>& p1, const std::pair<UINT32,UINT16>& p2) {
                    return p1.second < p2.second; });
            pos = x->first;
        }
        return pos;
    };

    std::vector<UINT32> runb;
    std::vector<UINT32> rune;
    std::vector<UINT16> runcnt;
    std::vector<UINT16> rununit;
    UINT32 numruns = 0;
    std::string dmodified = _draft.unpack();

    // pref
    for (UINT32 i = 0; i < _pre_arms.size(); ++i) {
        std::string modified = _pre_arms[i].unpack();
        std::string unit = "";
        auto trim_len = longest_periodic_pref(modified,unit,false);
        //std::cout << i << " patchP " << trim_len << unit << " " << modified << std::endl;
        if (trim_len>0) {
            // Find the corresponding run
            //std::cout << "Ends at " << _pends[i] << std::endl;
            bool run_found = false;
            for (auto j =0; j < numruns; ++j) {
                if (_pends[i] >runb[j] && _pends[i] <=rune[j] && unit.size()==rununit[j]) {
                    //std::cout << "Found in run " << j << " "<<runb[j]<<"-"<<rune[j]<<":"<<rununit[j]<<std::endl;
                    ++runcnt[j];
                    run_found=true;
                    break;
                }
            }
            if (!run_found) { // new run
                runb.emplace_back(extend_left(dmodified,unit,_pends[i]-1)); // _pends[i] can not be 0
                std::reverse(unit.begin(),unit.end());
                rune.emplace_back(extend_right(dmodified,unit,_pends[i]));
                runcnt.emplace_back(1);
                rununit.emplace_back(unit.size());
                ++numruns;
            }
        }
    }

    // suff
    for (UINT32 i = 0; i < _suf_arms.size(); ++i) {
        std::string modified = _suf_arms[i].unpack();
        std::string unit = "";
        auto trim_len = longest_periodic_pref(modified,unit,true);
        //std::cout << i << " patchS " << trim_len << unit << " " << modified << std::endl;
        if (trim_len>0) { 
            // Find the corresponding run
            //std::cout << "Begins at " << _sbegs[i] << std::endl;
            bool run_found = false;
            for (auto j =0; j < numruns; ++j) {
                if (_sbegs[i]>=runb[j] && _sbegs[i]<rune[j] &&  unit.size()==rununit[j]) {
                    //std::cout << "Found in run " << j << " "<<runb[j]<<"-"<<rune[j]<<":"<<rununit[j]<<std::endl;
                    ++runcnt[j];
                    run_found=true;
                    break;
                }
            }
            if (!run_found) { // new run
                rune.emplace_back(extend_right(dmodified,unit ,_sbegs[i]));
                std::reverse(unit.begin(),unit.end());
                runb.emplace_back(extend_left(dmodified,unit,_sbegs[i]-1)); // _sbegs[i] can not be 0
                runcnt.emplace_back(1);
                rununit.emplace_back(unit.size());
                ++numruns;
            }
         }
    }

    // find maxsupported if it is sufficiently long
    if (numruns==0)  {std::cout << "NO NEED OF PATCHING\n";return false;}
    auto result = std::max_element(runcnt.begin(), runcnt.end());
    auto max_count_id = std::distance(runcnt.begin(), result);
    INT run_len = rune[max_count_id]-runb[max_count_id];
    for (auto i=0; i < runcnt.size(); ++i) {
        std::cout << runb[i] <<"-"<<rune[i]<<":"<<rununit[i]<<":"<<runcnt[i]<<" ";
    }
    std::cout <<std::endl;
    if ((run_len/rununit[max_count_id]) >= RUNLEN_TH && (runcnt[max_count_id] >= PS_RUN_FRAC*(_num_pre+_num_suf))){
        patchb = runb[max_count_id];
        patche = rune[max_count_id] - (run_len%rununit[max_count_id]) -1;
        patch_unitlen = rununit[max_count_id];
        std::cout << "Max id: " << max_count_id<<" patch:" <<  patchb << "-"<<patche<< " (unit:"<< patch_unitlen<<") "<< dmodified.substr(patchb,patche-patchb+1)<<std::endl;
        return true;
    }
    else {std::cout << "NO NEED OF PATCHING\n";return false;}

}

void Window::generate_pscluster_consensus_short()
{    /////////////////////// PATCH /////////////////////////
    // hashtable: key: num_units of run; value; count supporting that num
    std::unordered_map<UINT32, UINT16> num_units;
    // patchb is the first base of the run; patche is the last base
    UINT32 patchb = 0;
    UINT32 patche = 0;
    UINT32 patch_unitlen = 0;
    bool need_patch = identify_patch_region(patchb,patche, patch_unitlen);
    std::string mod_draft = cHead+_draft.unpack()+cTail; 
    patchb += cMarkerLen;
    patche += cMarkerLen;
    // LAMBDA; Add b's run_units in the patch-region if b spans it; returns weather b sapns it or not
    auto update_run_support = [&](const AlnType at, const std::string& b, const UINT16 supp) -> bool {
        bool spans = false;
        UINT32 qb=0;
        UINT32 qe=0;
        if (at==AlnType::S) {
            _minimap_sw_engines[_engine_idx]->align_suffix_rev_cigar(mod_draft, b, _minimap_sw_results[_engine_idx].get());
            spans = get_query_pos_right(mod_draft, b, patchb, patche, qb, qe);
        }
        else if (at==AlnType::P) {
            _minimap_sw_engines[_engine_idx]->align_prefix(mod_draft, b, _minimap_sw_results[_engine_idx].get());
            spans = get_query_pos_left(mod_draft, b, patchb, patche, qb, qe);
        }
        else {
            _minimap_sw_engines[_engine_idx]->align(mod_draft, b, _minimap_sw_results[_engine_idx].get());
            spans = get_query_pos_left(mod_draft, b, patchb, patche, qb, qe);
        }
        if (spans && qb<qe) {
            // Check if periodic (no leftover)
            std::string unit;
            auto run_len = longest_periodic_pref(b.substr(qb,qe-qb+1),unit,true);
            std::cout << "qb:"<<qb<<" qe:"<<qe<<" u:"<<unit<<" runlen:"<<run_len<<" " << b.substr(qb,qe-qb+1) << " ";
            if (run_len==(qe-qb+1) && (run_len%unit.size()==0) && unit.size()==patch_unitlen) { // update support
                std::cout << "Diff: " << run_len - (patche-patchb) << std::endl;
                if (std::abs((int) run_len - int(patche-patchb)) < RELIABLE_DRAFT_STH) { // Trust only small indels
                    UINT32 nu = (UINT32) run_len/unit.size();
                    num_units[nu] += supp ;
                    std::cout << " supp:" << supp;
                }
            }
            std::cout << std::endl;
        }
        return spans;
    };

    // LAMBDA; Add b's run_units in the patch-region if b spans it including the clippped region; 
    auto update_cliprun_support = [&]() -> void {
        for (UINT i = 0; i < _pre_arms.size(); ++i)  {
            std::cout <<"PCLP: "<<i<<" pend:"<< _pends[i] << " pclip:"<<_pclips[i]<<std::endl;
            std::cout <<_pre_arms[i].unpack()<<std::endl;
            if (_pends[i]+cMarkerLen > patchb && _pends[i]+cMarkerLen <= patche && _pends[i]+cMarkerLen+_pclips[i]>patche) {
                std::string unit;
                auto run_len = longest_periodic_pref(_pre_arms[i].unpack(),unit,false);
                std::cout <<"unit: " <<unit<<" "<<run_len<<std::endl;
                if (unit.size()>0 && (run_len%unit.size()==0) && unit.size()==patch_unitlen) { // update support
                    UINT32 nu = (UINT32) run_len/unit.size();
                    num_units[nu] += _pre_awt[i] ;
                    std::cout << "Clp: P "<<i<<" "<<nu<<" "<<_pre_awt[i]<<std::endl;
                }
            }       
        }

        for (UINT i = 0; i < _suf_arms.size(); ++i)  {
            std::cout <<"SCLP: "<<i<<" sbeg:"<< _sbegs[i] << " sclip:"<<_sclips[i]<<std::endl;
            std::cout <<_suf_arms[i].unpack()<<std::endl;
            std::cout <<patchb<<"-"<<patche<<std::endl;
            std::cout<<(_sbegs[i] >= patchb)<<" "<<(_sbegs[i] < patche) <<""<<(_sbegs[i]<patchb+_sclips[i])<<std::endl;
            if (_sbegs[i]+cMarkerLen >= patchb && _sbegs[i]+cMarkerLen < patche && _sbegs[i]+cMarkerLen<patchb+_sclips[i]) {
                std::string unit;
                auto run_len = longest_periodic_pref(_suf_arms[i].unpack(),unit,true);
                std::cout <<"unit: " <<unit<<" "<<run_len<<std::endl;
                if (unit.size()>0 && (run_len%unit.size()==0) && unit.size()==patch_unitlen) { // update support
                    UINT32 nu = (UINT32) run_len/unit.size();
                    num_units[nu] += _suf_awt[i] ;
                    std::cout << "Clp: S "<<i<<" "<<nu<<" "<<_suf_awt[i]<<std::endl;
                }
            }       
        }
    };

    auto patch_it = [&](const std::string& temp_cons) -> std::string {
        UINT32 qb=0;
        UINT32 qe=0;
        std::string final_cons = temp_cons;
        std::cout << "Patching " << temp_cons<<std::endl;
        for (auto it=num_units.begin();it!=num_units.end();++it) {
            std::cout << it->first << ":\t"<< it->second<<std::endl;
        }
        _minimap_sw_engines[_engine_idx]->align(mod_draft, temp_cons, _minimap_sw_results[_engine_idx].get());
        bool spans = get_query_pos_left(mod_draft, temp_cons, patchb, patche, qb, qe);
        // Check if periodic (no leftover)
        if(qb>=qe) {return final_cons;}
        // use clippings to figure out patching if no arm spans the patch
        bool should_use_soft = (_num_internal==0 && num_units.empty());
        // Verify that no P/S spans the run
        if (should_use_soft) {
            for (auto p: _pends) {if (p+cMarkerLen > patche) {should_use_soft=false; std::cout << "P spans.no clip "<<p<<std::endl;break;}}
        }
        if (should_use_soft) {
            for (auto s: _sbegs) {if (s+cMarkerLen < patchb) {should_use_soft=false; std::cout << "S spans.no clip "<<s<<std::endl;break;}}
        }
        if (should_use_soft && patchb>cMarkerLen && patche<_draft.get_seq_size()+cMarkerLen) {
            std::cout << "Clip patching\n";
            update_cliprun_support();
            for (auto it=num_units.begin();it!=num_units.end();++it) {
                std::cout << it->first << ":\t"<< it->second<<std::endl;
            }
        }
        std::string unit;   
        auto run_len = longest_periodic_pref(temp_cons.substr(qb,qe-qb+1),unit,true);
        std::cout << "qb:"<<qb<<" qe:"<<qe<<" u:"<<unit<<" runlen:"<<run_len<<std::endl;
        if (run_len==(qe-qb+1) && (run_len%unit.size()==0) && unit.size()==patch_unitlen && !num_units.empty()) { // patch
            auto x = std::max_element(num_units.begin(), num_units.end(),
                [](const std::pair<UINT32,UINT16>& p1, const std::pair<UINT32,UINT16>& p2) {
                    return p1.second < p2.second; });
            auto chosen_nu = x->first;
            // Check if it is unique
            bool occ = false;
            for (auto it=num_units.begin(); it !=num_units.end();++it) {
                if (it->second==x->second){
                    if (occ) {std::cout << "NOT Patching as max count is not unique\n"; return final_cons;} //already an occ found; Not unique=>DO not patch
                    occ = true;
                }
            }
            std::cout << temp_cons.substr(0,qb) <<"  " << chosen_nu << "  "<<temp_cons.substr(qe+1)<<std::endl;
            final_cons = temp_cons.substr(0,qb);
            for (auto i = 0; i < chosen_nu; ++i) {
                final_cons += unit;
            }
            final_cons += temp_cons.substr(qe+1);
        }
        std::cout << final_cons<<std::endl;
        return final_cons;
    };

    auto patch_pre= [&]() -> void {
         if (_pre_arms.size()>0){
            // sort in decreasing order of lengths
            std::vector<uint32_t> rank;
            rank.reserve(_pre_arms.size());
            for (uint32_t i = 0; i < _pre_arms.size(); ++i) {
                rank.emplace_back(i);
            }
            std::stable_sort(rank.begin(), rank.end(), [&](uint32_t l, uint32_t r) {
                return _pre_arms[l].get_seq_size() > _pre_arms[r].get_seq_size(); });
            bool spans = true;
            for (UINT i = 0; i < _pre_arms.size(); ++i)   {
                auto ind1 = rank[i];
                std::string modified = cHead+_pre_arms[ind1].unpack();
                std::cout << i << " PatchSupp P"<<ind1<<std::endl;
                spans = update_run_support(AlnType::P, modified,_pre_awt[ind1]); 
                if (!spans) {break;}            
            }
        }
     };
     auto patch_suf= [&]() -> void {
         if (_suf_arms.size()>0){
            // sort in decreasing order of lengths
            std::vector<uint32_t> rank;
            rank.reserve(_suf_arms.size());
            for (uint32_t i = 0; i < _suf_arms.size(); ++i) {
                rank.emplace_back(i);
            }
            std::stable_sort(rank.begin(), rank.end(), [&](uint32_t l, uint32_t r) {
                return _suf_arms[l].get_seq_size() > _suf_arms[r].get_seq_size(); });
            bool spans = true;
            for (UINT i = 0; i < _suf_arms.size(); ++i)   {
                auto ind1 = rank[i];
                std::string modified = _suf_arms[ind1].unpack()+cTail;
                std::cout << i << " PatchSupp S"<<ind1<<std::endl;
                spans = update_run_support(AlnType::S, modified,_suf_awt[ind1]); 
                if (!spans) {break;}            
            }
        }
     };

    /////////////////////// Start with INTERNAL /////////////////////////
    // Check if any internal can be representative
    std::cout << "+++++++++++++++++++++++++++++++\t"<<_draft.unpack()<<std::endl;
    std::cout << _internal_arms.size() << " " << _pre_arms.size() << " "<< _suf_arms.size()<<std::endl;

    bool int_found = false;
    std::string irep = "";
    INT ichosen = -1;
    std::vector<std::string> iclusters;
    std::vector<std::vector<UINT16>> iwt;
    std::vector<std::vector<UINT16>> onlyiwt;
    std::vector<std::vector<UINT16>> icov;
    std::vector<UINT16>isupport;
    std::vector<UINT16>p_isupport;
    std::vector<UINT16>s_isupport;
    iclusters.reserve(_internal_arms.size());
    iwt.reserve(_internal_arms.size());
    onlyiwt.reserve(_internal_arms.size());
    isupport.reserve(_internal_arms.size());
    p_isupport.reserve(_internal_arms.size());
    s_isupport.reserve(_internal_arms.size());

    std::vector<uint32_t> ptrimmed_lens;
    std::vector<uint32_t> strimmed_lens;
    ptrimmed_lens.reserve(_pre_arms.size());
    strimmed_lens.reserve(_suf_arms.size());

    std::string unit;
    for (uint32_t i = 0; i < _pre_arms.size(); ++i) {
        auto tl = longest_periodic_pref(_pre_arms[i].unpack(),unit,false);
        ptrimmed_lens.emplace_back(_pre_arms[i].get_seq_size()-tl);
    }
    for (uint32_t i = 0; i < _suf_arms.size(); ++i) {
        auto tl = longest_periodic_pref(_suf_arms[i].unpack(),unit,true);
        strimmed_lens.emplace_back(_suf_arms[i].get_seq_size()-tl);
    }

    if (_internal_arms.size()>0) {
        for (UINT32 i = 0; i < _internal_arms.size(); ++i) {
            bool new_cluster = true;
            //std::cout << "INT arm " << i <<" " << _internal_arms[i].unpack()<<std::endl;
            for (UINT32 j = 0; j < iclusters.size(); ++j) {
                if (is_supported(AlnType::I,iclusters[j],cHead+_internal_arms[i].unpack()+cTail)) {
                    isupport[j]+=_internal_awt[i];
                    //std::cout <<"supported "<< j <<  " "<<iclusters[j] <<std::endl;
                    new_cluster = false;
                    for (auto k=0; k<iclusters[j].size(); ++k) {iwt[j][k]+=_internal_awt[i]; onlyiwt[j][k]+=_internal_awt[i];}
                    break;
                }
                for (auto k=0; k<iclusters[j].size(); ++k) {icov[j][k]+=_internal_awt[i];}
            }
            if (new_cluster) {
                //std::cout <<"New clutser\n";
                iclusters.push_back(cHead+_internal_arms[i].unpack()+cTail);
                iwt.emplace_back(std::vector<UINT16>(_internal_arms[i].get_seq_size()+(2*cMarkerLen),_internal_awt[i]));
                onlyiwt.emplace_back(std::vector<UINT16>(_internal_arms[i].get_seq_size()+(2*cMarkerLen),_internal_awt[i]));
                icov.emplace_back(std::vector<UINT16>(_internal_arms[i].get_seq_size()+(2*cMarkerLen),_internal_awt[i]));
                isupport.push_back(_internal_awt[i]);
                p_isupport.push_back(0);
                s_isupport.push_back(0);
            }
        }
        if (need_patch) { // Count the units in each internal arm
            for (UINT32 j = 0; j < iclusters.size(); ++j) {
                std::cout << "PatchSupp I " << j <<std::endl; 
                update_run_support(AlnType::I, iclusters[j],isupport[j]);                
            }
        }
        for (UINT32 i = 0; i < _pre_arms.size(); ++i) {
            std::string modified = cHead+_pre_arms[i].unpack();
            auto plen = cMarkerLen+ptrimmed_lens[i];
            for (UINT32 j = 0; j < iclusters.size(); ++j) {
                bool should_update = is_supported(AlnType::P,iclusters[j],modified);
                if (should_update) {
                    p_isupport[j]+=_pre_awt[i];isupport[j]+=_pre_awt[i];
                    for (INT k=0; k<plen; ++k) {iwt[j][k]+=_pre_awt[i];}
                }
                for (INT k=0; k<plen; ++k) {if (k < icov[j].size()){icov[j][k]+=_pre_awt[i];}}
            }
        }
        for (UINT32 i = 0; i < _suf_arms.size(); ++i) {
            std::string modified = _suf_arms[i].unpack()+cTail;
            auto slen = cMarkerLen+strimmed_lens[i];
            for (UINT32 j = 0; j < iclusters.size(); ++j) {
                bool should_update = is_supported(AlnType::S,iclusters[j],modified);
                if (should_update) {
                    s_isupport[j]+=_suf_awt[i];isupport[j]+=_suf_awt[i];
                    for (INT k=1; k<=slen; ++k) {iwt[j][iwt[j].size()-k]+=_suf_awt[i];}
                }
                for (INT k=1; k<=slen; ++k) {if (k <= icov[j].size()){icov[j][icov[j].size()-k]+=_suf_awt[i];}}
            }
        }
        
        // Find internal arms with max support
        auto result = std::max_element(isupport.begin(), isupport.end());
        auto max_supp_id = std::distance(isupport.begin(), result);

        std::cout << "***** IMax supp Id: " << max_supp_id << " " << isupport[max_supp_id] << " " << p_isupport[max_supp_id] << " " << s_isupport[max_supp_id]<<std::endl;
        std::cout << "I: " << iclusters[max_supp_id] << std::endl;
        for (auto i=0; i<iclusters[max_supp_id].size() ; ++i) {
            std::cout << iclusters[max_supp_id][i]<<":"<<iwt[max_supp_id][i]<< "/"<<icov[max_supp_id][i]<<" ";
        }
        std::cout << std::endl<< "*******"<<std::endl;
        for (auto i=0; i<iclusters.size() ; ++i){ 
            std::cout << iclusters[i] << std::endl; 
            for (auto j=0; j<iclusters[i].size() ; ++j) {
                std::cout << iclusters[i][j] <<":"<<iwt[i][j]<< "/"<<icov[i][j]<<" ";
            }
            std::cout << std::endl;
        }
        // If support is from more than half, it represents the arms in the window
        bool is_all_supp = true;
        for (auto k=0; k<iclusters[max_supp_id].size(); ++k) {if (2*iwt[max_supp_id][k]< icov[max_supp_id][k]) {std::cout <<"Low I supp at "<<k<<std::endl; is_all_supp=false;break;}}
        bool should_trust = !has_long_indel(_internal_cigs[max_supp_id]);
        //if (2*p_isupport[max_supp_id] > _num_pre && 2*s_isupport[max_supp_id] > _num_suf) {
        if (!is_sh_low() && is_all_supp && should_trust) {
            int_found = true;
            if (need_patch) {
                patch_pre;
                patch_suf;
                set_marked_consensus(patch_it(iclusters[max_supp_id]));
                std::cout << "IC\t" << iclusters[max_supp_id]<<std::endl;
            }
            else {
                set_marked_consensus(iclusters[max_supp_id]);
                std::cout << "IC\t" << iclusters[max_supp_id]<<std::endl;
            }
            
            info = "PSI";
        }
        //if (_internal_arms.size()==1 && (std::max(_longest_suf_len,_longest_pre_len) > _longest_int_len)) { // internal is not fit to be the backbone
        //if (std::abs((INT64)mod_draft.size()-(INT64)iclusters[max_supp_id].size()) >= RELIABLE_DRAFT_STH) {
        if (!should_trust || is_superlow_ps()) { // We should not trust internal arm with long del
            std::cout<<"Using draft\n";
            irep = cHead+_draft.unpack()+cTail;
        }
        else {
            irep = iclusters[max_supp_id];
            ichosen = max_supp_id;
        }   
    }
    std::vector<std::string> pclusters;
    std::vector<std::string> sclusters;
    std::vector<std::vector<UINT16>> pwt;
    std::vector<std::vector<UINT16>> swt;
    //std::vector<std::vector<UINT16>> pcov;
    //std::vector<std::vector<UINT16>> scov;
    //INT pchosen = -1;
    //INT schosen = -1;
    //std::vector<UINT16>psupport;
    //std::vector<UINT16>ssupport;
    pclusters.reserve(_pre_arms.size());
    //psupport.reserve(_pre_arms.size());
    sclusters.reserve(_suf_arms.size());
    //ssupport.reserve(_suf_arms.size());
     auto get_pre_clusters= [&]() -> void {
        if (_pre_arms.size()>0){
            // sort in decreasing order of lengths
            std::vector<uint32_t> rank;
            rank.reserve(_pre_arms.size());
            std::string unit;
            for (uint32_t i = 0; i < _pre_arms.size(); ++i) {
                rank.emplace_back(i);
            }
            std::stable_sort(rank.begin(), rank.end(), [&](uint32_t l, uint32_t r) {
                return ptrimmed_lens[l] > ptrimmed_lens[r]; });
            for (UINT i = 0; i < _pre_arms.size(); ++i)   {
                auto ind1 = rank[i];
                if (ptrimmed_lens[ind1]==0) {continue;}
                bool new_cluster = true;
                std::string modified = cHead+_pre_arms[ind1].unpack(0,ptrimmed_lens[ind1]);
                auto plen = cMarkerLen+ptrimmed_lens[ind1];
                for (UINT32 j = 0; j < pclusters.size(); ++j) {
                    if (is_supported(AlnType::P,pclusters[j],modified)) {
                        for (auto k=0; k<plen; ++k) {pwt[j][k]+=_pre_awt[ind1];}
                        new_cluster = false;
                        break;
                    }
                }
                if (new_cluster) {
                    pclusters.push_back(modified);
                    pwt.emplace_back(std::vector<UINT16>(plen,_pre_awt[ind1]));
                }               
            }
        }
        if (pclusters.size()>0) {
            std::cout << std::endl<< "*******"<<std::endl;
            for (auto i=0; i<pclusters.size() ; ++i){ 
                std::cout << pclusters[i] << std::endl; 
                for (auto j=0; j<pclusters[i].size() ; ++j) {
                    std::cout << j<<":"<<pclusters[i][j] <<":"<<pwt[i][j]<< " ";
                }
                std::cout << std::endl;
            }
        }
        else {std::cout << "No Untrimmed PREF\n";}
        
     };
     auto get_suf_clusters= [&]() -> void {
        if (_suf_arms.size()>0){
            // sort in decreasing order of lengths
            std::vector<uint32_t> rank;
            rank.reserve(_suf_arms.size());
            for (uint32_t i = 0; i < _suf_arms.size(); ++i) {
                rank.emplace_back(i);
            }
            std::stable_sort(rank.begin(), rank.end(), [&](uint32_t l, uint32_t r) {
                return strimmed_lens[l] > strimmed_lens[r]; });
            for (UINT i = 0; i < _suf_arms.size(); ++i)   {
                auto ind1 = rank[i];
                if (strimmed_lens[ind1]==0) {continue;}
                bool new_cluster = true;
                auto slen = cMarkerLen+strimmed_lens[ind1];
                std::string modified = _suf_arms[ind1].unpack(_suf_arms[ind1].get_seq_size()-strimmed_lens[ind1])+cTail;
                
                for (UINT32 j = 0; j < sclusters.size(); ++j) {
                    if (is_supported(AlnType::S,sclusters[j],modified)) {
                        new_cluster = false;
                        for (INT k=1; k<=slen; ++k) {swt[j][swt[j].size()-k]+=_suf_awt[ind1];}
                        break;
                    }
                }
                if (new_cluster) {
                    sclusters.push_back(modified);
                    swt.emplace_back(std::vector<UINT16>(slen,_suf_awt[ind1]));
                }           
            }
        }
        if (sclusters.size()>0) {
            std::cout << std::endl<< "*******"<<std::endl;
            for (auto i=0; i<sclusters.size() ; ++i){ 
                std::cout << "S"<<i<< " "<<sclusters.size()<<" "<<sclusters[i] << std::endl; 
                for (auto j=0; j<sclusters[i].size() ; ++j) {
                    std::cout << j<<":"<<sclusters[i][j] <<":"<<swt[i][j]<<" ";
                }
                std::cout << std::endl;
            }
        }
        else {std::cout << "No Untrimmed SUFF\n";}
        
     };
    if (is_sh_low() || !int_found) { // if low coverage window or good rep from internal was not found  

        // Find representative prefix
        if (_pre_arms.size()>0) {
           get_pre_clusters();
         }
        
        // Find representative suffix
        if (_suf_arms.size()>0) {
            get_suf_clusters();
        }


        // Combine
        const std::string& bb = (irep.empty()) ? (cHead+_draft.unpack()+cTail) : (irep);
        
        // Check whether all clusters should be used;
        std::string cons  = "";
        /*
        if (!is_sh_low() && is_s_good && is_p_good) {
            std::vector<UINT16> dummy;
            const std::vector<UINT16>& iw = (irep.empty()) ? (dummy) : (iwt[ichosen]);
            const std::vector<UINT16>& pw = (prep.empty()) ? (dummy) : (pwt[pchosen]);
            const std::vector<UINT16>& sw = (srep.empty()) ? (dummy) : (swt[schosen]);
            cons  = get_SW_pileup_consensus(bb, irep, iw, prep, pw, srep, sw); 
            info = "PSR";           
        }
        else { // Use all
        */
            std::vector<std::vector<UINT16>> dummy;
            //const std::vector<std::vector<UINT16>>& iw = (irep.empty()) ? (dummy) : (onlyiwt);
            //const std::vector<std::vector<UINT16>>& pw = (prep.empty()) ? (dummy) : (pwt);
            //const std::vector<std::vector<UINT16>>& sw = (srep.empty()) ? (dummy) : (swt);

            cons  = get_SW_pileup_consensus(bb, iclusters, onlyiwt, pclusters, pwt, sclusters, swt);
            info = "PSA";

        //}
        if (need_patch) {
            patch_pre();
            patch_suf();
            set_marked_consensus(patch_it(cons));

            std::cout << info<<"\t" << cons<<std::endl;
        }
        else {
            set_marked_consensus(cons);
            std::cout << info<<"\t" << cons<<std::endl;
        }
    }
}



void Window::generate_consensus_long()
{
    auto engine_idx = _engine_idx;
    std::string consensus = "";
    auto graph = spoa::createGraph();
    bool arms_added = false;
    auto num_non_empty_arms = _num_internal + _num_pre + _num_suf;
    UINT32 num_added_arms = 0;
    UINT32 max_clus_id = 0;
    auto dsize = INT64(_draft.get_seq_size());
    std::vector<UINT32> clus_ids(_internal_arms.size(),0); // index of the cluser towhich an armbelongs

///// Cluster according to length //////////////////////////////////////
    std::cout << "==========\t"<<_draft.unpack()<<std::endl;
    std::cout << _internal_arms.size() << " " << _pre_arms.size() << " "<< _suf_arms.size()<<std::endl;
    // sort in increasing order of the cigars of internal arms; 
    std::vector<UINT32> rank;
    rank.reserve(_internal_arms.size());
    for (uint32_t i = 0; i < _internal_arms.size(); ++i) {
        rank.emplace_back(i);
    }

    std::vector<UINT32> rank2;
    rank2.reserve(_internal_arms.size());
    for (uint32_t i = 0; i < _internal_arms.size(); ++i) {
        rank2.emplace_back(i);
    }
    
    INT64 chosen_cluslen = 0;
    if (_internal_arms.size()>0) { 
        /////////////////////////////////// DO CLUSTERING
        std::stable_sort(rank.begin(), rank.end(), [&](uint32_t l, uint32_t r) {
            return std::stoi(_internal_cigs[l]) < std::stoi(_internal_cigs[r]); });

        std::stable_sort(rank2.begin(), rank2.end(), [&](uint32_t l, uint32_t r) {
            return _areflen[l] >=_areflen[r]; });

        std::vector<INT64> cluslens;
        std::vector<UINT32> counts;  
        counts.emplace_back(0); //those for cigar=0; clus_id starts from 1
        cluslens.emplace_back(0);
        UINT32 cig_id = 0;
        UINT32 curr_centre = 0;
        UINT32 curr_min = 0;
        UINT32 curr_digitind = 0;
        // Find ins and del indices
        INT64 ins_bind = _internal_arms.size();
        INT64 del_eind = -1;
        for (UINT r = 0; r < _internal_arms.size(); ++r) {
            UINT i = rank[r];
            auto cig_len = std::stoi(_internal_cigs[i]);
            if (cig_len<0){++del_eind;}
            else if (cig_len>0) {ins_bind=r; break;} // first +, break;
            else {++counts[0];} // =0
        }
        std::cout<<"\n INSBIND:DELEIND " << ins_bind <<":"<<del_eind<<std::endl;
        auto cluster= [&](UINT r) -> void {
            UINT i = rank[r];
            auto cig_len = std::abs(std::stoi(_internal_cigs[i]));
            std::cout << r<<" "<<i<<" CL:"<<cig_len<<" radind:"<<curr_digitind << " curr_cen:"<<curr_centre <<" rad:"<<cRange[curr_digitind] << " range:"<<curr_centre+cRange[curr_digitind]<<std::endl;
            if (curr_centre+cRange[curr_digitind] >= cig_len) { // add to this cluster
                curr_centre = (curr_min+cig_len)/2;
                ++counts[cig_id];
                cluslens[cig_id] = curr_centre;
                std::cout << "Added to"<<cig_id<<std::endl;
            }
            else { // new cluster
                ++cig_id;
                curr_centre = cig_len;
                curr_min = cig_len;
                counts.emplace_back(1);
                cluslens.emplace_back(curr_centre); 
                std::cout << "New CLUSTER \n";   
            }                    
            clus_ids[i] = cig_id;    
            std::cout << r <<":"<<i<<":"<<cig_len<<":"<<_areflen[i]<<std::endl;                
            while (curr_centre>=cNumDigits[curr_digitind+1]) {++curr_digitind;std::cout<< "rad_ind becomes "<<curr_digitind<<std::endl;}   
        };
        // Process INS clusters
        for (INT64 r = ins_bind; r < _internal_arms.size(); ++r) {
            std::cout <<"INS\n";
            cluster(r);
        }
        curr_centre = 0;
        curr_min = 0;
        curr_digitind = 0;
         // Process DEL clusters
        for (INT64 r = del_eind; r >= 0; --r) {
            std::cout <<"DEL\n";
            cluster(r);
        }

        /////////////////////////////////// SELECT CLUSTER

        if (_long_del>0) { // this window has a long deletion; choosethe cluster  closest to length of long del
            UINT least_dist = 1000000; // a very large number
            for (UINT c=1; c<cluslens.size(); ++c) {
                if (std::abs(cluslens[c]-_long_del) < least_dist) {
                    least_dist = std::abs(cluslens[c]-_long_del);
                    max_clus_id = c;
                }
            }
            chosen_cluslen = cluslens[max_clus_id];
            num_added_arms = counts[max_clus_id];
            std::cout <<" LONG INDEL SELECTED "<< max_clus_id<<std::endl;
        }
        else{           
            // For each cluster, find the length of maximum aligned ref span amongst the reads belonging to that cluster
            std::vector<UINT64> maxalens(cluslens.size(),0);       
            for (UINT i = 0; i < _internal_arms.size(); ++i)
            {
                if (maxalens[clus_ids[i]] < _areflen[i]) {
                    maxalens[clus_ids[i]] = _areflen[i];
                }
            } 

            std::vector<UINT64> minalens(cluslens.size(),0);       
            for (UINT i = 0; i < _internal_arms.size(); ++i)
            {
                if (minalens[clus_ids[i]]==0 || minalens[clus_ids[i]] > _areflen[i]) {
                    minalens[clus_ids[i]] = _areflen[i];
                }
            } 



            for (UINT r = 0; r < _internal_arms.size(); ++r) {
                UINT i = rank2[r];
                std::cout << "alen "<<r <<":"<<i<<":"<<_areflen[i]<<":"<<clus_ids[i]<<std::endl;
            }

            // Find heaviest cluster   
            UINT max_cnt = 0;
            for (UINT c=0; c<cluslens.size(); ++c) {
                auto cnt = counts[c];
                if (cluslens[c]*10 > VERY_LONG_INDEL) {
                    // increase count by num of pre/suf (assuming their aln broke because of this indel)
                    cnt += (_num_pre+_num_suf);
                }
                if (cnt > max_cnt) {
                    max_cnt = cnt;
                    max_clus_id = c;
                }
                std::cout <<c << " CNT" << counts[c]<< " centre:" << cluslens[c]<<" ALEN:"<< maxalens[c]<<std::endl;
            }
            std::cout << "\nMax:"<<max_clus_id<<" Th:"<<std::ceil(max_cnt*LONG_CLUS_SECOND_TH)<<std::endl;
            // All clusters with strength about 40% of the heaviset clusters are considered; choose the one with max aligned len
            for (UINT c=0; c<cluslens.size(); ++c) {
                if (counts[c] >= std::ceil(max_cnt*LONG_CLUS_SECOND_TH)) { // consider this cluster
                    std::cout <<c << " cnt" << counts[c]<<" MAXalen:"<< maxalens[c]<<" MINalen:"<< minalens[c]<<" Avgalen:"<<(minalens[c]+ maxalens[c])/2<<std::endl;
                    if (maxalens[c] > maxalens[max_clus_id]) {
                        //max_clus_id = c;
                        //std::cout <<" Max becomes this\n";
                    }
                }
            }
            num_added_arms = counts[max_clus_id];
            chosen_cluslen = cluslens[max_clus_id];
            std::cout <<"SELECTED "<< max_clus_id<<std::endl;
            // See if empty should be considered as a cluster
            if(_num_empty > 0 && _num_empty > max_cnt) {
                //if (_num_empty >= std::ceil(max_cnt*LONG_CLUS_SECOND_TH)) { // consider this cluster
                    std::cout << "Empty alen:" << _empty_maxalen<<std::endl;
                    //if (_empty_maxalen > maxalens[max_clus_id]) {
                         set_consensus(""); // empty sequence is consensus
                         std::cout <<"Choosing empty\n";
                        return;
                    //}
                //}
            }

            info = std::to_string(num_added_arms) + "_"+ std::to_string(maxalens[max_clus_id]);
        }        
    }
    
    // If too few, use all alignment
    if (num_added_arms<3) {
        max_clus_id = 1;
        num_added_arms = _internal_arms.size();
        for (UINT i = 0; i < _internal_arms.size(); ++i) {
            clus_ids[i] = 1;
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////

    // internal
    _alignment_engines_long[engine_idx]->changeAlignType(spoa::AlignmentType::kNW);
    //  add draft
    std::string modified_str = cHead + _draft.unpack() + cTail;
    //std::string modified_str = _draft.unpack();
    auto alignment = _alignment_engines_long[engine_idx]->align(modified_str, graph);
    graph->add_alignment(alignment, modified_str);
    
    
    // Add from longer to shorter
    bool downsample = (_mode==Mode::LSA && num_added_arms>DOWNSAMPLE_TH);
    num_added_arms=0;
    int added=0;

    for (INT ind = _internal_arms.size()-1; ind >= 0; --ind) 
    {
        auto i = rank[ind];
        if (_internal_arms[i].get_seq_size() > 0 && clus_ids[i]==max_clus_id)
        {
            if(!downsample || added%2==0) {
                ++num_added_arms;
                std::string modified_str = cHead + _internal_arms[i].unpack() + cTail;
                //std::string modified_str = _internal_arms[i].unpack();
                arms_added = true;
                auto alignment = _alignment_engines_long[engine_idx]->align(modified_str, graph);
                if (_internal_qual.size() > 0 && _internal_qual[i].size() > 0)
                {
                    std::string q_str = HIGHEST_PHRED+_internal_qual[i]+HIGHEST_PHRED;
                    graph->add_alignment(alignment, modified_str, q_str);
                }
                else
                {
                    graph->add_alignment(alignment, modified_str);
                }
            }
            ++added;
        }
    }
    _alignment_engines_long[engine_idx]->changeAlignType(spoa::AlignmentType::kLOV);
    for (UINT i = 0; i < _pre_arms.size(); ++i)
    {
        if (_pre_arms[i].get_seq_size() > 0)
        {
            std::string modified_str = cHead + _pre_arms[i].unpack();
            //std::string modified_str = _pre_arms[i].unpack();
            arms_added = true;
            auto alignment = _alignment_engines_long[engine_idx]->align(modified_str, graph);
            if (_pre_qual.size() > 0 && _pre_qual[i].size() > 0)
            {
                std::string q_str = HIGHEST_PHRED+_pre_qual[i];
                graph->add_alignment(alignment, modified_str, q_str);
            }
            else
            {
                graph->add_alignment(alignment, modified_str);
            }
        }
    }
    //suffix
    _alignment_engines_long[engine_idx]->changeAlignType(spoa::AlignmentType::kROV);
    for (UINT i = 0; i < _suf_arms.size(); ++i)
    {
        if (_suf_arms[i].get_seq_size() > 0)
        {
            std::string modified_str = _suf_arms[i].unpack() + cTail;
            //std::string modified_str = _suf_arms[i].unpack();
            arms_added = true;
            auto alignment = _alignment_engines_long[engine_idx]->align(modified_str, graph);
            if (_suf_qual.size() > 0 && _suf_qual[i].size() > 0)
            {
                std::string q_str = _suf_qual[i]+HIGHEST_PHRED;
                graph->add_alignment(alignment, modified_str, q_str);
            }
            else
            {
                graph->add_alignment(alignment, modified_str);
            }
        }
    }
    if (arms_added)
    {
        std::vector<UINT32> dst;
        //consensus = graph->generate_consensus_custom(dst);
        consensus = graph->generate_consensus(dst);
        set_marked_consensus(curate_long(consensus, dst, num_added_arms));
    }
    else
    {
        set_consensus(_draft.unpack());
    }
}

std::string Window::curate_long(const std::string &con, const std::vector<UINT32> &dst, const UINT32 num_arms)
{
    //auto num_arms = _num_internal+_num_pre+_num_suf+_num_empty;
    auto conssize = con.size();
    std::string curated_consensus;
    curated_consensus.reserve(conssize);

    UINT cov_thres = (UINT)std::floor(num_arms * LONG_TH);
    INT beg = 0;
    INT lst = conssize - 1;
    for (UINT i = 0; i < conssize; ++i)
    {
        if (dst[i] >= cov_thres)
        {
            beg = i;
            break;
        }
    }
    for (INT i = conssize - 1; i >= 0; --i)
    {
        if (dst[i] >= cov_thres)
        {
            lst = i;
            break;
        }
    }
    curated_consensus = con.substr(beg, lst - beg + 1);
    curated_consensus.shrink_to_fit();
    return curated_consensus;
}

// Returns whether seq spans refe; Assumes SW has been called; refb < refe; refb is first base of run; refe is the last base after the run
// At beg: extra bases in query before beg are included; dels in query=> base after the del is the first
// At end: extra bases after the end not included; dels at end=> base before the del is the last
bool Window::get_query_pos_left(const std::string& bb, const std::string& seq, const UINT32 refb, const UINT32 refe, UINT32& qb, UINT32& qe) {   
    //_minimap_sw_engines[_engine_idx]->align(bb, seq, _minimap_sw_results[_engine_idx].get());
    assert(seq.size()>0);
    // we get the results now just process the CIGAR for updates
    UINT32 current_reference_index = 0;
    UINT32 current_query_index = 0;
    
    UINT32 current_cigar;
    UINT32 get_op;
    UINT32 get_oplen;
    
    UINT8 current_base;
    bool does_span = false;
    bool bset = false;
    
    std::string dbg = "";
    for(UINT32 i = 0; i < _minimap_sw_results[_engine_idx]->cigar.size(); i++) {        
        UINT32 current_cigar = _minimap_sw_results[_engine_idx]->cigar[i];
        get_op = bam_cigar_op(current_cigar);
        get_oplen = bam_cigar_oplen(current_cigar);
        if((bam_cigar_type(get_op) & 3) == 3) { // consume query and reference
            while(get_oplen > 0) {
                // update base
                if (!bset && current_reference_index==refb) {qb=current_query_index;bset=true;}
                if (current_reference_index==refe) {dbg+=seq[current_query_index];qe=current_query_index;does_span=true;break;}
                dbg+=seq[current_query_index];
                get_oplen--;
                current_reference_index++;
                current_query_index++;
            }
            if (does_span) {break;}
        } else if(bam_cigar_type(get_op) & 2) { // consume reference only: deletions
            while(get_oplen > 0) {
                if (!bset && current_reference_index==refb) {qb=current_query_index;bset=true;}
                if (current_reference_index==refe) {qe=(current_query_index>0)?(current_query_index-1):(0);;does_span=true;break;}
                get_oplen--;
                current_reference_index++;
                dbg+='-';
            }
            if (does_span) {break;}
        } else if(bam_cigar_type(get_op) & 1) { // consume query only: insertions
            if (!bset && current_reference_index==refb) {qb=current_query_index;bset=true;}
            // get substring of the inserted bases
            std::string inserted = "";
            while(get_oplen > 0) {
                inserted += seq[current_query_index++];
                get_oplen--;
            }
            // insert it to current_reference_index (right base after ins)
            if(current_reference_index >= 1) {
                dbg+='*';
            } 
            std::cout << "INS: " << current_reference_index<< " " << inserted << current_query_index<<std::endl;
        }
    }
    if (qe==seq.size()-1) {does_span=false;} // query ends in the last base of the run; can't be trusted
    std::cout <<bb<<std::endl;
    std::vector<UINT16> diff(DIST+1,0);
        sw_diff(bb, seq, diff);
    std::cout << "**"<<dbg<<std::endl;
    return does_span;
    
}

// Returns whether seq spans refb; Assumes SW has been called; refb < refe; refb is first base of run; refe is the last base after the run
// At beg: extra bases in query before beg are included; dels in query=> base after the del is the first
// At end: extra bases after the end not included; dels at end=> base before the del is the last
bool Window::get_query_pos_right(const std::string& bb, const std::string& seq, const UINT32 refb, const UINT32 refe, UINT32& qb, UINT32& qe) {   
    assert(seq.size()>0);
    //_minimap_sw_engines[_engine_idx]->align_suffix_rev_cigar(bb, seq, _minimap_sw_results[_engine_idx].get());
    
    // we get the results now just process the CIGAR for updates
    // cigar is reversed from end to start
    INT32 current_reference_index = bb.size() - 1;
    INT32 current_query_index = seq.size() - 1;
       
    UINT32 current_cigar;
    UINT32 get_op;
    UINT32 get_oplen;
    UINT32 orig_oplen;
    
    UINT8 current_base;
    bool does_span = false;
    bool bset = false;

   std::string dbg="";
    for(UINT32 i = 0; i < _minimap_sw_results[_engine_idx]->cigar.size(); i++) {
        current_cigar = _minimap_sw_results[_engine_idx]->cigar[i];
        get_op = bam_cigar_op(current_cigar);
        get_oplen = bam_cigar_oplen(current_cigar);
        if((bam_cigar_type(get_op) & 3) == 3) { // consume query and reference
            bool should_break = false;
            while(get_oplen > 0) {
                if ( current_reference_index>=refb) {qb=current_query_index;}
                if ( current_reference_index==refb) {does_span=true;} //can't break here; there may be INS in query before
                if ( current_reference_index<refb) {should_break = true;break;}
                if (!bset && current_reference_index<=refe) {qe=current_query_index;bset=true;}
                dbg+=seq[current_query_index];
                get_oplen--;
                current_reference_index--;
                current_query_index--;
            }
            if (should_break) {break;}
        } else if(bam_cigar_type(get_op) & 2) { // consume reference only: deletions
            orig_oplen = get_oplen;
            while(get_oplen > 0) {
                if ( current_reference_index<=refb) {does_span=true;break;} // qb is already set in the last match
                get_oplen--;
                current_reference_index--;
                dbg+='-';
            }
            if (does_span) {break;}
        } else if(bam_cigar_type(get_op) & 1) { // consume query only: insertions
            // get substring of the inserted bases
            std::string inserted = "";
            while(get_oplen > 0) {
                inserted += seq[current_query_index--];
                get_oplen--;
            }
            //std::reverse(inserted.begin(), inserted.end());
            //std::cout << "INS: " << current_reference_index<< " " << inserted << current_query_index<<std::endl;
            if ( current_reference_index<refb) {qb=current_query_index+1;break;}            
        }
    }
    if (qb==0) {does_span=false;} // query begins in the first base of the run; can't be trusted

    std::reverse(dbg.begin(), dbg.end());
    std::string bl(bb.size()-dbg.size(),'-');
    std::cout << bb<<std::endl;
    std::vector<UINT16> cov(bb.size(),0);
    //std::vector<UINT16> num_dels(backbone.size(),0);
    std::vector<std::unordered_map<UINT16, UINT16>> dels(bb.size(), std::unordered_map<UINT16, UINT16>());
    std::vector<std::unordered_map<std::string,UINT16>> ins(bb.size(),std::unordered_map<std::string,UINT16>());
    std::vector<UINT16>sw(seq.size(),1);
    std::vector<std::vector<UINT16>> bases_cover(bb.size(), std::vector<UINT16>(cCode.size(), 0));
    sw_align_suffix(bb, seq, sw, cov, bases_cover, dels, ins);
    std::cout << "**"<<bl+dbg << std::endl;
    return does_span;

}


void Window::sw_diff(const std::string& bb, const std::string& seq, std::vector<UINT16>& diff) {   
    _minimap_sw_engines[_engine_idx]->align(bb, seq, _minimap_sw_results[_engine_idx].get());
    
    // we get the results now just process the CIGAR for updates
    UINT32 current_reference_index = 0;
    UINT32 current_query_index = 0;
    
    UINT32 current_cigar;
    UINT32 get_op;
    UINT32 get_oplen;
    
    UINT8 current_base;
    
    std::string dbg = "";
    for(UINT32 i = 0; i < _minimap_sw_results[_engine_idx]->cigar.size(); i++) {        
        UINT32 current_cigar = _minimap_sw_results[_engine_idx]->cigar[i];
        get_op = bam_cigar_op(current_cigar);
        get_oplen = bam_cigar_oplen(current_cigar);
        if((bam_cigar_type(get_op) & 3) == 3) { // consume query and reference
            while(get_oplen > 0) {
                // update base
                (bb[current_reference_index]==seq[current_query_index]) ? (++diff[M]) : (++diff[X]);
                dbg+=seq[current_query_index];
                get_oplen--;
                current_reference_index++;
                current_query_index++;
            }
        } else if(bam_cigar_type(get_op) & 2) { // consume reference only: deletions
            ++diff[D];
            diff[DIST] += get_oplen;
            while(get_oplen > 0) {
                get_oplen--;
                current_reference_index++;
                dbg+='-';
            }
        } else if(bam_cigar_type(get_op) & 1) { // consume query only: insertions
            ++diff[I];
            diff[DIST] += get_oplen;
            // get substring of the inserted bases
            std::string inserted = "";
            while(get_oplen > 0) {
                inserted += seq[current_query_index++];
                get_oplen--;
            }
            // insert it to current_reference_index (right base after ins)
            if(current_reference_index >= 1) {
                dbg+='*';
            } 
        }
    }
    std::cout << dbg<<std::endl;
    
}


void Window::sw_align_internal(const std::string& bb, const std::string& seq, const std::vector<UINT16>& w, std::vector<UINT16>& cov, std::vector<std::vector<UINT16>>& bases, std::vector<std::unordered_map<UINT16, UINT16>> & num_dels, std::vector<std::unordered_map<std::string,UINT16>>& ins) {   
    _minimap_sw_engines[_engine_idx]->align(bb, seq, _minimap_sw_results[_engine_idx].get());
    
    // we get the results now just process the CIGAR for updates
    UINT32 current_reference_index = 0;
    UINT32 current_query_index = 0;
    
    UINT32 current_cigar;
    UINT32 get_op;
    UINT32 get_oplen;
    
    UINT8 current_base;
    
    std::string dbg = "";
    for(UINT32 i = 0; i < _minimap_sw_results[_engine_idx]->cigar.size(); i++) {        
        UINT32 current_cigar = _minimap_sw_results[_engine_idx]->cigar[i];
        get_op = bam_cigar_op(current_cigar);
        get_oplen = bam_cigar_oplen(current_cigar);
        if((bam_cigar_type(get_op) & 3) == 3) { // consume query and reference
            while(get_oplen > 0) {
                cov[current_reference_index]+=w[current_query_index]; // increase coverage
                // update base
                current_base = cNt4Table[seq[current_query_index]];
                bases[current_reference_index][current_base]+=w[current_query_index];
                dbg+=seq[current_query_index];
                get_oplen--;
                current_reference_index++;
                current_query_index++;
            }
        } else if(bam_cigar_type(get_op) & 2) { // consume reference only: deletions
            // deletions after the previous current_reference_index
            num_dels[current_reference_index][get_oplen]+=w[current_query_index]; // wt of del is that of the following base
            while(get_oplen > 0) {
                cov[current_reference_index]+=w[current_query_index];; // increase coverage
                get_oplen--;
                current_reference_index++;
                dbg+='-';
            }
        } else if(bam_cigar_type(get_op) & 1) { // consume query only: insertions
            // get substring of the inserted bases
            std::string inserted = "";
            while(get_oplen > 0) {
                inserted += seq[current_query_index++];
                get_oplen--;
            }
            //std::cout << "INS: " << current_reference_index<< " " << inserted << current_query_index<<std::endl;
            // insert it to current_reference_index (right base after ins)
            if(current_reference_index >= 1) {
                ins[current_reference_index-1][inserted]+=w[current_query_index]; // wt of ins is that of the following base
                dbg+='*';
                //std::cout << "Inserting an insertion of " << inserted << " to index " << current_reference_index - 1 << std::endl;
            } else {
                // should never happen since "J" will match
                //std::cout << "ERROR: J not matching for internal" << std::endl;
            }
        }
    }
    std::cout << dbg<<std::endl;
    
}

//void Window::sw_align_prefix(const UINT32 engine_idx, const std::string& bb, const std::string& seq, std::vector<UINT16>& cov, std::vector<UINT16> & num_dels, std::vector<std::unordered_map<std::string,UINT16>>& ins) {
void Window::sw_align_prefix(const std::string& bb, const std::string& seq, const std::vector<UINT16>& w, std::vector<UINT16>& cov, std::vector<std::vector<UINT16>>& bases, std::vector<std::unordered_map<UINT16, UINT16>> & num_dels, std::vector<std::unordered_map<std::string,UINT16>>& ins) {
    _minimap_sw_engines[_engine_idx]->align_prefix(bb, seq, _minimap_sw_results[_engine_idx].get());
    
    // we get the results now just process the CIGAR for updates
    UINT32 current_reference_index = 0;
    UINT32 current_query_index = 0;
    
    UINT32 current_cigar;
    UINT32 get_op;
    UINT32 get_oplen;
    
    UINT8 current_base;

    std::string dbg = "";
    for(UINT32 i = 0; i < _minimap_sw_results[_engine_idx]->cigar.size(); i++) {
        current_cigar = _minimap_sw_results[_engine_idx]->cigar[i];
        get_op = bam_cigar_op(current_cigar);
        get_oplen = bam_cigar_oplen(current_cigar);
        if((bam_cigar_type(get_op) & 3) == 3) { // consume query and reference
            while(get_oplen > 0) {
                // update base
                current_base = cNt4Table[seq[current_query_index]];
                //if (current_query_index<strim_pos) { // bases in non-trimmed part are reliable
                    bases[current_reference_index][current_base]+=w[current_query_index];
                //}
                dbg+=seq[current_query_index];
                get_oplen--;
                //if (current_query_index<strim_pos) {
                    cov[current_reference_index]+=w[current_query_index];;
                    //}
                //if (current_query_index>=strim_pos) {++trimcov[current_reference_index];}
                current_reference_index++;
                current_query_index++;
            }
        } else if(bam_cigar_type(get_op) & 2) { // consume reference only: deletions
            // deletions starting at current_reference_index; 
            // Wt of del is that of the next/last base in the non-trimmed/trimmed part resp.
            //auto chosenw = (current_query_index<strim_pos) ? (w[current_query_index]) : (w[seq.size()-1]);
            auto chosenw = (current_query_index<seq.size()) ? (w[current_query_index]) : (w[seq.size()-1]);
            num_dels[current_reference_index][get_oplen]+=chosenw;
            while(get_oplen > 0) {
                //if (current_query_index<strim_pos) {
                    cov[current_reference_index]+=chosenw;
                //}
                //if (current_query_index>=strim_pos) {++trimcov[current_reference_index];}
                get_oplen--;
                current_reference_index++;
                dbg+='-';
            }
        } else if(bam_cigar_type(get_op) & 1) { // consume query only: insertions
            // get substring of the inserted bases
            std::string inserted = "";
            while(get_oplen > 0) {
                inserted += seq[current_query_index++];
                get_oplen--;
            }
            // insert it to current_reference_index (left base before ins)
            // Wt of ins is that of the next/last base in the non-trimmed/trimmed part resp.
            auto chosenw = (current_query_index<seq.size()) ? (w[current_query_index]) : (w[seq.size()-1]);
            if(current_reference_index >= 1) {
                ins[current_reference_index-1][inserted]+=chosenw;
                dbg+='*';
                //std::cout << "Inserting an insertion of " << inserted << " to index " << current_reference_index - 1 << std::endl;
            } else {
                // should never happen since "J" will match
                //std::cout << "ERROR: J not matching for prefix" << std::endl;
            }
        }
    }
    std::cout <<dbg<<std::endl;
}

void Window::sw_align_suffix(const std::string& bb, const std::string& seq, const std::vector<UINT16>& w, std::vector<UINT16>& cov, std::vector<std::vector<UINT16>>& bases, std::vector<std::unordered_map<UINT16, UINT16>> & num_dels, std::vector<std::unordered_map<std::string,UINT16>>& ins) {
    _minimap_sw_engines[_engine_idx]->align_suffix_rev_cigar(bb, seq, _minimap_sw_results[_engine_idx].get());
    
    // we get the results now just process the CIGAR for updates
    // cigar is reversed from end to start
    INT32 current_reference_index = bb.size() - 1;
    INT32 current_query_index = seq.size() - 1;
       
    UINT32 current_cigar;
    UINT32 get_op;
    UINT32 get_oplen;
    UINT32 orig_oplen;
    
    UINT8 current_base;

   std::string dbg="";
    for(UINT32 i = 0; i < _minimap_sw_results[_engine_idx]->cigar.size(); i++) {
        current_cigar = _minimap_sw_results[_engine_idx]->cigar[i];
        get_op = bam_cigar_op(current_cigar);
        get_oplen = bam_cigar_oplen(current_cigar);
        if((bam_cigar_type(get_op) & 3) == 3) { // consume query and reference
            while(get_oplen > 0) {
                //if (current_query_index>=strim_pos) {
                    cov[current_reference_index]+=w[current_query_index];
                //}
                //if (current_query_index<strim_pos) {++trimcov[current_reference_index];}
                // update base
                current_base = cNt4Table[seq[current_query_index]];
                // bases in non-trimmed part are reliable
                //if (current_query_index>=strim_pos) {
                    bases[current_reference_index][current_base]+=w[current_query_index];
                //}
                dbg+=seq[current_query_index];
                get_oplen--;
                current_reference_index--;
                current_query_index--;
            }
        } else if(bam_cigar_type(get_op) & 2) { // consume reference only: deletions
            orig_oplen = get_oplen;
            auto chosenw = (current_query_index>0) ? (w[current_query_index]) : (w[0]);
            while(get_oplen > 0) {
                if (current_query_index>=0) {
                    cov[current_reference_index]+=chosenw;
                }
                //if (current_query_index<strim_pos) {++trimcov[current_reference_index];}
                get_oplen--;
                current_reference_index--;
                dbg+='-';
            }
            // deletions are after current_reference_index
            // Wt of del is that of the pvs base (in trimmed first base is the pvs as del is left-aligned)
            num_dels[current_reference_index+1][orig_oplen]+=chosenw;
        } else if(bam_cigar_type(get_op) & 1) { // consume query only: insertions
            // get substring of the inserted bases
            std::string inserted = "";
            while(get_oplen > 0) {
                inserted += seq[current_query_index--];
                get_oplen--;
            }
            std::reverse(inserted.begin(), inserted.end());
            //std::cout << "INS: " << current_reference_index<< " " << inserted << current_query_index<<std::endl;
            // insert it to the base before
            // Wt of ins is that of the pvs base (in trimmed, first base is the pvs as ins is left-aligned)
            if(current_reference_index >= 0) {
                auto chosenw = (current_query_index>0) ? (w[current_query_index]) : (w[0]);
                ins[current_reference_index][inserted]+=chosenw;
            } else {
                //std::cout << "ERROR: Insertion before J should not happen here." << std::endl;
            }
            
        }
    }
    
    std::reverse(dbg.begin(), dbg.end());
    std::string bl(bb.size()-dbg.size(),'-');
    std::cout << bl+dbg << std::endl;

}


} // namespace hypo


