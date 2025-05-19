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
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include "Window.hpp"

namespace hypo
{
    const float Window::_cThresh = 0.4;
    std::vector<std::shared_ptr<spoa::AlignmentEngine>> Window::_alignment_engines;
    std::vector<std::shared_ptr<spoa::AlignmentEngine>> Window::_alignment_engines_long;
    void Window::prepare_for_poa(const ScoreParams& sp, const UINT32 num_threads) {
        for (UINT32 i = 0; i < num_threads; ++i){
            _alignment_engines.emplace_back(spoa::createAlignmentEngine(
                spoa::AlignmentType::kNW, sp.sr_match_score,
                sp.sr_misMatch_score, sp.sr_gap_penalty));
            _alignment_engines_long.emplace_back(spoa::createAlignmentEngine(
                spoa::AlignmentType::kNW, sp.lr_match_score,
                sp.lr_misMatch_score, sp.lr_gap_penalty));
        }
        _alignment_engines.shrink_to_fit();
        _alignment_engines_long.shrink_to_fit();
    }
    
bool Window::compare_arm_sequence(const PackedSeq<2> & left, const PackedSeq<2> & right) {
    if(left.get_seq_size() != right.get_seq_size()) {
        return left.get_seq_size() > right.get_seq_size();
    }
    for(int i = 0; i < left.get_seq_size(); i++) {
        if(left.base_at(i) != right.base_at(i)) {
            return left.base_at(i) < right.base_at(i);
        }
    }
    return false;
}
    
void Window::generate_consensus(const UINT32 engine_idx) {
    std::stable_sort(_internal_arms.begin(), _internal_arms.end(), compare_arm_sequence);
    std::stable_sort(_pre_arms.begin(), _pre_arms.end(), compare_arm_sequence);
    std::stable_sort(_suf_arms.begin(), _suf_arms.end(), compare_arm_sequence);
    
    auto wt = _wtype;
    auto num_non_empty_arms = _num_internal+_num_pre+_num_suf;
    if (_num_empty > num_non_empty_arms) {
        set_consensus(""); // empty sequence is consensus
        set_consensus_2(""); // empty sequence is consensus
    }
    else if (num_non_empty_arms >= 2) { // 1 is draft; there should be at least 2 arms to do poa
        if (wt == WindowType::SHORT) {
            generate_consensus_short(engine_idx);
        }
        else {
            generate_consensus_long(engine_idx);
        }
    }
    else {
        set_consensus(_draft.unpack());
        set_consensus_2(_draft.unpack());
    }  
}

std::ostream &operator<<(std::ostream &os, const Window &wnd) {
        // Write numbers
        os << wnd._num_internal << "\t" << wnd._num_pre << "\t" << wnd._num_suf << "\t" << wnd._num_empty << std::endl;
        // Write draft sequence
        os << "++\t" << wnd._draft.unpack() << std::endl;

        // Write consensus sequence
        os << "++\t" << wnd._consensus << std::endl;
        // Write internal arms
        for (UINT32 i=0; i < wnd._num_internal; ++i) {
            os << wnd._internal_arms[i].unpack() << std::endl;
        }
        // Write pref arms
        for (UINT32 i=0; i < wnd._num_pre; ++i) {
            os << wnd._pre_arms[i].unpack() << std::endl;
        }
        // Write suff arms
        for (UINT32 i=0; i < wnd._num_suf; ++i) {
            os << wnd._suf_arms[i].unpack() << std::endl;
        }
        return os;
    }


void Window::generate_consensus_short(const UINT32 engine_idx) {
    std::string consensus = "";
    auto graph = spoa::createGraph();
    bool arms_added = false;
    
    // internal
    // Avoid draft if an internal arm exists
    _alignment_engines[engine_idx]->changeAlignType(spoa::AlignmentType::kNW);
    if (_internal_arms.size()==0) { // use POA only when no internal arm
        // draft will never have zero size
        //std::string unpacked_str(std::move(_internal_arms[i].unpack()));
        std::string modified_str = cHead + _draft.unpack() + cTail;
        auto alignment = _alignment_engines[engine_idx]->align(modified_str, graph);
        graph->add_alignment(alignment, modified_str);
        // prefix (add in reverse order (Since bam is sorted. last one should be the longest)
        _alignment_engines[engine_idx]->changeAlignType(spoa::AlignmentType::kLOV);
        for (auto it = _pre_arms.rbegin(); it!=_pre_arms.rend(); ++it) {
            if (it->get_seq_size() > 0) {
            //std::string unpacked_str(std::move(it.unpack()));
            std::string modified_str = cHead + it->unpack();
            arms_added = true;
            auto alignment = _alignment_engines[engine_idx]->align(modified_str, graph);
            graph->add_alignment(alignment, modified_str);
            }
        }
        //suffix
        _alignment_engines[engine_idx]->changeAlignType(spoa::AlignmentType::kROV);
        for (auto it = _suf_arms.begin(); it!=_suf_arms.end(); ++it) {
            if (it->get_seq_size() > 0) {
            //std::string unpacked_str(std::move(it.unpack()));
            std::string modified_str = it->unpack() + cTail;
            arms_added = true;
            auto alignment = _alignment_engines[engine_idx]->align(modified_str, graph);
            graph->add_alignment(alignment, modified_str);
            }
        }
        ///*************** FOR DEBUGGING ***************
        #ifdef DEBUG3
        std::vector<std::string> msa;
        graph->generate_multiple_sequence_alignment(msa, true);
        auto seq_num = 0;
        std::cout << "================================="<< seg.get_id() <<std::endl;
        for (const auto &it : msa){
            std::cout << seq_num<< "\t" << it.c_str() << std::endl;
            ++seq_num;
        }
        #endif
        //*/ // /********************************************/
        if (arms_added)
        { 
            consensus = graph->generate_consensus();
            set_marked_consensus(consensus);  
            set_marked_consensus_2(consensus);  
        }
        else {
            set_consensus(_draft.unpack());
            set_consensus_2(_draft.unpack());
        }
    } else {
        std::unordered_map<std::string, int> counter;
        
        int max_counter = 0;
        int max_counter_2 = 0;
        std::string max_str;
        std::string max_str_2;
        
        for (UINT i = 0; i <  _internal_arms.size(); ++i) {
            std::string current_str = _internal_arms[i].unpack();
            
            if(counter.find(current_str) == counter.end()) {
                counter[current_str] = 1;
            } else {
                counter[current_str] = counter[current_str] + 1;
            }
            
            if(counter[current_str] > max_counter) {
                max_counter_2 = max_counter;
                max_str_2 = max_str;
                
                max_counter = counter[current_str];
                max_str = current_str;
            } else if(counter[current_str] > max_counter_2) {
                max_counter_2 = counter[current_str];
                max_str_2 = current_str;
            }
        }
        
        set_consensus(max_str);        
        if(max_counter_2 > 0) set_consensus_2(max_str_2);
        else set_consensus_2(max_str);
    }
}

void Window::generate_consensus_long(const UINT32 engine_idx) {
    std::string consensus = "";
    auto graph = spoa::createGraph();
    bool arms_added = false;
    auto num_non_empty_arms = _num_internal+_num_pre+_num_suf;
    
    std::vector<std::vector<UINT32>> counts;
    counts.reserve(num_non_empty_arms);
    
    // for consistency, use only the internal arms.
    _alignment_engines[engine_idx]->changeAlignType(spoa::AlignmentType::kNW);
    for ( UINT i=0; i <  _internal_arms.size(); ++i) {
        std::string unpacked_str(std::move(_internal_arms[i].unpack()));
        if (unpacked_str.size() > 0) {
            arms_added = true;
            auto alignment = _alignment_engines_long[engine_idx]->align(unpacked_str, graph);
            graph->add_alignment(alignment, unpacked_str);
        }
    }
    
    //*/ // /********************************************/
    if (arms_added)
    {
        std::vector<UINT32> dst;
        consensus = graph->generate_consensus_custom(dst);
        set_consensus(curate(consensus, dst));
        set_consensus_2(consensus);  
    }
    else {
        set_consensus(_draft.unpack());
        set_consensus_2(_draft.unpack());
    } 
}


std::string Window::curate(const std::string &con, const std::vector<UINT32> &dst) {
    //auto num_arms = _num_internal+_num_pre+_num_suf+_num_empty;
    auto conssize = con.size();
    std::string curated_consensus;
    curated_consensus.reserve(conssize);

    UINT cov_thres = (UINT)std::floor(_num_internal * _cThresh);

    for (UINT i = 0; i < conssize; ++i) {
        if (dst[i] >= cov_thres) {
            curated_consensus.push_back(con[i]);
        }
    }
    curated_consensus.shrink_to_fit();
    return curated_consensus;
}
    
} // namespace hypo
