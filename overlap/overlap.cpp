#include <zlib.h> // for reading compressed .fq file
#include <iostream>
#include <ctime>
#include <chrono>
#include <unordered_set>
#include <unordered_map>
#include <tuple>
#include <vector>
#include <algorithm>
#include <htslib/kseq.h>
#include <htslib/kstring.h>
#include <omp.h>
// #include <sdsl/bit_vectors.hpp>
// #include <sdsl/util.hpp>
#include "CustomBitvector.hpp"
using namespace std;
KSEQ_INIT(gzFile, gzread)

unsigned char seq_nt4_table[256] = {
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

tuple<int, int, int, int, int, int> find_overlap(const vector<unordered_map<uint64_t, vector<tuple<uint8_t, uint64_t> > > > & contig_solids, const vector<int64_t> & contig_lens, uint64_t i, uint64_t j, uint64_t i_from, uint64_t i_to, uint64_t j_from, uint64_t j_to) {
    vector<tuple<uint64_t, uint64_t> > solid_matches_forward;
    vector<tuple<uint64_t, uint64_t> > solid_matches_reverse;
    for(auto & content : contig_solids[i]) {
        uint64_t current_key = content.first;
        auto current_associated = contig_solids[j].find(current_key);
        if(current_associated != contig_solids[j].end()) {
            for(auto & pos : content.second) {
                if(i_from <= get<1>(pos) && get<1>(pos) <= i_to) { 
                    for(auto & pos2 : current_associated->second) {
                        if(j_from <= get<1>(pos2) && get<1>(pos2) <= j_to) {
                            if(get<0>(pos) == get<0>(pos2)) solid_matches_forward.push_back(make_tuple(get<1>(pos), get<1>(pos2)));
                            else solid_matches_reverse.push_back(make_tuple(get<1>(pos), get<1>(pos2)));
                        }
                    }
                }
            }
        }
    }
    
    // cout << "Found matches:" << solid_matches_forward.size() << " " << solid_matches_reverse.size() << endl;
    // find best chaining score
    sort(solid_matches_forward.begin(), solid_matches_forward.end());
    sort(solid_matches_reverse.begin(), solid_matches_reverse.end());
    
    int max_score_forward = -1;
    int max_score_reverse = -1;
    
    int forward_max_index, reverse_max_index;
    vector<int> score_forward, score_reverse;
    vector<int> pred_forward, pred_reverse;
    
    int best_contig_match = -1;
    int best_score = -1;
    int best_contig_s1 = -1;
    int best_contig_s2 = -1;
    int best_contig_e1 = -1;
    int best_contig_e2 = -1;
    int best_orientation = -1;
    
    if(solid_matches_forward.size() > 0) {
        int max_size = solid_matches_forward.size();
        score_forward.resize(max_size, -1);
        pred_forward.resize(max_size, -1);
        score_forward[0] = 1;
        pred_forward[0] = 0;
        
        int64_t distance1, distance2;
        
        //cout << "Find forward" << endl;
        
        forward_max_index = 0;
        for(int it = 1; it < max_size; it++) {
            // can't start here if too far from start / end
            if(get<0>(solid_matches_forward[it]) > 0.2 * contig_lens[i] || get<0>(solid_matches_forward[it]) < 0.8 * contig_lens[i]) {
                score_forward[it] = -1000;
                pred_forward[it] = it;
            }
            score_forward[it] = 1;
            pred_forward[it] = it;
            
            for(int it2 = it-1; it2 >= 0; it2--) {
                distance1 = get<0>(solid_matches_forward[it]) - get<0>(solid_matches_forward[it2]);
                distance2 = get<1>(solid_matches_forward[it]) - get<1>(solid_matches_forward[it2]);
                if(distance1 > 0 && distance2 > 0 && abs(distance1 - distance2) <= 0.1 * max(distance1, distance2) && abs(distance1 - distance2) <= 50000) { // a valid extension
                    if(score_forward[it2] + 1 > score_forward[it]) {
                        score_forward[it] = score_forward[it2] + 1;
                        pred_forward[it] = it2;
                    }
                }
            }
            
            if(score_forward[it] > score_forward[forward_max_index]) forward_max_index = it;
        }
        if(score_forward[forward_max_index] > best_score) {
            best_score = score_forward[forward_max_index];
            best_contig_match = j;
            best_contig_e1 = get<0>(solid_matches_forward[forward_max_index]);
            best_contig_e2 = get<1>(solid_matches_forward[forward_max_index]);
            
            int start_from = forward_max_index;
            while(pred_forward[start_from] != start_from) start_from = pred_forward[start_from];
            best_contig_s1 = get<0>(solid_matches_forward[start_from]);
            best_contig_s2 = get<1>(solid_matches_forward[start_from]);
            
            best_orientation = 1;
        }
        
        max_score_forward = score_forward[forward_max_index];
    }
    
    if(solid_matches_reverse.size() > 0) {
    
        //cout << "Find reverse" << endl;
        
        int max_size = solid_matches_reverse.size();
        score_reverse.resize(max_size, -1);
        pred_reverse.resize(max_size, -1);
        score_reverse[0] = 1;
        pred_reverse[0] = 0;
        
        int distance1, distance2;
        
        int reverse_max_index = 0;
        for(int it = 1; it < max_size; it++) {
            // can't start here if too far from start / end
            if(get<0>(solid_matches_reverse[it]) > 0.2 * contig_lens[i] || get<0>(solid_matches_reverse[it]) < 0.8 * contig_lens[i]) {
                score_reverse[it] = -1000;
                pred_reverse[it] = it;
            }
            else {
                score_reverse[it] = 1;
                pred_reverse[it] = it;
            }
            
            for(int it2 = it-1; it2 >= 0; it2--) {
                distance1 = get<0>(solid_matches_reverse[it]) - get<0>(solid_matches_reverse[it2]);
                distance2 = get<1>(solid_matches_reverse[it]) - get<1>(solid_matches_reverse[it2]);
                if(((distance1 > 0 && distance2 < 0) || (distance1 < 0 && distance2 > 0)) && (abs(abs(distance1) - abs(distance2)) <= 0.1 * max(abs(distance1), abs(distance2))) && (abs(abs(distance1) - abs(distance2)) <= 50000)) { // a valid extension
                    if(score_reverse[it2] + 1 > score_reverse[it]) {
                        score_reverse[it] = score_reverse[it2] + 1;
                        pred_reverse[it] = it2;
                    }
                }
            }
            
            if(score_reverse[it] > score_reverse[reverse_max_index]) reverse_max_index = it;
        }
        
        if(score_reverse[reverse_max_index] > best_score) {
            best_score = score_reverse[reverse_max_index];
            best_contig_match = j;
            best_contig_e1 = get<0>(solid_matches_reverse[reverse_max_index]);
            best_contig_e2 = get<1>(solid_matches_reverse[reverse_max_index]);
            
            int start_from = reverse_max_index;
            while(pred_reverse[start_from] != start_from) start_from = pred_reverse[start_from];
            best_contig_s1 = get<0>(solid_matches_reverse[start_from]);
            best_contig_s2 = get<1>(solid_matches_reverse[start_from]);
            
            best_orientation = -1;
        }
        
        max_score_reverse = score_reverse[reverse_max_index];
    }
    
    return make_tuple(best_contig_s1, best_contig_e1, best_contig_s2, best_contig_e2, best_score, best_orientation);
}

int main(int argc, char* argv[]) {
    if(argc < 4) {
        cerr << "Usage: " << argv[0] << " <k> <kmers bitvector> <contigs fasta> <optional: thread count (default=1)> <optional: remove >2 occurence solid kmers (1/0) (default=0)>" << endl;
        return 1;
    }
    cerr << "Command: ";
    for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
    cerr << endl;
    
    
    int thread_count = 1;
    try {
        thread_count = stoi(argv[4]);
    } catch(exception const &e) {
        cerr << "Thread count need to be int. Defaults to 1" << endl;
    }
    
    int remove_dupe = 0;
    try {
        remove_dupe = stoi(argv[5]);
    } catch(exception const &e) {
        cerr << "Remove > 2 solids defaults to 0 (no)." << endl;
    }
    
    cerr << "Using " << thread_count << " threads." << endl;
    
    int k = 0;
    try {
        k = stoi(argv[1]);
    } catch(exception const &e) {
        cerr << "k must be int." << endl;
        return 1;
    }
    cerr << "k = " << k << endl;
    
    unordered_set<uint64_t> solid_kmers;
    
    cerr << "Reading kmers from " << argv[2] << "." << endl;
    //sdsl::bit_vector bv_sd;
    //sdsl::load_from_file(bv_sd, argv[2]);
    //auto get_ones = sdsl::bit_vector::rank_1_type(&bv_sd)(bv_sd.size());
    CustomBitvector bv_sd(1ULL<<(2*k));
    bv_sd.load(argv[2]);
    auto get_ones = bv_sd.count();
    
    cerr << "Got " << get_ones << " solid kmers." << endl;
    
    uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
    
    vector<string> contig_names;
    vector<int64_t> contig_lens;
    vector<unordered_map<uint64_t, vector<tuple<uint8_t, uint64_t> > > > contig_solids;
    
    unordered_map<uint64_t, uint64_t> kmer_counter;
    
    cerr << "Processing contigs from " << argv[3] << "." << endl;
    gzFile fp = gzopen(argv[3], "r");
    kseq_t *seq;
    seq = kseq_init(fp);
    int l;
    while((l = kseq_read(seq)) >= 0) {
        string read_name = seq->name.s;
        cerr << "Processing contig " << read_name << endl;
        int count_not_N = 0;
        
        contig_names.push_back(read_name);
        contig_lens.push_back(seq->seq.l);
        unordered_map<uint64_t, vector<tuple<uint8_t, uint64_t> > > current_contig_solids;
        for(size_t i = 0; i < seq->seq.l; i++) {
            int c = seq_nt4_table[(uint8_t)seq->seq.s[i]];
            if(c < 4) {
                kmer[0] = (kmer[0] << 2ull | c) & mask;           // forward k-mer
                kmer[1] = (kmer[1] >> 2ull) | (3ULL^c) << shift1; // reverse k-mer
                count_not_N++;
                int z = kmer[0] < kmer[1] ? 0 : 1;
                if(count_not_N >= k) {
                    if(bv_sd[kmer[z]]) {
                        //found a solid kmer, do something here
                        current_contig_solids[kmer[z]].push_back(make_tuple(z, i));
                        kmer_counter[kmer[z]]++;
                    }
                }
            } else {
                count_not_N = 0;
            }
        }
        
        contig_solids.push_back(current_contig_solids);
        
        cerr << "Found " << current_contig_solids.size() << " solids." << endl;
    }
    
    if(remove_dupe) {
        unordered_set<uint64_t> multi_solid;
        for(auto & content : kmer_counter) {
            if(content.second > 2) {
                multi_solid.insert(content.first);
            }
        }
        
        cerr << "Remove " << multi_solid.size() << " kmers." << endl;
        
        for(int i = 0; i < contig_names.size(); i++) {
            for(auto & remove : multi_solid) {
                contig_solids[i].erase(remove);
            }
        }
    }
    
    cerr << "Finding overlaps" << endl;
    
    
    #pragma omp parallel for num_threads(thread_count)
    for(int i = 0; i < contig_names.size(); i++) {
        int best_contig_match = -1;
        int best_score = -1;
        int best_contig_s1 = -1;
        int best_contig_s2 = -1;
        int best_contig_e1 = -1;
        int best_contig_e2 = -1;
        int best_orientation = -1;
        
        /*
        for(int j = 0; j < i - 1; j++) {
            if(get<0>(previous_scores[i][j]) > best_score) {
                best_contig_match = j;
                best_score = get<0>(previous_scores[i][j]);
                best_contig_s1 = get<3>(previous_scores[i][j]);
                best_contig_e1 = get<4>(previous_scores[i][j]);
                best_contig_s2 = get<1>(previous_scores[i][j]);
                best_contig_e2 = get<2>(previous_scores[i][j]);
            }
        }
        */
        
        for(int j = 0; j < contig_names.size(); j++) {
            if(i == j) continue;
            // find solid kmer matches by orientation + position
            vector<tuple<uint64_t, uint64_t> > solid_matches_forward;
            vector<tuple<uint64_t, uint64_t> > solid_matches_reverse;
            for(auto & content : contig_solids[i]) {
                uint64_t current_key = content.first;
                auto current_associated = contig_solids[j].find(current_key);
                if(current_associated != contig_solids[j].end()) {
                    for(auto & pos : content.second) {
                        for(auto & pos2 : current_associated->second) {
                            if(get<0>(pos) == get<0>(pos2)) solid_matches_forward.push_back(make_tuple(get<1>(pos), get<1>(pos2)));
                            else solid_matches_reverse.push_back(make_tuple(get<1>(pos), get<1>(pos2)));
                        }
                    }
                }
            }
            
            // cout << "Found matches:" << solid_matches_forward.size() << " " << solid_matches_reverse.size() << endl;
            // find best chaining score
            sort(solid_matches_forward.begin(), solid_matches_forward.end());
            sort(solid_matches_reverse.begin(), solid_matches_reverse.end());
            
            int max_score_forward = -1;
            int max_score_reverse = -1;
            
            int forward_max_index, reverse_max_index;
            vector<int> score_forward, score_reverse;
            vector<int> pred_forward, pred_reverse;
            
            if(solid_matches_forward.size() > 0) {
                int max_size = solid_matches_forward.size();
                score_forward.resize(max_size, -1);
                pred_forward.resize(max_size, -1);
                score_forward[0] = 1;
                pred_forward[0] = 0;
                
                int64_t distance1, distance2;
                
                //cout << "Find forward" << endl;
                
                forward_max_index = 0;
                for(int it = 1; it < max_size; it++) {
                    // can't start here if too far from start / end
                    if(get<0>(solid_matches_forward[it]) > 0.2 * contig_lens[i] || get<0>(solid_matches_forward[it]) < 0.8 * contig_lens[i]) {
                        score_forward[it] = -1000;
                        pred_forward[it] = it;
                    }
                    score_forward[it] = 1;
                    pred_forward[it] = it;
                    
                    for(int it2 = it-1; it2 >= 0; it2--) {
                        distance1 = get<0>(solid_matches_forward[it]) - get<0>(solid_matches_forward[it2]);
                        distance2 = get<1>(solid_matches_forward[it]) - get<1>(solid_matches_forward[it2]);
                        if(distance1 > 0 && distance2 > 0 && abs(distance1 - distance2) <= 0.1 * max(distance1, distance2) && abs(distance1 - distance2) <= 50000) { // a valid extension
                            if(score_forward[it2] + 1 > score_forward[it]) {
                                score_forward[it] = score_forward[it2] + 1;
                                pred_forward[it] = it2;
                            }
                        }
                    }
                    
                    if(score_forward[it] > score_forward[forward_max_index]) forward_max_index = it;
                }
                if(score_forward[forward_max_index] > best_score) {
                    best_score = score_forward[forward_max_index];
                    best_contig_match = j;
                    best_contig_e1 = get<0>(solid_matches_forward[forward_max_index]);
                    best_contig_e2 = get<1>(solid_matches_forward[forward_max_index]);
                    
                    int start_from = forward_max_index;
                    while(pred_forward[start_from] != start_from) start_from = pred_forward[start_from];
                    best_contig_s1 = get<0>(solid_matches_forward[start_from]);
                    best_contig_s2 = get<1>(solid_matches_forward[start_from]);
                    
                    best_orientation = 1;
                }
                
                max_score_forward = score_forward[forward_max_index];
            }
            
            if(solid_matches_reverse.size() > 0) {
            
                //cout << "Find reverse" << endl;
                
                int max_size = solid_matches_reverse.size();
                score_reverse.resize(max_size, -1);
                pred_reverse.resize(max_size, -1);
                score_reverse[0] = 1;
                pred_reverse[0] = 0;
                
                int distance1, distance2;
                
                int reverse_max_index = 0;
                for(int it = 1; it < max_size; it++) {
                    // can't start here if too far from start / end
                    if(get<0>(solid_matches_reverse[it]) > 0.2 * contig_lens[i] || get<0>(solid_matches_reverse[it]) < 0.8 * contig_lens[i]) {
                        score_reverse[it] = -1000;
                        pred_reverse[it] = it;
                    }
                    else {
                        score_reverse[it] = 1;
                        pred_reverse[it] = it;
                    }
                    
                    for(int it2 = it-1; it2 >= 0; it2--) {
                        distance1 = get<0>(solid_matches_reverse[it]) - get<0>(solid_matches_reverse[it2]);
                        distance2 = get<1>(solid_matches_reverse[it]) - get<1>(solid_matches_reverse[it2]);
                        if(((distance1 > 0 && distance2 < 0) || (distance1 < 0 && distance2 > 0)) && abs(abs(distance1) - abs(distance2)) <= 0.1 * max(abs(distance1), abs(distance2)) && abs(abs(distance1) - abs(distance2)) <= 50000) { // a valid extension
                            if(score_reverse[it2] + 1 > score_reverse[it]) {
                                score_reverse[it] = score_reverse[it2] + 1;
                                pred_reverse[it] = it2;
                            }
                        }
                    }
                    
                    if(score_reverse[it] > score_reverse[reverse_max_index]) reverse_max_index = it;
                }
                
                if(score_reverse[reverse_max_index] > best_score) {
                    best_score = score_reverse[reverse_max_index];
                    best_contig_match = j;
                    best_contig_e1 = get<0>(solid_matches_reverse[reverse_max_index]);
                    best_contig_e2 = get<1>(solid_matches_reverse[reverse_max_index]);
                    
                    int start_from = reverse_max_index;
                    while(pred_reverse[start_from] != start_from) start_from = pred_reverse[start_from];
                    best_contig_s1 = get<0>(solid_matches_reverse[start_from]);
                    best_contig_s2 = get<1>(solid_matches_reverse[start_from]);
                    
                    best_orientation = -1;
                }
                
                max_score_reverse = score_reverse[reverse_max_index];
            }
            
            /*
            if(max_score_forward > max_score_reverse && max_score_forward > -1) {
                int start_from = forward_max_index;
                while(pred_forward[start_from] != start_from) start_from = pred_forward[start_from];
                previous_scores[j].push_back(make_tuple(score_forward[forward_max_index], get<0>(solid_matches_forward[start_from]), get<0>(solid_matches_forward[forward_max_index]), get<1>(solid_matches_forward[start_from]), get<1>(solid_matches_forward[forward_max_index])));                
            } else if(max_score_reverse > -1) {
                int start_from = reverse_max_index;
                while(pred_reverse[start_from] != start_from) start_from = pred_reverse[start_from];
                previous_scores[j].push_back(make_tuple(score_reverse[reverse_max_index], get<0>(solid_matches_reverse[start_from]), get<0>(solid_matches_reverse[reverse_max_index]), get<1>(solid_matches_reverse[start_from]), get<1>(solid_matches_reverse[reverse_max_index])));
            } else {
                previous_scores[j].push_back(make_tuple(-1, -1, -1, -1, -1));
            }
            */
            //cout << "With contig " << j << " fwd: " << max_score_forward << " rev: " << max_score_reverse << endl;
        }
        
        // found best contig, now find the overlap parts
        
        vector<tuple<int, int, int, int, int, int> > results;
        if(best_score > -1) {
            // recursively find overlapping parts
            
            // first overlapping part is from s_1 e_1 to s_2 e_2
            
            
            results.push_back(make_tuple(best_contig_s1, best_contig_e1, best_contig_s2, best_contig_e2, best_score, best_orientation));
            
            if(best_orientation == 1) {
                tuple<int, int, int, int, int, int> temp = find_overlap(contig_solids, contig_lens, i, best_contig_match, 0, best_contig_s1, 0, best_contig_s2);
                if(get<4>(temp) > best_score / 4) results.push_back(temp);
                temp = find_overlap(contig_solids, contig_lens, i, best_contig_match, best_contig_e1, contig_lens[i], best_contig_e2, contig_lens[best_contig_match]);
                if(get<4>(temp) > best_score / 4) results.push_back(temp);
            }
            
            if(best_orientation == -1) {
                tuple<int, int, int, int, int, int> temp = find_overlap(contig_solids, contig_lens, i, best_contig_match, 0, best_contig_s1, best_contig_e2, contig_lens[best_contig_match]);
                if(get<4>(temp) > best_score / 4) results.push_back(temp);
                temp = find_overlap(contig_solids, contig_lens, i, best_contig_match, best_contig_e1, contig_lens[i], 0, best_contig_s2);
                if(get<4>(temp) > best_score / 4) results.push_back(temp);
            }
        }
        
        
        #pragma omp critical 
        {
            if(best_score > -1) {
                tuple<int, int, int, int, int, int> result;
                for(int zz = 0; zz < results.size(); zz++) {
                    result = results[zz];
                    cout << contig_names[i] << "\t" << get<0>(result) << "\t" << get<1>(result) << "\t" << contig_names[best_contig_match] << "\t" << get<2>(result) << "\t" << get<3>(result) << "\t" << get<4>(result) << "\t" << get<5>(result) << "\t";
                    if(zz == 0) cout << "P" << endl;
                    else cout << "EXT" << endl;
                } 
            }
            else cout << contig_names[i] << "\t-1\t-1\tNA\t-1\t-1\t-1\t-1" << endl;
        }
    }
}
