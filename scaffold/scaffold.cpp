#include <zlib.h> // for reading compressed .fq file
#include <iostream>
#include <ctime>
#include <chrono>
#include <unordered_set>
#include <unordered_map>
#include <map>
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

int main(int argc, char* argv[]) {
    if(argc < 8) {
        cerr << "Usage: " << argv[0] << " <k> <kmers input> <contigs 1 fasta> <contigs 2 fasta> <long reads> <output 1> <output 2> <optional: thread count (default=1)> <optional: remove >2 occurence solid kmers (1/0) (default=0)> <optional: debug_output>" << endl;
        return 1;
    }
    cerr << "Command: ";
    for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
    cerr << endl;
    
    bool is_debug = false;
    ofstream out_debug;
    
    if(argc > 10) {
        is_debug = true;
    }
    
    if(is_debug) {
        try {
            out_debug.open(argv[10]);
        } catch(exception const &e) {
            cerr << "Can't open file" << endl;
            cerr << e.what() << endl;
            cerr << strerror(errno) << endl;
            return -1;
        }
    }
    
    int thread_count = 1;
    try {
        thread_count = stoi(argv[8]);
    } catch(exception const &e) {
        cerr << "Thread count need to be int. Defaults to 1" << endl;
    }
    
    int remove_dupe = 0;
    try {
        remove_dupe = stoi(argv[9]);
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
    CustomBitvector bv_sd(1ULL<<(2*k));
    bv_sd.load(argv[2]);
    auto get_ones = bv_sd.count();
    
    cerr << "Got " << get_ones << " solid kmers." << endl;
    cerr << "Bitvector size: " << bv_sd.size() << endl;
    
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
    
    int total_kmers = 0;
    while((l = kseq_read(seq)) >= 0) {
        string read_name = seq->name.s;
        cerr << "Processing contig " << read_name << endl;
        int count_not_N = 0;
        int count_kmers = 0;
        
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
                    count_kmers++;
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
        
        cerr << "Found " << current_contig_solids.size() << " solids out of " << count_kmers << endl;
        total_kmers += current_contig_solids.size();
    }
    
    cerr << "Total: " << total_kmers << " " << kmer_counter.size() << endl;
    
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
    
    cerr << "Processing reads from " << argv[5] << "." << endl;
    fp = gzopen(argv[5], "r");
    seq = kseq_init(fp);
    
    map<tuple<uint64_t, uint64_t>, uint64_t> matching_solids;
    uint64_t minkmer, maxkmer, tempkmer1, tempkmer2;
    map<tuple<uint64_t, uint64_t>, vector<string> > sequences_between;
    map<tuple<uint64_t, uint64_t>, uint32_t > counter_between;
    
    uint64_t total_reads = 0;
    uint64_t read_total_kmers = 0;
    while((l = kseq_read(seq)) >= 0) {
        string read_name = seq->name.s;
        int count_not_N = 0;
        int count_kmer = 0;
        int previous_solid = -1;
        uint64_t previous_kmer;
        vector<uint64_t> current_solids;
        for(size_t i = 0; i < seq->seq.l; i++) {
            int c = seq_nt4_table[(uint8_t)seq->seq.s[i]];
            if(c < 4) {
                kmer[0] = (kmer[0] << 2ull | c) & mask;           // forward k-mer
                kmer[1] = (kmer[1] >> 2ull) | (3ULL^c) << shift1; // reverse k-mer
                count_not_N++;
                int z = kmer[0] < kmer[1] ? 0 : 1;
                if(count_not_N >= k) {
                    count_kmer += 1;
                    if(bv_sd[kmer[z]]) {
                        read_total_kmers++;
                        current_solids.push_back(kmer[z]);
                        
                        if(previous_solid != -1) {
                            std::string get_sequence = string(seq->seq.s, seq->seq.s+i-previous_solid);
                            if(is_debug) sequences_between[make_tuple(previous_kmer, kmer[z])].push_back(get_sequence);
                            counter_between[make_tuple(previous_kmer, kmer[z])]++;
                        }
                        
                        previous_kmer = kmer[z];
                        previous_solid = i;
                    }
                }
            } else {
                count_not_N = 0;
            }
        }
        for(size_t i = 0; i < current_solids.size(); i++) {
            for(size_t j = i + 1; j < current_solids.size(); j++) {
                minkmer = current_solids[i] < current_solids[j] ? current_solids[i] : current_solids[j];
                maxkmer = current_solids[i] < current_solids[j] ? current_solids[j] : current_solids[i];
                matching_solids[make_tuple(minkmer, maxkmer)]++;
            }
        }
        total_reads++;
    }
    
    cerr << "Found " << read_total_kmers << " out of " << total_reads << endl;
    
    
    cerr << "Finding scaffolds" << endl;
    
    ofstream output1(argv[6]);
    
    // initialize ordering based on the solid kmers count
    vector<size_t> ordering;
    for(int i = 0; i < contig_names.size(); i++) ordering.push_back(i);
    sort(ordering.begin(), ordering.end(), [&](size_t i1, size_t i2) {return contig_solids[i1].size() > contig_solids[i2].size();});
    
    #pragma omp parallel for num_threads(thread_count)
    for(int i_idx = 0; i_idx < contig_names.size(); i_idx++) {
        int i = ordering[i_idx];
        
        int best_contig_match = -1;
        int best_score = -1;
        int best_contig_s1 = -1;
        int best_contig_s2 = -1;
        int best_contig_e1 = -1;
        int best_contig_e2 = -1;
        int best_orientation = -1;
        
        // vector of all matches instead of keeping one
        vector<vector<tuple<uint64_t, uint64_t, uint64_t> > > match_fwd;
        vector<vector<tuple<uint64_t, uint64_t, uint64_t> > > match_rev;
        
        
        for(int j_idx = 0; j_idx < contig_names.size(); j_idx++) {
            int j = ordering[j_idx];
            
            vector<tuple<uint64_t, uint64_t, uint64_t> > solid_matches_forward;
            vector<tuple<uint64_t, uint64_t, uint64_t> > solid_matches_reverse;
            
            
            if(i == j) {
                // just put empty as it will be discarded
                match_fwd.push_back(solid_matches_forward);
                match_rev.push_back(solid_matches_reverse);
                continue;
            }
            // find solid kmer matches by orientation + position
            
            for(auto & content : contig_solids[i]) {
                uint64_t current_key = content.first;
                auto current_associated = contig_solids[j].find(current_key);
                if(current_associated != contig_solids[j].end()) {
                    for(auto & pos : content.second) {
                        for(auto & pos2 : current_associated->second) {
                            if(get<0>(pos) == get<0>(pos2)) solid_matches_forward.push_back(make_tuple(get<1>(pos), get<1>(pos2), current_key));
                            else solid_matches_reverse.push_back(make_tuple(get<1>(pos), get<1>(pos2), current_key));
                        }
                    }
                }
            }
            
            // cout << "Found matches:" << solid_matches_forward.size() << " " << solid_matches_reverse.size() << endl;
            // find best chaining score
            sort(solid_matches_forward.begin(), solid_matches_forward.end());
            sort(solid_matches_reverse.begin(), solid_matches_reverse.end());
            
            match_fwd.push_back(solid_matches_forward);
            match_rev.push_back(solid_matches_reverse);
        }
        
        // after all the matches are found, process them in order by most matches
        vector<size_t> new_ordering;
        for(int i = 0; i < contig_names.size(); i++) new_ordering.push_back(i);
        sort(new_ordering.begin(), new_ordering.end(), [&](size_t i1, size_t i2) {return (match_fwd[i1].size() + match_rev[i1].size()) > (match_fwd[i2].size() + match_rev[i2].size());});
        
        vector<tuple<int, int, int, int, int, int> > all_matches;
        
        for(int j_idx = 0; j_idx < contig_names.size(); j_idx++) {
            int j = new_ordering[j_idx];
            
            vector<tuple<uint64_t, uint64_t, uint64_t> > solid_matches_forward = match_fwd[j];
            vector<tuple<uint64_t, uint64_t, uint64_t> > solid_matches_reverse = match_rev[j];
            
            #pragma omp critical 
            {
                if(is_debug) {
                    out_debug << contig_names[i] << "\t" << contig_names[j] << "\n";
                    
                    for(int it = 0; it < solid_matches_forward.size(); it++) out_debug << "F\t" << get<0>(solid_matches_forward[it]) << "\t" << get<1>(solid_matches_forward[it]) << "\n";
                    for(int it = 0; it < solid_matches_reverse.size(); it++) out_debug << "R\t" << get<0>(solid_matches_reverse[it]) << "\t" << get<1>(solid_matches_reverse[it]) << "\n";
                }
            }
            
            int max_score_forward = -1;
            int max_score_reverse = -1;
            
            int forward_max_index, reverse_max_index;
            vector<int> score_forward, score_reverse;
            vector<int> pred_forward, pred_reverse;
            
            // find overlap parts from all matches
            
            if(solid_matches_forward.size() > 0) {
                int max_size = solid_matches_forward.size();
                score_forward.resize(max_size, -1);
                pred_forward.resize(max_size, -1);
                score_forward[0] = 1;
                pred_forward[0] = 0;
                
                int64_t distance1, distance2;
                
                int current_best_contig_match = -1;
                int current_best_score = -1;
                int current_best_contig_s1 = -1;
                int current_best_contig_s2 = -1;
                int current_best_contig_e1 = -1;
                int current_best_contig_e2 = -1;
                int current_best_orientation = -1;
                
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
                
                if(score_forward[forward_max_index] > current_best_score) {
                    current_best_score = score_forward[forward_max_index];
                    current_best_contig_match = j;
                    current_best_contig_e1 = get<0>(solid_matches_forward[forward_max_index]);
                    current_best_contig_e2 = get<1>(solid_matches_forward[forward_max_index]);
                    
                    int start_from = forward_max_index;
                    while(pred_forward[start_from] != start_from) start_from = pred_forward[start_from];
                    current_best_contig_s1 = get<0>(solid_matches_forward[start_from]);
                    current_best_contig_s2 = get<1>(solid_matches_forward[start_from]);
                    
                    current_best_orientation = 1;
                }
                
                max_score_forward = score_forward[forward_max_index];
                
                all_matches.push_back(make_tuple(current_best_contig_s1, current_best_contig_e1, current_best_contig_s2, current_best_contig_e2, current_best_score, current_best_orientation));
                
                if(current_best_score > best_score) {
                    best_score = current_best_score;
                    best_contig_match = current_best_contig_match;
                    best_contig_e1 = current_best_contig_e1;
                    best_contig_e2 = current_best_contig_e2;
                    best_contig_s1 = current_best_contig_s1;
                    best_contig_s2 = current_best_contig_s2;
                    best_orientation = current_best_orientation;
                }
            }
            
            if(solid_matches_reverse.size() > 0) {
            
                //cout << "Find reverse" << endl;
                
                int max_size = solid_matches_reverse.size();
                score_reverse.resize(max_size, -1);
                pred_reverse.resize(max_size, -1);
                score_reverse[0] = 1;
                pred_reverse[0] = 0;
                
                int distance1, distance2;
                
                int current_best_contig_match = -1;
                int current_best_score = -1;
                int current_best_contig_s1 = -1;
                int current_best_contig_s2 = -1;
                int current_best_contig_e1 = -1;
                int current_best_contig_e2 = -1;
                int current_best_orientation = -1;
                
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
                                
                                // adjust score by read solid kmer support
                                
                                tempkmer1 = get<2>(solid_matches_reverse[it]);
                                tempkmer2 = get<2>(solid_matches_reverse[it2]);
                                
                                minkmer = tempkmer1 < tempkmer2 ? tempkmer1 : tempkmer2;
                                maxkmer = tempkmer1 < tempkmer2 ? tempkmer2 : tempkmer1;
                                
                                // this match is not supported by reads!
                                if(matching_solids.find(make_tuple(minkmer, maxkmer)) == matching_solids.end()) {
                                    score_reverse[it] = 1;
                                } 
                                else {
                                    score_reverse[it] += matching_solids[make_tuple(minkmer, maxkmer)];
                                }
                            }
                        }
                    }
                    
                    if(score_reverse[it] > score_reverse[reverse_max_index]) reverse_max_index = it;
                }
                
                if(score_reverse[reverse_max_index] > current_best_score) {
                    current_best_score = score_reverse[reverse_max_index];
                    current_best_contig_match = j;
                    current_best_contig_e1 = get<0>(solid_matches_reverse[reverse_max_index]);
                    current_best_contig_e2 = get<1>(solid_matches_reverse[reverse_max_index]);
                    
                    int start_from = reverse_max_index;
                    while(pred_reverse[start_from] != start_from) start_from = pred_reverse[start_from];
                    current_best_contig_s1 = get<0>(solid_matches_reverse[start_from]);
                    current_best_contig_s2 = get<1>(solid_matches_reverse[start_from]);
                    
                    current_best_orientation = -1;
                }
                
                max_score_reverse = score_reverse[reverse_max_index];
                
                
                all_matches.push_back(make_tuple(current_best_contig_s1, current_best_contig_e1, current_best_contig_s2, current_best_contig_e2, current_best_score, current_best_orientation));
                
                if(current_best_score > best_score) {
                    best_score = current_best_score;
                    best_contig_match = current_best_contig_match;
                    best_contig_e1 = current_best_contig_e1;
                    best_contig_e2 = current_best_contig_e2;
                    best_contig_s1 = current_best_contig_s1;
                    best_contig_s2 = current_best_contig_s2;
                    best_orientation = current_best_orientation;
                }
            }
            
        }
        
        vector<tuple<int, int, int, int, int, int> > results;
        if(best_score > -1) {
            
            #pragma omp critical
            {
                out_debug << best_contig_s1 << " " << best_contig_e1 << " " << best_contig_s2 << " " << " " << best_contig_e2 << " "  << best_score << " " << best_orientation << endl;
            }
            
            results.push_back(make_tuple(best_contig_s1, best_contig_e1, best_contig_s2, best_contig_e2, best_score, best_orientation));
            
            for(int am = 0; am < all_matches.size(); am++) if(get<4>(all_matches[am]) > best_score / 10) results.push_back(all_matches[am]);
        }
        
        
        #pragma omp critical 
        {
            if(best_score > -1) {
                tuple<int, int, int, int, int, int> result;
                for(int zz = 0; zz < results.size(); zz++) {
                    result = results[zz];
                    output1 << contig_names[i] << "\t" << get<0>(result) << "\t" << get<1>(result) << "\t" << contig_names[best_contig_match] << "\t" << get<2>(result) << "\t" << get<3>(result) << "\t" << get<4>(result) << "\t" << get<5>(result) << "\t";
                    if(zz == 0) cout << "P" << endl;
                    else cout << "EXT" << endl;
                } 
            }
            else output1 << contig_names[i] << "\t-1\t-1\tNA\t-1\t-1\t-1\t-1" << endl;
        }
    }
    output1.close();
    
    cerr << "Processing contigs from " << argv[4] << "." << endl;
    fp = gzopen(argv[4], "r");
    seq = kseq_init(fp);
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
    
    cerr << "Finding scaffolds" << endl;
    
    ofstream output2(argv[7]);
    
    // initialize ordering based on the solid kmers count
    ordering.clear();
    for(int i = 0; i < contig_names.size(); i++) ordering.push_back(i);
    sort(ordering.begin(), ordering.end(), [&](size_t i1, size_t i2) {return contig_solids[i1].size() > contig_solids[i2].size();});
    
    #pragma omp parallel for num_threads(thread_count)
    for(int i_idx = 0; i_idx < contig_names.size(); i_idx++) {
        int i = ordering[i_idx];
        
        int best_contig_match = -1;
        int best_score = -1;
        int best_contig_s1 = -1;
        int best_contig_s2 = -1;
        int best_contig_e1 = -1;
        int best_contig_e2 = -1;
        int best_orientation = -1;
        
        // vector of all matches instead of keeping one
        vector<vector<tuple<uint64_t, uint64_t, uint64_t> > > match_fwd;
        vector<vector<tuple<uint64_t, uint64_t, uint64_t> > > match_rev;
        
        
        for(int j_idx = 0; j_idx < contig_names.size(); j_idx++) {
            int j = ordering[j_idx];
            
            vector<tuple<uint64_t, uint64_t, uint64_t> > solid_matches_forward;
            vector<tuple<uint64_t, uint64_t, uint64_t> > solid_matches_reverse;
            
            
            if(i == j) {
                // just put empty as it will be discarded
                match_fwd.push_back(solid_matches_forward);
                match_rev.push_back(solid_matches_reverse);
                continue;
            }
            // find solid kmer matches by orientation + position
            
            for(auto & content : contig_solids[i]) {
                uint64_t current_key = content.first;
                auto current_associated = contig_solids[j].find(current_key);
                if(current_associated != contig_solids[j].end()) {
                    for(auto & pos : content.second) {
                        for(auto & pos2 : current_associated->second) {
                            if(get<0>(pos) == get<0>(pos2)) solid_matches_forward.push_back(make_tuple(get<1>(pos), get<1>(pos2), current_key));
                            else solid_matches_reverse.push_back(make_tuple(get<1>(pos), get<1>(pos2), current_key));
                        }
                    }
                }
            }
            
            // cout << "Found matches:" << solid_matches_forward.size() << " " << solid_matches_reverse.size() << endl;
            // find best chaining score
            sort(solid_matches_forward.begin(), solid_matches_forward.end());
            sort(solid_matches_reverse.begin(), solid_matches_reverse.end());
            
            
            match_fwd.push_back(solid_matches_forward);
            match_rev.push_back(solid_matches_reverse);
        }
        
        // after all the matches are found, process them in order by most matches
        vector<size_t> new_ordering;
        for(int i = 0; i < contig_names.size(); i++) new_ordering.push_back(i);
        sort(new_ordering.begin(), new_ordering.end(), [&](size_t i1, size_t i2) {return (match_fwd[i1].size() + match_rev[i1].size()) > (match_fwd[i2].size() + match_rev[i2].size());});
        
        vector<tuple<int, int, int, int, int, int> > all_matches;
        
        for(int j_idx = 0; j_idx < contig_names.size(); j_idx++) {
            int j = new_ordering[j_idx];
            
            vector<tuple<uint64_t, uint64_t, uint64_t> > solid_matches_forward = match_fwd[j];
            vector<tuple<uint64_t, uint64_t, uint64_t> > solid_matches_reverse = match_rev[j];
            
            #pragma omp critical 
            {
                if(is_debug) {
                    out_debug << contig_names[i] << "\t" << contig_names[j] << "\n";
                    
                    for(int it = 0; it < solid_matches_forward.size(); it++) out_debug << "F\t" << get<0>(solid_matches_forward[it]) << "\t" << get<1>(solid_matches_forward[it]) << "\n";
                    for(int it = 0; it < solid_matches_reverse.size(); it++) out_debug << "R\t" << get<0>(solid_matches_reverse[it]) << "\t" << get<1>(solid_matches_reverse[it]) << "\n";
                }
            }
            
            int max_score_forward = -1;
            int max_score_reverse = -1;
            
            int forward_max_index, reverse_max_index;
            vector<int> score_forward, score_reverse;
            vector<int> pred_forward, pred_reverse;
            
            // find overlap parts from all matches
            
            if(solid_matches_forward.size() > 0) {
                int max_size = solid_matches_forward.size();
                score_forward.resize(max_size, -1);
                pred_forward.resize(max_size, -1);
                score_forward[0] = 1;
                pred_forward[0] = 0;
                
                int64_t distance1, distance2;
                
                int current_best_contig_match = -1;
                int current_best_score = -1;
                int current_best_contig_s1 = -1;
                int current_best_contig_s2 = -1;
                int current_best_contig_e1 = -1;
                int current_best_contig_e2 = -1;
                int current_best_orientation = -1;
                
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
                
                if(score_forward[forward_max_index] > current_best_score) {
                    current_best_score = score_forward[forward_max_index];
                    current_best_contig_match = j;
                    current_best_contig_e1 = get<0>(solid_matches_forward[forward_max_index]);
                    current_best_contig_e2 = get<1>(solid_matches_forward[forward_max_index]);
                    
                    int start_from = forward_max_index;
                    while(pred_forward[start_from] != start_from) start_from = pred_forward[start_from];
                    current_best_contig_s1 = get<0>(solid_matches_forward[start_from]);
                    current_best_contig_s2 = get<1>(solid_matches_forward[start_from]);
                    
                    current_best_orientation = 1;
                }
                
                max_score_forward = score_forward[forward_max_index];
                
                all_matches.push_back(make_tuple(current_best_contig_s1, current_best_contig_e1, current_best_contig_s2, current_best_contig_e2, current_best_score, current_best_orientation));
                
                if(current_best_score > best_score) {
                    best_score = current_best_score;
                    best_contig_match = current_best_contig_match;
                    best_contig_e1 = current_best_contig_e1;
                    best_contig_e2 = current_best_contig_e2;
                    best_contig_s1 = current_best_contig_s1;
                    best_contig_s2 = current_best_contig_s2;
                    best_orientation = current_best_orientation;
                }
            }
            
            if(solid_matches_reverse.size() > 0) {
            
                //cout << "Find reverse" << endl;
                
                int max_size = solid_matches_reverse.size();
                score_reverse.resize(max_size, -1);
                pred_reverse.resize(max_size, -1);
                score_reverse[0] = 1;
                pred_reverse[0] = 0;
                
                int distance1, distance2;
                
                int current_best_contig_match = -1;
                int current_best_score = -1;
                int current_best_contig_s1 = -1;
                int current_best_contig_s2 = -1;
                int current_best_contig_e1 = -1;
                int current_best_contig_e2 = -1;
                int current_best_orientation = -1;
                
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
                                
                                // adjust score by read solid kmer support
                                
                                tempkmer1 = get<2>(solid_matches_reverse[it]);
                                tempkmer2 = get<2>(solid_matches_reverse[it2]);
                                
                                minkmer = tempkmer1 < tempkmer2 ? tempkmer1 : tempkmer2;
                                maxkmer = tempkmer1 < tempkmer2 ? tempkmer2 : tempkmer1;
                                
                                if(matching_solids.find(make_tuple(minkmer, maxkmer)) == matching_solids.end()) {
                                    score_reverse[it] = 1;
                                } 
                                else {
                                    score_reverse[it] += matching_solids[make_tuple(minkmer, maxkmer)];
                                }
                            }
                        }
                    }
                    
                    if(score_reverse[it] > score_reverse[reverse_max_index]) reverse_max_index = it;
                }
                
                if(score_reverse[reverse_max_index] > current_best_score) {
                    current_best_score = score_reverse[reverse_max_index];
                    current_best_contig_match = j;
                    current_best_contig_e1 = get<0>(solid_matches_reverse[reverse_max_index]);
                    current_best_contig_e2 = get<1>(solid_matches_reverse[reverse_max_index]);
                    
                    int start_from = reverse_max_index;
                    while(pred_reverse[start_from] != start_from) start_from = pred_reverse[start_from];
                    current_best_contig_s1 = get<0>(solid_matches_reverse[start_from]);
                    current_best_contig_s2 = get<1>(solid_matches_reverse[start_from]);
                    
                    current_best_orientation = -1;
                }
                
                max_score_reverse = score_reverse[reverse_max_index];
                
                
                all_matches.push_back(make_tuple(current_best_contig_s1, current_best_contig_e1, current_best_contig_s2, current_best_contig_e2, current_best_score, current_best_orientation));
                
                if(current_best_score > best_score) {
                    best_score = current_best_score;
                    best_contig_match = current_best_contig_match;
                    best_contig_e1 = current_best_contig_e1;
                    best_contig_e2 = current_best_contig_e2;
                    best_contig_s1 = current_best_contig_s1;
                    best_contig_s2 = current_best_contig_s2;
                    best_orientation = current_best_orientation;
                }
            }
            
        }
        
        vector<tuple<int, int, int, int, int, int> > results;
        if(best_score > -1) {
            
            #pragma omp critical
            {
                out_debug << best_contig_s1 << " " << best_contig_e1 << " " << best_contig_s2 << " " << " " << best_contig_e2 << " "  << best_score << " " << best_orientation << endl;
            }
            
            results.push_back(make_tuple(best_contig_s1, best_contig_e1, best_contig_s2, best_contig_e2, best_score, best_orientation));
            
            for(int am = 0; am < all_matches.size(); am++) if(get<4>(all_matches[am]) > best_score / 10) results.push_back(all_matches[am]);
        }
        
        
        #pragma omp critical 
        {
            if(best_score > -1) {
                tuple<int, int, int, int, int, int> result;
                for(int zz = 0; zz < results.size(); zz++) {
                    result = results[zz];
                    output1 << contig_names[i] << "\t" << get<0>(result) << "\t" << get<1>(result) << "\t" << contig_names[best_contig_match] << "\t" << get<2>(result) << "\t" << get<3>(result) << "\t" << get<4>(result) << "\t" << get<5>(result) << "\t";
                    if(zz == 0) cout << "P" << endl;
                    else cout << "EXT" << endl;
                } 
            }
            else output2 << contig_names[i] << "\t-1\t-1\tNA\t-1\t-1\t-1\t-1" << endl;
        }
    }
    output2.close();
    
    if(is_debug) out_debug.close();
}
