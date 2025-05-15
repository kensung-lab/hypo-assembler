#include <zlib.h> // for reading compressed .fq file
#include <iostream>
#include <ctime>
#include <chrono>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <set>
#include <tuple>
#include <vector>
#include <algorithm>
#include <htslib/kseq.h>
#include <htslib/kstring.h>
#include <omp.h>
#include <sys/stat.h>
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

typedef tuple<uint32_t, uint8_t, uint32_t, uint8_t> hash_key;

struct key_hash : public std::unary_function<hash_key, std::size_t> {
    std::size_t operator()(const hash_key & k) const {
        size_t seed = hash<uint32_t>()(get<0>(k));
        seed ^= hash<uint32_t>()(get<2>(k)) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        seed ^= hash<uint8_t>()(get<1>(k)) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        seed ^= hash<uint8_t>()(get<3>(k)) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        return seed;
    }
};

int main(int argc, char* argv[]) {
    size_t THRESHOLD_END_OF_CONTIG = 10000;
    
    auto start_time = chrono::high_resolution_clock::now();
    
    if(argc < 8) {
        cerr << "Usage: " << argv[0] << " <k> <kmers input> <contigs 1 fasta> <contigs 2 fasta> <long reads> <output 1> <output 2> <reads output> <optional: thread count (default=1)> <optional: remove >2 occurence solid kmers (1/0) (default=0)> <optional: debug_output>" << endl;
        return 1;
    }
    cerr << "Command: ";
    for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
    cerr << endl;
    
    bool is_debug = false;
    ofstream out_debug;
    
    ofstream out_reads(argv[8]);;
    
    if(argc > 11) {
        is_debug = true;
    }
    
    if(is_debug) {
        try {
            out_debug.open(argv[11]);
        } catch(exception const &e) {
            cerr << "Can't open file" << endl;
            cerr << e.what() << endl;
            cerr << strerror(errno) << endl;
            return -1;
        }
    }
    
    int thread_count = 1;
    try {
        thread_count = stoi(argv[9]);
    } catch(exception const &e) {
        cerr << "Thread count need to be int. Defaults to 1" << endl;
    }
    
    int remove_dupe = 0;
    try {
        remove_dupe = stoi(argv[10]);
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
    
    cerr << "Debug mode: " << is_debug << endl;
    
    uint32_t SCAFFOLD_THRESHOLD = 5;
    
    unordered_set<uint64_t> solid_kmers;
    
    cerr << "Reading kmers from " << argv[2] << "." << endl;
    struct stat buffer;
    if(stat(argv[2], &buffer) != 0) {
        cerr << "Error, " << string(argv[2]) << " does not exist!" << endl;
        return 1;
    }
    
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
    
    // kmer to contig_id, pos, orientation
    unordered_map<uint64_t, tuple<uint32_t, uint32_t, uint8_t> > kmer_locations;
    
    auto begin_time = chrono::high_resolution_clock::now();
    cerr << "Processing contigs from " << argv[3] << "." << endl;
    if(stat(argv[3], &buffer) != 0) {
        cerr << "Error, " << string(argv[3]) << " does not exist!" << endl;
        return 1;
    }
    
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
        
        unordered_map<uint64_t, vector<tuple<uint8_t, uint64_t> > > current_contig_solids;
        for(uint32_t i = 0; i < seq->seq.l; i++) {
            int c = seq_nt4_table[(uint8_t)seq->seq.s[i]];
            if(c < 4) {
                kmer[0] = (kmer[0] << 2ull | c) & mask;           // forward k-mer
                kmer[1] = (kmer[1] >> 2ull) | (3ULL^c) << shift1; // reverse k-mer
                count_not_N++;
                int z = kmer[0] < kmer[1] ? 0 : 1;
                if(count_not_N >= k) {
                    count_kmers++;
                    if(bv_sd[kmer[z]] && (i <= THRESHOLD_END_OF_CONTIG || i >= seq->seq.l - THRESHOLD_END_OF_CONTIG) ) {
                        //found a solid kmer, do something here
                        current_contig_solids[kmer[z]].push_back(make_tuple(z, i));
                        kmer_counter[kmer[z]]++;
                        kmer_locations[kmer[z]] = make_tuple(contig_names.size(), i, z);
                    }
                }
            } else {
                count_not_N = 0;
            }
        }
        
        contig_names.push_back(read_name);
        contig_lens.push_back(seq->seq.l);
        
        contig_solids.push_back(current_contig_solids);
        
        cerr << "Found " << current_contig_solids.size() << " solids out of " << count_kmers << endl;
        total_kmers += current_contig_solids.size();
    }
    
    if(is_debug) {
        for(auto & x : contig_solids) {
            for(auto & y : x) {
                out_debug << "DG\t" << y.first;
                for(auto & z : y.second) {
                    out_debug << "\t" << get<0>(z) << "\t" << get<1>(z);
                }
                out_debug << "\n";
            }
        }
    }
    
    auto current_time = chrono::high_resolution_clock::now();
    
    chrono::duration<double, milli> time_diff = current_time - begin_time;
    chrono::duration<double, milli> time_from_start = current_time - start_time;
    chrono::duration<double, milli> total_time_scan, total_time_paired, total_time_hash_table;
    
    cerr << "Processing contigs done: " << kmer_counter.size() << endl;
    cerr << "Processing contigs time: " << time_diff.count() << " , time from start: " << time_from_start.count() << endl;
    
    
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
    
    if(is_debug) {
        for(auto & x : contig_solids) {
            for(auto & y : x) {
                out_debug << "DG\t" << y.first;
                for(auto & z : y.second) {
                    out_debug << "\t" << get<0>(z) << "\t" << get<1>(z);
                }
                out_debug << "\n";
            }
        }
    }
    
    cerr << "Processing reads from " << argv[5] << "." << endl;
    if(stat(argv[5], &buffer) != 0) {
        cerr << "Error, " << string(argv[5]) << " does not exist!" << endl;
        return 1;
    }
    
    fp = gzopen(argv[5], "r");
    seq = kseq_init(fp);
    
    // contig 1, orientation 1, contig 2, orientation 2
    unordered_map<tuple<uint32_t, uint8_t, uint32_t, uint8_t>, uint32_t, key_hash> matching_contigs;
    uint64_t minkmer, maxkmer, tempkmer1, tempkmer2;
    
    uint64_t total_reads = 0;
    
    // kmer, orientation
    
    begin_time = chrono::high_resolution_clock::now();
    
    vector<string> bundle_process;
    vector<string> read_names;
    
    while((l = kseq_read(seq)) >= 0) {
        string read_name = seq->name.s;
        auto start_set_time = chrono::high_resolution_clock::now();
        
        bundle_process.push_back(seq->seq.s);
        read_names.push_back(read_name);
        
        
        if(bundle_process.size() == 1000) {
            #pragma omp parallel for num_threads(thread_count)
            for(size_t j = 0; j < bundle_process.size(); j++) {
                int count_not_N = 0;
                int count_kmer = 0;
                int previous_solid = -1;
                uint64_t previous_kmer;
                uint32_t count_solid = 0;
                uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
                
                string get_seq = bundle_process[j];
                vector<tuple<uint64_t, uint8_t> > current_contigs;
                for(uint32_t i = 0; i < get_seq.size(); i++) {
                    int c = seq_nt4_table[(uint8_t)get_seq[i]];
                    if(c < 4) {
                        kmer[0] = (kmer[0] << 2ull | c) & mask;           // forward k-mer
                        kmer[1] = (kmer[1] >> 2ull) | (3ULL^c) << shift1; // reverse k-mer
                        count_not_N++;
                        int z = kmer[0] < kmer[1] ? 0 : 1;
                        if(count_not_N >= k) {
                            count_kmer += 1;
                            if(bv_sd[kmer[z]] && kmer_counter.find(kmer[z]) != kmer_counter.end()) {
                                current_contigs.push_back(make_tuple(get<0>(kmer_locations[kmer[z]]), z ^ get<2>(kmer_locations[kmer[z]])));
                                count_solid += 1;
                            }
                        }
                    } else {
                        count_not_N = 0;
                    }
                }
                
                unordered_set<tuple<uint32_t, uint8_t, uint32_t, uint8_t>, key_hash > to_add;
                for(int i = 0; i < current_contigs.size(); i++) {
                    for(int j = i + 1; j < current_contigs.size(); j++) {
                        uint32_t contig_i = get<0>(current_contigs[i]);
                        uint32_t contig_j = get<0>(current_contigs[j]);
                        
                        if(contig_i != contig_j) {
                            uint8_t orientation_i = get<1>(current_contigs[i]);
                            uint8_t orientation_j = get<1>(current_contigs[j]);
                            
                            to_add.insert(make_tuple(contig_i, orientation_i, contig_j, orientation_j));
                        }
                    }
                }
                
                #pragma omp critical
                {
                    for(auto & x : to_add) matching_contigs[x]++;
                    if(count_solid >= 10) {
                        out_reads << ">" << read_names[j] << "\n" << bundle_process[j] << "\n";
                    }
                }
            }
            bundle_process.clear();
            read_names.clear();
            
            total_reads += 1000;
            current_time = chrono::high_resolution_clock::now();
            time_diff = current_time - begin_time;
            time_from_start = current_time - start_time;
            
            cerr << "Processed " << total_reads << " reads. Size of matching contigs: " << matching_contigs.size() << endl;
            cerr << "Time: " << time_diff.count() << " , time from start: " << time_from_start.count() << endl;
        }
    }
    // go through the reads again to get the sequences, check if can get consensus!
    out_reads.close();
    
    cerr << "Finding scaffolds" << endl;
    
    ofstream output1(argv[6]);
    for(auto & x : matching_contigs) {
        if(x.second >= SCAFFOLD_THRESHOLD) {
            output1 << "R\t" << contig_names[get<0>(x.first)] << "\t" << unsigned(get<1>(x.first)) << "\t" << contig_names[get<2>(x.first)] << "\t" << unsigned(get<3>(x.first)) << "\t" << x.second << "\n";
        }
    }
    
    // initialize ordering based on the solid kmers count
    vector<size_t> ordering;
    for(int i = 0; i < contig_names.size(); i++) ordering.push_back(i);
    sort(ordering.begin(), ordering.end(), [&](size_t i1, size_t i2) { return contig_solids[i1].size() > contig_solids[i2].size(); });
    output1.close();    
    
    contig_names.clear();
    contig_lens.clear();
    contig_solids.clear();
    kmer_counter.clear();
    kmer_locations.clear();
    
    begin_time = chrono::high_resolution_clock::now();
    cerr << "Processing contigs from " << argv[4] << "." << endl;
    if(stat(argv[4], &buffer) != 0) {
        cerr << "Error, " << string(argv[4]) << " does not exist!" << endl;
        return 1;
    }
    
    fp = gzopen(argv[4], "r");
    seq = kseq_init(fp);
    
    total_kmers = 0;
    while((l = kseq_read(seq)) >= 0) {
        string read_name = seq->name.s;
        cerr << "Processing contig " << read_name << endl;
        int count_not_N = 0;
        int count_kmers = 0;
        
        unordered_map<uint64_t, vector<tuple<uint8_t, uint64_t> > > current_contig_solids;
        for(uint32_t i = 0; i < seq->seq.l; i++) {
            int c = seq_nt4_table[(uint8_t)seq->seq.s[i]];
            if(c < 4) {
                kmer[0] = (kmer[0] << 2ull | c) & mask;           // forward k-mer
                kmer[1] = (kmer[1] >> 2ull) | (3ULL^c) << shift1; // reverse k-mer
                count_not_N++;
                int z = kmer[0] < kmer[1] ? 0 : 1;
                if(count_not_N >= k) {
                    count_kmers++;
                    if(bv_sd[kmer[z]] && (i <= THRESHOLD_END_OF_CONTIG || i >= seq->seq.l - THRESHOLD_END_OF_CONTIG) ) {
                        //found a solid kmer, do something here
                        current_contig_solids[kmer[z]].push_back(make_tuple(z, i));
                        kmer_counter[kmer[z]]++;
                        kmer_locations[kmer[z]] = make_tuple(contig_names.size(), i, z);
                    }
                }
            } else {
                count_not_N = 0;
            }
        }
        
        contig_names.push_back(read_name);
        contig_lens.push_back(seq->seq.l);
        
        contig_solids.push_back(current_contig_solids);
        
        cerr << "Found " << current_contig_solids.size() << " solids out of " << count_kmers << endl;
        total_kmers += current_contig_solids.size();
    }
    
    current_time = chrono::high_resolution_clock::now();
    
    time_diff = current_time - begin_time;
    time_from_start = current_time - start_time;
    
    cerr << "Processing contigs done: " << kmer_counter.size() << endl;
    cerr << "Processing contigs time: " << time_diff.count() << " , time from start: " << time_from_start.count() << endl;
    
    
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
    if(stat(argv[5], &buffer) != 0) {
        cerr << "Error, " << string(argv[5]) << " does not exist!" << endl;
        return 1;
    }
    fp = gzopen(argv[5], "r");
    seq = kseq_init(fp);
    
    // contig 1, orientation 1, contig 2, orientation 2
    matching_contigs.clear();
    
    total_reads = 0;
    
    // kmer, orientation
    
    begin_time = chrono::high_resolution_clock::now();
    
    while((l = kseq_read(seq)) >= 0) {
        string read_name = seq->name.s;
        auto start_set_time = chrono::high_resolution_clock::now();
        
        bundle_process.push_back(seq->seq.s);
        read_names.push_back(read_name);
        
        
        if(bundle_process.size() == 1000) {
            #pragma omp parallel for num_threads(thread_count)
            for(size_t j = 0; j < bundle_process.size(); j++) {
                int count_not_N = 0;
                int count_kmer = 0;
                int previous_solid = -1;
                uint64_t previous_kmer;
                uint32_t count_solid = 0;
                uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
                
                string get_seq = bundle_process[j];
                vector<tuple<uint64_t, uint8_t> > current_contigs;
                for(uint32_t i = 0; i < get_seq.size(); i++) {
                    int c = seq_nt4_table[(uint8_t)get_seq[i]];
                    if(c < 4) {
                        kmer[0] = (kmer[0] << 2ull | c) & mask;           // forward k-mer
                        kmer[1] = (kmer[1] >> 2ull) | (3ULL^c) << shift1; // reverse k-mer
                        count_not_N++;
                        int z = kmer[0] < kmer[1] ? 0 : 1;
                        if(count_not_N >= k) {
                            count_kmer += 1;
                            if(bv_sd[kmer[z]] && kmer_counter.find(kmer[z]) != kmer_counter.end()) {
                                current_contigs.push_back(make_tuple(get<0>(kmer_locations[kmer[z]]), z ^ get<2>(kmer_locations[kmer[z]])));
                                count_solid += 1;
                            }
                        }
                    } else {
                        count_not_N = 0;
                    }
                }
                
                unordered_set<tuple<uint32_t, uint8_t, uint32_t, uint8_t>, key_hash > to_add;
                for(int i = 0; i < current_contigs.size(); i++) {
                    for(int j = i + 1; j < current_contigs.size(); j++) {
                        uint32_t contig_i = get<0>(current_contigs[i]);
                        uint32_t contig_j = get<0>(current_contigs[j]);
                        
                        if(contig_i != contig_j) {
                            uint8_t orientation_i = get<1>(current_contigs[i]);
                            uint8_t orientation_j = get<1>(current_contigs[j]);
                            
                            to_add.insert(make_tuple(contig_i, orientation_i, contig_j, orientation_j));
                        }
                    }
                }
                
                #pragma omp critical
                {
                    for(auto & x : to_add) matching_contigs[x]++;
                }
            }
            bundle_process.clear();
            read_names.clear();
            
            total_reads += 1000;
            current_time = chrono::high_resolution_clock::now();
            time_diff = current_time - begin_time;
            time_from_start = current_time - start_time;
            
            cerr << "Processed " << total_reads << " reads. Size of matching contigs: " << matching_contigs.size() << endl;
            cerr << "Time: " << time_diff.count() << " , time from start: " << time_from_start.count() << endl;
        }
    }
    // go through the reads again to get the sequences, check if can get consensus!
    
    cerr << "Finding scaffolds" << endl;
    
    ofstream output2(argv[7]);
    for(auto & x : matching_contigs) {
        if(x.second >= SCAFFOLD_THRESHOLD) {
            output2 << "R\t" << contig_names[get<0>(x.first)] << "\t" << unsigned(get<1>(x.first)) << "\t" << contig_names[get<2>(x.first)] << "\t" << unsigned(get<3>(x.first)) << "\t" << x.second << "\n";
        }
    }
    
    // initialize ordering based on the solid kmers count
    ordering.clear();
    for(int i = 0; i < contig_names.size(); i++) ordering.push_back(i);
    sort(ordering.begin(), ordering.end(), [&](size_t i1, size_t i2) { return contig_solids[i1].size() > contig_solids[i2].size(); });
    output2.close();
}
