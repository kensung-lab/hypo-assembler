#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iterator>
#include <cstdint>

class CustomBitvector {
    public:
        inline CustomBitvector(const uint64_t input_size) { data = std::vector<uint64_t>((input_size + (uint64_t)63) / (uint64_t)64);}
        ~CustomBitvector() = default;
        
        inline uint64_t real_index(uint64_t i) const { return (i >> (uint64_t)6); }
        inline uint64_t mask(uint64_t i) const { return (uint64_t)1 << (i & (uint64_t)63); }
        
        inline void set(uint64_t i) { data[real_index(i)] |= mask(i); }
        inline bool at(uint64_t i) const { return 0 != (data[real_index(i)] & mask(i)); }
        
        inline bool operator[] (uint64_t i) const { return at(i); }
        
        inline uint64_t count() const { uint64_t total = 0; for(uint64_t i = 0; i < data.size(); i++) total += __builtin_popcountll(data[i]); return total; }
        inline uint64_t size() const { return data.size() * (uint64_t)64; }
        
        inline void write(std::string out) const { std::ofstream f(out); for(uint64_t i = 0; i < data.size(); i++) f << data[i] << "\n"; f.close();  }
        inline void load(std::string in) { std::ifstream f(in); uint64_t read; data.clear(); while(f >> read) data.push_back(read); f.close(); }
        
    private:
        std::vector<uint64_t> data;        
};
