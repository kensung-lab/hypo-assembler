#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iterator>
#include <cstdint>

class CustomBitvector {
    public:
        CustomBitvector(const uint64_t size) : _bitvector(size) {}
        ~CustomBitvector() = default;
        
        inline bool at(uint64_t i) const { return _bitvector[i]; }
        inline void assign(uint64_t i, bool b) { _bitvector[i] = b; }
        inline uint64_t count() const { return std::count(_bitvector.begin(), _bitvector.end(), true); }
        inline uint64_t size() const { return _bitvector.size(); }
        
        inline void write(std::string out) const { std::ofstream f(out); std::copy(_bitvector.begin(), _bitvector.end(), std::ostream_iterator<bool>(f, " ")); }
        inline void load(std::string in) { std::ifstream f(in); std::copy(std::istream_iterator<bool>(f), {}, std::back_inserter(_bitvector)); }
        
        inline bool operator[] (uint64_t i) { return _bitvector[i]; }
        inline bool operator[] (uint64_t i) const { return _bitvector[i]; }
        
    private:
        std::vector<bool> _bitvector;
        
};
