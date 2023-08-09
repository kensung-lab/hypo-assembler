#pragma once
#include <unordered_map>
#include <thread>
#include <memory>
#include <chrono>

#include "slog/Monitor.hpp"

#include "globalDefs.hpp"
#include "Contig.hpp"
#include "Alignment.hpp"

namespace hypo {
    class Objects {
        public:
            std::vector<std::string> contigs;
            std::vector<std::string> contig_name;
            
            std::vector<std::string> contigs_1;
            std::vector<std::string> contig_name_1;
            
            std::vector<std::string> contigs_2;
            std::vector<std::string> contig_name_2;
            
            std::vector<Alignment> short_alignments;
            std::vector<Alignment> long_alignments;
            //SolidKmers solids;
            
            slog::Monitor monitor;
            
    };
}
