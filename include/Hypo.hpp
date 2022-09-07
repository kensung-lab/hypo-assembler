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
/** Class Hypo.
 * It is the master class carrying out the polishing.
 */
#pragma once
#include <unordered_map>
#include <thread>
#include <memory>
#include <chrono>

#include "slog/Monitor.hpp"

#include "globalDefs.hpp"
#include "Contig.hpp"
#include "Alignment.hpp"

namespace hypo
{
    const UINT16 HTSIZE = 20; // About 1M bp
    const int MIN_RP_SUPP = 5;
    const double LI_READ_PERC1 = 0.6;
    const double LI_READ_PERC2 = 0.8;
    const double LI_READ_PERC3 = 0.25;
class SolidKmers;

class Hypo {
public:
    Hypo(const InputFlags&);
    ~Hypo() = default;
    void polish(const Mode mode, const FileNames& filenames);
private:
    const InputFlags& _cFlags;
    std::ofstream _gStagefile;
    std::ofstream _remapfile;
    bool _sk_done;
    UINT32 _contig_batch_size;
    samFile *_sf_short;
    bam_hdr_t *_sam_header_short;
    bam1_t *_hts_align_short;
    samFile *_sf_long;
    bam_hdr_t *_sam_header_long;
    bam1_t *_hts_align_long;
    slog::Monitor _monitor;
    
    std::unordered_map<std::string, UINT32> _cname_to_id;
    std::vector<std::unique_ptr<Contig>> _contigs;
    std::vector<std::vector<std::unique_ptr<Alignment>>> _lalignment_store;
    std::vector<std::vector<std::unique_ptr<Alignment>>> _salignment_store;

    void polish_so(const UINT32 batch_id, const UINT32 initial_cid, const UINT32 final_cid, const Mode mode);
    void polish_lo(const UINT32 batch_id, const UINT32 initial_cid, const UINT32 final_cid, const Mode mode);

    void create_alignments(const bool is_sr, const Mode mode, const UINT32 batch_id); 
    void create_alignments_old(const bool is_sr, const Mode mode, const UINT32 batch_id); 
    std::unique_ptr<suk::SolidKmers> get_solid_kmers(const FileNames& filenames);
    void write_read(const UINT32 cid, std::ofstream& remapfile);
    void long_insert(UINT32 cid);
    void long_indel (const UINT32 cid, const std::vector<UINT32>& bp, const std::vector<UINT32>& ep);
    void gap_long_ins (const UINT32 cid, const std::vector<UINT32>& bp, const std::vector<UINT32>& ep);


    inline void add_read_to_realign(const std::string name) {
        // Remove multiplenull characters at the end
        INT e = name.size();
        while (name[e-1]=='\0') {--e;}
        std::string trimmed(name,0,e);
        _remapfile << trimmed <<std::endl;
    }
    inline void print_time() {
        std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        std::string s(30, '\0');
        std::strftime(&s[0], s.size(), "%Y-%m-%d %H:%M:%S", std::localtime(&now));
        std::cout << " " << s << " ";
    }
}; // Hypo
} // namespace hypo
