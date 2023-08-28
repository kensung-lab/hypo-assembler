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
/** Class PackedSequence.
 * It represents a sequence with a DNA base with 4-bits letter.
 * Either 2-bits or 4-bits packing
 * 4-bits Valid letters: A,C, G, T,N
 * 2-bits valid letter
 */
#pragma once
#ifndef UTILS_HPP
#define UTILS_HPP
#include "globalDefs.hpp"
#include "Hypo.hpp"

#define MAX_LEN_LIMIT 0xffffffffu


namespace hypo
{
    namespace utils {
        int load_contigs(hypo::Objects & objects, const std::string load_path);
        int write_contigs(const hypo::Objects & objects, const std::string write_path);
        int write_contigs_1(const hypo::Objects & objects, const std::string write_path);
        int write_contigs_2(const hypo::Objects & objects, const std::string write_path);
        
        int load_short_alignments(hypo::Objects & objects, const std::string load_path);
        int load_long_alignments(hypo::Objects & objects, const std::string load_path);
        int initialize_solid_kmers(hypo::Objects & objects, const int k, const std::string load_path_1, const std::string load_path_2);
        
        int realign_long_reads(hypo::Objects & objects, const hypo::InputFlags & flags, const std::string contig_path, const std::string write_path);
        int realign_short_reads(hypo::Objects & objects, const hypo::InputFlags & flags, const std::string contig_path, const std::string write_path);
    }
} // namespace hypo
#endif
