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

#pragma once

#include <cstdlib>
#include <cassert>
#include <fstream>
#include <cstdint>
#include <stdio.h>
#include <iterator>
//#include <tuple>
//#include <utility>
#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <stdexcept>
#include <sstream> 
#include <zlib.h> // for reading compressed .fq file
#include <htslib/kseq.h>
KSEQ_INIT(gzFile, gzread)

#include <htslib/sam.h>

namespace hypo{
#define DEBUG

  // Types
  using UINT = unsigned int;
  using INT = int;
  using INT8 = int8_t;
  using UINT8 = uint8_t;
  using BYTE = uint8_t;
  using INT16 = int16_t;
  using UINT16 = uint16_t;
  using INT32 = int32_t;
  using UINT32 = uint32_t;
  using INT64 = int64_t;
  using UINT64 = uint64_t;

  #define VERSION 2.0

  enum Mode{
  LSA, // Long read polishing with short read alignment on the contigs given
  LO, // Only long reads polishing
  SO, // Only short reads polishing
  SECOND // Second round in LSA
  };

  using InputFlags = struct SInputFlags{
    UINT8 map_qual_th; // mapping_qual_th
    UINT8 norm_edit_th; // normalised edit distance threshold for long reads
    UINT32 threads;
    UINT32 processing_batch_size;
    UINT16 k;
    UINT16 sz_in_gb;
    UINT done_stage;
    Mode mode;
    // whether intermediate results should be used and stored to disk; 
    // if false, neither intermed files will be written nor used (even if exist)
    bool intermed; 
    std::string wdir; //working directory
    
    std::string output_directory;
    std::string short_initial_mapping_path;
    std::string long_initial_mapping_path;
    std::string initial_assembly_path;
    std::string short_path_1;
    std::string short_path_2;
    std::string hic_path_1;
    std::string hic_path_2;
    std::string long_path;
    
    bool run_hic;
    
    std::string flye_path;
    std::string samtools_path;
    std::string minimap2_path;
    
    uint32_t nano_type;
    
    uint64_t genome_size;
    
    std::string run_mode;
    
    uint32_t samtools_threads;
    std::string samtools_memory;
    std::string samtools_temp;
    
    std::string initial_contigs;

  };

using FileNames = struct SFileNames{
    std::vector<std::string> sr_filenames;
    std::string sr_bam_filename;
    std::string lr_bam_filename;
    std::string draft_filename;
    std::string output_filename;
    };


  // Stages
  #define STAGE_BEG 0u
  #define STAGE_SK 3u
  #define STAGE_FIRST 1u
  #define STAGE_REMAP 2u

  // Region types (Needed by Contig and Alignment Classes)
  enum class RegionType: UINT8 {
      SWS,
      SW,
      WS,
      MWM,
      MW,
      WM,
      SWM,
      MWS,
      OTHER,
      LONG,
      SR,
      MSR, // Minimser is same as SR (no polishing)
      PSEUDO,
      NOPOL, // NoPolishnig in this window; consensus same as the draft
      INVALID, // No polishing in this window;consensus is empty string
      LI // Consensus by prefix-suffix
  };
  
  

  // Folders and files
  #define AUX_DIR  "aux/"
  #define SR_DIR "aux/SR/"
  #define SKFILE "aux/solid_kmers.bvsd"
  #define STAGEFILE "aux/stage.txt"
  #define INSPECTFILEPREF "aux/inspect_"
  #define BEDFILE "aux/regions.bed"
  #define SBAMFILE "aux/mapped-sr2.bam"




  // SR related settings
  using SRSettings = struct SSRSettings{
    UINT cov_th;
    double supp_frac;
  };
  extern const SRSettings Sr_settings;
    
  // Minimiser based window cutting thresholds/constants
  using MinimizerSettings = struct SMinimizerSettings{
    UINT k; // Change Poly-base macros, if this k gets changed.; Should be <= 16 (minimiser is 32 bits)
    UINT w;
    UINT cov_th;
    double supp_frac;
    UINT polyA;
    UINT polyC;
    UINT polyG;
    UINT polyT;
  };
  extern const MinimizerSettings Minimizer_settings;
    
  // Window-size related settings (set according to param)
  using WindowSettings = struct SWindowSettings{
    UINT ideal_swind_size;  // Ideal short window size
    UINT ideal_lwind_size;  // Ideal long window size
    UINT wind_size_search_th; //This should always be <= ideal_wind_size
  };
  extern WindowSettings Window_settings;

  // Window-coverage related settings (set according to param)
  using WindowCovSettings = struct SWindowCovSettings{
    UINT16 mean_cov;  // Mean coverage of short reads
    double high_frac; // fraction for deciding upper-cutoff to classify window as high-coverage [high_th = high_frac*mean_cov]
    double low_frac; // fraction for deciding lower-cutoff to classify window as lower-coverage [low_th = low_frac*mean_cov]
    UINT16 high_th; // if cov > (high_th)*mean_cov, it is high
    UINT16 low_th; // if cov < (low_th)*mean_cov, it is low
  };
  extern WindowCovSettings SWindow_cov_settings; // Short window
  extern WindowCovSettings LWindow_cov_settings; // Long window
  
  using FilteringSettings = struct SFilteringSettings{
      double clip_ratio; // filter out reads only if it has (total clip / aligned length) > this threshold
      UINT32 contig_end_leniency; // clips around contig_end_leniency is not considered
  };
  extern FilteringSettings Filtering_settings;

  // Arm-filling realted settings
  using ArmsSettings = struct SArmsSettings{
    UINT min_short_num; // Minimum num of short-reads arms for a window to be declared short
    UINT min_internal_num1; // Minimum num of internal arms for a short-arms window to discard pref/suff arms
    UINT min_internal_num2; // Minimum num of internal arms for a spcl (SW,WS,SWS) short-arms window to discard pref/suff arms
    UINT min_internal_num3; // Minimum num of internal arms for a long-arms window to discard pref/suff arms
    UINT min_contrib; // Minimum num of arms for a short-arms window to be considered for discarding pref/suff arms
    double min_internal_contrib; // Minimum fraction of internal arms (wrt total arms) for a short-arms window to discard pref/suff arms
    UINT short_arm_coef; // Short arm len should be >= window_len/coefficient
    UINT long_arm_min_len; // Minimum Long arm len
    UINT base_qual_th; // threshold for base quality to be considered as high quality
    UINT allowed_num_lq_bases; // number of bases which are allowed to be of low qualilty in an arm (used for Short armonly)
    double lq_frac_th; // threshold of what fraction of internal arms containing low qual bases are tolerated
  };
  extern const ArmsSettings Arms_settings;

  const std::vector<char> cCode = {'A','C','G','T','N'};
// A: 0, C: 1, G:2, T: 3; everything else is 4.
// Unpacking will map everything other ACGT to N.
const BYTE cNt4Table[256] = {
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
  
enum Stage{
  SK, // solid kmer
  SP // solid pos

    };

} // end namespace


