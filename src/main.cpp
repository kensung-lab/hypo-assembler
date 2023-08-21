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

/** Module containing main() method.
 */

#include <sys/time.h>
#include <getopt.h>
#include <string>
#include <cctype>
#include <math.h> 
#include <algorithm>
#include <sys/stat.h>
#include <sys/wait.h>
#include <cstdlib>
#include "slog/Monitor.hpp"

#include "globalDefs.hpp"
#include "utils.hpp"
#include "Hypo.hpp"
#include "Polish.hpp"

/** Module containing main() method for reading and processing arguments.
 */

namespace hypo{
void usage (void);
void decodeFlags(int argc, char* argv [], InputFlags& flags, FileNames& filenames);
std::vector<std::string> split(const std::string &str, char delimiter);
void print_byte_to_strings(int nb);
UINT get_kmer_len(const std::string& arg);
UINT64 get_given_size(const std::string& arg);
UINT16 get_expected_file_sz(const std::string& given_size, UINT16 cov);
void set_kind(const std::string& kind);
    
static struct option long_options[] = {
    {"short-read-1", required_argument, NULL, '1'},
    {"short-read-2", required_argument, NULL, '2'},
    {"hic-read-1", required_argument, NULL, '3'},
    {"hic-read-2", required_argument, NULL, '4'},
    {"long-read", required_argument, NULL, 'l'},
    {"size-ref", required_argument, NULL, 's'},
    {"coverage-short", required_argument, NULL, 'c'},
    {"coverage-long", required_argument, NULL, 'C'},
    {"working-dir", required_argument, NULL, 'w'},
    {"threads", required_argument, NULL, 't'},
    {"initial-contigs", required_argument, NULL, 'i'},
    {"flye-path", required_argument, NULL, 'F'},
    {"samtools-path", required_argument, NULL, 'S'},
    {"minimap2-path", required_argument, NULL, 'M'},
    {"nanopore-type", required_argument, NULL, 'n'},
    {"samtools-thread", required_argument, NULL, '@'},
    {"samtools-memory", required_argument, NULL, 'Z'},
    {"samtools-temp", required_argument, NULL, 'X'},
    {"debug-mode", no_argument, NULL, 'I'},   
    {NULL, 0, NULL, 0}};

inline bool file_exists (const std::string& name) {
  struct stat st;   
  return (stat (name.c_str(), &st) == 0); 
}

inline bool create_dir (const std::string& name) {
  struct stat st; 
  int status = 0;  
  if (stat(name.c_str(), &st) == -1) {
    status = mkdir(name.c_str(), 0777);
  }
  return (status==0);
}

/** Define various settings as globals
  */
const SRSettings Sr_settings = {5u,0.8};
const MinimizerSettings Minimizer_settings = {10u,10u,5u,0.8,0x000000u,0x055555u,0x0aaaaau,0x0fffffu};
WindowSettings Window_settings = {100u,500u,80u};
//WindowCovSettings SWindow_cov_settings = {50u,1.4,70};
WindowCovSettings SWindow_cov_settings = {50u,1.2,0.3,70,30};
WindowCovSettings LWindow_cov_settings = {50u,1.4,0.1,70,30};
FilteringSettings Filtering_settings = {0.1, 1000u};
//const ArmsSettings Arms_settings = {3u,20u,5u,10u,10u,0.4,10u,20u,10u,2u,0.4};
const ArmsSettings Arms_settings = {3u,10u,5u,10u,10u,0.3,10u,20u,20u,2u,0.4};


/** Decode the input flags
   */
void decodeFlags(int argc, char *argv[], InputFlags &flags, FileNames& filenames)
{
    int opt;
    std::string infile;

    std::string kind="sr";
    filenames.lr_bam_filename = "";
    filenames.sr_bam_filename = "";
    filenames.output_filename = "";

    flags.map_qual_th = 2;
    flags.norm_edit_th = 30;
    flags.threads = 1;
    flags.processing_batch_size = 1;
    flags.intermed = false;
    flags.sz_in_gb = 12;
    flags.wdir = "./";

    flags.output_directory = "hypo_wd/";
    flags.flye_path = "flye";
    flags.samtools_path = "samtools";
    flags.minimap2_path = "minimap2";
    
    flags.samtools_threads = 1;
    flags.samtools_memory = "";
    flags.samtools_temp = "";
    
    flags.initial_contigs = "";

    flags.nano_type = 0;

    flags.threads = 1;

    flags.run_mode = "full";
    flags.genome_size = 3000000000;
    
    bool short_1_check = false;
    bool short_2_check = false;
    
    bool hic_1_check = false;
    bool hic_2_check = false;
    
    bool long_check = false;
    
    bool size_check = false;
    bool cov_check = false;
    bool cov2_check = false;
    bool thread_check = false;
    bool workdir_check = false;
    
    std::string given_sz;
    std::string cmd = "hypo ";
    std::string err_string = "";
    std::string cv = "";
    int cov = 0;

    /* static struct option long_options[] = {
    {"short-read-1", required_argument, NULL, '1'},
    {"short-read-2", required_argument, NULL, '2'},
    {"hic-read-1", required_argument, NULL, '3'},
    {"hic-read-2", required_argument, NULL, '4'},
    {"long-read", required_argument, NULL, 'l'},
    {"size-ref", required_argument, NULL, 's'},
    {"coverage-short", required_argument, NULL, 'c'},
    {"coverage-long", required_argument, NULL, 'C'},
    {"working-dir", required_argument, NULL, 'w'},
    {"threads", required_argument, NULL, 't'},
    {"flye-path", required_argument, NULL, 'F'},
    {"samtools-path", required_argument, NULL, 'S'},
    {"minimap2-path", required_argument, NULL, 'M'},
    {"nanopore-type", required_argument, NULL, 'n'},    
    {NULL, 0, NULL, 0}}; */
    
    std::string name;
    
    /* initialisation */
    while ((opt = getopt_long(argc, argv, "1:2:3:4:l:s:c:C:w:t:F:S:M:n:@:Z:X:i:I", long_options,
                        nullptr)) != -1) {
        switch (opt) {
            case '1':
                name = std::string(optarg);
                fprintf(stderr, "[Hypo] Short read file 1 is %s\n", name.c_str());
                if(!file_exists(name)) {
                    fprintf(stderr, "[Hypo::utils] Error: File Error: Reads file does not exist %s!\n",name.c_str());
                    exit(1);
                }
                flags.short_path_1 = name;
                short_1_check = true;
                break;
            case '2':
                name = std::string(optarg);
                fprintf(stderr, "[Hypo] Short read file 2 is %s\n", name.c_str());
                if(!file_exists(name)) {
                    fprintf(stderr, "[Hypo::utils] Error: File Error: Reads file does not exist %s!\n",name.c_str());
                    exit(1);
                }
                flags.short_path_2 = name;
                short_2_check = true;
                break;
            case '3':
                name = std::string(optarg);
                fprintf(stderr, "[Hypo] Hi-C read file 1 is %s\n", name.c_str());
                if(!file_exists(name)) {
                    fprintf(stderr, "[Hypo::utils] Error: File Error: Reads file does not exist %s!\n",name.c_str());
                    exit(1);
                }
                flags.hic_path_1 = name;
                hic_1_check = true;
                break;
            case '4':
                name = std::string(optarg);
                fprintf(stderr, "[Hypo] Hi-C read file 2 is %s\n", name.c_str());
                if(!file_exists(name)) {
                    fprintf(stderr, "[Hypo::utils] Error: File Error: Reads file does not exist %s!\n",name.c_str());
                    exit(1);
                }
                flags.hic_path_2 = name;
                hic_2_check = true;
                break;
            case 'l':
                name = std::string(optarg);
                fprintf(stderr, "[Hypo] Long reads file is %s\n", name.c_str());
                if(!file_exists(name)) {
                    fprintf(stderr, "[Hypo::utils] Error: File Error: Reads file does not exist %s!\n",name.c_str());
                    exit(1);
                }
                flags.long_path = name;
                long_check = true;
                break;
            case 's':  
                flags.k = std::max(2U,get_kmer_len(std::string(optarg)));
                flags.genome_size = get_given_size(std::string(optarg));
                cmd += (" -s " + std::string(optarg));
                given_sz = std::string(optarg);
                UINT power;
                fprintf(stderr, "[Hypo] Estimated genome size is %s\n", given_sz.c_str());
                size_check = true;
                break;
            case 'c': 
                if (atoi(optarg) <= 0 ) {
                    fprintf(stderr, "[Hypo::utils] Error: Arg Error: Coverage should be positive %d!\n",atoi(optarg));
                    exit(1);
                }   
                cov = atoi(optarg);
                SWindow_cov_settings.mean_cov = cov;
                SWindow_cov_settings.high_th = UINT16(std::ceil(SWindow_cov_settings.high_frac * cov));
                SWindow_cov_settings.low_th = UINT16(std::ceil(SWindow_cov_settings.low_frac * cov));
                cmd += (" -c " + std::string(optarg));
                cv += (" SHORT: "+ std::to_string(SWindow_cov_settings.high_th)+" and " + std::to_string(SWindow_cov_settings.low_th));
                
                fprintf(stderr, "[Hypo] Estimated short read coverage is %d\n", cov);
                
                cov_check = true;
                break;
            case 'C': 
                if (atoi(optarg) <= 0 ) {
                    fprintf(stderr, "[Hypo::utils] Error: Arg Error: Coverage should be positive %d!\n",atoi(optarg));
                    exit(1);
                }   
                cov = atoi(optarg);
                LWindow_cov_settings.mean_cov = cov;
                LWindow_cov_settings.high_th = UINT16(std::ceil(LWindow_cov_settings.high_frac * cov));
                LWindow_cov_settings.low_th = UINT16(std::ceil(LWindow_cov_settings.low_frac * cov));
                cmd += (" -C " + std::string(optarg));
                cv += (" LONG: "+ std::to_string(LWindow_cov_settings.high_th)+" and " + std::to_string(LWindow_cov_settings.low_th));
                
                
                fprintf(stderr, "[Hypo] Estimated long read coverage is %d\n", cov);
                
                cov2_check = true;
                break;
            case 'k':
                kind = std::string(optarg);
                cmd += (" -k " + kind);
                break;
            case 't':
                if (atoi(optarg) <= 0) {
                    fprintf(stderr, "[Hypo::utils] Error: Arg Error: Number of threads (t) must be positive %d!\n",atoi(optarg));
                    exit(1);
                }
                //flags.threads = std::max((UINT32)atoi(optarg)-1,(UINT32)1);
                flags.threads = std::max((UINT32)atoi(optarg),(UINT32)1);
                cmd += (" -t " + std::string(optarg));
                
                fprintf(stderr, "[Hypo] Using %d threads.", flags.threads);
                
                thread_check = true;
                break;
            case 'w':
                flags.output_directory = (std::string(optarg)+"/");
                
                cmd += (" -w " + std::string(optarg));
                
                workdir_check = true;
                break;
            case 'F':
                flags.flye_path = std::string(optarg);
                break;
            case 'S':
                flags.samtools_path = std::string(optarg);
                break;
            case 'M':
                flags.minimap2_path = std::string(optarg);
                break;
            case 'Z':
                flags.samtools_memory = std::string(optarg);
                break;
            case '@':
                if (atoi(optarg) <= 0) {
                    fprintf(stderr, "[Hypo::utils] Error: Arg Error: Number of samtools threads (@) must be positive %d!\n",atoi(optarg));
                    exit(1);
                }
                flags.samtools_threads = std::max((UINT32)atoi(optarg),(UINT32)1);
                break;
            case 'X':
                flags.samtools_temp = std::string(optarg);
                break;
            case 'i':
                name = std::string(optarg);
                fprintf(stderr, "[Hypo] Initial contigs file is %s\n", name.c_str());
                if(!file_exists(name)) {
                    fprintf(stderr, "[Hypo::utils] Error: File Error: Initial contigs file does not exist %s!\n",name.c_str());
                    exit(1);
                }
                flags.initial_contigs = name;
                break;
            case 'I':
                break;
            default:
                usage();
                exit(0);
        }
    }
    
    // arguments checker
    if(!short_1_check || !short_2_check) {
        fprintf(stderr, "[Hypo] Short reads file must be provided\n");
        exit(1);
    }
    
    if(!hic_1_check || !hic_2_check) {
        fprintf(stderr, "[Hypo] Hi-C reads file must both be provided. Running without Hi-C.\n");
        flags.run_hic = false;
    } else {
        flags.run_hic = true;
    }
    
    if(!long_check) {
        fprintf(stderr, "[Hypo] Long reads file must be provided.\n");
        exit(1);
    }
    
    if(!size_check) {
        fprintf(stderr, "Estimated size not provided, defaults to 3G for human genome.\n");
    }
    
    fprintf(stderr, "[Hypo] Running on %d threads.\n", flags.threads);
        
    if(!workdir_check) {
        fprintf(stderr, "Work directory not provided, defaults to hypo_wd on current working directory.\n");
    }
    
    if (!file_exists(flags.output_directory)) {
        fprintf(stderr, "[Hypo] Working directory %s does not exist. Creating.\n",flags.output_directory.c_str());
        int mkwdir = mkdir(flags.output_directory.c_str(), 0755);
        if(mkwdir == -1) {
            fprintf(stderr, "[Hypo] ERROR: Failed to create directory.\n");
            exit(1);
        }
    }
}

/*
   * Usage of the tool
   */
void usage(void)
{
    std::cout << "Param error" << std::endl;
    /*
  std::cout << "\n Hypo Version:"<<VERSION<< "\n\n";
  std::cout << "\n Usage: ./hypo <args>\n\n";
  std::cout << " ****** Mandatory args:\n";
  std::cout << "\t-d, --draft <str>\n"
            << "\tInput file name containing the draft contigs (in fasta/fastq format; can be compressed). \n\n";
  std::cout << "\t-s, --size-ref <str>\n"
            << "\tApproximate size of the genome (a number; could be followed by units k/m/g; e.g. 10m, 2.3g). \n\n\n";

  std::cout << " ****** Mode dependent (mandatory) args (Mode will be selected based on the given arguments):\n";
  std::cout << " ** Long Only Mode (LO):\n";
  std::cout << "\t-B, --bam-lr <str>\n"
            << "\tInput file name containing the alignments of long reads against the draft (in bam/sam format; must have CIGAR information). \n";

  std::cout << "\t-C, --coverage-long <int>\n"
            << "\tApproximate mean coverage of the long reads. \n\n";

  std::cout << " ** Short Only Mode (SO):\n";
  std::cout << "\t-r, --reads-short <str>\n"
            << "\tInput file name containing reads (in fasta/fastq format; can be compressed). "
            << "\t\tA list of files containing file names in each line can be passed with @ prefix.\n"
            << "\t\tIn the list of files containing file names (passed with @ prefix), first line contains R1 of paired-end reads and the second line contina the name of R2.\n";
  std::cout << "\t-b, --bam-sr <str>\n"
            << "\tInput file name containing the alignments of short reads against the draft (in bam/sam format; must have CIGAR information). \n";
  std::cout << "\t-c, --coverage-short <int>\n"
            << "\tApproximate mean coverage of the short reads. \n\n";

  std::cout << " ** LONG + SHORT Mode (LSA):\n";
  std::cout << "\t-B, --bam-lr <str>\n"
            << "\tInput file name containing the alignments of long reads against the draft (in bam/sam format; must have CIGAR information). \n";
  std::cout << "\t-C, --coverage-long <int>\n"
            << "\tApproximate mean coverage of the long reads. \n";
  std::cout << "\t-r, --reads-short <str>\n"
            << "\tInput file name containing reads (in fasta/fastq format; can be compressed).\n "
            << "\t\tA list of files containing file names in each line can be passed with @ prefix.\n"
            << "\t\tIn the list of files containing file names (passed with @ prefix), first line contains R1 of paired-end reads and the second line contina the name of R2.\n";
  std::cout << "\t-b, --bam-sr <str>\n"
            << "\tInput file name containing the alignments of short reads against the draft (in bam/sam format; must have CIGAR information). \n";
  std::cout << "\t-c, --coverage-short <int>\n"
            << "\tApproximate mean coverage of the short reads. \n\n";

  std::cout << " ****** Optional args:\n";

  std::cout << " ** Short reads related (used in SO, LSA modes):\n";
  std::cout << "\t-k, --kind-sr <str>\n"
            << "\tKind of the short reads. \n"
            << "\t[Valid Values] \n"
            << "\t\tsr\t(Corresponding to NGS reads like Illumina reads) \n"
            << "\t\tccs\t(Corresponding to HiFi reads like PacBio CCS reads) \n"
            << "\t[Default] sr.\n\n ";

  std::cout << " ** Output related:\n";
  std::cout << "\t-o, --output <str>\n"
            << "\tOutput file name. \n"
            << "\t[Default] hypo_<draft_file_name>.fasta in the working directory.\n ";
  std::cout << "\t-i, --intermed\n"
            << "\tStore or use (if already exist) the intermediate files in the working directory. \n"
            << "\t[Currently, only Solid kmers are stored as an intermediate file.].\n "
            << "\t[Always Set for LSA mode; Not useful for LO mode; Setting/unsetting only affects SO mode.].\n\n ";
  std::cout << "\t-w, --wdir <str>\n"
            << "\tPath to working directory. \n"
            << "\t[Default] Present Working Directory (.).\n\n ";
  
  std::cout << " ** Processing related:\n";
  std::cout << "\t-t, --threads <int>\n"
            << "\tNumber of threads. \n"
            << "\t[Default] 1.\n\n ";
  std::cout << "\t-p, --processing-size <int>\n"
            << "\tNumber of contigs to be processed in one batch. Lower value means less memory usage. \n"
            << "\t[Default] One contig.\n\n ";
    
  std::cout << "\t-h, --help\n"
            << "\tPrint the usage. \n\n"; */
}


std::vector<std::string> split(const std::string &str, char delimiter) {
  std::string token;
  std::istringstream ss{str};
  std::vector<std::string> tokens;
  while (std::getline(ss, token, delimiter))
  {
    tokens.push_back(token);
  }
  return tokens;
}

void print_byte_to_strings(int nb) {
    // 0x3 for 2-bits and 0x7 for 4-bits
    BYTE base_mask = (BYTE(1) << nb)-1;
    
    // for 2-bits, every BYTE is valid
    if (nb==2) {
        for (UINT i=0; i < 256; ++i) {
            auto num = BYTE(i);
            if (i%16==0) {
              std::cout << "\n\t";
            }
            std::string token(4,'A');
            for (UINT j=4; j > 0; --j) {
                token[j-1] = cCode[(num & base_mask)];
                num = num >> nb;
            }
            if (i%4==0 && i%16!=0) {
              std::cout << " ";
            }
            std::cout << "\"" << token << "\"" << ", ";
        }
    }
    else { // for 4-bits, not every UINT8 is valid
        std::string token("NN");
        for (UINT i=0; i < 256; ++i) {
            auto num = BYTE(i);
            if (i%16==0) {
              std::cout << "\n\t";
            }
            BYTE second_qbit = num & base_mask;
            if (second_qbit < cCode.size() ) { // A,C,G,T,N
                token[1] = cCode[second_qbit];
            }
            num = num >> nb;
            BYTE first_qbit = num & base_mask;
            if (first_qbit < cCode.size()) { // A,C,G,T,N
                token[0] = cCode[first_qbit];
            }
            if (i%4==0 && i%16!=0) {
              std::cout << " ";
            }
            std::cout << "\"" << token << "\"" << ", ";
        }
    }       
}

uint64_t get_given_size(const std::string& given_size) {
    // Find minimum k such that it is odd and  4^k >= genome_size => 2**2k >= genome_size => 2k = log_2(genome_size)
  UINT power;
  size_t ind;
  auto val = std::stof(given_size,&ind);
  if (ind >= given_size.size()) { // should be absolute number
    if (floor(val)!=ceil(val)) {
      fprintf(stderr, "[Hypo::Utils] Error: Wrong format for genome-size: Genome-size with no units (K,M,G etc,) should be absolute number!\n");
      exit(1);
    }
    power = 0;
  }
  else {
    char unit = std::toupper(given_size[ind]);
    switch (unit)
    {
    case 'K':
      power = 10;
      break;
    case 'M':
      power = 20;
      break;
    case 'G':
      power = 30;
      break;
    case 'T':
      power = 40;
      break;
    default:
      fprintf(stderr, "[Hypo::Utils] Error: Wrong format for genome-size: Allowed units for Genome-size are K (10^3),M (10^6),G (10^9),T (10^12)!\n");
      exit(1);
    }
  }
  return val;
}

UINT get_kmer_len(const std::string& given_size) {
  // Find minimum k such that it is odd and  4^k >= genome_size => 2**2k >= genome_size => 2k = log_2(genome_size)
  UINT power;
  size_t ind;
  auto val = std::stof(given_size,&ind);
  if (ind >= given_size.size()) { // should be absolute number
    if (floor(val)!=ceil(val)) {
      fprintf(stderr, "[Hypo::Utils] Error: Wrong format for genome-size: Genome-size with no units (K,M,G etc,) should be absolute number!\n");
      exit(1);
    }
    power = 0;
  }
  else {
    char unit = std::toupper(given_size[ind]);
    switch (unit)
    {
    case 'K':
      power = 10;
      break;
    case 'M':
      power = 20;
      break;
    case 'G':
      power = 30;
      break;
    case 'T':
      power = 40;
      break;
    default:
      fprintf(stderr, "[Hypo::Utils] Error: Wrong format for genome-size: Allowed units for Genome-size are K (10^3),M (10^6),G (10^9),T (10^12)!\n");
      exit(1);
    }
  }
  UINT kmer_len = (power) + ceil(log2(val));
  kmer_len = ceil(kmer_len/2); // divide by 2
  if (kmer_len%2==0) {++kmer_len;}
  fprintf(stdout, "[Hypo::Utils] Info: Value of K chosen for the given genome size (%s): %u\n",given_size.c_str(),kmer_len);
  return kmer_len;
}

UINT16 get_expected_file_sz(const std::string& given_size, UINT16 cov) {
  // Find minimum k such that it is odd and  4^k >= genome_size => 2**2k >= genome_size => 2k = log_2(genome_size)
  UINT64 sz_in_gb=1024;
  size_t ind;
  auto val = std::stof(given_size,&ind);
  val = 2 * cov * val;
  if (ind >= given_size.size()) { // should be absolute number
    if (floor(val)!=ceil(val)) {
      fprintf(stderr, "[Hypo::Utils] Error: Wrong format for genome-size: Genome-size with no units (K,M,G etc,) should be absolute number!\n");
      exit(1);
    }
    else {
      sz_in_gb = val/1000000000;
    }
  }
  else {
    char unit = std::toupper(given_size[ind]);
    switch (unit)
    {
    case 'K':
      sz_in_gb = val/1000000;
      break;
    case 'M':
      sz_in_gb = val/1000;
      break;
    case 'G':
      sz_in_gb = val;
      break;
    case 'T':
      sz_in_gb = 1024;
      break;
    default:
      fprintf(stderr, "[Hypo::Utils] Error: Wrong format for genome-size: Allowed units for Genome-size are K (10^3),M (10^6),G (10^9),T (10^12)!\n");
      exit(1);
    }
  }
  if (sz_in_gb < 12) {sz_in_gb=12;}
  else if (sz_in_gb > 1024) {sz_in_gb=1024;}
  fprintf(stdout, "[Hypo::Utils] Info: File size expected for the given genome size (%s) and cov (%u): %luG\n",given_size.c_str(),cov,sz_in_gb); 
  return UINT16(sz_in_gb);
}

void set_kind(const std::string& kind) {
  if (kind == "sr") {
    Window_settings.ideal_swind_size = 100;
    Window_settings.wind_size_search_th = 80;
  }
  else if (kind == "ccs") {
    Window_settings.ideal_swind_size = 500;
    Window_settings.wind_size_search_th = 400;
  }
  else {
    fprintf(stderr, "[Hypo::Utils] Error: Wrong value for kind-sr: Allowed values are <sr>, <ccs>!\n");
    exit(1);
  }
}


} // end namespace

int main(int argc, char **argv) {
    
    hypo::InputFlags flags;
    hypo::FileNames filenames;
    hypo::decodeFlags(argc, argv, flags, filenames);
    
    hypo::Objects main_objects;
    
    slog::Monitor monitor;
    std::string tm="";
    
    
    // full pipeline: start with assembly and mapping of reads
    if(flags.run_mode == "full") {
        
        /*
         * First step: run FlyE
         */
         
         if(flags.initial_contigs.size() == 0) {
                std::string flye_command = flags.flye_path + " ";
                
                if(flags.nano_type == 0)  // two different modes of FlyE for old and new nanopore reads
                    flye_command += "--nano-raw ";
                else 
                    flye_command += "--nano-hq ";
                
                flye_command += flags.long_path;
                
                flye_command += " --genome-size " + std::to_string(flags.genome_size);
                flye_command += " --threads " + std::to_string(flags.threads) ;
                flye_command += " --out-dir " + flags.output_directory + "/flye/";
                flye_command += " -i 0";
                
                
                monitor.start();
                
                std::cout << flye_command << std::endl;
                system(flye_command.c_str());
                
                tm = monitor.stop("[Hypo:main]: Run FlyE.");
                fprintf(stdout, "[Hypo::main] //////////////////\n Running FlyE done. \n [Hypo::Hypo] ////////////////// \n%s\n", tm.c_str());
        
            flags.initial_assembly_path = flags.output_directory + "/flye/assembly.fasta";
        } else {
            flags.initial_assembly_path = flags.initial_contigs;
        }
        
        // included in flye run: map short and long reads
        std::string minimap2_short_command = flags.minimap2_path + " -ax sr -t " + std::to_string(flags.threads)  + " ";
        minimap2_short_command += flags.initial_assembly_path + " " + flags.short_path_1 + " " + flags.short_path_2;
        minimap2_short_command += " | " + flags.samtools_path + " view -bS | ";
        minimap2_short_command += flags.samtools_path + " sort -@ " + std::to_string(flags.samtools_threads);
        if(flags.samtools_memory.size() > 0) minimap2_short_command += " -m " + flags.samtools_memory;
        if(flags.samtools_temp.size() > 0) minimap2_short_command += " -T " + flags.samtools_temp;
        minimap2_short_command += " -o " + flags.output_directory + "/short_read_initial.bam";
        
        monitor.start();
        
        std::cout << minimap2_short_command << std::endl;
        system(minimap2_short_command.c_str());
        tm = monitor.stop("[Hypo:main]: Minimap2: Mapping short reads to initial assembly.");
        fprintf(stdout, "[Hypo::main] //////////////////\n Minimap2 - short reads done. \n [Hypo::Hypo] ////////////////// \n%s\n", tm.c_str());
        
        flags.short_initial_mapping_path = flags.output_directory + "/short_read_initial.bam";
        
        std::string minimap2_long_command = flags.minimap2_path + " -ax map-ont -t " + std::to_string(flags.threads) + " ";
        minimap2_long_command += flags.initial_assembly_path + " " + flags.long_path;
        minimap2_long_command += " | " + flags.samtools_path + " view -bS | ";
        minimap2_long_command += flags.samtools_path + " sort -@ " + std::to_string(flags.samtools_threads);
        if(flags.samtools_memory.size() > 0) minimap2_long_command += " -m " + flags.samtools_memory;
        if(flags.samtools_temp.size() > 0) minimap2_long_command += " -T " + flags.samtools_temp;
        minimap2_long_command += " -o " + flags.output_directory + "/long_read_initial.bam";
        
        monitor.start();
        std::cout << minimap2_long_command << std::endl;
        
        system(minimap2_long_command.c_str());
        tm = monitor.stop("[Hypo:main]: Minimap2: Mapping long reads to initial assembly.");
        fprintf(stdout, "[Hypo::main] //////////////////\n Minimap2 - long reads done. \n [Hypo::Hypo] ////////////////// \n%s\n", tm.c_str());
        
        flags.long_initial_mapping_path = flags.output_directory + "/long_read_initial.bam";
    }
    
    monitor.start();
    hypo::utils::initialize_solid_kmers(main_objects, flags.k, flags.short_path_1, flags.short_path_2);
    tm = monitor.stop("[Hypo:main]: Extracted solid kmers.");
    fprintf(stdout, "[Hypo::main] //////////////////\n Solid kmers initialization done. \n [Hypo::Hypo] ////////////////// \n%s\n", tm.c_str());
    
    // Assembly enchancement or full run: Run the overlap detection
    if(flags.run_mode == "full" || flags.run_mode == "enhance_assembly") { 
        hypo::overlap_detection_main(main_objects, flags);
    }
    else if(flags.run_mode == "polish") {
        hypo::polish_main(main_objects, flags);
    } // else if(flags.mode == "scaffold") {
        // scaffold_main(flags);
    // }
    return 0;
}
