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
#include "Hypo.hpp"

/** Module containing main() method for reading and processing arguments.
 */

namespace hypo{
void usage (void);
void decodeFlags(int argc, char* argv [], InputFlags& flags, FileNames& filenames);
std::vector<std::string> split(const std::string &str, char delimiter);
void print_byte_to_strings(int nb);
UINT get_kmer_len(const std::string& arg);
UINT16 get_expected_file_sz(const std::string& given_size, UINT16 cov);
void set_kind(const std::string& kind);

static struct option long_options[] = {
    {"reads-short", required_argument, NULL, 'r'},
    {"draft", required_argument, NULL, 'd'},
    {"size-ref", required_argument, NULL, 's'},
    {"coverage-short", required_argument, NULL, 'c'},
    {"coverage-long", required_argument, NULL, 'C'},
    {"bam-sr", required_argument, NULL, 'b'},
    {"bam-lr", required_argument, NULL, 'B'},
    {"output", required_argument, NULL, 'o'},
    {"threads", required_argument, NULL, 't'},
    {"processing-size", required_argument, NULL, 'p'},
    {"kind-sr", required_argument, NULL, 'k'},
    {"wdir", required_argument, NULL, 'w'},
    {"intermed", no_argument, NULL, 'i'},
    {"help", no_argument, NULL, 'h'},
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
  int args = 0;
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

  bool is_sr = false;
  bool is_draft = false;
  bool is_size = false;
  bool is_scov = false;
  bool is_lcov = false;
  bool is_bamsr = false;
  bool is_bamlr = false;
  std::string given_sz;
  std::string cmd = "hypo ";
  std::string err_string = "";
  std::string cv = "";
  int cov = 0;
  /* initialisation */
  while ((opt = getopt_long(argc, argv, "r:d:s:c:C:b:B:o:t:p:k:w:ihv", long_options,
                            nullptr)) != -1)
  {
    switch (opt)
    {
    case 'r':
      infile = optarg;
      cmd += (" -r " + std::string(optarg));
      if(infile[0]=='@') {
          std::string name;
          std::ifstream sr_list(infile.substr(1,infile.size()).c_str());
          if (!sr_list.good()) {
              fprintf(stderr, "[Hypo::utils] Error: File Error: Could not open the file %s!\n",infile.substr(1,infile.size()).c_str());
              exit(1); 
          }
          while(std::getline(sr_list,name)) {
            if (name.size()>0) {
              filenames.sr_filenames.push_back(name);
              if (!file_exists(name)) {
                fprintf(stderr, "[Hypo::utils] Error: File Error: Reads file does not exist %s!\n",name.c_str());
                exit(1);
              }
            }
        }
      }
      else {
          filenames.sr_filenames.push_back(infile);
          if (!file_exists(infile)) {
            fprintf(stderr, "[Hypo::utils] Error: File Error: Reads file does not exist %s!\n",infile.c_str());
            exit(1);
          }
      }
      is_sr = true;
      args++;
      break;
    
    case 'd':
      filenames.draft_filename = std::string(optarg);
      if (!file_exists(filenames.draft_filename)) {
        fprintf(stderr, "[Hypo::utils] Error: File Error: Draft file does not exist %s!\n",filenames.draft_filename.c_str());
        exit(1);
      }
      cmd += (" -d " + std::string(optarg));
      args++;
      is_draft = true;
      break;

    case 's':    
      flags.k = std::max(2U,get_kmer_len(std::string(optarg)));
      cmd += (" -s " + std::string(optarg));
      given_sz = std::string(optarg);
      is_size = true;
      args++;
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
      is_scov = true;
      args++;
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
      is_lcov = true;
      args++;
      break;

    case 'b':
      filenames.sr_bam_filename = std::string(optarg);
      if (!file_exists(filenames.sr_bam_filename)) {
        fprintf(stderr, "[Hypo::utils] Error: File Error: Short reads BAM file does not exist %s!\n",filenames.sr_bam_filename.c_str());
        exit(1);
      }
      cmd += (" -b " + std::string(optarg));
      is_bamsr = true;
      args++;
      break;

    case 'B':
      filenames.lr_bam_filename = std::string(optarg);
      if (!file_exists(filenames.lr_bam_filename)) {
        fprintf(stderr, "[Hypo::utils] Error: File Error: Long reads BAM file does not exist %s!\n",filenames.lr_bam_filename.c_str());
        exit(1);
      }
      cmd += (" -B " + std::string(optarg));
      is_bamlr = true;
      args++;
      break;

    case 'o':
      filenames.output_filename = std::string(optarg);
      cmd += (" -o " + std::string(optarg));
      args++;
      break;

    case 'k':
      kind = std::string(optarg);
      cmd += (" -k " + kind);
      args++;
      break;
    case 't':
      if (atoi(optarg) <= 0) {
        fprintf(stderr, "[Hypo::utils] Error: Arg Error: Number of threads (t) must be positive %d!\n",atoi(optarg));
        exit(1);
      }
      //flags.threads = std::max((UINT32)atoi(optarg)-1,(UINT32)1);
      flags.threads = std::max((UINT32)atoi(optarg),(UINT32)1);
      cmd += (" -t " + std::string(optarg));
      args++;
      break;
    case 'p':
      if (atoi(optarg) < 0) {
        fprintf(stderr, "[Hypo::utils] Error: Arg Error: Processing-size, i.e. number of contigs processed in a batch, (p) must NOT be negative %d!\n",atoi(optarg));
        exit(1);
      }
      flags.processing_batch_size = std::max((UINT32)atoi(optarg),(UINT32)0);
      cmd += (" -p " + std::string(optarg));
      args++;
      break;
    case 'i':
      flags.intermed=true;
      cmd += (" -i ");
      break;
    case 'w':
      flags.wdir = (std::string(optarg)+"/");
      if (!file_exists(flags.wdir)) {
        fprintf(stderr, "[Hypo::utils] Error: Directory Error: Working directory does not exist %s!\n",flags.wdir.c_str());
        exit(1);
      }
      cmd += (" -w " + std::string(optarg));
      break;
    default:
      usage();
      exit(0);
    }
  }
  // TODO: Check if the int args conform to the assumpions
  bool is_complete = false;
  std::string mode = "";
  if (is_draft && is_size) {
      if (is_bamlr && is_lcov) { // Long reads alignment
        is_complete = true;
        if (is_sr && is_scov) { // Long + short 
          if (is_bamsr) {
            flags.mode = Mode::LSA;
            mode = "LONG + SHORT (Using initial alignments)";
            flags.intermed = true;
          }
          else {
            flags.mode = Mode::LO;
            mode = "LONG ONLY";
          }
        }
        else { // Long read only          
          flags.mode = Mode::LO;
          mode = "LONG ONLY";
        }
    }
    else { // short reads only
      if (is_sr && is_scov && is_bamsr) {
        is_complete = true;
        flags.mode = Mode::SO;
        mode = "SHORT ONLY";
      }
    }
  }
  if (is_complete) {
    // Set output name
    if (filenames.output_filename == "") {
      size_t ind = filenames.draft_filename.find_last_of("(/\\");
      std::string dfullname (filenames.draft_filename,ind + 1);
      ind = dfullname.find_last_of(".");
      std::string dname(dfullname,0,ind);
      filenames.output_filename = "hypo_" + dname + ".fasta";
    }
    fprintf(stdout, "[Hypo::Utils] Info: Given Command: %s.\n",cmd.c_str()); 
    fprintf(stdout, "[Hypo::Utils] Info: Operating Mode: %s.\n",mode.c_str());
    fprintf(stdout, "[Hypo::Utils] Info: High and Low Coverage thresholds: %s.\n",cv.c_str()); 
    // Set stage to start from  
    create_dir(flags.wdir+AUX_DIR);
    create_dir(flags.wdir+SR_DIR);
    if (flags.intermed) {  
      if (file_exists(flags.wdir+STAGEFILE)) {
        std::ifstream ifs(flags.wdir+STAGEFILE);
        if (!ifs.is_open()) {
          fprintf(stderr, "[Hypo::Utils] Error: File open error: Stage File (%s) exists but could not be opened!\n",STAGEFILE);
          exit(1);
        }
        std::string dummy1,dummy2,dummy3;
        UINT stage_num=STAGE_BEG;
        while (ifs >> dummy1 >> dummy2 >> dummy3 >> stage_num){}
        flags.done_stage = stage_num;
      }
      else {
        flags.done_stage = STAGE_BEG;
      }
      fprintf(stdout, "[Hypo::Utils] Info: Intermediate Files will be stored.\n"); 
    }
    else {
      flags.done_stage = STAGE_BEG;
      fprintf(stdout, "[Hypo::Utils] Info: Intermediate Files will NOT be stored.\n"); 
    }
    fprintf(stdout, "[Hypo::Utils] Info: Beginning from stage: %u\n",flags.done_stage); 
  }
  else
  {
    fprintf(stderr, "[Hypo::] Error: Invalid command: Too few arguments!\n");
    usage();
    exit(1);
  }
  
  if (flags.mode==Mode::LSA || flags.mode==Mode::SO)
  {
    // Set window settings
    void set_kind(const std::string& kind);
    // Set expected short reads file size
    flags.sz_in_gb = get_expected_file_sz(given_sz,SWindow_cov_settings.mean_cov); 
  }
}

/*
   * Usage of the tool
   */
void usage(void)
{
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
            << "\tPrint the usage. \n\n";
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
    
  slog::Monitor monitor;
  /* Decode arguments */
  hypo::InputFlags flags;
  hypo::FileNames filenames;
  hypo::decodeFlags(argc, argv, flags, filenames);

  hypo::Mode mode = flags.mode; 
  hypo::Hypo hypo(flags);
  monitor.start();
  std::string tm="";
  if (mode==hypo::Mode::LSA) {
    // Call first round of polishing using short reads/aln
    std::string orignal_out_name = filenames.output_filename;
    std::string first_out_name = (orignal_out_name+".0.fasta");
    filenames.output_filename = first_out_name;
    
    if (flags.done_stage < STAGE_FIRST) {
      hypo.polish(mode, filenames);
      flags.done_stage = STAGE_FIRST;
    }
    
    tm = monitor.stop("[Hypo:main]: First round. ");
    fprintf(stdout, "[Hypo::main] //////////////////\n FIRST ROUND POLISHING DONE\n [Hypo::Hypo] ////////////////// \n%s\n", tm.c_str());

    // Map short reads to long read polished draft.
    // IT IS SAFE TO RUN SYSTEM AS USER-GIVEN ARGS ARE ALREADY SANITISED
    monitor.start();
    std::string sf_name (flags.wdir+SBAMFILE);
    std::string numth = std::to_string(flags.threads);

    std::string rn_file = "";
    if (filenames.sr_filenames.size()==2) { // R1 nd R2  
      fprintf(stdout, "[Hypo::main] Assuming Paired end reads with R1 as %s and R2 as %s \n", filenames.sr_filenames[0].c_str(), filenames.sr_filenames[1].c_str());
      rn_file = filenames.sr_filenames[0]+" "+filenames.sr_filenames[1];
    }
    else {
      rn_file = filenames.sr_filenames[0];
    }

    std::string aln_cmd = "minimap2 --MD -ax sr --secondary=no -t ";
    aln_cmd += (numth + " " + first_out_name) + " " + rn_file;

    std::string sam_cmd = ("samtools view -@ " + numth + " -Sb - | samtools sort -@" + numth + " -m 2g -o " + sf_name + " -");
    
    if (flags.done_stage < STAGE_REMAP) {
      std::string final_cmd = aln_cmd + " | " + sam_cmd;
      fprintf(stdout, "[Hypo::main] Following re-aligns the reads \n%s \n", final_cmd.c_str());
      auto status = system(final_cmd.c_str());
      if(!WIFEXITED(status)) {
          fprintf(stderr, "[Hypo::main] Error: Failed to run minimap2/samtools.\n");
          exit(1);
      }
    }
    if (!hypo::file_exists(sf_name)) {
        fprintf(stderr, "[Hypo::main] Error: File Error: (Intermediate) Short reads BAM file does not exist %s!\n",sf_name.c_str());
        exit(1);
    }
    tm = monitor.stop("[Hypo:main]: Mapping. ");
    fprintf(stdout, "[Hypo::main] //////////////////\n MAPPING SELECTIVE SHORT READS ONTO FIRST ROUND POLISHED DRAFT DONE\n [Hypo::Hypo] ////////////////// \n%s\n", tm.c_str());

    // Call second round of polishing using short reads only
    monitor.start();
    filenames.sr_bam_filename = sf_name;
    filenames.lr_bam_filename = "";
    filenames.draft_filename = first_out_name;
    filenames.output_filename = orignal_out_name;
    hypo.polish(hypo::Mode::SECOND,filenames);
    tm = monitor.stop("[Hypo:main]: Short round. ");
    fprintf(stdout, "[Hypo::Hypo] //////////////////\n FINAL POLISHING DONE\n [Hypo::Hypo] ////////////////// \n%s\n", tm.c_str());
  }
  else {
    hypo::Hypo hypo(flags);
    hypo.polish(mode,filenames);
  }
  monitor.total("[Hypo:main]: TOTAL (Combined). ");
  return 0;
}

