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
/** Defines the class PackedSequence.
 * It represents a sequence with a DNA base with 4-bits letter.
 * Either 2-bits or 4-bits packing
 * 4-bits Valid letters: A,C, G, T,N
 * 2-bits valid letter
 */
#include "utils.hpp"

namespace hypo {    
    namespace utils {
        int load_contigs(hypo::Objects & objects, const std::string load_path) {
            gzFile fp = gzopen(load_path.c_str(), "r");
            kseq_t *seq;
            seq = kseq_init(fp);
            int l;
            while((l = kseq_read(seq)) >= 0) {
                std::string get_sequence = "";
                std::string get_id(seq->name.s, seq->name.l);
                
                for(int i = 0; i < seq->seq.l; i++) {
                    int c = seq->seq.s[i];
                    if(c % 4 == 0) get_sequence += "A";
                    if(c % 4 == 1) get_sequence += "C";
                    if(c % 4 == 2) get_sequence += "G";
                    if(c % 4 == 3) get_sequence += "T";
                }
                
                objects.contigs.push_back(get_sequence);
                objects.contig_name.push_back(get_id);
            }
            return 0;
        }
        
        int write_contigs(const hypo::Objects & objects, const std::string write_path) {
            std::ofstream ofile(write_path);
            if (!ofile.is_open()) {
                fprintf(stderr, "[Hypo::Hypo] Error: File open error: Output File (%s) could not be opened!\n",write_path.c_str());
                exit(1);
            }
            for(int i = 0; i < objects.contigs.size(); i++) {
                ofile << ">" << objects.contig_name[i] << "\n";
                ofile << objects.contigs[i] << "\n";
            }
            ofile.close();
            return 0;
        }
        
        int write_contigs_1(const hypo::Objects & objects, const std::string write_path) {
            std::ofstream ofile(write_path);
            if (!ofile.is_open()) {
                fprintf(stderr, "[Hypo::Hypo] Error: File open error: Output File (%s) could not be opened!\n",write_path.c_str());
                exit(1);
            }
            for(int i = 0; i < objects.contigs_1.size(); i++) {
                ofile << ">" << objects.contig_name_1[i] << "\n";
                ofile << objects.contigs_1[i] << "\n";
            }
            ofile.close();
            return 0;
        }
        
        int write_contigs_2(const hypo::Objects & objects, const std::string write_path) {
            std::ofstream ofile(write_path);
            if (!ofile.is_open()) {
                fprintf(stderr, "[Hypo::Hypo] Error: File open error: Output File (%s) could not be opened!\n",write_path.c_str());
                exit(1);
            }
            for(int i = 0; i < objects.contigs_2.size(); i++) {
                ofile << ">" << objects.contig_name_2[i] << "\n";
                ofile << objects.contigs_2[i] << "\n";
            }
            ofile.close();
            return 0;
        }
        
        int load_short_alignments(hypo::Objects & objects, const std::string load_path) {
            return 0;
        }
        
        int load_long_alignments(hypo::Objects & objects, const std::string load_path) {
            auto sam_file = sam_open(load_path.c_str(), "r");
            auto sam_header = sam_hdr_read(sam_file);
            auto current_align = bam_init1();
            
            while(sam_read1(sam_file, sam_header, current_align)>=0) {
                // ignore unmapped reads
                if (current_align->core.flag & (BAM_FUNMAP)) {
                    continue;
                }
                
                // ignore secondary and failed reads
                if (current_align->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) { 
                    continue;
                }
                
                // filter out quality value (long reads)
                if (current_align->core.qual < 2) {
                    continue;
                }
                
                std::string contig_name(sam_hdr_tid2name(sam_header, current_align->core.tid));
                
                // int contig_id = objects.find_contig(contig_name);
                // if(contig_id == -1) {
                //    fprintf(stderr, "[Hypo::Hypo] Error: Alignment File error: Contig-reference (%s) does not exist in the draft!\n",cname.c_str());
                //    exit(1);
                //}
                
                // objects.add_alignment(contig_id, Alignment(current_align));
            }
        }
        
        int initialize_solid_kmers(hypo::Objects & objects, const int k, const std::string load_path_1, const std::string load_path_2) {
            return 0;
        }
        
    }
} // namespace hypo
