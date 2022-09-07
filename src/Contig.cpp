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
/** Defines the class Contig.
 * It represents a contig and contains the functionality for different processing that will be carried out on a contig.
 */

#include <algorithm>
#include <numeric>
#include "Contig.hpp"
#include "MinimizerDeque.hpp"

namespace hypo
{
    Mode Contig::_mode = Mode::LSA;
    Contig::Contig(const UINT32 id, const std::string& name, const std::string seq): 
    _id(id), _name(name), _len(seq.size()), _pseq(seq), _solid_pos(seq.size(),0), _reg_pos(seq.size()+1,0), 
    _numSR(0), _lenSR(0), _num_wind(0) {_reg_pos[0]=1; _reg_pos[_len]=1;}
    
    Contig::Contig(const UINT32 id, const kseq_t * ks): 
    _id(id), _name(std::string(ks->name.s,ks->name.l)), _len(ks->seq.l), _pseq(ks) , _solid_pos(ks->seq.l,0), 
    _reg_pos(ks->seq.l+1,0), _numSR(0), _lenSR(0), _num_wind(0) {_reg_pos[0]=1; _reg_pos[_len]=1;}
    

    void Contig::find_solid_pos(const std::unique_ptr<suk::SolidKmers>& pSK) {
        const UINT16 cv = Sr_settings.cov_th+1;
        const UINT k = pSK->get_k();
        UINT64 kmer=0;
        UINT kmer_len=0;
        const UINT64 kmask = (1ULL<<2*k) - 1;
        size_t num_bases = _len;
        for (size_t i=0; i < num_bases; ++i) {
            BYTE b = _pseq.enc_base_at(i);
            if (b < 4) { //ACGT
                kmer = (kmer << 2 | b) & kmask;
                if (kmer_len<k) { ++kmer_len;}                
            }
            else { // N; reset
                kmer_len = 0;
                kmer=0;
            }
            if (kmer_len==k && pSK->is_solid(kmer)) {
                // check if starts/ends in HP
                bool should_add = true;
                if ((i < num_bases-1) && (_pseq.enc_base_at(i+1) == b) ) { // next base is same as this one
                    should_add = false;
                }
                auto beg_pos = i+1-k;
                if ((beg_pos > 0) && (_pseq.enc_base_at(beg_pos-1) == _pseq.enc_base_at(beg_pos)) ) { // pvs base is same as the first one ogf kmer
                    should_add = false;
                }
                if (should_add) {
                    _solid_pos[beg_pos]=1;
                    if (Contig::_mode==Mode::LSA) { // Artificially fill in coverage and support of solid kmers so that SO/SECOND code can be used to find SR
                        _kmerinfo.emplace_back(std::make_unique<KmerInfo>(kmer,cv,cv));
                    }
                    else {
                        _kmerinfo.emplace_back(std::make_unique<KmerInfo>(kmer));
                    }
                }                
            }
        }
        sdsl::util::init_support(_Rsolid_pos,&_solid_pos);
        sdsl::util::init_support(_Ssolid_pos,&_solid_pos);
    }

    void Contig::prepare_for_division(const UINT k, const std::string wdir) {        
        ///////////////////////////////////////////////////////
        /* Find SR */
        std::vector<UINT32> sr_pos;
        std::vector<UINT32> sr_len;
        
        auto num_sk = _kmerinfo.size();
        _anchor_kmers.reserve(2*num_sk);
        // 0th index has dummy
        _anchor_kmers.push_back(0);
        sr_pos.reserve(num_sk);
        sr_len.reserve(num_sk);
        std::vector<UINT32> kinds;
        std::vector<UINT32> sr_poss;
        UINT32 last_sr_pos = 0;
        bool in_sr = false;
        UINT i = 0;
        for (UINT pos=0; pos < _solid_pos.size(); ++pos) {
            if (_solid_pos[pos]==1) { // kmer
                // Check if valid                
                bool is_valid =false;
                if (_kmerinfo[i]->coverage >=Sr_settings.cov_th) {
                    UINT supp_th = UINT(Sr_settings.supp_frac * _kmerinfo[i]->coverage);
                    if  (_kmerinfo[i]->support >= (supp_th)) { // this kmer has >=80% support; supported by both haplotypes
                        is_valid = true;
                    }
                }		
                if (is_valid) { // valid will either be the first or ctd of sr                    
                    if (!in_sr) { // first kmer of sr
                        kinds.clear();
                        sr_poss.clear();
                        in_sr = true;
                    }
                    kinds.emplace_back(i);
                    sr_poss.emplace_back(pos);
                    last_sr_pos = pos+k; // 1 past last pos covered
                }
                ++i;
            }

            // Add if reached the end of sr; Add until the second last (effectively removes the last kmer)
            if (in_sr && pos==last_sr_pos) {
                auto first_sr_pos = sr_poss[0];
                if (Contig::_mode!=Mode::LSA || last_sr_pos>(first_sr_pos+MIN_SR_LEN_LSA)) {
                    sr_pos.emplace_back(first_sr_pos);
                    sr_len.emplace_back(last_sr_pos-first_sr_pos);
                    _anchor_kmers.emplace_back((_kmerinfo[kinds[0]]->kid));
                    _anchor_kmers.emplace_back((_kmerinfo[kinds[kinds.size()-1]]->kid));                        
                }
                in_sr=false;             
            }                                
        }
        if (in_sr) { // add the last SR (if not already added)
            auto first_sr_pos = sr_poss[0];
            if (Contig::_mode!=Mode::LSA || last_sr_pos>(first_sr_pos+MIN_SR_LEN_LSA)) {
                sr_pos.emplace_back(first_sr_pos);
                sr_len.emplace_back(last_sr_pos-first_sr_pos);
                _anchor_kmers.emplace_back((_kmerinfo[kinds[0]]->kid));
                _anchor_kmers.emplace_back((_kmerinfo[kinds[kinds.size()-1]]->kid));                         
            }
        }
        // Free solid_pos and supporting DSas they are not needed any longer
        _kmerinfo.clear();
        _kmerinfo.shrink_to_fit();
        sdsl::util::clear(_solid_pos);
        sdsl::util::clear(_Rsolid_pos);
        sdsl::util::clear(_Ssolid_pos);
        _anchor_kmers.shrink_to_fit();
    
        // Set the number and the length of SRs 
        _numSR = sr_pos.size();
        _lenSR = std::accumulate(sr_len.begin(), sr_len.end(), 0);
        ///////////////////////////////////////////////////////
        /* Divide contig into SR and MegaWindows */
        _is_win_even = !(_numSR > 0 && sr_pos[0]==0);    
        _reg_pos[0]=1;
        // dummy
        auto contig_len = _len;
        UINT32 dummy_sr_pos = UINT32(contig_len);
        sr_pos.emplace_back(dummy_sr_pos);
        _reg_pos[dummy_sr_pos]=1;

        // Handle window after each SR
        for (UINT32 ind=0; ind < _numSR; ++ind) {
            _reg_pos[sr_pos[ind]]=1;
            UINT32 mw_start = sr_pos[ind]+sr_len[ind];
            _reg_pos[mw_start]=1;           
        }
        //////////////////////////////////////////////////////
        /* Minimiser set-up */
        if (Contig::_mode==Mode::SECOND|| Contig::_mode==Mode::SO) {
            _minimserinfo.reserve(_numSR+1);  
            // handle 0th window
            UINT32 windex = 0;
            if (_is_win_even) {            
                UINT32 mw_len = sr_pos[0];
                _minimserinfo.emplace_back(std::make_unique<MWMinimiserInfo>());
                if (mw_len>Window_settings.ideal_swind_size) {
                    initialise_minimserinfo(_pseq.unpack(0,sr_pos[0]), 0);
                }
                ++windex;    
            }
            for (UINT32 ind=0; ind < _numSR; ++ind,++windex) {               
                _minimserinfo.emplace_back(std::make_unique<MWMinimiserInfo>());
                UINT32 mw_start = sr_pos[ind]+sr_len[ind];
                UINT32 mw_len = sr_pos[ind+1]-mw_start;
                if (mw_len>Window_settings.ideal_swind_size) {
                    initialise_minimserinfo(_pseq.unpack(mw_start,sr_pos[ind+1]), windex);
                }            
            }
            sdsl::util::init_support(_RMreg_pos,&_reg_pos);
            sdsl::util::init_support(_SMreg_pos,&_reg_pos);
        }
    }

    

    void Contig::divide_into_regions() {
        auto contig_len = _len;
        auto mode = Contig::_mode;
        UINT32 sr_rank = 1; // rank of the first SR in contigs (i.e. 1) + number of window-type so that 0th SR
        // process region at the beginning of the next one (works because dummy at the end)
        UINT32 reg_start = 0;
        UINT32 reg_ind = 0;
        UINT32 li_ind = 0;
        bool look_for_gap = ((mode==Mode::LSA || mode==Mode::LO) && _gap_li.size()>0) ;
        UINT32 gb = 0;
        UINT32 ge = 0;
        if (look_for_gap) {
            gb = std::get<0>(_gap_li[0]);
            ge = std::get<1>(_gap_li[0]);
        }
        const UINT32 cToo_large = 2*Window_settings.ideal_swind_size;
        for (size_t i=1; i < contig_len+1; ++i) {
            if (_reg_pos[i]==1) {
                UINT32 reg_end = i;
                // MW
                if ((_is_win_even && reg_ind%2==0) || (!_is_win_even && reg_ind%2==1)) {
                    if (mode==Mode::LO || mode==Mode::LSA) {
                        // This MW may have one or more LI
                        while (look_for_gap && gb>=reg_start && ge<=reg_end) {
                            if (gb>reg_start) {fixed_divide(reg_start,gb);}
                            //gap-li
                            std::cout << _name << " Invalidated LI "<< gb<<"-"<<ge<<std::endl;
                            _reg_pos[gb] =1; //it's redundant;                 
                            _reg_info.emplace_back(0);
                            _reg_type.emplace_back(RegionType::INVALID);
                            reg_start = ge;
                            ++li_ind;
                            if (li_ind < _gap_li.size()) {
                                gb = std::get<0>(_gap_li[li_ind]);
                                ge = std::get<1>(_gap_li[li_ind]);
                            }
                            else {look_for_gap=false;}
                        }                      
                        
                        // This function fills in the pos, info, type
                        fixed_divide(reg_start,reg_end);
                    }
                    else { // SO/SECOND
                        // next will be sr if this is not the last; pvs will be sr if this is not 0
                        char pvs = (reg_ind==0) ? ('n') : ('s');
                        char nxt = (i==contig_len) ? ('n') : ('s');
                        // This function fills in the pos, info, type
                        divide(reg_ind,reg_start,reg_end,pvs,nxt);
                    }
                }
                else { // SR
                    _reg_info.emplace_back(sr_rank);
                    ++sr_rank;
                    _reg_type.emplace_back(RegionType::SR);
                }
                ++reg_ind;
                reg_start = reg_end;
                /*
                if (_reg_type.size()!=sdsl::bit_vector::rank_1_type(&_reg_pos)(reg_end)) {
                    std::cout << "ERROR!!!! " << i << " " <<_reg_type.size() << " " << sdsl::bit_vector::rank_1_type(&_reg_pos)(reg_end) << std::endl;
                    exit(1);
                }
                */
            }
        }

        // Include dummy
        _reg_type.emplace_back(RegionType::SR);

        if (mode==Mode::LO|| mode==Mode::LSA) {
            auto num_reg = _reg_info.size()+1; // 0th is nullptr
            _lcovinfo.reserve(num_reg);
            _lcovinfo.emplace_back(std::unique_ptr<LCovInfo>());
            _plwindows.reserve(num_reg); 
            _plwindows.emplace_back(std::unique_ptr<Window>());
        }
        else {// SO/SECOND
            // Free memory
            _minimserinfo.clear();
            _minimserinfo.shrink_to_fit();
            sdsl::util::clear(_RMreg_pos);
            sdsl::util::clear(_SMreg_pos);

            // Window pointers
            _pswindows.reserve(_num_wind+1); // 0th is nullptr
            _pswindows.emplace_back(std::unique_ptr<Window>());
            _scovinfo.reserve(_num_wind+1);
            _scovinfo.emplace_back(std::unique_ptr<SCovInfo>());
        }
        

        // Initialise Rank-Select
        sdsl::util::init_support(_Rreg_pos,&_reg_pos);
        sdsl::util::init_support(_Sreg_pos,&_reg_pos);
        // Assign Window pointers
        for (size_t i=0; i < _reg_type.size();++i) {
            if (_reg_type[i] != RegionType::SR && _reg_type[i] != RegionType::MSR && _reg_type[i] != RegionType::INVALID) { // a valid window
                if (mode==Mode::LO || mode==Mode::LSA) {
                    _plwindows.emplace_back(std::make_unique<Window>(_pseq,_Sreg_pos(i+1),_Sreg_pos(i+2),WindowType::LONG));
                    _lcovinfo.emplace_back(std::make_unique<LCovInfo>());
                }
                else { //SO/SECOND                   
                    _pswindows.emplace_back(std::make_unique<Window>(_pseq,_Sreg_pos(i+1),_Sreg_pos(i+2),WindowType::SHORT));
                    _scovinfo.emplace_back(std::make_unique<SCovInfo>());
                }
            }
        }
        _reg_type.shrink_to_fit();
        _reg_info.shrink_to_fit();
    }

    
    // This will destroy alignments after use. 
    void Contig::prune_short_windows() {
        ///////////////////////////////////////////////////////
        // Free memeory
        _anchor_kmers.clear();
        _anchor_kmers.shrink_to_fit();

        ///////////////////////////////////////////////////////
        /* Prune Short Windows */
        UINT MIN_SHORT_NUM = Arms_settings.min_short_num;
        UINT MIN_CONTRIB = Arms_settings.min_contrib;
        double MIN_INTERNAL_CONTRIB = Arms_settings.min_internal_contrib;
        UINT NUM1 = Arms_settings.min_internal_num1;
        UINT NUM2 = Arms_settings.min_internal_num2;

        for (size_t i = 0; i < _reg_type.size(); ++i) {   
            if (Contig::_mode==Mode::SO || Contig::_mode==Mode::SECOND) { // prune all windows in SO/Second   ; redundant; will be calledonly while filling short
                if (_reg_type[i] !=  RegionType::SR && _reg_type[i] !=  RegionType::MSR) { // a short window  
                    UINT32 pi = _reg_info[i];
                    UINT64 internal_contrib = _pswindows[pi]->get_num_internal();
                    UINT64 contrib = _pswindows[pi]->get_num_total(); 
                    bool cond0 = internal_contrib > NUM1;
                    bool cond1 = contrib >=MIN_CONTRIB && (internal_contrib >= std::floor(MIN_INTERNAL_CONTRIB*contrib));
                    if (!_pswindows[pi]->has_low_qual_internal()) { // Use only internal arms wherever possible;   
                        bool cond2 = (_reg_type[i]==RegionType::SWS || _reg_type[i]==RegionType::SW || _reg_type[i]==RegionType::WS || _reg_type[i]==RegionType::MWS || _reg_type[i]==RegionType::SWM) && (internal_contrib >= NUM2);
                        if (cond0 || cond1 || cond2) { // get rid of pre/suf arms
                            _pswindows[pi]->set_internal_only();
                        }                        
                    }
                }
            }
        }        
    }

    std::ostream &operator<<(std::ostream &os, const Contig &ctg) {
        // Write name of the contig
        os << ">" << ctg._name << std::endl;
        // Write sequence
        auto num_reg = ctg._reg_type.size()-1; //excluding the dummy
        auto curr = ctg._Sreg_pos(1);
        for (UINT32 i=0; i < num_reg; ++i) {
            auto next = ctg._Sreg_pos(i+2);
            if (ctg._reg_type[i] == RegionType::SR || ctg._reg_type[i] == RegionType::MSR) { // an SR
                os << ctg._pseq.unpack(curr,next);
            }
            else if (ctg._reg_type[i] == RegionType::INVALID) { // window can be invalid only in LSA mode
                ; // Do nothing                
            }
            else if (ctg._reg_type[i] == RegionType::LONG || ctg._reg_type[i] == RegionType::PSEUDO) { // a long window
                auto pi = ctg._reg_info[i];
                os << (ctg._plwindows[pi])->get_consensus();
            }
            else if (ctg._reg_type[i] == RegionType::LI) { // a long ins window
                auto pi = ctg._reg_info[i];
                os << (ctg._pliwindows[pi])->get_consensus();
            }
            else { // a valid short window
                auto pi = ctg._reg_info[i];
                os << (ctg._pswindows[pi])->get_consensus();
            }
            curr = next;
        }
        os << std::endl;
        return os;
    }

    void Contig::generate_inspect_file(std::string wdir, std::ofstream& bedfile) {
        std::string suf = ".txt";
        if (Contig::_mode==Mode::LSA){suf=".0.txt";}
        std::string fn(wdir+INSPECTFILEPREF+_name+suf);
        std::ofstream ofs(fn);
        if (!ofs.is_open()) {
            fprintf(stderr, "[Hypo::Contig] Error: File open error: Inspection File (%s) could not be opened!\n",fn.c_str());
            exit(1);
        }
        auto num_reg = _reg_type.size()-1; //excluding the dummy
        // Write name of the contig
        ofs << ">" << _name << std::endl;
        // Write number of regions
        ofs << "#" << num_reg << std::endl;

        // lambda
        auto get_reg_type = [&](RegionType j)->std::string {
            std::string short_type = "";
            switch (j) {
                case RegionType::SWS: 
                    short_type = "SWS";
                    break;
                case RegionType::WS: 
                    short_type = "WS";
                    break;
                case RegionType::SW: 
                    short_type = "SW";
                    break;
                case RegionType::MWM: 
                    short_type = "MWM";
                    break;
                case RegionType::WM: 
                    short_type = "WM";
                    break;
                case RegionType::MW: 
                    short_type = "MW";
                    break;
                case RegionType::SWM: 
                    short_type = "SWM";
                    break;
                case RegionType::MWS: 
                    short_type = "MWS";
                    break;
                case RegionType::OTHER: 
                    short_type = "OTH";
                    break;
                case RegionType::INVALID: 
                    short_type = "INV";
                    break;
                case RegionType::NOPOL: 
                    short_type = "NPL";
                    break;
                case RegionType::LONG: 
                    short_type = "LNG";
                    break;
                 case RegionType::PSEUDO: 
                    short_type = "PSD";
                    break;
                case RegionType::LI: 
                    short_type = "LI";
                    break;
                case RegionType::SR: 
                    short_type = "SR";
                    break;
                case RegionType::MSR: 
                    short_type = "MSR";
                    break;
            };
            return short_type;
        };
        auto curr = _Sreg_pos(1);
        // Write regions        
        for (UINT32 i=0; i < num_reg; ++i) {
            auto next = _Sreg_pos(i+2);
            if (_reg_type[i] == RegionType::SR || _reg_type[i] == RegionType::MSR) { // an SR                
                ofs << "==========(" << curr << "-" <<  next-1 << ")\t";
                ofs << get_reg_type(_reg_type[i])<< "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << std::endl;
                std::string draft(_pseq.unpack(curr,next));
                ofs << "++\t" <<  draft <<std::endl;
                ofs << "++\t" <<  draft <<std::endl;
                bedfile << _name << "\t" << curr << "\t" << next << "\t" << get_reg_type(_reg_type[i]) << std::endl;
            }
            else if (_reg_type[i] == RegionType::INVALID) { // no-mapping short window (nopol can be only if long reads not used)
                    ofs << "==========(" << curr << "-" <<  next-1 << ")\t";
                    ofs << get_reg_type(_reg_type[i]) << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << std::endl;
                    std::string draft(_pseq.unpack(curr,next));
                    ofs << "++\t" << draft << std::endl;
                    ofs << "++\t" <<  std::endl;
                    bedfile << _name << "\t" << curr << "\t" <<curr+1 << "\t" <<  get_reg_type(_reg_type[i]) << std::endl;
            }
            else if (_reg_type[i] == RegionType::LONG || _reg_type[i] == RegionType::PSEUDO) { // a long window
                auto pi = _reg_info[i];
                ofs << "==========(" << curr << "-" <<  curr+_plwindows[pi]->get_window_len()-1 << ")\t";
                ofs << get_reg_type(_reg_type[i])<<"\t";
                ofs << *(_plwindows[pi]);
                bedfile << _name << "\t" << curr << "\t" << curr+1 << "\t" <<  get_reg_type(_reg_type[i]) << std::endl;
            }
            else if (_reg_type[i] == RegionType::LI) { // a long window
                auto pi = _reg_info[i];
                ofs << "==========(" << curr << "-" <<  curr+_pliwindows[pi]->get_window_len()-1 << ")\t";
                ofs << get_reg_type(_reg_type[i])<<"\t";
                ofs << *(_pliwindows[pi]);
                bedfile << _name << "\t" << curr << "\t" << curr+1 << "\t" <<  get_reg_type(_reg_type[i]) << std::endl;
            }
            else { // a valid short window
                auto pi = _reg_info[i];
                ofs << "==========(" << curr << "-" <<  curr+_pswindows[pi]->get_window_len()-1 << ")\t";
                ofs << get_reg_type(_reg_type[i])<<"\t";
                ofs << *(_pswindows[pi]);
                bedfile << _name << "\t" << curr << "\t" << curr+1 << "\t" <<  get_reg_type(_reg_type[i]) << std::endl;
            }
            curr = next;
        }
        ofs.close();
    }

    void Contig::initialise_minimserinfo(const std::string& draft_seq, UINT32 minfoind) {
    
        UINT32 last_found_position = draft_seq.size() + 1; //a unique identifier for 'first minimizer'
        UINT32 MINIMIZER_K = Minimizer_settings.k;
        UINT32 MINIMIZER_W = Minimizer_settings.w;
        //UINT32 shift = 2 * (Minimizer_settings.k - 1);
        UINT32 mask = (1ULL<<2*MINIMIZER_K) - 1;
        UINT32 kmer[2] = {0,0};
        MinimizerDeque<UINT32> minimizer_window(MINIMIZER_W + 1);
        UINT32 count_not_N = 0;
        UINT32 processed_kmer = 0;
        UINT32 current_start_position = 0;
        
        //we have to make a vector of found minimizers and a map to validity, to handle duplicates
        //is there a better solution?
        std::vector<UINT32> found_minimizers;
        std::vector<UINT32> found_minimizers_positions;
        std::unordered_map<UINT32, UINT8> counter;
        
        for(size_t i = 0; i < draft_seq.size(); ++i) {
            BYTE c = cNt4Table[(BYTE)draft_seq[i]];
            if(c < 4) {
                ++count_not_N;
                
                kmer[0] = (kmer[0] << 2ull | c) & mask;           // forward k-mer
                //kmer[1] = (kmer[1] >> 2ull) | (3ULL^c) << shift; // reverse k-mer
                //int z = kmer[0] < kmer[1] ? 0 : 1;
                int z=0;
                if(count_not_N >= MINIMIZER_K) {
                    while(!minimizer_window.empty() && std::get<0>(minimizer_window.back()) > kmer[z]) minimizer_window.pop_back();
                    minimizer_window.push_back(std::make_tuple(kmer[z], i));
                    while(std::get<1>(minimizer_window.front()) + MINIMIZER_W <= i) minimizer_window.pop_front();
                    ++processed_kmer;
                    if(processed_kmer >= MINIMIZER_W) {
                        current_start_position = std::get<1>(minimizer_window.front()) - MINIMIZER_K + 1;
                        if(current_start_position != last_found_position) { //first minimizer
                            found_minimizers_positions.push_back(current_start_position);
                            found_minimizers.push_back(std::get<0>(minimizer_window.front()));
                            ++counter[std::get<0>(minimizer_window.front())];
                        }
                        last_found_position = current_start_position;
                    }
                    
                }
            } else {
                count_not_N = 0;
            }
        }

        // initialise
        last_found_position = 0;
        _minimserinfo[minfoind]->minimisers.reserve(found_minimizers.size());
        _minimserinfo[minfoind]->rel_pos.reserve(found_minimizers.size());
        for(size_t i = 0; i < found_minimizers.size(); ++i) {
            if(counter[found_minimizers[i]] == 1) { //unique minimizer
                UINT32 p = found_minimizers_positions[i];
                //if (draft_seq[p]!=draft_seq[p+1] && draft_seq[p+MINIMIZER_K-1]!=draft_seq[p+MINIMIZER_K-2]) { // doesn't have HP at terminals
                //if (found_minimizers[i]!=Minimizer_settings.polyA && found_minimizers[i]!=Minimizer_settings.polyC && found_minimizers[i]!=Minimizer_settings.polyG && found_minimizers[i]!=Minimizer_settings.polyT) {
                if (!is_periodic(draft_seq.substr(p,MINIMIZER_K))){    
                    _minimserinfo[minfoind]->minimisers.push_back(found_minimizers[i]);
                    _minimserinfo[minfoind]->rel_pos.push_back(p - last_found_position);                    
                    last_found_position = found_minimizers_positions[i];
                }
            }
        }
        _minimserinfo[minfoind]->minimisers.shrink_to_fit();
        _minimserinfo[minfoind]->rel_pos.shrink_to_fit();
        auto num_minimisers = _minimserinfo[minfoind]->rel_pos.size();
        _minimserinfo[minfoind]->support.resize(num_minimisers,0);
        _minimserinfo[minfoind]->coverage.resize(num_minimisers,0);
    }
        
    void Contig::divide(const UINT32 reg_index, const UINT32 beg, const UINT32 end, char pvs, char nxt) { 
        const UINT32 IDEAL_WIND_SIZE =  Window_settings.ideal_swind_size;
        const UINT32 MINIMIZER_K = Minimizer_settings.k;
        const UINT32 cToo_large = 2*IDEAL_WIND_SIZE;     
        //// Collect supported minimisers 
        
        UINT32 minfoidx = (_is_win_even) ? (reg_index/2) : ((reg_index-1)/2);
        auto num_minimsers = _minimserinfo[minfoidx]->rel_pos.size();    
        UINT32 minimiser_pos = beg;
        std::vector<UINT32> supp_pos;
        supp_pos.reserve(num_minimsers);
        std::vector<UINT32> supp_minimisers;
        supp_minimisers.reserve(num_minimsers);

        for (UINT32 mi=0; mi< num_minimsers; ++mi) {
            minimiser_pos += _minimserinfo[minfoidx]->rel_pos[mi]; // absolute pos
            if (_minimserinfo[minfoidx]->coverage[mi]>=Minimizer_settings.cov_th) { 
                UINT supp_th = UINT(Minimizer_settings.supp_frac*_minimserinfo[minfoidx]->coverage[mi]);
                if (_minimserinfo[minfoidx]->support[mi]>=supp_th && (minimiser_pos+MINIMIZER_K<end)) { // supported minimiser found which is not adjacent to next sr
                    supp_pos.push_back(minimiser_pos);
                    supp_minimisers.push_back(_minimserinfo[minfoidx]->minimisers[mi]);
                } 
            }
        }
        ////// Find Minimser-based window cutting pos
        UINT remaining_size = end-beg;
        auto start = beg;
        std::vector<UINT32> cut_minimiser_index; // store those minimsers which will be used for cutting
        for (UINT32 mi=0; mi< supp_pos.size() && remaining_size>IDEAL_WIND_SIZE; ++mi) {
            // if this is the last minimiser or the next one makes the window go bigger than ideal, break at this one
            bool should_break = (mi==supp_pos.size()-1) ? (true) :(supp_pos[mi+1]>IDEAL_WIND_SIZE+start);
            if (should_break && supp_pos[mi]>start) { // break at this minimiser
                cut_minimiser_index.push_back(mi);
                start = supp_pos[mi]+MINIMIZER_K;
                remaining_size = end-start;
            }
        }

        /////////// Create windows
        UINT32 num_minimiser_windows=cut_minimiser_index.size();
        if (num_minimiser_windows==0) { // Only window
            if (end > beg + cToo_large) { // force break
                force_divide(beg,end,pvs,nxt);
            }
            else {
                ++_num_wind;
                _reg_pos[beg] =1; //it's redundant;                 
                _reg_info.emplace_back(_num_wind);
                if (pvs=='s' && nxt=='s') {_reg_type.emplace_back(RegionType::SWS);}
                else if (pvs=='s') {_reg_type.emplace_back(RegionType::SW);}
                else if (nxt=='s') {_reg_type.emplace_back(RegionType::WS);}
                else {_reg_type.emplace_back(RegionType::OTHER);}
            }
        }
        else { // at least one minimser
            // First window
            UINT32 win_end = supp_pos[cut_minimiser_index[0]];
            if (win_end > beg + cToo_large) { // force break
                force_divide(beg,win_end,pvs,'m');
            }
            else {
                ++_num_wind;
                _reg_pos[beg] =1; //it's redundant; 
                _reg_info.emplace_back(_num_wind);
                if (pvs=='s') {_reg_type.emplace_back(RegionType::SWM);}
                else {_reg_type.emplace_back(RegionType::WM);}
            }
            // internal windows (pvs window set at the beginning of a minimser)
            for (UINT32 cmi=1; cmi< num_minimiser_windows; ++cmi) {
                UINT32 pvs_mi = cut_minimiser_index[cmi-1];
                // Mark minimser window;
                _reg_pos[supp_pos[pvs_mi]] = 1;
                _reg_info.emplace_back(supp_minimisers[pvs_mi]);
                _reg_type.emplace_back(RegionType::MSR);
                // Window following the MSR
                UINT32 win_start = supp_pos[pvs_mi]+MINIMIZER_K;
                UINT32 win_end = supp_pos[cut_minimiser_index[cmi]];
                if (win_end > cToo_large+win_start) { // this window is too large
                    force_divide(win_start,win_end,'m','m');
                }
                else {
                    ++_num_wind;
                    _reg_pos[win_start] = 1;
                    _reg_info.emplace_back(_num_wind);
                    _reg_type.emplace_back(RegionType::MWM);
                }
            }
            // Last window
            UINT32 pvs_mi = cut_minimiser_index[num_minimiser_windows-1];
            // Mark minimser window;
            _reg_pos[supp_pos[pvs_mi]] = 1;
            _reg_info.emplace_back(supp_minimisers[pvs_mi]);
            _reg_type.emplace_back(RegionType::MSR);
            // Window following the MSR
            UINT32 win_start = supp_pos[pvs_mi]+MINIMIZER_K;
            if (end > cToo_large+win_start) { // this window is too large
                force_divide(win_start,end,'m',nxt);
            }
            else {
                ++_num_wind;
                _reg_pos[start] = 1;
                _reg_info.emplace_back(_num_wind);
                if (nxt=='s') {_reg_type.emplace_back(RegionType::MWS);}
                else {_reg_type.emplace_back(RegionType::MW);}
            }
        }
    }
    
    void Contig::force_divide(const UINT32 beg, const UINT32 end, char pvs, char nxt) {
        const UINT32 IDEAL_WIND_SIZE =  Window_settings.ideal_swind_size;
        const UINT32 cToo_large = 2*IDEAL_WIND_SIZE;
        // Find cutting pos
        int start = beg;
        int remaining_size = end-start;
        std::vector<UINT32> cut_pos;
        while (remaining_size >  IDEAL_WIND_SIZE) {
            UINT32 search_ind = start + Window_settings.wind_size_search_th;
            
            /*  # Breaking points(bp) of a window conceptually is the first and last bases of a window. Either should not fall into a HP.
                # => bp is a base that may surround a HP but can not be a part of it.
                # search_ind (end_bp) will be the index of the last base of this window. 
                # It should not be part of HP (Ensured by cond1 and 2 below)
                # cond1: end_bp and its pvs base are diffrent; cond2: end_bp and its next base (beg_bp for next wondow) are different
                # cond3: beg_bp of next win and its next base should be different (cond 2 and 3 ensure next wind start-bp surrounds a HP).
                // See slides for justification.
                // In short, this is how breaking should be done:  ----AAAB || CDDDD---- where A!=B and B!=C and C!=D
            */ 
            while (search_ind < end) {
                BYTE base = _pseq.enc_base_at(search_ind);
                if (base == _pseq.enc_base_at(search_ind-1)) { // cond1: ends in hp
                    search_ind += 1;
                }
                else if (search_ind+1 < end && base == _pseq.enc_base_at(search_ind+1)){// cond2: ends in hp
                    search_ind += 2; // next one will fail cond 1; thus increase by 2
                }                 
                else if (search_ind+2<end && _pseq.enc_base_at(search_ind+2) == _pseq.enc_base_at(search_ind+1)){  // see cond3
                    search_ind += 3; // next one will fail cond 2; next to next one will fail cond 1; thus increase by 3
                }
                else { //(all cond pass=> good bp)
                    break;
                }
                if (search_ind > start+cToo_large) { // Window has become too large
                    break;
                }
            }
            if (search_ind < end) {
                cut_pos.emplace_back(start);
                start = search_ind+1;
                remaining_size = end-start;
            }
            else { //exhausted
                break;
            }
        }
        if (start < end) { // Add last
            cut_pos.emplace_back(start);
        }
        
        /////// Create windows
        UINT32 num_windows=cut_pos.size();
        if (num_windows==1) { // only window
            ++_num_wind;
            _reg_pos[beg] =1; //it's redundant; 
            _reg_info.emplace_back(_num_wind);
            if (pvs=='s' && nxt=='s') {_reg_type.emplace_back(RegionType::SWS);}
            else if (pvs=='s' && nxt=='m') {_reg_type.emplace_back(RegionType::SWM);}
            else if (pvs=='s' && nxt=='n') {_reg_type.emplace_back(RegionType::SW);}
            else if (pvs=='m' && nxt=='s') {_reg_type.emplace_back(RegionType::MWS);}
            else if (pvs=='m' && nxt=='m') {_reg_type.emplace_back(RegionType::MWM);}
            else if (pvs=='m' && nxt=='n') {_reg_type.emplace_back(RegionType::MW);}
            else if (pvs=='n' && nxt=='s') {_reg_type.emplace_back(RegionType::WS);}
            else if (nxt=='n' && nxt=='m') {_reg_type.emplace_back(RegionType::WM);}
            else {_reg_type.emplace_back(RegionType::OTHER);}
        }
        else { // at least 2 windows
            // First window
            ++_num_wind;
            _reg_pos[beg] =1; //it's redundant; 
            _reg_info.emplace_back(_num_wind);
            if (pvs=='s') {_reg_type.emplace_back(RegionType::SW);}
            else if (pvs=='m') {_reg_type.emplace_back(RegionType::MW);}
            else {_reg_type.emplace_back(RegionType::OTHER);}

            // Internal windows
            for (UINT32 i=1; i<num_windows-1; ++i) {
                ++_num_wind;
                _reg_pos[cut_pos[i]] =1; 
                _reg_info.emplace_back(_num_wind);
                _reg_type.emplace_back(RegionType::OTHER);
            }
            // Last
            ++_num_wind;
            _reg_pos[cut_pos[num_windows-1]] =1; 
            _reg_info.emplace_back(_num_wind);
            if (nxt=='s') {_reg_type.emplace_back(RegionType::WS);}
            else if (nxt=='m') {_reg_type.emplace_back(RegionType::WM);}
            else {_reg_type.emplace_back(RegionType::OTHER);}
        }
    }  
    // Does nothing if beg and end are the same
    void Contig::fixed_divide(const UINT32 beg, const UINT32 end) {
        // Find cutting pos
        UINT start = beg;
        UINT remaining_size = end-start;
        std::vector<UINT32> cut_pos;
        auto offset = Window_settings.ideal_lwind_size-1;
        while (remaining_size >  Window_settings.ideal_lwind_size) {
            UINT32 search_ind = start + offset;
            if (search_ind < end) {
                cut_pos.emplace_back(start);
                start = search_ind+1;
                remaining_size = end-start;
            }
            else { //exhausted
                break;
            }
        }
        if (start < end) { // Add last
            cut_pos.emplace_back(start);
        }
        
        /////// Create windows
        UINT32 num_windows=cut_pos.size();
        // Internal windows
        for (UINT32 i=0; i<num_windows; ++i) {
            ++_num_wind;
            _reg_pos[cut_pos[i]] =1; 
            _reg_info.emplace_back(_num_wind);
            _reg_type.emplace_back(RegionType::LONG);
        }        
    }  

    // Called only in LSA
    void Contig::handle_ovl_li() {
        INT64 li_ind = -1; // index ofthe window that will be labelled Li and will store the result
        auto get_seq= [&](UINT32 i, UINT32 curr, UINT32 next) -> std::string {
            std::string sq="";
             if (_reg_type[i] == RegionType::SR) { // an SR
                sq =_pseq.unpack(curr,next);
            }
            else if (_reg_type[i] == RegionType::INVALID) { // window can be invalid only in LSA mode
                ; // Do nothing                
            }
            else if (_reg_type[i] == RegionType::LONG || _reg_type[i] == RegionType::PSEUDO) { // a long window
                auto pi = _reg_info[i];
                sq = (_plwindows[pi])->get_consensus();
                if (li_ind==-1) {li_ind=i;}
            }
            else { // a valid short window
                auto pi = _reg_info[i];
                sq =(_pswindows[pi])->get_consensus();
                if (li_ind==-1) {li_ind=i;}
            }
            return sq;
        };
        for (auto ov: _ovl_li) {
            auto ob = std::get<0>(ov);
            auto ol = std::get<1>(ov);
            std::cout << _name <<" Rectified LI "<< ob << ":"<<ol<<std::endl;
            UINT32 lft = (ob > ol) ? (ob-ol) : (0);
            UINT32 rt = (ob + ol < _len) ? (ob+ol) : (_len);
            auto lind = _Rreg_pos(lft);
            if (_reg_pos[lft]==0) { // starts before
                --lind;
            }
            // rind will point to 1 past the last region falling in the li (i.e. starting pos of the next reg)
            auto rind = _Rreg_pos(rt);
            auto mind = _Rreg_pos(ob);
            if (_reg_pos[ob]==0) { // starts before
                --mind;
            }
            std::string left_s="";
            std::string rt_s="";
            auto curr = _Sreg_pos(lind+1);
            for (UINT32 i=lind; i <= mind; ++i) {                
                auto next = _Sreg_pos(i+2);
                left_s+=get_seq(i,curr,next);
                curr = next;
            }
            curr = _Sreg_pos(mind+1);
            for (UINT32 i=mind; i < rind; ++i) {                
                auto next = _Sreg_pos(i+2);
                rt_s+=get_seq(i,curr,next);
                curr = next;
            }
            // Reset to INV, except li_ind which gets set to li
            for (UINT32 i=lind; i < rind; ++i) {_reg_type[i] = RegionType::INVALID; _reg_info[i] = 0;}
            _reg_type[li_ind] = RegionType::LI;
            _reg_info[li_ind] = _pliwindows.size();
            // draft here will not be used; for inspect file, better touse the original draft
            _pliwindows.emplace_back(std::make_unique<Window>(_pseq,_Sreg_pos(li_ind+1),_Sreg_pos(li_ind+2),WindowType::LONG));

            // Reset consensus of li_ind
            auto pi = _reg_info[li_ind];
            (_pliwindows[pi])->reset_consensus(left_s,rt_s);
        }
    } 

} // namespace hypo

