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

namespace hypo {
    int scaffold_main(hypo::Flags flags);

} // end namespace

int scaffold_main(hypo::Flags flags) {
    
    // initialize the objects that were initialized on the previous modules
    hypo::utils::load_contigs(flags.initial_assembly_contigs);
    hypo::utils::load_long_alignments(flags.long_initial_mapping_path);
    
    hypo::filter_alignments();
    
    hypo::scaffold(flags);
    
}

