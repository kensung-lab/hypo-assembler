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
#include "Polish.hpp"
#include "utils.hpp"

namespace hypo {
    int polish_main(hypo::Objects & objects, hypo::InputFlags & flags) {
        
        // initialize the objects that were initialized on the previous modules
        hypo::utils::load_contigs(objects, flags.initial_assembly_path);
        hypo::utils::load_short_alignments(objects, flags.short_initial_mapping_path);
        hypo::utils::load_long_alignments(objects, flags.long_initial_mapping_path);
        
        hypo::polish(objects, flags);
        
        return 0;
    }
}
