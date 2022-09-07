#include <string>
#include <iostream>
#include "parser.hpp"

int main() {
    std::string a = "AAACGT";
    std::string b = "AAAGGT";
    std::string cigar;
    
    SWParser::Aligner x(10, 3, 4, 7);
    
    std::cout << x.align(a, b, cigar) << std::endl;
    std::cout << cigar << std::endl;
}
