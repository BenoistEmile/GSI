
#include "Fonction.h"

#include <cctype>
#include <math.h>

bool Is_Valid_Sequence(std::string sequence) {
    std::string::iterator iter = sequence.begin();
    if (iter == sequence.end()) {
        return false;
    }
    bool found = false;
    while (iter != sequence.end() && !found) {
        found = !isalpha(*iter) || *iter == 'X';
        ++iter;
    }
    return !found;
}

double Round_Precision(double value, unsigned int precision) {
    return round(value * pow(10, precision)) / (pow(10, precision));
}