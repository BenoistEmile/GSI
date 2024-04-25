#include "Spectrum.h"
#include "Fonction.h"

#include<set>
#include <iterator>
#include <algorithm>
#include <cmath>

//__________________________________________________________________________________________________________

Spectrum::Spectrum(const std::size_t id, const std::vector<Pic*>* pics) : id(id), pics(pics), origin(nullptr) ,is_copy(false) {
    if (pics->size() > 1) {
        std::size_t i = 0;
        for (std::size_t j = 1; j < pics->size(); ++j) {
            if ((*pics)[i]->mass > (*pics)[j]->mass) {
                throw std::string("Error : pics are not sorted.");
            }
            ++i;
        }
    }
}

Spectrum::Spectrum(const std::size_t id, const std::vector<Pic*>* pics, const Origin* origin ,const bool is_copy) : id(id), pics(pics), origin(origin) ,is_copy(is_copy) {
    if (pics->size() > 1) {
        std::size_t i = 0;
        for (std::size_t j = 1; j < pics->size(); ++j) {
            if ((*pics)[i]->mass > (*pics)[j]->mass) {
                throw std::string("Error : pics are not sorted.");
            }
            ++i;
        }
    }
}

Spectrum::~Spectrum() {
    if (!is_copy) {
        delete origin;
        for (std::size_t i = 0; i < (*pics).size(); ++i) {
            delete (*pics)[i];
        }
        delete pics;
    }
}

//__________________________________________________________________________________________________________

const std::size_t Spectrum::Get_Id() const {
    return id;
}

const std::vector<Pic*>* Spectrum::Get_Pics() const{
    return pics;
}

const Origin* Spectrum::Get_Origin() const {
    return origin;
}

const bool Spectrum::Is_Simulated() const {
    return origin != nullptr;
}

//__________________________________________________________________________________________________________

std::vector<Pic*>* Spectrum::Filter1(unsigned int minimum_number_of_masses, unsigned int maximum_number_of_masses, int accuracy, unsigned int number_of_copies) const {
    std::vector<std::size_t> sorted_pics;
    sorted_pics.reserve(pics->size());
    std::size_t i = 0;
    for (; i < pics->size(); ++i) {
        sorted_pics.push_back(i);
    }
    std::sort(sorted_pics.begin(), sorted_pics.end(), [this](std::size_t pic1, std::size_t pic2)->bool {return pics->at(pic1)->intensity > pics->at(pic2)->intensity; });
    if (Is_Simulated()){
        std::vector<Pic*>* result = new std::vector<Pic*>;
        if (pics->size() >= minimum_number_of_masses) {
            for (i = 0; (i < maximum_number_of_masses || maximum_number_of_masses == 0) && i < pics->size(); ++i) {
                result->push_back(new Pic(Round_Precision(pics->at(i)->mass, accuracy), pics->at(i)->intensity));
            }
        }
        return result;
    }
    else {
        std::set<std::size_t> selected;
        i = 0;
        std::size_t j;
        const double min_accuracy = pow(10, -accuracy);
        bool found;
        while (i < pics->size() && (selected.size() <= maximum_number_of_masses || maximum_number_of_masses == 0)) {
            found = false;
            j = sorted_pics[i] + 1;
            while (!found && j < pics->size() && pics->at(j)->mass - pics->at(sorted_pics[i])->mass <= 2.0 * number_of_copies * min_accuracy) {
                found = selected.contains(j);
                ++j;
            }
            j = sorted_pics[i];
            while (!found && j > 0 && pics->at(sorted_pics[i])->mass - pics->at(j - 1)->mass <= 2.0 * number_of_copies * min_accuracy) {
                found = selected.contains(j - 1);
                --j;
            }
            if (!found) {
                selected.insert(i);
            }
            ++i;
        }
        std::vector<Pic*>* result = new std::vector<Pic*>;
        if (selected.size() >= minimum_number_of_masses) {
            result->reserve(selected.size());
            for (std::size_t pic : selected) {
                result->push_back(new Pic(Round_Precision(pics->at(pic)->mass, accuracy), pics->at(pic)->intensity));
            }
        }
        return result;
    }
}

//__________________________________________________________________________________________________________

std::ostream& operator<<(std::ostream& os, const Spectrum& spectrum)
{
    os << "Spectrum " << spectrum.id;
    if (spectrum.origin != nullptr) {
        os << "(origin : protein " << spectrum.origin->protein << " ,peptide " << spectrum.origin->peptide << " ,position " << spectrum.origin->position << ")";
    }
    os << " :\n    pics (mass ,intensity) :";
    for (auto iter = spectrum.pics->begin(); iter != spectrum.pics->end(); ++iter) {
        os << "\n         (" << (*iter)->mass << " ," << (*iter)->intensity << ")";
    }
    return os;
}