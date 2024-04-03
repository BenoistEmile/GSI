#pragma once

#include <iostream>
#include <vector>

//__________________________________________________________________________________________________________

struct Pic {
    const double mass;
    const double intensity;

    Pic(double mass, double intensity) : mass(mass), intensity(intensity) {}

    friend std::ostream& operator<<(std::ostream& os, const Pic& pic) {
        os << "(" << pic.mass << " ," << pic.intensity << ")";
        return os;
    };
};

struct Origin {
    const std::size_t peptide;
    const std::size_t protein;
    const std::size_t position;

    Origin(const Origin& origin) : peptide(origin.peptide), protein(origin.protein), position(origin.position) {};
    Origin(const std::size_t peptide, const std::size_t protein, const std::size_t position) : peptide(peptide), protein(protein), position(position)  {}
};

class Spectrum {
private:
    const std::size_t id;
    const std::vector<Pic*>* pics;
    const Origin* origin;
    const bool is_copy;
public:
    Spectrum(const std::size_t id, const std::vector<Pic*>* pics);
    Spectrum(const std::size_t id, const std::vector<Pic*>* pics, const Origin* origin, const bool is_copy);
    ~Spectrum();

    const std::size_t Get_Id() const;
    const std::vector<Pic*>* Get_Pics() const;
    const Origin* Get_Origin() const;

    const bool Is_Simulated() const;
    std::vector<Pic*>* Filter1(unsigned int minimum_number_of_masses, unsigned int maximum_number_of_masses, int accuracy, unsigned int number_of_copies) const;

    friend std::ostream& operator<<(std::ostream& os, const Spectrum& spectrum);
};

