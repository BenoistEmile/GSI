#pragma once

#include <iostream>
#include <string>
#include <vector>

class Protein {
private:
    const std::size_t id;
    const std::string sequence;
    std::vector<std::size_t> peptides;
public:
    Protein(const std::size_t id, const std::string sequence);
    ~Protein();

    friend std::ostream& operator<<(std::ostream& os, const Protein& protein);

    const std::size_t Get_Id() const;
    const std::string& Get_Sequence() const;
    const std::vector<std::size_t>& Get_Peptides() const;

    bool Is_Digested() const;

    void Add_Peptide(std::size_t peptide);
};