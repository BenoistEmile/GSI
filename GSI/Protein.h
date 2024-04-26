#pragma once

#include <iostream>
#include <string>
#include <vector>

class Protein {
private:
    const std::size_t id;
    const std::string sequence, accession;
    std::vector<std::size_t> peptides;
public:
    Protein(const std::size_t id, const std::string sequence);
    Protein(const std::size_t id, const std::string sequence, const std::string accession);
    ~Protein();

    friend std::ostream& operator<<(std::ostream& os, const Protein& protein);

    const std::size_t Get_Id() const;
    const std::string& Get_Sequence() const;
    const std::vector<std::size_t>& Get_Peptides() const;
    const std::string& Get_Accession() const;

    bool Is_Digested() const;

    void Add_Peptide(std::size_t peptide);
};