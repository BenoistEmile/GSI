#pragma once

#include <iostream>
#include <string>
#include <vector>

class Protein {
private:
    const std::size_t id;
    const std::string sequence, accession;
    std::vector<std::size_t> peptides;
    std::vector<std::size_t> deactivated_peptides;
public:
    Protein(const std::size_t id, const std::string sequence);
    Protein(const std::size_t id, const std::string sequence, const std::string accession);
    ~Protein();

    friend std::ostream& operator<<(std::ostream& os, const Protein& protein);

    const std::size_t Get_Id() const;
    const std::string& Get_Sequence() const;
    const std::vector<std::size_t>& Get_Peptides() const;
    const std::string& Get_Accession() const;
    const std::size_t Get_Removed_Edges() const;
    const bool Get_Peptide_Activation(std::size_t peptide) const;

    bool Is_Digested() const;

    void Add_Peptide(std::size_t peptide);
    void Remove_Peptide(std::size_t peptide);
};