#include "Model.h"

//__________________________________________________________________________________________________________

void Solution::Print(const std::vector<Peptide*> peptides , std::vector<Spectrum*> spectra) const {
    std::cout << "\nSelected proteins :" << std::endl;
    for (auto& couple : abundances) {
        std::cout << "   - " << couple.first << " : " << couple.second << std::endl;
    }
    for (const Identification* identification : identifications) {
        std::cout << *identification << std::endl;
    }
    if (spectra.size() && spectra[0]->Is_Simulated()) {
        std::cout << "\nNumber of wrong selected edges : ";
        unsigned int compteur = 0;
        for (const Identification* identification : identifications) {
            if (spectra[identification->spectrum]->Get_Origin()->peptide != identification->peptide) {
                compteur++;
            }
        }
        std::cout << compteur << std::endl;
    }
    /*
    std::cout << "\nSelected edges :" << std::endl;
    for (const Identification* identification : identifications) {
        std::cout << "   - (p" << identification->peptide << " ,s" << identification->spectrum << ")" << std::endl;
    }
    */
}

void Solution::Save(std::vector<Spectrum*> spectra, const std::string file_name, bool overwrite, bool save_proteins, bool save_ident) const {
    std::ofstream output_file;
    std::filesystem::path file_path = std::filesystem::current_path() / "solution" / file_name;
    file_path += ".csv";
    if (!fileExists(file_path) || overwrite) {
        output_file.open(file_path);
        output_file << "protein_id,abundance" << std::endl;
        if (save_proteins) {
            output_file << "Selected proteins :" << std::endl;
            for (auto& couple : abundances) {
                output_file << couple.first << "," << couple.second << std::endl;
            }
        }
        if (save_ident) {
            for (const Identification* identification : identifications) {
                output_file << *identification << std::endl;
            }
        }
        if (spectra.size() && spectra[0]->Is_Simulated()) {
            output_file << "\nNumber of wrong selected edges : ";
            unsigned int compteur = 0;
            for (const Identification* identification : identifications) {
                if (spectra[identification->spectrum]->Get_Origin()->peptide != identification->peptide) {
                    compteur++;
                }
            }
            output_file << compteur << std::endl;
        }
        output_file.close();
        std::cout << "Saved solution to " << file_name << std::endl;
    }
    else {
        std::cout << "ERROR : There already is a file named : " << file_name << std::endl;
    }
}

void Solution::Save(std::vector<Spectrum*> spectra, const std::string file_name, const std::vector<Protein*>& proteins, bool overwrite, bool save_proteins, bool save_ident) const {
    std::ofstream output_file;
    std::filesystem::path file_path = std::filesystem::current_path() / "solution" / file_name;
    file_path += ".csv";
    if (!fileExists(file_path) || overwrite) {
        output_file.open(file_path);
        output_file << "protein_id,abundance" << std::endl;
        if (save_proteins) {
            output_file << "Selected proteins :" << std::endl;
            for (auto& couple : abundances) {
                output_file << proteins.at(couple.first)->Get_Accession() << "," << couple.second << std::endl;
            }
        }
        if (save_ident) {
            for (const Identification* identification : identifications) {
                output_file << *identification << std::endl;
            }
        }
        if (spectra.size() && spectra[0]->Is_Simulated()) {
            output_file << "\nNumber of wrong selected edges : ";
            unsigned int compteur = 0;
            for (const Identification* identification : identifications) {
                if (spectra[identification->spectrum]->Get_Origin()->peptide != identification->peptide) {
                    compteur++;
                }
            }
            output_file << compteur << std::endl;
        }
        output_file.close();
        std::cout << "Saved solution to " << file_name << std::endl;
    }
    else {
        std::cout << "ERROR : There already is a file named : " << file_name << std::endl;
    }
}

void Solution::Save(std::vector<Spectrum*> spectra, std::ofstream& output_file, bool save_proteins, bool save_ident) const {
    output_file << "SOLUTION :" << std::endl;
    output_file << "Selected " << abundances.size() << " proteins." << std::endl;
    if (save_proteins) {
        output_file << "Selected proteins :" << std::endl;
        for (auto& couple : abundances) {
            output_file << "   - " << couple.first << " : " << couple.second << std::endl;
        }
    }
    if (save_ident) {
        for (const Identification* identification : identifications) {
            output_file << *identification << std::endl;
        }
    }
    if (spectra.size() && spectra[0]->Is_Simulated()) {
        output_file << "\nNumber of wrong selected edges : ";
        unsigned int compteur = 0;
        for (const Identification* identification : identifications) {
            if (spectra[identification->spectrum]->Get_Origin()->peptide != identification->peptide) {
                compteur++;
            }
        }
        output_file << compteur << std::endl;
    }
    output_file << std::endl << std::endl;
}

void Solution::Save(std::vector<Spectrum*> spectra, std::ofstream& output_file, const std::vector<Protein*>& proteins, bool save_proteins, bool save_ident) const {
    output_file << "SOLUTION :" << std::endl;
    output_file << "Selected " << abundances.size() << " proteins." << std::endl;
    if (save_proteins) {
        output_file << "Selected proteins :" << std::endl;
        for (auto& couple : abundances) {
            output_file << "   - " << proteins.at(couple.first)->Get_Accession() << " : " << couple.second << std::endl;
        }
    }
    if (save_ident) {
        for (const Identification* identification : identifications) {
            output_file << *identification << std::endl;
        }
    }
    if (spectra.size() && spectra[0]->Is_Simulated()) {
        output_file << "\nNumber of wrong selected edges : ";
        unsigned int compteur = 0;
        for (const Identification* identification : identifications) {
            if (spectra[identification->spectrum]->Get_Origin()->peptide != identification->peptide) {
                compteur++;
            }
        }
        output_file << compteur << std::endl;
    }
    output_file << std::endl << std::endl;
}