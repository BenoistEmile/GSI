#include "Model.h"

//__________________________________________________________________________________________________________

static std::size_t Next_Cut(const std::string& sequence, std::size_t start) {
    std::size_t position_R = sequence.find_first_of('R', start);
    std::size_t position_K = sequence.find_first_of('K', start);
    std::size_t position = std::min(position_R, position_K);
    if (position == std::string::npos) {
        return sequence.size() - 1;
    }
    else if (position == sequence.size() - 1 || sequence.at(position + 1) != 'P') {
        return position;
    }
    else {
        return Next_Cut(sequence, position + 1);
    }
}

//__________________________________________________________________________________________________________

void Model::In_Silico_Digestion(int minimum_number_of_amino_acids ,int maximum_number_of_amino_acids) {
    std::unordered_map<std::string, std::size_t>::const_iterator position;
    std::size_t cut;
    std::size_t previous_cut;
    std::vector<std::string> sequences;
    std::string sequence;
    for (std::size_t i = 0; i < proteins.size(); ++i) {
        if (!(*proteins[i]).Is_Digested()) {
            sequences = {};
            sequence = (*proteins[i]).Get_Sequence();
            cut = 0;
            previous_cut = 0;
            while (cut < sequence.size()) {
                cut = Next_Cut(sequence, cut) + 1;
                sequences.push_back(sequence.substr(previous_cut, cut - previous_cut));
                previous_cut = cut;
            }
            for (std::size_t j = 0; j < sequences.size(); ++j) {
                if (sequences[j].size() >= minimum_number_of_amino_acids && (maximum_number_of_amino_acids == 0 || (sequences[j].size() <= maximum_number_of_amino_acids))) {
                    position = peptides_sequences.find(sequences[j]);
                    if (position == peptides_sequences.end()) {
                        peptides_sequences[sequences[j]] = peptides.size();
                        (*proteins[i]).Add_Peptide(peptides.size());
                        peptides.push_back(new Peptide(peptides.size(), sequences[j]));
                        peptides.back()->Add_Protein(i);
                    }
                    else {
                        (*proteins[i]).Add_Peptide(position->second);
                        peptides[position->second]->Add_Protein(i);
                    }
                }
            }
        }
    }
}

void Model::In_Silico_Digestion(std::vector<std::string>(digest_function)(const std::string& sequence) ,int minimum_number_of_amino_acids, int maximum_number_of_amino_acids) {
    std::unordered_map<std::string, std::size_t>::const_iterator position;
    std::vector<std::string> sequences;
    for (unsigned int i = 0; i < proteins.size(); ++i) {
        if (!(*proteins[i]).Is_Digested()) {
            sequences = digest_function((*proteins[i]).Get_Sequence());
            for (std::size_t j = 0; j < sequences.size(); ++j) {
                if (sequences[j].size() >= minimum_number_of_amino_acids && (maximum_number_of_amino_acids == 0 || (sequences[j].size() <= maximum_number_of_amino_acids))) {
                    position = peptides_sequences.find(sequences[j]);
                    if (position == peptides_sequences.end()) {
                        peptides_sequences[sequences[j]] = peptides.size();
                        (*proteins[i]).Add_Peptide(peptides.size());
                        peptides.push_back(new Peptide(peptides.size(), sequences[j]));
                        peptides.back()->Add_Protein(i);
                    }
                    else {
                        (*proteins[i]).Add_Peptide(position->second);
                        peptides[position->second]->Add_Protein(i);
                    }
                }
            }
        }
    }
}

void Model::In_Silico_Digestion(std::string file_name, int minimum_number_of_amino_acids ,int maximum_number_of_amino_acids) {
    std::unordered_map<std::string, std::size_t>::const_iterator position;
    std::size_t cut;
    std::size_t previous_cut;
    std::vector<std::string> sequences;
    std::string sequence;
    std::vector<std::string> n_terms;
    std::vector<std::string> c_terms;
    std::string n_term;
    std::string c_term;
    std::ofstream output_file;
    output_file.open(std::filesystem::current_path().generic_string() + "/data/digestion/" + file_name + ".csv");
    output_file << "peptide,nterm,cterm,miss1,miss2,label,accession,protein_id,peptide_id" << std::endl;
    for (std::size_t i = 0; i < proteins.size(); ++i) {
        if (!(*proteins[i]).Is_Digested()) {
            sequences = {};
            n_terms = {};
            c_terms = {};
            sequence = (*proteins[i]).Get_Sequence();
            cut = 0;
            previous_cut = 0;
            while (cut < sequence.size()) {
                cut = Next_Cut(sequence, cut) + 1;
                sequences.push_back(sequence.substr(previous_cut, cut - previous_cut));
                if (previous_cut < 8) {
                    n_term = std::string(8 - previous_cut, 'Z') + sequence.substr(0, 7 + previous_cut);
                }
                else {
                    n_term = sequence.substr(previous_cut - 8, 15);
                }
                if (n_term.length() < 15) {
                    n_term += std::string(15 - n_term.length(), 'Z');
                }
                n_terms.push_back(n_term);
                if (cut < 8) {
                    c_term = std::string(7 - cut, 'Z') + sequence.substr(cut, 8 + cut);
                }
                else {
                    c_term = sequence.substr(cut - 8, 15);
                }
                if (c_term.length() < 15) {
                    c_term += std::string(15 - c_term.length(), 'Z');
                }
                c_terms.push_back(c_term);
                previous_cut = cut;
            }
            for (std::size_t j = 0; j < sequences.size(); ++j) {
                if (sequences[j].size() >= minimum_number_of_amino_acids && (maximum_number_of_amino_acids == 0 || (sequences[j].size() <= maximum_number_of_amino_acids))) {
                    if (sequences[j].find("U") != std::string::npos || n_terms[j].find("U") != std::string::npos || c_terms[j].find("U") != std::string::npos) {
                        continue;
                    }
                    if (sequences[j].find("B") != std::string::npos || n_terms[j].find("B") != std::string::npos || c_terms[j].find("B") != std::string::npos) {
                        continue;
                    }
                    if (sequences[j].find("Z") != std::string::npos || n_terms[j].find("Z") != std::string::npos || c_terms[j].find("Z") != std::string::npos) {
                        continue;
                    }
                    position = peptides_sequences.find(sequences[j]);
                    if (position == peptides_sequences.end()) {
                        peptides_sequences[sequences[j]] = peptides.size();
                        (*proteins[i]).Add_Peptide(peptides.size());
                        peptides.push_back(new Peptide(peptides.size(), sequences[j]));
                        peptides.back()->Add_Protein(i);
                        output_file << sequences[j] << ',' << n_terms[j] << ',' << c_terms[j] << ",ZZZZZZZZZZZZZZZ,ZZZZZZZZZZZZZZZ,NA," << proteins.at(i)->Get_Accession() << "," << i << "," << peptides.back()->Get_Id() << std::endl;
                    }
                    else {
                        (*proteins[i]).Add_Peptide(position->second);
                        peptides[position->second]->Add_Protein(i);
                        output_file << sequences[j] << ',' << n_terms[j] << ',' << c_terms[j] << ",ZZZZZZZZZZZZZZZ,ZZZZZZZZZZZZZZZ,NA," << proteins.at(i)->Get_Accession() << "," << i << "," << peptides[position->second]->Get_Id() << std::endl;
                    }
                }
            }
        }
    }
    output_file.close();
}

void Model::In_Silico_Digestion(std::string file_name, std::ofstream& log_file, int minimum_number_of_amino_acids ,int maximum_number_of_amino_acids) {
    std::unordered_map<std::string, std::size_t>::const_iterator position;
    std::size_t cut;
    std::size_t previous_cut;
    std::vector<std::string> sequences;
    std::string sequence;
    std::vector<std::string> n_terms;
    std::vector<std::string> c_terms;
    std::string n_term;
    std::string c_term;
    std::ofstream output_file;
    output_file.open(std::filesystem::current_path().generic_string() + "/data/digestion/" + file_name + ".csv");
    output_file << "peptide,nterm,cterm,miss1,miss2,label,accession,protein_id,peptide_id" << std::endl;
    unsigned int count = 0;
    for (std::size_t i = 0; i < proteins.size(); ++i) {
        if (!(*proteins[i]).Is_Digested()) {
            sequences = {};
            n_terms = {};
            c_terms = {};
            sequence = (*proteins[i]).Get_Sequence();
            cut = 0;
            previous_cut = 0;
            while (cut < sequence.size()) {
                cut = Next_Cut(sequence, cut) + 1;
                sequences.push_back(sequence.substr(previous_cut, cut - previous_cut));
                if (previous_cut < 8) {
                    n_term = std::string(8 - previous_cut, 'Z') + sequence.substr(0, 7 + previous_cut);
                }
                else {
                    n_term = sequence.substr(previous_cut - 8, 15);
                }
                if (n_term.length() < 15) {
                    n_term += std::string(15 - n_term.length(), 'Z');
                }
                n_terms.push_back(n_term);
                if (cut < 8) {
                    c_term = std::string(7 - cut, 'Z') + sequence.substr(cut, 8 + cut);
                }
                else {
                    c_term = sequence.substr(cut - 8, 15);
                }
                if (c_term.length() < 15) {
                    c_term += std::string(15 - c_term.length(), 'Z');
                }
                c_terms.push_back(c_term);
                previous_cut = cut;
            }
            for (std::size_t j = 0; j < sequences.size(); ++j) {
                if (sequences[j].size() >= minimum_number_of_amino_acids && (maximum_number_of_amino_acids == 0 || (sequences[j].size() <= maximum_number_of_amino_acids))) {
                    if (sequences[j].find("U") != std::string::npos || n_terms[j].find("U") != std::string::npos || c_terms[j].find("U") != std::string::npos) {
                        continue;
                    }
                    if (sequences[j].find("B") != std::string::npos || n_terms[j].find("B") != std::string::npos || c_terms[j].find("B") != std::string::npos) {
                        continue;
                    }
                    if (sequences[j].find("Z") != std::string::npos || n_terms[j].find("Z") != std::string::npos || c_terms[j].find("Z") != std::string::npos) {
                        continue;
                    }
                    count++;
                    position = peptides_sequences.find(sequences[j]);
                    if (position == peptides_sequences.end()) {
                        peptides_sequences[sequences[j]] = peptides.size();
                        (*proteins[i]).Add_Peptide(peptides.size());
                        peptides.push_back(new Peptide(peptides.size(), sequences[j]));
                        peptides.back()->Add_Protein(i);
                        output_file << sequences[j] << ',' << n_terms[j] << ',' << c_terms[j] << ",ZZZZZZZZZZZZZZZ,ZZZZZZZZZZZZZZZ,NA," << proteins.at(i)->Get_Accession() << "," << i << "," << peptides.back()->Get_Id() << std::endl;
                    }
                    else {
                        (*proteins[i]).Add_Peptide(position->second);
                        peptides[position->second]->Add_Protein(i);
                        output_file << sequences[j] << ',' << n_terms[j] << ',' << c_terms[j] << ",ZZZZZZZZZZZZZZZ,ZZZZZZZZZZZZZZZ,NA," << proteins.at(i)->Get_Accession() << "," << i << "," << peptides[position->second]->Get_Id() << std::endl;
                    }
                }
            }
        }
    }
    log_file << "Digested " << proteins.size() << " proteins, obtained " << peptides.size() << " individual peptides, " << count << " in total" << std::endl << std::endl;
    output_file.close();
}

void Model::In_Silico_Digestion_2(std:: string file_name, int minimum_number_of_amino_acids, int maximum_number_of_amino_acids) {
    std::unordered_map<std::string, std::size_t>::const_iterator position;
    std::size_t cut;
    std::size_t previous_cut;
    std::vector<std::string> sequences;
    std::string sequence;
    std::vector<std::string> n_terms;
    std::vector<std::string> c_terms;
    std::string n_term;
    std::string c_term;
    std::ofstream output_file;
    std::size_t pep_count = 0;
    output_file.open(std::filesystem::current_path().generic_string() + "/data/digestion/" + file_name + ".csv");
    output_file << "peptide,nterm,cterm,miss1,miss2,label,accession,protein_id,peptide_id" << std::endl;
    for (std::size_t i = 0; i < proteins.size(); ++i) {
        if (!(*proteins[i]).Is_Digested()) {
            sequences = {};
            n_terms = {};
            c_terms = {};
            sequence = (*proteins[i]).Get_Sequence();
            cut = 0;
            previous_cut = 0;
            while (cut < sequence.size()) {
                cut = Next_Cut(sequence, cut) + 1;
                sequences.push_back(sequence.substr(previous_cut, cut - previous_cut));
                if (previous_cut < 8) {
                    n_term = std::string(8 - previous_cut, 'Z') + sequence.substr(0, 7 + previous_cut);
                }
                else {
                    n_term = sequence.substr(previous_cut - 8, 15);
                }
                if (n_term.length() < 15) {
                    n_term += std::string(15 - n_term.length(), 'Z');
                }
                n_terms.push_back(n_term);
                if (cut < 8) {
                    c_term = std::string(7 - cut, 'Z') + sequence.substr(cut, 8 + cut);
                }
                else {
                    c_term = sequence.substr(cut - 8, 15);
                }
                if (c_term.length() < 15) {
                    c_term += std::string(15 - c_term.length(), 'Z');
                }
                c_terms.push_back(c_term);
                previous_cut = cut;
            }
            for (std::size_t j = 0; j < sequences.size(); ++j) {
                if (sequences[j].size() >= minimum_number_of_amino_acids && (maximum_number_of_amino_acids == 0 || (sequences[j].size() <= maximum_number_of_amino_acids))) {
                    if (sequences[j].find("U") != std::string::npos || n_terms[j].find("U") != std::string::npos || c_terms[j].find("U") != std::string::npos) {
                        continue;
                    }
                    if (sequences[j].find("B") != std::string::npos || n_terms[j].find("B") != std::string::npos || c_terms[j].find("B") != std::string::npos) {
                        continue;
                    }
                    if (sequences[j].find("Z") != std::string::npos || n_terms[j].find("Z") != std::string::npos || c_terms[j].find("Z") != std::string::npos) {
                        continue;
                    }
                    position = peptides_sequences.find(sequences[j]);
                    if (position == peptides_sequences.end()) {
                        peptides_sequences[sequences[j]] = pep_count;
                        output_file << sequences[j] << ',' << n_terms[j] << ',' << c_terms[j] << ",ZZZZZZZZZZZZZZZ,ZZZZZZZZZZZZZZZ,NA," << proteins.at(i)->Get_Accession() << "," << i << "," << pep_count << std::endl;
                        pep_count++;
                    }
                    else {
                        output_file << sequences[j] << ',' << n_terms[j] << ',' << c_terms[j] << ",ZZZZZZZZZZZZZZZ,ZZZZZZZZZZZZZZZ,NA," << proteins.at(i)->Get_Accession() << "," << i << "," << position->second << std::endl;
                    }
                }
            }
        }
    }
    output_file.close();
}