#pragma once

#include <iostream>
#include <string>
#include <unordered_map>

/*
* Ce fichier permet la gestion des acides aminés.
* Il contient une classe AminoAcid définissant un acide aminé par un ensemble d'attributs.
* Il contient également une map associant une lettre à un acide aminé.
*/

class AminoAcid {
private:
    const std::string name;
    const std::string short_name;
    const char letter;
    const double mono_isotopic_mass;
    const double average_mass;
    const double pI;
    const double hydrophobicity;
    const double occurance;
    const double pK;
public:
    AminoAcid(std::string name, std::string short_name, char letter, double mono_isotopic_mass, double average_mass, double pI, double hydrophobicity, double occurance, double pK = 0.0) :
        name(name), short_name(short_name), letter(letter), mono_isotopic_mass(mono_isotopic_mass), average_mass(average_mass), pI(pI), hydrophobicity(hydrophobicity), occurance(occurance), pK(pK) {};
    ~AminoAcid() {};

    const std::string Get_Name() const { return name; }
    const std::string Get_Short_Name() const { return short_name; }
    const char Get_Letter() const { return letter; }
    const double Get_Mono_Isotopic_Mass() const { return mono_isotopic_mass; }
    const double Get_Average_Mass() const { return average_mass; }
    const double Get_PI() const { return pI; }
    const double Get_Hydrophobicity() const { return hydrophobicity; }
    const double Get_Occurance() const { return occurance; }
    const double Get_PK() const { return pK; }

    static double Get_B_Mass() { return 1.007276466879; } // H+
    static double Get_Y_Mass() { return 15.99491461956 + 1.007825032241 + 1.007825032241 + 1.007276466879; } // CT(OH) + NT + H+
};

std::unordered_map<char, AminoAcid> Amino_Acids = {
    {'A' ,AminoAcid("alaline" ,"Ala" ,'A' ,71.037113785565 ,71.08 ,6.01 ,1.8 ,0.078755)},
    {'R' , AminoAcid("arginine" ,"Arg" ,'R' ,156.101111025652 ,156.19 ,10.76 ,-4.5 ,0.054255 ,12.0)},
    {'N' , AminoAcid("asparagine" ,"Asn" ,'N' ,114.042927442166 ,114.10 ,5.41 ,-3.5 ,0.041355)},
    {'D' , AminoAcid("aspartic acid" ,"Asp" ,'D' ,115.026943024685 ,115.09 ,2.77 ,-3.5 ,0.053455 ,4.4)},
    {'C' , AminoAcid("cysteine" ,"Cys" ,'C' ,103.009184785565 ,103.15 ,5.07 ,2.5 ,0.0150 ,8.5)},
    {'Q' , AminoAcid("glutamine" ,"Gln" ,'Q' ,128.058577506648 ,128.13 ,5.65 ,-3.5 ,0.039655)},
    {'E' , AminoAcid("glutamic acid" ,"Glu" ,'E' ,129.042593089167 ,129.12 ,3.22 ,-3.5 ,0.066655 ,4.4)},
    {'G' , AminoAcid("glycine" ,"Gly" ,'G' ,57.021463721083 ,57.05 ,5.97 ,-0.4 ,0.069555)},
    {'H' , AminoAcid("histidine" ,"His" ,'H' ,137.058911859647 ,137.14 ,7.59 ,-3.2 ,0.022955 ,6.5)},
    {'I' , AminoAcid("isoleucine" ,"Ile" ,'I' ,113.047678469607 ,113.16 ,6.02 ,4.5 ,0.059155)},
    {'L' , AminoAcid("leucine" ,"Leu" ,'L' ,113.084063979011 ,113.16 ,5.98 ,3.8 ,0.096555)},
    {'K' , AminoAcid("lysine" ,"Lys" ,'K' ,128.094963016052 ,128.17 ,9.74 ,-3.9 ,0.059255 ,10.0)},
    {'M' , AminoAcid("methionine" ,"Met" ,'M' ,131.040484914529 ,131.20 ,5.74 ,1.9 ,0.023955)},
    {'F' , AminoAcid("phenylalanine" ,"Phe" ,'F' ,147.068413914529 ,147.18 ,5.48 ,2.8 ,0.039555)},
    {'P' , AminoAcid("proline" ,"Pro" ,'P' ,97.052763850047 ,97.12 ,6.48 ,-1.6 ,0.048255)},
    {'S' , AminoAcid("serine" ,"Ser" ,'S' ,87.032028405125,87.08,5.68,-0.8,0.068455)},
    {'T' , AminoAcid("threonine" ,"Thr" ,'T' ,101.047678469607 ,101.11 ,5.87 ,-0.7 ,0.054155)},
    {'W' , AminoAcid("tryptophan" ,"Trp" ,'W' ,186.07931295157 ,186.21 ,5.89 ,-0.9 ,0.011355)},
    {'Y' , AminoAcid("tyrosine" ,"Tyr" ,'Y' ,163.063328534089 ,163.18 ,5.66 ,-1.3 ,0.030255 ,10.0)},
    {'V' , AminoAcid("valine" ,"Val" ,'V' ,99.068413914529 ,99.13 ,5.97 ,4.2 ,0.067355)},
    {'U' , AminoAcid("selenocysteine" ,"Sel" ,'U' ,168.964198469607 ,168.00 ,100.0 ,100.0 ,100.0)}, // A améliorer
    {'O' , AminoAcid("pyrrolysine" ,"Pyr" ,'O' ,255.158291550141 ,255.00 ,100.0 ,100.0 ,100.0)} // A améliorer
};