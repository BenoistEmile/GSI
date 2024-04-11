#pragma once

#include <iostream>
#include <string>
#include <vector>
//#include <unordered_map>
#include <sys/stat.h>

/*
* Ce fichier poss�de un ensemble de fonctions utilitaires
*/


/*
* Renvoie vrai si et seulement si la s�quence d'acides amin�s correspond � une s�quence valide, c'est � dire avec des acides amin�s valides.
*/
bool Is_Valid_Sequence(std::string sequence);

/*
* Effectue un arrondi sur la valeur en param�tre avec un nombre de chiffres apr�s la virgule donn�.
*/
double Round_Precision(double value, unsigned int precision);

/*
* Affiche les �l�ments d'un vecteur
*/
template<typename T>
void PrintVector(std::vector<T> vector) {
    if (vector.size()) {
        std::cout << "[";
        for (auto iter = vector.begin(); iter != vector.end() - 1; ++iter) {
            //Print<Elem>(*iter, false);
            std::cout << (*iter) << " ,";
        }
        std::cout << vector.back() << "]";
    }
    else {
        std::cout << "[]";
    }
    std::cout << std::endl;
};

/*
* Affiche les �l�ments situ�s aux adresses pr�sent dans un vecteur
*/
template<typename T>
void PrintVector2(std::vector<T> vector) {
    if (vector.size()) {
        std::cout << "[";
        for (auto iter = vector.begin(); iter != vector.end() - 1; ++iter) {
            //Print<Elem>(*iter, false);
            std::cout << *(*iter) << " ,";
        }
        std::cout << *vector.back() << "]";
    }
    else {
        std::cout << "[]";
    }
    std::cout << std::endl;
};

/*
* Renvoie vrai si le fichier indiqué existe.
*/
bool fileExists(const std::string& filename);