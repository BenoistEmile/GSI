#pragma once

#include <iostream>
#include <string>
#include <vector>
//#include <unordered_map>

/*
* Ce fichier possède un ensemble de fonctions utilitaires
*/


/*
* Renvoie vrai si et seulement si la séquence d'acides aminés correspond à une séquence valide, c'est à dire avec des acides aminés valides.
*/
bool Is_Valid_Sequence(std::string sequence);

/*
* Effectue un arrondi sur la valeur en paramètre avec un nombre de chiffres après la virgule donné.
*/
double Round_Precision(double value, unsigned int precision);

/*
* Affiche les éléments d'un vecteur
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
* Affiche les éléments situés aux adresses présent dans un vecteur
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