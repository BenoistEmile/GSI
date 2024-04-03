#pragma once

#include <iostream>
#include <string>
#include <vector>
//#include <unordered_map>

bool Is_Valid_Sequence(std::string sequence);

double Round_Precision(double value, unsigned int precision);

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
template<typename type>
void Print(type elem, bool final = true) {
    std::cout << type_info(elem);
    if (final) {
        std::cout << std::endl;
    }
};

template<typename Elem>
void Print(std::vector<Elem> table ,bool final = true) {
    if (table.size()) {
        std::cout << "[";
        for (auto iter = table.begin(); iter != table.end() - 1; ++iter) {
            Print<Elem>(*iter ,false);
            std::cout << " ,";
            if (new_line) {
                std::cout << "\n";
            }
        }
        std::cout << table.back() << "]";
    }
    else {
        std::cout << "[]";
    }
    if (final) {
        std::cout << std::endl;
    }
};
template<typename Elem>
void Print_Vector(std::vector<Elem> table, bool new_line = false) {
    if (table.size()) {
        std::cout << "[";
        for (auto iter = table.begin(); iter != table.end() - 1; ++iter) {
            std::cout << (*iter) << " ,";
            if (new_line) {
                std::cout << "\n";
            }
        }
        std::cout << table.back() << "]" << std::endl;
    }
    else {
        std::cout << "[]" << std::endl;
    }
};

template<typename Key ,typename Elem>
void Print_Unordered_Map(std::unordered_map<Key, Elem> map ,bool new_line = false) {
    if (map.empty()) {
        std::cout << "{}" << std::endl;
    }
    else {
        std::cout << "{";
        for (auto iter = map.begin(); iter != --map.end(); ++iter) {
            std::cout << iter->first << " : " << iter->second << ", ";
            if (new_line) {
                std::cout << "\n";
            }
        }
        std::cout << (--map.end())->first << " : " << (--map.end())->second << "}" << std::endl;
    }
};
*/