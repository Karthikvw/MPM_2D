// Some Output formates for used STL Containers
#pragma once

#include <fstream>
#include <iostream>
#include <map>
#include <array>
#include <vector>

inline std::ostream& operator<<(std::ostream &os, std::map<int, double> &m)
{
    for (std::map<int, double>::iterator it = m.begin(); it != m.end(); ++it) 
    { 
    os << "\t[ " << it->first << " ]\t[ " << it->second << " ] \n"; 
    } 
    return os;
};

inline std::ostream& operator<<(std::ostream &os, std::map<int, std::array<double, 3>> &m)
{
    for (std::map<int, std::array<double, 3>>::iterator it = m.begin(); it != m.end(); ++it) 
    { 
    os << "\t[ " << it->first << " ]\t[ " << it->second[0] << ", " << it->second[1] << ", " << it->second[2] << " ] \n"; 
    } 
    return os;
};

template<typename T, size_t size>
inline std::ostream& operator<<(std::ostream &os, std::array<T, size> v)
{
    os << "[ ";
    for (auto value : v) os << value << " ";
    os << " ] \n" ;
    return os;
};

template<typename T>
inline std::ostream& operator<<(std::ostream &os, std::vector<T> v)
{
    os << "[ ";
    for (auto value : v) os << value << " ";
    os << " ]" ;
    return os;
};
