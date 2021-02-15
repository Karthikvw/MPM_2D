#include <iostream>
#include <vector>
#include <string>
#include <MPM_Grid.hpp>

//using namespace std;

int main(int argc, char **argv)
{

    std::vector<std::string> msg {"Hello", "C++", "World", "from", "VS Code!"};

    for (const std::string& word : msg)
    {
        std::cout << word << " ";
    } 
    std::cout << std::endl;

    int four = 4;
    double six = 6e0;
    double res = four * six;
    std::cout << "Number: " << res << std::endl;

}

