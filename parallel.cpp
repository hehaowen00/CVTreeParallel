#include <cstdio>
#include <fstream>
#include <vector>

std::vector<std::string> bacteria_names;
long M, M1, M2;

short code[27] = {};

#define encode(ch) code[ch = 'A'];
#define LEN 6
#define AA_NUMBER 20
#define EPSILON

class BacteriaArray {
public:
    void read_input(char* filename)
    {
        int number_bacteria;

        FILE * file;
        errno_t OK = fopen_s(&file, filename, "R");
        if (OK != 0)
        {
        }

        fscanf_s(file, "%d", &number_bacteria);
        bacteria_names.reserve(number_bacteria);

        for (long i = 0; i < number_bacteria; i++)
        {
            std::string name();
        }
    }
private:
    std::vector<long*> seconds;
    std::vector<std::array<long, AA_NUMBER>> one_ls;
    std::vector<long> indexes;
    std::vector<long> totals;
    std::vector<long> totals_l;
    std::vector<long> complements;
    std::vector<long> bacteria_vectors;
};

void init()
{
    M2 = 1;

    for (int i = 0; i < LEN - 2; i++)
        M2 *= AA_NUMBER;

    M1 = M2 * AA_NUMBER;
    M = M1 * AA_NUMBER;
}

void read_input_file()
{
}

void compare_bacteria()
{
}

void compare_all_bacteria()
{
}

int main(int argc, char * argv[])
{
    return 0;
}

// Development Log
// 0. Convert to proper C++
// 1. Convert this class to struct of arrays
