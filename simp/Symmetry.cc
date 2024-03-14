#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <cassert>
#include <cstring>
#include "Symmetry.h"
#include "core/SAGA.h"

using namespace SAGA;

//

namespace SAGA
{
    int num_of_clauses;
    int num_of_vars;
    std::vector<int> clause_lit_count;
    std::vector<int> clause_delete;
    std::vector<int> var_lit_count;
    std::vector<int> fix;
    std::vector<int> fixed_vars;
}
std::vector<std::vector<SAGA::lit>> SAGA::parse_dimacs_formula(const std::string &filename)
{
    std::ifstream input_file(filename);
    if (!input_file)
    {
        std::cerr << "Failed to open input file" << std::endl;
        return {};
    }

    std::vector<std::vector<lit>> clauses;
    std::string line;
    std::istringstream iss(line);
    char c;
    int clause_num = 0;

    std::vector<lit> clause;
    std::cout << "Parsing dimacs formula" << std::endl;

    while (std::getline(input_file, line))
    {
        iss.clear();
        iss.str(line);
        // iss >> c;

        if (iss.peek() == 'c')
        {
            iss >> c;
            // Ignore comment line
            continue;
        }
        else if (iss.peek() == 'p')
        {
            iss >> c;
            // Parse problem line
            std::string problem_type;
            int num_vars, num_clauses;
            iss >> problem_type >> num_vars >> num_clauses;

            std::cout << "Problem type: " << problem_type << std::endl;
            std::cout << "Number of variables: " << num_vars << std::endl;
            std::cout << "Number of clauses: " << num_clauses << std::endl;

            // Set number of variables and clauses in SAGA class
            SAGA::num_of_vars = num_vars;
            SAGA::num_of_clauses = num_clauses;

            // Initialize clauses vector with given size
            clauses.resize(SAGA::num_of_clauses);

            SAGA::clause_lit_count.resize(SAGA::num_of_clauses);
            SAGA::clause_delete.resize(SAGA::num_of_clauses);

            for (int c = 0; c < SAGA::num_of_clauses; c++)
            {
                SAGA::clause_lit_count[c] = 0;
                SAGA::clause_delete[c] = 0;
            }

            SAGA::var_lit_count.resize(SAGA::num_of_vars + 1);
            SAGA::fix.resize(SAGA::num_of_vars + 1);

            for (int v = 1; v <= SAGA::num_of_vars; ++v)
            {
                SAGA::var_lit_count[v] = 0;
                SAGA::fix[v] = 0;
            }
        }
        else
        {
            // Parse clause line
            int literal;

            clause.clear();
            // literal = c - '0';
            while (iss >> literal)
            {
                if (literal == 0)
                {
                    // End of clause
                    break;
                }
                else if (abs(literal) > num_of_vars)
                {
                    throw std::runtime_error("Invalid variable number");
                }
                else
                {
                    // Literal in clause
                    SAGA::clause_lit_count[clause_num]++;
                    SAGA::var_lit_count[abs(literal)]++;
                    clause.push_back(lit(clause_num + 1, abs(literal), (literal > 0) ? 1 : 0));
                }
            }

            SAGA::fixed_vars.resize(num_of_vars + 1);

            if (SAGA::clause_lit_count[clause_num] == 1)
            {
                SAGA::clause_delete[clause_num] = 1;
                SAGA::fix[abs(clause[0].var_num)] = 1;
                SAGA::fixed_vars[abs(clause[0].var_num)] = clause[0].sense;
            }

            if (clause.empty())
            {
                throw std::runtime_error("Empty clause");
            }

            if (clause_num > num_of_clauses)
            {
                throw std::runtime_error("Too many clauses");
            }

            // Assign clause to clauses vector at given index
            clauses.at(clause_num) = clause;
            clause_num++;
        }
    }

    return clauses;
}

// A function that prints a vector of clauses
void SAGA::print_clauses(const std::vector<std::vector<lit>> &clauses)
{
    for (const auto &clause : clauses)
    {
        std::cout << "(";
        for (const auto &literal : clause)
        {
            std::cout << literal.var_num * (literal.sense ? 1 : -1) << " ";
        }
        std::cout << ") ";
    }
    std::cout << "\n";
}

// == == == == == == == == == == ==     //
//     == Parse symmetries ==           //
//     == == == == == == == == == == == //

std::vector<Generator> SAGA::parse_bliss(const std::string &file_name, const int &n_vars, int *sym_variables, int &num_sym_vars)
{
    std::vector<Generator> generators;
    std::ifstream fin(file_name);
    std::string line;
    char _;
    std::fill(sym_variables, sym_variables + n_vars + 1, 0);

    while (std::getline(fin, line))
    {
        if (line == "[")
        {
            continue;
        }
        else if (line == "]")
        {
            break;
        }
        else
        {
            SAGA::Generator g;
            std::istringstream iss(line);
            int a, b;

            // read (a,b)
            while (iss >> _ >> a >> _ >> b >> _)
            {
                if (sym_variables[var_lit(unbliss(a, n_vars))] == 0)
                {
                    sym_variables[var_lit(unbliss(a, n_vars))] = 1;
                    // cout << var(unbliss(a, n_vars)) << " ";
                    num_sym_vars++;
                    // cout << num_sym_vars << endl;
                }
                else
                {
                    sym_variables[var_lit(unbliss(a, n_vars))] += 1;
                }

                if (sym_variables[var_lit(unbliss(b, n_vars))] == 0)
                {
                    sym_variables[var_lit(unbliss(b, n_vars))] = 1;
                    // cout << var(unbliss(b, n_vars)) << " ";
                    num_sym_vars++;
                    // cout << num_sym_vars << endl;
                }
                else
                {
                    sym_variables[var_lit(unbliss(b, n_vars))] += 1;
                }

                g.add(lit(unbliss(a, n_vars)), lit(unbliss(b, n_vars)));
            }

            generators.push_back(g);
        }
    }

    return generators;
}