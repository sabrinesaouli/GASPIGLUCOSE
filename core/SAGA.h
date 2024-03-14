#ifndef _SAGA_H_
#define _SAGA_H_

#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <chrono>
#include <unordered_map>
#include <sstream>
#include <math.h>
#include <utility>
#include <random>
#include "core/SolverTypes.h"
#include "core/Solver.h"
#include "simp/Symmetry.h"

#include <boost/compute/detail/lru_cache.hpp>

namespace SAGA
{

    //  int num_of_clauses;
    //  int num_of_vars;
    //  std::vector<int> clause_lit_count;
    //  std::vector<int> clause_delete;
    //  std::vector<int> var_lit_count;
    //  std::vector<int> fix;
    //  std::vector<int> fixed_vars;

    class Formula
    {
    public:
        Formula() = default;
        Formula(int numVariables, int numClauses)
            : numVariables(numVariables), numClauses(numClauses),
              //   clauses(numClauses, std::vector<lit>(numVariables)),
              //   clause_lit_count(numClauses, 0),
              //   clause_delete(numClauses, 0),
              //   var_lit_count(numVariables + 1, 0),
              fix(numVariables + 1, 0),
              fixed_vars(numVariables + 1, 0)
        {
        }

        Formula(const Formula &other) : numVariables(other.numVariables), numClauses(other.numClauses),
                                        // clauses(other.clauses),
                                        // clause_lit_count(other.clause_lit_count),
                                        // clause_delete(other.clause_delete),
                                        // var_lit_count(other.var_lit_count),
                                        fix(other.fix),
                                        fixed_vars(other.fixed_vars),
                                        generators(other.generators),
                                        sym_variables(other.sym_variables),
                                        num_sym_vars(other.num_sym_vars),
                                        degree_centrality_variables(other.degree_centrality_variables) //,
                                                                                                       // centrality_computation_time(other.centrality_computation_time)

        {
        }
        Formula(int numVariables, int numClauses, const std::vector<std::vector<lit>> &clauses)
            : numVariables(numVariables), numClauses(numClauses),
              //   clauses(clauses),
              //   clause_lit_count(numClauses, 0),
              //   clause_delete(numClauses, 0),
              //   var_lit_count(numVariables + 1, 0),
              fix(numVariables + 1, false),
              fixed_vars(numVariables + 1, false)
        {
        }
        ~Formula() = default;

        Formula &operator=(const Formula &other)
        {
            numVariables = other.numVariables;
            numClauses = other.numClauses;
            // clauses = other.clauses;
            // clause_lit_count = other.clause_lit_count;
            // clause_delete = other.clause_delete;
            // var_lit_count = other.var_lit_count;
            fix = other.fix;
            fixed_vars = other.fixed_vars;
            degree_centrality_variables = other.degree_centrality_variables;
            generators = other.generators;
            sym_variables = other.sym_variables;
            num_sym_vars = other.num_sym_vars;
            // centrality_computation_time = other.centrality_computation_time;
            return *this;
        }

        // void addClause(const std::vector<lit> &clause, int clause_num)
        // {
        //     clauses.at(clause_num) = clause;
        // }

        int getNumVariables() const
        {
            return numVariables;
        }

        int getNumClauses() const
        {
            return numClauses;
        }

        // const std::vector<std::vector<lit>> &getClauses() const
        // {
        //     return clauses;
        // }

        // const std::vector<lit> &operator[](std::size_t i) const
        // {
        //     return clauses[i];
        // }

        // std::vector<lit> &operator[](std::size_t i)
        // {
        //     return clauses[i];
        // }

        // std::vector<int> &getClauseLitCount()
        // {
        //     return clause_lit_count;
        // }

        // std::vector<bool> &getClauseDelete()
        // {
        //     return clause_delete;
        // }

        // std::vector<int> &getVarLitCount()
        // {
        //     return var_lit_count;
        // }

        std::vector<bool> &getFix()
        {
            return fix;
        }

        std::vector<bool> &getFixedVars()
        {
            return fixed_vars;
        }

        void setNumVariables(int numVariables)
        {
            this->numVariables = numVariables;
            // this->var_lit_count.resize(numVariables + 1);
            this->fix.resize(numVariables + 1);
            this->fixed_vars.resize(numVariables + 1);
        }

        void setNumClauses(int numClauses)
        {
            this->numClauses = numClauses;
            // this->clauses.resize(numClauses);
            // this->clause_lit_count.resize(numClauses);
            // this->clause_delete.resize(numClauses);
        }

        // void setClauses(const std::vector<std::vector<lit>> &clauses)
        // {
        //     this->clauses = clauses;
        // }

        void parse_bliss(std::string &in)
        {
            std::ifstream fin(in);
            std::string line;
            char _;
            sym_variables.resize(numVariables + 1, 0);
            generators.clear();

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
                        if (sym_variables[var_lit(unbliss(a, numVariables))] == 0)
                        {
                            sym_variables[var_lit(unbliss(a, numVariables))] = 1;
                            // cout << var(unbliss(a, numVariables)) << " ";
                            num_sym_vars++;
                            // cout << num_sym_vars << endl;
                        }
                        else
                        {
                            sym_variables[var_lit(unbliss(a, numVariables))] += 1;
                        }

                        if (sym_variables[var_lit(unbliss(b, numVariables))] == 0)
                        {
                            sym_variables[var_lit(unbliss(b, numVariables))] = 1;
                            // cout << var(unbliss(b, numVariables)) << " ";
                            num_sym_vars++;
                            // cout << num_sym_vars << endl;
                        }
                        else
                        {
                            sym_variables[var_lit(unbliss(b, numVariables))] += 1;
                        }

                        g.add(lit(unbliss(a, numVariables)), lit(unbliss(b, numVariables)));
                    }

                    generators.push_back(g);
                    // std::cout << "Generator: " << generators.size() << std::endl;
                }
            }
        }

        std::vector<unsigned> &get_degree_centrality_variables()
        {
            return degree_centrality_variables;
        }
        void set_degree_centrality_variables(std::vector<unsigned> &degree_centrality_variables)
        {
            this->degree_centrality_variables = degree_centrality_variables;
        }

        // void read_important_variables(const std::string &in)
        // {
        //     std::ifstream fin(in);
        //     std::string line;
        //     int var;
        //     // std::cout << "Reading centrality variables:" << std::endl;
        //     // read computation time
        //     std::getline(fin, line);
        //     std::istringstream iss(line);
        //     iss >> centrality_computation_time;
        //     while (std::getline(fin, line))
        //     {
        //         std::istringstream iss(line);

        //         if (iss >> var)
        //         {
        //             // std::cout << var << std::endl;
        //             degree_centrality_variables.push_back(var);
        //         }
        //     }
        // }

        // double getCentralityTime()
        // {
        //     return centrality_computation_time;
        // }
        std::vector<SAGA::Generator> generators; // tab of all generators
        std::vector<int> sym_variables;          // tab of all variables present in generators (one hot encoding)
        int num_sym_vars = 0;                    // number of variables present in generators

        // std::vector<int> clause_lit_count; // Number of literals in each clause.
        // std::vector<bool> clause_delete;   // True if clause is deleted.
        // std::vector<int> var_lit_count;    // Number of times a variable is used as a literal.
        std::vector<bool> fix;        // fixed variable.
        std::vector<bool> fixed_vars; // Values of  Fixed variables.

    private:
        std::vector<unsigned> degree_centrality_variables;
        int numVariables;
        int numClauses;
        // float centrality_computation_time = 0;
    };

    class Solution
    {
    public:
        Solution(const Solution &other) : solution(other.solution), fitness(other.fitness), unsatisfying_variables(other.unsatisfying_variables)
        {
        }
        Solution(std::vector<int> solution_, int fitness_) : solution(solution_), fitness(fitness_), unsatisfying_variables()
        {
        }
        Solution(size_t solution_size, int nclauses)
            : solution(solution_size), fitness(nclauses), unsatisfying_variables()
        {
            solution.reserve(solution_size);
        }

        Solution() : solution(), fitness(1000000) {}

        ~Solution() {}

        // Operators overloading
        Solution &operator=(const Solution &other)
        {
            if (this != &other)
            {
                solution = other.solution;
                fitness = other.fitness;
                unsatisfying_variables = other.unsatisfying_variables;
            }

            return *this;
        }

        bool operator==(const Solution &s) const
        {
            return solution == s.solution && fitness == s.fitness;
        }

        bool operator!=(const Solution &s) const
        {
            return !(*this == s);
        }

        const int &operator[](std::size_t i) const { return solution[i]; }
        int &operator[](std::size_t i) { return solution[i]; }

        std::string toString() const
        {
            std::ostringstream oss;
            for (int i = 1; i < size(); ++i)
                oss << (solution[i] ? "" : "-") << i << " ";
            oss << std::endl;
            return oss.str();
        }

        // Getters and setters
        int getFitness() const { return fitness; }
        void setFitness(int fitness_) { fitness = fitness_; }
        std::size_t size() const { return solution.size(); }
        void resize(std::size_t size_) { solution.resize(size_); }
        std::vector<int> getSolution() const { return solution; }
        void setSolution(std::vector<int> solution_) { solution = solution_; }
        std::vector<int> get_unsatisfying_variables() const { return unsatisfying_variables; }

        bool no_unsatisfying_variables() const
        {
            return unsatisfying_variables.empty();
        }

        void add_unsatisfying_variable(int var)
        {
            // Check if the variable already exists in the vector
            if (std::find(unsatisfying_variables.begin(), unsatisfying_variables.end(), var) == unsatisfying_variables.end())
            {
                unsatisfying_variables.push_back(var);
            }
        }

        // int getRandomUnsatisfyingVariable()
        // {
        //     std::random_device rd;
        //     std::mt19937 rng(rd());
        //     // std::cout << "unsat var size: " << unsatisfying_variables.size() << std::endl;
        //     std::uniform_int_distribution<int> dis(0, unsatisfying_variables.size() - 1);
        //     return unsatisfying_variables[dis(rng)].first;
        // }

        // int getMaxUnsatisfyingVariable()
        // {

        //     int max = 0;
        //     int max_var = 0;
        //     for (auto var : unsatisfying_variables)
        //     {
        //         if (var.second > max)
        //         {
        //             max = var.second;
        //             max_var = var.first;
        //         }
        //     }
        //     return max_var;
        // }

        bool is_symmetric(Solution s, std::vector<SAGA::Generator> &generators)
        {

            if (fitness != s.getFitness())
                return false;

            for (auto g : generators)
            {
                // std::cout << "is_symmetric =" << g.are_symmetrical(solution, s.solution, num_vars+1) << std::endl;
                if (g.are_symmetrical(solution, s.getSolution(), solution.size()))
                    return true;
            }
            return false;
        }

    private:
        std::vector<int> solution;
        std::vector<int> unsatisfying_variables; // Vector that maps unsatisfying variables to the number of clauses they do not satisfy
        int fitness;
    };

    template <typename T>
    void hash_combine(std::size_t &seed, const T &value)
    {
        std::hash<T> hasher;
        seed ^= hasher(value) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }

    // custom hash function for the Solution class
    struct SolutionHash
    {
        size_t operator()(const Solution &solution) const
        {
            size_t seed = 0;
            // use the hash_combine function to combine hashes of individual members
            hash_combine(seed, solution.size());
            for (int i = 0; i <= solution.size(); ++i)
            {
                hash_combine(seed, solution[i]);
            }
            return seed;
        }
    };

    class Population
    {
    public:
        Population(int population_size_) : population_size(population_size_)
        {
            population.reserve(population_size_);
        }
        Population(const Population &other) : population(other.population), population_size(other.population_size) {}
        Population(int population_size_, int solution_size, int nclauses)
            : population(population_size_), population_size(population_size_)
        {
            population.reserve(population_size_);
            for (int i = 0; i < population_size_; ++i)
            {
                population.emplace_back(solution_size, nclauses);
            }
        }
        ~Population() {}

        // Operators overloading
        Population &operator=(const Population &other)
        {
            if (this != &other)
            {
                population = other.population;
                population_size = other.population_size;
            }

            return *this;
        }

        bool operator==(const Population &p) const
        {
            return population == p.population;
        }

        bool operator!=(const Population &p) const
        {
            return !(*this == p);
        }

        const Solution &operator[](std::size_t i) const { return population[i]; }
        Solution &operator[](std::size_t i) { return population[i]; }

        // Getters and setters
        std::size_t size() const { return population.size(); }
        void resize(std::size_t size_) { population.resize(size_); }
        void setPopulation(std::vector<Solution> population_)
        {
            population.clear();
            population.shrink_to_fit();
            population = population_;
        }
        std::vector<Solution> getPopulation() const { return population; }
        void push_back(Solution solution) { population.push_back(std::move(solution)); }

        void sort()
        {
            std::sort(population.begin(), population.end(), [](const Solution &a, const Solution &b)
                      { return a.getFitness() < b.getFitness(); });
        }
        std::string toString() const
        {
            std::ostringstream oss;
            for (const auto &sol : population)
            {
                oss << sol.toString() << "\t fitness = " << sol.getFitness() << std::endl;
            }
            return oss.str();
        }

    private:
        int population_size;              // Number of solutions per generation
        std::vector<Solution> population; // Population of solutions
    };

    class GeneticAlgorithm
    {
    public:
        GeneticAlgorithm(int population_size, int solution_size, int max_iterations, float mutation_rate,
                         float crossover_rate, Formula &formula, Glucose::Solver &solver)
            : population_size_(population_size),
              solution_size_(solution_size + 1),
              max_iterations_(max_iterations),
              mutation_rate_(mutation_rate),
              crossover_rate_(crossover_rate),
              population_(population_size),
              formula_(formula),
              solver_(solver)
        {
            population_ = Population(population_size_);
        }
        ~GeneticAlgorithm() {}

        Solution solve()
        {
            // Create a random number generator with a fixed seed for reproducibility
            std::random_device rd;
            std::mt19937 rng(rd());
            std::cout << "c |  Initializing the population ..." << std::endl;
            initialize_population(rng);

            evaluate_fitness();

            for (int iteration = 0; iteration < max_iterations_; ++iteration)
            {
                std::vector<Solution> parents = select_parents_tournament(rng);

                std::vector<Solution> offspring = create_offspring(parents, rng);

                evaluate_fitness(offspring);
                select_survivors_ellitist(offspring);

                // Clear the parents and offspring vectors
                parents.clear();
                offspring.clear();
                parents.shrink_to_fit();
                offspring.shrink_to_fit();

                // Check if a solution has been found
                if (solution_found())
                    break;

                // if (iteration % 10 == 0)
                // {
                //     std::cout << "c |  Iteration " << iteration << std::endl
                //               << "c |  Best fitness: " << population_[0].getFitness() << std::endl
                //               << "c |  Worst fitness: " << population_.getPopulation().back().getFitness() << std::endl;
                //     std::cout << "c |  \t ==================================== " << std::endl;
                // }
            }

            return population_[0];
        }

        boost::compute::detail::lru_cache<std::vector<int>, int> memo{500}; // cache memory to store the results of previous calls to fitness_unsat
        // std::unordered_map<Solution, int, SolutionHash> cache_memo{};       // cache memory to store the results of previous calls to fitness_unsat

        Solution &getBestSolution();
        Solution &getWorstSolution();
        bool is_symmetric(const Solution &s1, const Solution &s2)
        {
            if (s1.getFitness() != s2.getFitness())
                return false;
            for (auto g : formula_.generators)
            {

                if (g.are_symmetrical(s1.getSolution(), s2.getSolution(), formula_.getNumVariables() + 1))
                    return true;
            }
            return false;
        }

    private:
        size_t population_size_;
        size_t solution_size_;
        int max_iterations_;
        float mutation_rate_;
        float crossover_rate_;
        Population population_;
        Formula formula_;
        Glucose::Solver &solver_;

        void initialize_population(std::mt19937 rng);
        void evaluate_fitness(std::vector<Solution> &offspring);
        void evaluate_fitness();
        int fitness(Solution &solution);
        Solution select_parent(std::mt19937 rng);
        std::vector<Solution> create_offspring_two_points(std::mt19937 rng);
        std::vector<Solution> create_offspring_three_points(std::mt19937 rng);
        std::vector<Solution> create_offspring(std::mt19937 rng);
        void mutate(Solution &solution, std::mt19937 rng);
        std::vector<Solution> select_parents_tournament(std::mt19937 rng);
        std::vector<Solution> select_parents_random(std::mt19937 rng);
        std::vector<Solution> create_offspring(const std::vector<Solution> &parents, std::mt19937 rng);
        std::vector<Solution> create_offspring_two_points(const std::vector<Solution> &parents, std::mt19937 rng);
        std::vector<Solution> create_offspring_three_points(const std::vector<Solution> &parents, std::mt19937 rng);
        void select_survivors(const std::vector<Solution> &offspring);
        void select_survivors_ellitist(const std::vector<Solution> &offspring);
        bool solution_found();
    };

    void initialize_polarity(Solution &solution, Glucose::Solver &solver);

}

#endif // __SAGA_H__
