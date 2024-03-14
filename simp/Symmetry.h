#ifndef __SYMMETRY_H__
#define __SYMMETRY_H__

#include <vector>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <sstream>
#include <math.h>

namespace SAGA
{
#define var_lit(x) (((x) > 0) ? (x) : (-x))
#define sign_lit(x) (((x) > 0) ? (1) : (0))
#define unbliss(a, n) ((a <= n) ? (a) : (n - a))

	// int sym_variables[1000000];

	// Define a data structure for a literal in the SAT problem.
	typedef struct lit
	{
	public:
		int clause_num; // clause num, begin with 0
		int var_num;	// variable num, begin with 1
		int sense;		// is 1 for true literals, 0 for false literals.

		lit(const int &clause_num_,
			const unsigned int &var_num_,
			const unsigned int &sense_) : clause_num(clause_num_),
										  var_num(var_num_),
										  sense(sense_) {}

		lit(const unsigned int &var_num_,
			const unsigned int &sense_) : lit(-1, var_num_, sense_) {}

		lit(const unsigned int &var_num_) : lit(var_num_, 1) {}

		lit(const int &var_num_) : lit(var_lit(var_num_), sign_lit(var_num_)) {}

		lit()
		{
			clause_num = -1;
			var_num = -1;
			sense = -1;
		}

		~lit() {}

		inline operator int() const { return var_num * (sense ? 1 : -1); }
		inline bool operator==(const lit &l) const { return var_num == l.var_num && sense == l.sense; }
		inline bool operator!=(const lit &l) const { return var_num != l.var_num || sense != l.sense; }
		inline lit operator=(const lit &l)
		{
			clause_num = l.clause_num;
			var_num = l.var_num;
			sense = l.sense;
			return *this;
		}
	} lit;

	typedef struct Generator
	{

		std::unordered_map<unsigned int, int> images;

		~Generator() {}

		void add(const lit &a, const lit &b)
		{

			if (a.var_num == b.var_num)
				return;

			images[a.var_num] = a.sense ? ((int)b) : (-(int)b);
			images[b.var_num] = b.sense ? ((int)a) : (-(int)a);
		}

		lit get_image(const lit &a)
		{

			if (images.find(a.var_num) != images.end())
				return a.sense ? lit(images[a.var_num]) : lit(-images[a.var_num]);
			else
				return a;
		}

		lit get_image(const int &a)
		{

			if (images.find(var_lit(a)) != images.end())
				return a > 0 ? lit(images[a]) : lit(-images[-a]);
			else
				return lit(a);
		}

		int *get_symmetrical_assignment(int *assignment, int len)
		{

			int *sym_assignment = new int[len];
			sym_assignment[0] = -1;
			for (int i = 1; i < len; i++)
			{
				lit sym_i = get_image(i);
				sym_assignment[i] = sym_i.sense ? assignment[sym_i.var_num] : 1 - assignment[sym_i.var_num];
			}

			return sym_assignment;
		}

		bool are_symmetrical(const std::vector<int> &a, const std::vector<int> &b, int len)
		{

			for (int i = 1; i < len; i++)
			{
				if (a[i] == b[i])
					continue;
				lit sym_i = get_image(i);

				if (sym_i.sense)
				{
					if (a[i] != a[sym_i.var_num])
					{

						if (a[i] != b[sym_i.var_num])
						{
							// cout << "a[" << i<<"] != a[" << sym_i.var_num <<"] and a[" <<i <<"] != b[" << sym_i.var_num <<"] with sy_de_i.sense != 0" << endl;
							return false;
						}
					}
					else
					{

						if (a[i] == b[sym_i.var_num])
						{
							// cout << "a[" << i<<"] == a[" << sym_i.var_num <<"] and a[" <<i <<"] == b[" << sym_i.var_num <<"] with sy_de_i.sense != 0" << endl;
							return false;
						}
					}
				}
				else
				{
					if (a[i] != a[sym_i.var_num])
					{
						if (a[i] == b[sym_i.var_num])
						{
							// cout << "a[" << i<<"] != a[" << sym_i.var_num <<"] and a[" <<i <<"] == b[" << sym_i.var_num <<"] with sy_de_i.sense == 0" << endl;
							return false;
						}
					}
					else
					{
						if (a[i] != b[sym_i.var_num])
						{
							// cout << "a[" << i<<"] == a[" << sym_i.var_num <<"] and a[" <<i <<"] !	= b[" << sym_i.var_num <<"] with sy_de_i.sense == 0" << endl;
							return false;
						}
					}
				}
			}

			return true;
		}

		std::string to_string()
		{
			std::ostringstream oss;
			for (auto entry : images)
				oss << "(" << entry.first << "," << entry.second << ")";

			return oss.str();
		}
	} Generator;

	// Symmetry generators
	// std::vector<Generator> symmetry_generators;

	std::vector<Generator> parse_bliss(const std::string &file_name, const int &n_vars, int *sym_variables, int &num_sym_vars);
	std::vector<std::vector<lit>> parse_dimacs_formula(const std::string &filename);
	void print_clauses(const std::vector<std::vector<lit>> &clauses);
}

#endif // __SYMMETRY_H__