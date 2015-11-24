//*******************************************************************************************/
//  Collect.h
//  Created by Daniel Falster, 2006-2010
//*******************************************************************************************/
#ifndef __COLLECT_H_
#define __COLLECT_H_

#include <vector>
#include <string>

class Collect {
public:
	Collect();
	Collect(int cols);
	~Collect();
	void set(int Cols);
	void clear_all(void);
	void feed(int col, double val);
	void feed_col(int col, double* vals, int n);
	void feed_row(double* vals, int n);
	void print_file(std::string name, int format, int append);
	double** return_as_matrix(int& Rows, int& Cols);
	void print_last_row(std::ofstream* OutFile, bool end_line);

private:
	int cols;
	std::vector<std::vector<double> > out;
	int isSetup;
};

#endif
