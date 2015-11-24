#include "Collect.h"
#include <iostream>
#include <fstream>
#include <iomanip> // for formatting
#include <algorithm>
#include "MatrixCPP.h"

Collect::Collect()
     {cols=0;
     isSetup=0;}

Collect::Collect(int Cols)
     {cols=0;
      set(Cols);}

Collect::~Collect()
    {clear_all();}

void Collect::set(int Cols)
     {
     if(cols!=0) clear_all();
     cols=Cols;
     out.resize(cols);
     isSetup=1;}

void Collect::clear_all(void)
     {for(int i =0; i<cols; i++) out[i].clear();
     out.clear();
     cols=0;
     }

void Collect::feed(int col, double val)
     {
     if(col>=cols) {std::cerr << "column number too high in collect::feed" << std::endl; system("pause"); exit(EXIT_FAILURE);}
     out[col].push_back(val);}

// Feed an array of number. If the container is a vector, then need to feed it &(*temp.begin())
void Collect::feed_col(int col, double* vals, int n)
     {for(int i =0; i<n; i++) feed(col, vals[i]);}


void Collect::feed_row(double* vals, int n){
     if(n!=cols) {std::cerr << "Incorrect vector length in collect::feed_row" << std::endl; system("pause"); exit(EXIT_FAILURE);}
     for(int i =0; i<n; i++) feed(i, vals[i]);}


/*Prints last row corresponding to data in column1 to file. Columns w missing data are printed as blank*/
void Collect::print_last_row(std::ofstream* OutFile, bool end_line)
	{
	if(!(*OutFile)) {std::cerr << "File not opened" << std::endl; system("pause"); exit(EXIT_FAILURE);}

    int row=out[0].size();
    for(int i =0; i<cols; i++)
	{if((int) out[i].size() < row)
		(*OutFile) << " \t";
	 else
		(*OutFile)<<out[i][row-1] << "\t";
         }
    if(end_line)
	(*OutFile)<<std::endl;
    }


// Creates a matrix of type double** based on values
double** Collect::return_as_matrix(int& Rows, int& Cols)
     {
     std::vector<int> rows(cols, 0);
     int i, j, max_rows=0;
     for(j =0; j<cols; j++)
           {rows[j]=out[j].size();
            max_rows=std::max(max_rows, rows[j]);}
     Rows =max_rows; Cols = cols;
     double** Matrix = M2Dd_alloc(Rows, Cols);
     for(i =0; i<max_rows; i++)
         {for(j =0; j<cols; j++)
                {
                if(rows[j]>i) Matrix[i][j] = out[j][i];
                else          Matrix[i][j] = 0;
                }
          }
     return Matrix;
     }

/*Output contents of vector to Matrix. append specifies whether to try and append
details to previous file, if it exists */
void Collect::print_file(std::string name, int format, int append)
      {
     std::ofstream OutFile;
     if(append)
               OutFile.open(name.c_str(), std::ios::app);
     else
               OutFile.open(name.c_str(), std::ios::out);
     if(!OutFile) {std::cerr << "File not opened" << std::endl; system("pause"); exit(EXIT_FAILURE);}
     if(format ==0) // Integers
               OutFile<<std::setprecision(0);
     else if(format ==1) // Floating
               OutFile<<std::setiosflags(std::ios::scientific|std::ios::floatfield);

     std::vector<int> rows(cols, 0);
     int i, j, max_rows=0;
     for(j =0; j<cols; j++)
           {rows[j]=out[j].size();
            max_rows=std::max(max_rows, rows[j]);
            }
     for(i =0; i<max_rows; i++)
         {for(j =0; j<cols; j++)
                {if(rows[j]>i) OutFile<<out[j][i] << "\t";
                else OutFile << "\t";}
          OutFile<<std::endl;
         }
     OutFile.close();
     }

