#include "Utils.h"
#include "MatrixCPP.h"


// First derivative calculated numerically with h2 truncation error
double dfdx(double f_xhigh, double f_xlow, double h)
     {return((f_xhigh - f_xlow)/2.0/h);}
// Second derivative calculated numerically with h2 truncation error
double df2d2x(double f_x_l, double f_x, double f_x_h, double h)
     {return((f_x_l - 2.0*f_x + f_x_h)/h/h);}

bool IsDirectory(const char* FileName)
     {struct stat my_stat;
     if (stat(FileName, &my_stat) != 0) return false;
     return ((my_stat.st_mode & S_IFDIR) != 0);}

bool FileExist(const char* FileName)
     {struct stat my_stat;
     return (stat(FileName, &my_stat) == 0);}

bool equals(double d1, double d2, double precision)
	{
	double eps1 = fabs(d1), eps2 = fabs(d2), eps;
	eps = (eps1 > eps2) ? eps1 : eps2;
	if (eps == 0.0)
	return true; // Both numbers are 0.0
	// Eps hold the minimum distance between the values that will be considered as the numbers are equal
	// Considering the magnitude of the numbers
	eps *= precision;
	return (fabs(d1 - d2) < eps);
	}

// Check sign of a double number
double sign(double temp){
       if(temp>= 0) return 1.0;
       else return -1.0;}

// Returns number with max number of decimal places set by precis (without rounding)
double set_decimal(double T, int precis)
	{
	if(T==0) return 0.0;
	int X = int(log10(T));
	if(X<0) X = -X+precis;
	else X =precis;
	return floor(T*pow(10, X)) /pow(10, X);
	}

void print_file(std::string name, double** out, int rows, int cols)
     {
     std::ofstream OutFile;
     OutFile.open(name.c_str()); if(!OutFile) {std::cerr << "File not opened" << std::endl; system("pause"); exit(EXIT_FAILURE);}
     OutFile<<std::setiosflags(std::ios::scientific|std::ios::floatfield);
     int i, j;
     for(i =0; i<rows; i++)
         {for(j =0; j<cols; j++) OutFile<<out[i][j] << "\t";
          OutFile<<std::endl;
         }
     OutFile.close();
     }

double** vec2D_to_matrix(std::vector<std::vector<double> >& Vec, int &rows, int &cols)
    {
    double **a=NULL;
    int i, j;
    rows = Vec.size();
    if(rows==0) return a; // Safety

    // Find length of second vector dimension
    cols= 0;
    for(i=0; i< rows; i++)   cols=std::max(cols, (int) Vec[i].size());
    // Alocate memory
    a = M2Dd_alloc(rows, cols);
    // Copy data
    for(i = 0; i< rows; i++)
            {
            for(j=0; j< (int) Vec[i].size(); j++)    a[i][j] = Vec[i][j];
            for(j=  (int) Vec[i].size(); j<cols; j++)  a[i][j] = 0;
            }
    return a;
    }

void hold(void){
	char g;
	std::cout << "press return to continue... ";
	std::cin.get(g);
	}

void go_to_dir(std::string dir)
	{
	// New version - uses sheel command mkdir
	 if(!IsDirectory(dir.c_str()))
		system(("mkdir -p " + dir).c_str());

	chdir(dir.c_str());
	std::cout << "Moving to directory " << dir<<std::endl;
	// OLD VERSION - can only create one dir and does recognise home ditrectory symbol ~
	 // if(!IsDirectory(dir.c_str()))
	 //	mkdir(dir.c_str());
	 //     chdir(dir.c_str());
    }

