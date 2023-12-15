#pragma once
#include<vector>
using namespace std;
class MyMatrix2d {
private:
	int rows, cols;
	double** data;
	void initialize();
public:
	MyMatrix2d(int, int);
	MyMatrix2d(int, int, double);
	MyMatrix2d(int, int, double*);
	MyMatrix2d(const MyMatrix2d&);
	MyMatrix2d(const vector<vector<double> >);
	virtual ~MyMatrix2d();
	double at(int y, int x);
	MyMatrix2d getPatch(int ys, int ye, int xs, int xe);
	void Mat2MEXpr(double*);
	void show();
	int size(int dimension);
	MyMatrix2d operator=(const MyMatrix2d&);
	double** getPr();
};