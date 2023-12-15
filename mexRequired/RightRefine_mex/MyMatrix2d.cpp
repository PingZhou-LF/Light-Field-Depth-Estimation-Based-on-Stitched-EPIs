#include "MyMatrix2d.h"
#include<vector>
#include<iostream>
using namespace std;

void MyMatrix2d::initialize() {
	data = new double* [rows];
	for (int i = 0; i < rows; i++) data[i] = new double[cols]();
}

MyMatrix2d::MyMatrix2d(int m, int n) {
	rows = m;
	cols = n;
	initialize();
}

MyMatrix2d::MyMatrix2d(int m, int n, double value) {
	rows = m;
	cols = n;
	initialize();
	for (int y = 0; y < rows; y++) {
		for (int x = 0; x < cols; x++) {
			data[y][x] = value;
		}
	}
}

MyMatrix2d::MyMatrix2d(int m, int n, double* pr_m) {
	rows = m;
	cols = n;
	initialize();
	for (int y = 0; y < rows; y++) {
		for (int x = 0; x < cols; x++) {
			data[y][x] = pr_m[x * m + y];
		}
	}
}

MyMatrix2d::MyMatrix2d(const MyMatrix2d& source) {
	rows = source.rows;
	cols = source.cols;
	initialize();
	for (int y = 0; y < rows; y++) {
		for (int x = 0; x < cols; x++) {
			data[y][x] = source.data[y][x];
		}
	}
}

MyMatrix2d::MyMatrix2d(const vector<vector<double> > source) {
	rows = source.size();
	cols = source[0].size();
	initialize();
	for (int y = 0; y < rows; y++) {
		for (int x = 0; x < cols; x++) {
			data[y][x] = source[y][x];
		}
	}
}

MyMatrix2d::~MyMatrix2d() {
	for (int i = 0; i < rows; i++) delete[]data[i];
	delete[]data;
}

double MyMatrix2d::at(int _row, int _col) {
	//if (_row < 0 || _row >= rows || _col < 0 || _col >= cols) return -1;
	return data[_row][_col];
}

MyMatrix2d MyMatrix2d::getPatch(int ys, int ye, int xs, int xe) {
	int height = ye - ys + 1;
	int length = xe - xs + 1;
	if (height <= 0 || length <= 0 || ys < 0 || ye > rows - 1 || xs < 0 || xe > cols - 1) {
		MyMatrix2d ans = MyMatrix2d(rows, cols);
		cout << "ERROR COORDINATES!" << endl;
		return ans; //坐标错误则返回一个与原矩阵大小一致的全0矩阵
	}
	MyMatrix2d ans = MyMatrix2d(height, length);
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < length; x++) {
			ans.data[y][x] = this->data[y + ys][x + xs];
		}
	}
	return ans;
}

void MyMatrix2d::Mat2MEXpr(double* pr_m) {
	for (int y = 0; y < rows; y++) {
		for (int x = 0; x < cols; x++) {
			pr_m[x * rows + y] = this->data[y][x];
		}
	}
}

void MyMatrix2d::show() {
	for (int y = 0; y < rows; y++) {
		for (int x = 0; x < cols; x++) {
			cout << this->data[y][x] << ' ';
		}
		cout << endl;
	}
	cout << endl;
}

int MyMatrix2d::size(int dimension) {
	if (dimension == 0) return rows;
	else return cols;
}

MyMatrix2d MyMatrix2d::operator=(const MyMatrix2d& source) {
	if (rows != source.rows || cols != source.cols) {
		for (int i = 0; i < rows; i++) delete[]data[i];
		delete[]data;
		rows = source.rows;
		cols = source.cols;
		initialize();
	}
	for (int y = 0; y < rows; y++) {
		for (int x = 0; x < cols; x++) {
			data[y][x] = source.data[y][x];
		}
	}
	return *this;
}

double** MyMatrix2d::getPr() {
	return data;
}