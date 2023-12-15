#include "mex.h"
#include "MyMatrix2d.h"
#include<queue>
#include<unordered_set>
#include<cmath>
#include<algorithm>

using namespace std;

#define w_m       prhs[0]
#define dis_m     prhs[1]
#define info_m    prhs[2]
#define hs_m      prhs[3]
#define hr_m      prhs[4]
#define deno_m    prhs[5]
#define occ_x_m   prhs[6]
#define occ_y_m   prhs[7]
#define key_m     plhs[0]
//#define key1      plhs[1]

struct coordinate {
    int y; //行
    int x; //列
    coordinate(int i, int j) {
        y = i;
        x = j;
    }
};
struct arguments {
    int Nt, Ns, Nx, Ny, cc, t0, s0, res;
    double alpha1, alpha2, dpth, angle_min, angle_max, step, dis_max, dis_min;
    double* csai;
};

template<typename T> T threshold(T input, T max, T min) {
    if (input > max) input = max;
    if (input < min) input = min;
    return input;
}

void setStructInfo(arguments& info, const mxArray* prhs[]) {
    mxArray* temp = mxGetField(info_m, 0, "Nt");
    info.Nt = mxGetScalar(temp);
    temp = mxGetField(info_m, 0, "Ns");
    info.Ns = mxGetScalar(temp);
    temp = mxGetField(info_m, 0, "Nx");
    info.Nx = mxGetScalar(temp);
    temp = mxGetField(info_m, 0, "Ny");
    info.Ny = mxGetScalar(temp);
    temp = mxGetField(info_m, 0, "s0");
    info.s0 = mxGetScalar(temp);
    temp = mxGetField(info_m, 0, "t0");
    info.t0 = mxGetScalar(temp);
    temp = mxGetField(info_m, 0, "dis_max");
    info.dis_max = mxGetScalar(temp);
    temp = mxGetField(info_m, 0, "dis_min");
    info.dis_min = mxGetScalar(temp);
    temp = mxGetField(info_m, 0, "res");
    info.res = mxGetScalar(temp);
    temp = mxGetField(info_m, 0, "alpha1");
    info.alpha1 = mxGetScalar(temp);
    temp = mxGetField(info_m, 0, "alpha2");
    info.alpha2 = mxGetScalar(temp);
    temp = mxGetField(info_m, 0, "dpth");
    info.dpth = mxGetScalar(temp);
    temp = mxGetField(info_m, 0, "cc");
    info.cc = mxGetScalar(temp);
    temp = mxGetField(info_m, 0, "angle_min");
    info.angle_min = mxGetScalar(temp);
    temp = mxGetField(info_m, 0, "angle_max");
    info.angle_max = mxGetScalar(temp);
    temp = mxGetField(info_m, 0, "step");
    info.step = mxGetScalar(temp);
    temp = mxGetField(info_m, 0, "csai");
    mxArray* temp1 = mxDuplicateArray(temp);
    info.csai = mxGetPr(temp1);
}

void setInput1d(double* input, double* arr, int length) {
    //if (length != arr.size()) {
    //    mexPrintf("ERROR DIMENSION !");
    //    //cout << "ERROR DIMENSION !" << endl;
    //    return;
    //}
    for (int i = 0; i < length; i++) {
        arr[i] = input[i];
    }
}

int mybwlabel(MyMatrix2d& im, MyMatrix2d& label) {
    if (im.size(0) != label.size(0) || label.size(1) != im.size(1)) {
        //cout << "ERROR DIMENSION" << endl;
        return -1;
    }
    int m = im.size(0), n = im.size(1);
    int label_num = 1;
    int neighbour[4][2] = { -1, 0, 0, -1, 0, 1, 1, 0 };
    queue<coordinate> queu;
    double index[2] = { 0 };
    double cur_index[2] = { 0 };
    MyMatrix2d im1(m + 2, n + 2);
    double** im1_pr = im1.getPr();
    MyMatrix2d label1(m + 2, n + 2);
    double** label1_pr = label1.getPr();
    for (int y = 0; y < m; y++) {
        for (int x = 0; x < n; x++) {
            im1_pr[y + 1][x + 1] = im.at(y, x);
        }
    }
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            if (im1_pr[i][j] == 1.0 && label1_pr[i][j] == 0) {
                //queue<coordinate> queu;
                label1_pr[i][j] = label_num;
                coordinate temp(i, j);
                queu.push(temp);
                while (!queu.empty()) {
                    index[0] = queu.front().y;
                    index[1] = queu.front().x;
                    for (int k = 0; k < 4; k++) {
                        cur_index[0] = index[0] + neighbour[k][0];
                        cur_index[1] = index[1] + neighbour[k][1];
                        if (cur_index[0] >= 1 && cur_index[0] <= m && cur_index[1] >= 1 && cur_index[1] <= n) {
                            if (im1_pr[(int)cur_index[0]][(int)cur_index[1]] == 1 && label1_pr[(int)cur_index[0]][(int)cur_index[1]] == 0) {
                                label1_pr[(int)cur_index[0]][(int)cur_index[1]] = label_num;
                                coordinate temp1((int)cur_index[0], (int)cur_index[1]);
                                queu.push(temp1);
                            }
                        }
                    }
                    queu.pop();
                }
                label_num++;
                if (label_num > m * n) {
                    mexPrintf("Error!");
                    //cout << "Error!" << endl;
                    return -1;
                }
            }
        }
    }
    double** label_pr = label.getPr();
    for (int y = 0; y < m; y++) {
        for (int x = 0; x < n; x++) {
            label_pr[y][x] = label1_pr[y + 1][x + 1];
        }
    }
    return label_num - 1;
}

void meanshsegm(MyMatrix2d& im, MyMatrix2d& label, double hs, double hr) {
    if (im.size(0) != label.size(0) || im.size(1) != label.size(1)) {
        mexPrintf("ERROR DIMENSION !");
        //cout << "ERROR DIMENSION !" << endl;
        return;
    }
    double** label_pr = label.getPr();
    int height = im.size(0);
    int length = im.size(1);
    //vector<vector<vector<double> > > clusterCenter(height, vector<vector<double> >(length, vector<double>(3, 0))); //7*7*3
    double clusterCenter[7][7][3] = { 0 };
    double h[3] = { 0 };
    h[0] = hs; h[1] = hs; h[2] = hr;
    for (int ixm = 1; ixm <= length; ixm++) {
        for (int iym = 1; iym <= height; iym++) {
            int xmin = max((double)(ixm - hs), 1.0);
            int xmax = min((double)(ixm + hs), (double)length);
            int xlen = xmax - xmin + 1;
            int ymin = max((double)(iym - hs), 1.0);
            int ymax = min((double)(iym + hs), (double)height);
            int ylen = ymax - ymin + 1;
            int len = xlen * ylen;
            //vector<int> iw(len, 0);
            int* iw = new int[len];
            for (int i = 0; i < len; i++) iw[i] = i;

            MyMatrix2d fw(len, 3);
            double** fw_pr = fw.getPr();
            MyMatrix2d im_patch = im.getPatch(ymin - 1, ymax - 1, xmin - 1, xmax - 1);
            double** im_patch_pr = im_patch.getPr();
            for (int i = 0; i < len; i++) {
                double temp = iw[i] / ylen + xmin;
                fw_pr[i][0] = (int)temp;
                fw_pr[i][1] = fmod(iw[i], ylen) + ymin;
                fw_pr[i][2] = im_patch_pr[i % ylen][(int)(i / ylen)];
            }
            //vector<double> center(3, 0);
            double center[3] = { 0 };
            double** im_pr = im.getPr();
            center[0] = ixm; center[1] = iym; center[2] = im_pr[iym - 1][ixm - 1];
            //vector<double> temp(3, 0);
            double temp[3] = { 0 };
            //vector<vector<double> > dis(len, vector<double>(3, 0));
            MyMatrix2d dis(len, 3);
            double** dis_pr = dis.getPr();
            //vector<double> dis1(len, 0);
            double* dis1 = new double[len]();
            while (1) {
                vector<int> index;
                for (int i = 0; i < len; i++) {
                    for (int j = 0; j < 3; j++) {
                        dis_pr[i][j] = (fw_pr[i][j] - center[j]) / h[j];
                    }
                    dis1[i] = pow(dis_pr[i][0], 2.0) + pow(dis_pr[i][1], 2.0) + pow(dis_pr[i][2], 2.0);
                    if (dis1[i] < 1.0) index.push_back(i);
                }
                temp[0] = center[0]; temp[1] = center[1]; temp[2] = center[2];
                center[0] = 0; center[1] = 0; center[2] = 0;
                for (int i = 0; i < index.size(); i++) {
                    center[0] += fw_pr[index[i]][0];
                    center[1] += fw_pr[index[i]][1];
                    center[2] += fw_pr[index[i]][2];
                }
                center[0] /= index.size(); center[1] /= index.size(); center[2] /= index.size();
                double residual = 0;
                for (int i = 0; i < 3; i++) {
                    residual += pow(center[i] - temp[i], 2.0);
                }
                if (sqrt(residual) < 1e-3 * max(hs, hr)) break;
            }
            for (int i = 0; i < 3; i++) clusterCenter[iym - 1][ixm - 1][i] = center[i];
            delete[]iw;
            delete[]dis1;
        }
    }
    //vector<vector<double> > s(2 * height + 1, vector<double>(2 * length + 1, 1));
    MyMatrix2d s(2 * height + 1, 2 * length + 1);
    double** s_pr = s.getPr();
    for (int y = 0; y < 2 * height + 1; y += 2) {
        for (int x = 0; x < 2 * length + 1; x++) {
            s_pr[y][x] = 0;
        }
    }
    for (int y = 0; y < 2 * height + 1; y++) {
        for (int x = 0; x < 2 * length + 1; x += 2) {
            s_pr[y][x] = 0;
        }
    }
    MyMatrix2d cat11(height, length - 1);
    MyMatrix2d cat12(height, length - 1);
    MyMatrix2d cat13(height, length - 1);
    MyMatrix2d cat21(height - 1, length);
    MyMatrix2d cat22(height - 1, length);
    MyMatrix2d cat23(height - 1, length);
    MyMatrix2d all1(height, length - 1);
    MyMatrix2d all2(height - 1, length);
    double** cat11_pr = cat11.getPr(); double** cat12_pr = cat12.getPr();
    double** cat13_pr = cat13.getPr(); double** cat21_pr = cat21.getPr();
    double** cat22_pr = cat22.getPr(); double** cat23_pr = cat23.getPr();
    double** all1_pr = all1.getPr(); double** all2_pr = all2.getPr();
    for (int y = 0; y < height; y++) {
        for (int x = 1; x < length; x++) {
            if (abs(clusterCenter[y][x][0] - clusterCenter[y][x - 1][0]) < hs)
                cat11_pr[y][x - 1] = 1;
            if (abs(clusterCenter[y][x][1] - clusterCenter[y][x - 1][1]) < hs)
                cat12_pr[y][x - 1] = 1;
            if (abs(clusterCenter[y][x][2] - clusterCenter[y][x - 1][2]) < hr)
                cat13_pr[y][x - 1] = 1;
            if (cat11_pr[y][x - 1] + cat12_pr[y][x - 1] + cat13_pr[y][x - 1] == 3)
                all1_pr[y][x - 1] = 1;
        }
    }
    for (int y = 1; y < height; y++) {
        for (int x = 0; x < length; x++) {
            if (abs(clusterCenter[y][x][0] - clusterCenter[y - 1][x][0]) < hs)
                cat21_pr[y - 1][x] = 1;
            if (abs(clusterCenter[y][x][1] - clusterCenter[y - 1][x][1]) < hs)
                cat22_pr[y - 1][x] = 1;
            if (abs(clusterCenter[y][x][2] - clusterCenter[y - 1][x][2]) < hr)
                cat23_pr[y - 1][x] = 1;
            if (cat21_pr[y - 1][x] + cat22_pr[y - 1][x] + cat23_pr[y - 1][x] == 3)
                all2_pr[y - 1][x] = 1;
        }
    }
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < length - 1; x++) {
            s_pr[1 + y * 2][2 + x * 2] = all1_pr[y][x];
        }
    }
    for (int y = 0; y < height - 1; y++) {
        for (int x = 0; x < length; x++) {
            s_pr[2 + 2 * y][1 + 2 * x] = all2_pr[y][x];
        }
    }
    //vector<vector<double> > label_all(2 * height + 1, vector<double>(2 * length + 1, 0));
    MyMatrix2d label_all(2 * height + 1, 2 * length + 1);
    double** label_all_pr = label_all.getPr();
    int nlabel = mybwlabel(s, label_all);
    if (nlabel == -1) {
        mexPrintf("ERROR IN BWLABEL!");
        //cout << "ERROR IN BWLABEL!" << endl;
        return;
    }
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < length; x++) {
            label_pr[y][x] = label_all_pr[1 + 2 * y][1 + 2 * x];
        }
    }
}

void tabulate(MyMatrix2d& im, vector<vector<double> >& result) { //result参数只声明，不要初始化
    unordered_set<double> elements;
    vector<double> element(2, 0);
    double** im_pr = im.getPr();
    //vector<double> ele_list;
    for (int y = 0; y < im.size(0); y++) {
        for (int x = 0; x < im.size(1); x++) {
            if (elements.count(im_pr[y][x]) == 0) {
                elements.insert(im_pr[y][x]);
                element[0] = im_pr[y][x];
                element[1] = 1;
                result.push_back(element);
                //ele_list.push_back(im[y][x]);
            }
            else {
                for (int i = 0; i < result.size(); i++) {
                    if (result[i][0] == im_pr[y][x]) result[i][1] += 1;
                }
            }
        }
    }
}

int occRefine(int y, int x, int w, MyMatrix2d& dis, arguments info, double hs, double hr, double deno) { //坐标为C++坐标，传入MatLab坐标时注意处理
    int key = 0;
    //vector<vector<double> > patch(2 * w + 1, vector<double>(2 * w + 1, 0));
    MyMatrix2d patch = dis.getPatch(y - w, y + w, x - w, x + w);
    double** patch_pr = patch.getPr();
    //getPatch(dis, patch, x - w, x + w, y - w, y + w);
    //vector<vector<double> > v(2 * w + 1, vector<double>(2 * w + 1, 0));
    MyMatrix2d v(2 * w + 1, 2 * w + 1);
    double** v_pr = v.getPr();
    for (int y = 0; y < 2 * w + 1; y++) {
        for (int x = 0; x < 2 * w + 1; x++) {
            v_pr[y][x] = 255 / (info.dis_max - info.dis_min) * patch_pr[y][x] + 255 * info.dis_min / (info.dis_min - info.dis_max);
        }
    }
    //vector<vector<double> > l(2 * w + 1, vector<double>(2 * w + 1, 0));
    MyMatrix2d l(2 * w + 1, 2 * w + 1);
    meanshsegm(v, l, hs, hr);
    vector<vector<double> > sort_initial;
    tabulate(l, sort_initial);
    if (sort_initial.size() == 1) {
        key = 0;
    }
    else {
        vector<int> sort1;
        for (int i = 0; i < sort_initial.size(); i++) sort1.push_back(sort_initial[i][1]);
        sort(sort1.begin(), sort1.end());
        if (sort1[sort1.size() - 2] > v.size(0) * v.size(1) / deno) key = 1;
        else key = 0;
    }
    return key;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    arguments info;
    setStructInfo(info, prhs);
    int dis_rows = (int)mxGetM(dis_m);
    int dis_cols = (int)mxGetN(dis_m);
    int occ_x_rows = (int)mxGetM(occ_x_m);
    int occ_y_rows = (int)mxGetM(occ_y_m);
    int w = mxGetScalar(w_m);
    double hs = mxGetScalar(hs_m);
    double hr = mxGetScalar(hr_m);
    double deno = mxGetScalar(deno_m);
    double* dis_pr = mxGetPr(dis_m);
    MyMatrix2d dis(dis_rows, dis_cols, dis_pr);
    key_m = mxCreateDoubleMatrix(info.Ny, info.Nx, mxREAL);
    double* key_m_pr = mxGetPr(key_m);
    MyMatrix2d key_mat(info.Ny, info.Nx);
    double** key_mat_pr = key_mat.getPr();
    double* key_pr = mxGetPr(key_m);
    double* occ_x_pr = mxGetPr(occ_x_m);
    double* occ_y_pr = mxGetPr(occ_y_m);
    for (int i = 0; i < occ_x_rows; i++) {
        int x_m = threshold(occ_x_pr[i], (double)(info.Nx - w), (double)(1 + w));
        int y_m = threshold(occ_y_pr[i], (double)(info.Ny - w), (double)(1 + w));
        key_mat_pr[y_m - 1][x_m - 1] = (double)occRefine(y_m - 1, x_m - 1, w, dis, info, hs, hr, deno);
    }
    key_mat.Mat2MEXpr(key_m_pr);
}
