#include "mex.h"
#include "MyMatrix2d.h"
#include<vector>
#include<math.h>
#include<algorithm>
using namespace std;

#define Nt_m  prhs[0]        //int
#define Ns_m  prhs[1]        //int
#define dis_m prhs[2]        //512*512
#define info_m  prhs[3]        //struct
#define maskRight_m  prhs[4]  //512*512
#define kRes_m prhs[5]       //int
#define sai_m  prhs[6]       //9*9*512*512

#define disRefine_m plhs[0]    //512*512

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

void calRange(double k, double b, int n, vector<double>& nqpsk1k2) {
    vector<int> xd(n + 1, 0);
    vector<double> yd(n + 1, 0);
    vector<int> cy(n + 1, 0);
    vector<int> Chain(n, 0);
    for (int i = 0; i <= n; i++) {
        xd[i] = i;
        yd[i] = k * xd[i] + b;
        cy[i] = floor(yd[i]);
    }
    for (int i = 0; i < n; i++) Chain[i] = cy[i + 1] - cy[i];
    int q = 0;
    for (q = 1; q <= n; q++) {
        int flag = 1;
        for (int i = 0; i < n - q; i++) {
            if (Chain[i + q] != Chain[i]) {
                flag = 0;
                break;
            }
        }
        if (flag) break;
    }
    q = min(q, n);
    int p = 0;
    for (int i = 0; i < q; i++) {
        p += Chain[i];
    }
    int s = 0;
    for (s = 0; s <= q - 1; s++) {
        int flag = 1;
        for (int i = 1; i <= q; i++) {
            double temp = floor(p * (i - s) / q) - floor(p * (i - s - 1) / q);
            if (Chain[i - 1] != temp) flag = 0;
        }
        if (flag) break;
    }
    if (s == q) s--;
    int l = 0;
    for (l = 0; l <= q - 1; l++) {
        if (fmod(l * p + 1, q) == 0) break;
    }
    int Fs = s;
    int Ls = s + floor((n - s) / q) * q;
    int Fs1 = s + l - floor((s + l) / q) * q;
    int Ls1 = s + l + floor((n - s - l) / q) * q;
    double pq = (p * 1.0) / (q * 1.0);
    double k1 = pq - 1.0 / (q * (Ls - Fs1));
    double k2 = pq + 1.0 / (q * (Ls1 - Fs));
    nqpsk1k2[0] = n; nqpsk1k2[1] = q; nqpsk1k2[2] = p; nqpsk1k2[3] = s;
    nqpsk1k2[4] = k1; nqpsk1k2[5] = k2;
}

void calRangek(double k0, double b0, int n, double& kmin, double& kmax) {
    vector<double> nqpsk1k2(6, 0);
    //double kmin = 0, kmax = 0;
    if (k0 >= 1) {
        double k = 1 / k0;
        double b = -b0 / k0;
        calRange(k, b, n, nqpsk1k2);
        kmin = 1.0 / nqpsk1k2[5];
        kmax = 1.0 / nqpsk1k2[4];
    }
    else {
        if (k0 > 0) {
            calRange(k0, b0, n, nqpsk1k2);
            kmin = nqpsk1k2[4];
            kmax = nqpsk1k2[5];
        }
        else {
            if (k0 >= -1) {
                calRange(-k0, b0, n, nqpsk1k2);
                kmin = -1 * nqpsk1k2[5];
                kmax = -1 * nqpsk1k2[4];
            }
            else {
                double k = -1 / k0;
                double b = b0 / k0;
                calRange(k, b, n, nqpsk1k2);
                kmin = -1 / nqpsk1k2[4];
                kmax = -1 / nqpsk1k2[5];
            }
        }
    }
    if (isinf(kmin)) kmin = k0;
    if (isinf(kmax)) kmax = k0;
}

double indexCalUp(double k, int x0, int y0, arguments info, const mxArray* prhs[]) {
    const int s0 = 5;
    const int Nt = 9;
    int radius = info.t0 - 1;
    double* index_x = new double[s0 * Nt];
    ::memset(index_x, 0, sizeof(index_x));
    double* index_y = new double[s0 * Nt];
    ::memset(index_y, 0, sizeof(index_y));
    double* cost = new double[s0 * Nt];
    ::memset(cost, 0, sizeof(cost));
    double* index_x_temp = new double[s0];
    ::memset(index_x_temp, 0, sizeof(index_x_temp));
    double* index_y_temp = new double[s0];
    ::memset(index_y_temp, 0, sizeof(index_y_temp));
    double* center = new double[3];
    // mxArray* csai_copy = mxDuplicateArray(csai);
    double* csai_d = info.csai;
    for (int i = 0; i < 3; i++) {
        center[i] = csai_d[i * info.Nx * info.Ny + x0 * info.Ny + y0];
    }
    double* temp = new double[3];
    ::memset(temp, 0, sizeof(temp));
    double* temp1 = new double[3];
    ::memset(temp1, 0, sizeof(temp1));
    double* temp2 = new double[3];
    ::memset(temp2, 0, sizeof(temp2));
    double* temp3 = new double[3];
    ::memset(temp3, 0, sizeof(temp3));
    double* temp4 = new double[3];
    ::memset(temp4, 0, sizeof(temp4));

    for (int delta_t = -radius; delta_t <= radius; delta_t++) {
        double index_x0 = x0 + (double)delta_t * s0 / k;
        double index_y_ = y0 + (double)delta_t / k;
        for (int i = 0; i < s0; i++) {
            index_x_temp[i] = index_x0 + (double)i / k;
            index_y_temp[i] = index_y_;
        }
        for (int i = 0; i < s0; i++) {
            index_x[i + s0 * (delta_t + radius)] = index_x_temp[i] + 1;
            index_y[i + s0 * (delta_t + radius)] = index_y_temp[i] + 1;
        }
    }
    for (int t_ind = 0; t_ind < Nt; t_ind++) {
        for (int s_ind = 0; s_ind < s0; s_ind++) {
            double x_ind = index_x[t_ind * s0 + s_ind] - (double)(t_ind + 1 - info.t0) * s0 / k;
            double y_ind = index_y[t_ind * s0 + s_ind];
            int x_floor = (int)x_ind;
            x_floor = threshold(x_floor, info.Nx, 1);
            int x_ceil = x_floor < x_ind ? x_floor + 1 : x_ind;
            x_ceil = threshold(x_ceil, info.Nx, 1);
            int y_floor = (int)y_ind;
            y_floor = threshold(y_floor, info.Ny, 1);
            int y_ceil = y_floor < y_ind ? y_floor + 1 : y_ind;
            y_ceil = threshold(y_ceil, info.Ny, 1);
            double dis_x_c = x_ceil - x_ind;
            double dis_x_f = 1.0 - dis_x_c;
            double dis_y_c = y_ceil - y_ind;
            double dis_y_f = 1.0 - dis_y_c;
            mxArray* sai1 = mxGetCell(sai_m, (-t_ind + Nt - 1) + (-s_ind + s0 - 1) * Nt);
            double* sai_ = (double*)mxGetPr(sai1);
            double cost_ = 0;
            for (int i = 0; i < 3; i++) {
                temp1[i] = sai_[i * info.Nx * info.Ny + (x_floor - 1) * info.Ny + y_floor - 1] * dis_y_c * dis_x_c;
                temp2[i] = sai_[i * info.Nx * info.Ny + (x_ceil - 1) * info.Ny + y_floor - 1] * dis_y_c * dis_x_f;
                temp3[i] = sai_[i * info.Nx * info.Ny + (x_floor - 1) * info.Ny + y_ceil - 1] * dis_y_f * dis_x_c;
                temp4[i] = sai_[i * info.Nx * info.Ny + (x_ceil - 1) * info.Ny + y_ceil - 1] * dis_y_f * dis_x_f;
                temp[i] = temp1[i] + temp2[i] + temp3[i] + temp4[i];
                cost_ += abs(temp[i] - center[i]) / 3.0;
            }
            cost[t_ind * s0 + s_ind] = cost_ * cost_;
        }
    }
    double cost_sum = 0;
    for (int i = 0; i < s0 * Nt; i++) {
        cost_sum += cost[i];
    }
    delete[]index_x;
    delete[]index_y;
    delete[]index_x_temp;
    delete[]index_y_temp;
    delete[]cost;
    delete[]temp;
    delete[]temp1;
    delete[]temp2;
    delete[]temp3;
    delete[]temp4;
    return cost_sum;
}

double indexCalDown(double k, int x0, int y0, arguments info, const mxArray* prhs[]) {
    const int s0 = 5;
    const int Nt = 9;
    int radius = info.t0 - 1;
    double* index_x = new double[s0 * Nt];
    memset(index_x, 0, sizeof(index_x));
    double* index_y = new double[s0 * Nt];
    memset(index_y, 0, sizeof(index_y));
    double* cost = new double[s0 * Nt];
    memset(cost, 0, sizeof(cost));
    double* index_x_temp = new double[s0];
    memset(index_x_temp, 0, sizeof(index_x_temp));
    double* index_y_temp = new double[s0];
    memset(index_y_temp, 0, sizeof(index_y_temp));
    double* center = new double[3];
    // mxArray* csai_copy = mxDuplicateArray(csai);
    double* csai_d = info.csai;
    for (int i = 0; i < 3; i++) {
        center[i] = csai_d[i * info.Nx * info.Ny + x0 * info.Ny + y0];
    }
    double* temp = new double[3];
    memset(temp, 0, sizeof(temp));
    double* temp1 = new double[3];
    memset(temp1, 0, sizeof(temp1));
    double* temp2 = new double[3];
    memset(temp2, 0, sizeof(temp2));
    double* temp3 = new double[3];
    memset(temp3, 0, sizeof(temp3));
    double* temp4 = new double[3];
    memset(temp4, 0, sizeof(temp4));
    for (int delta_t = -radius; delta_t <= radius; delta_t++) {
        double index_x0 = x0 + (double)delta_t * s0 / k;
        double index_y_ = y0 + (double)delta_t / k;
        for (int i = -s0 + 1; i <= 0; i++) {
            index_x_temp[i + s0 - 1] = index_x0 + (double)i / k;
            index_y_temp[i + s0 - 1] = index_y_;
        }
        for (int i = 0; i < s0; i++) {
            index_x[i + s0 * (delta_t + radius)] = index_x_temp[i] + 1;
            index_y[i + s0 * (delta_t + radius)] = index_y_temp[i] + 1;
        }
    }
    for (int t_ind = 0; t_ind < Nt; t_ind++) {
        for (int s_ind = 0; s_ind < s0; s_ind++) {
            double x_ind = index_x[t_ind * s0 + s_ind] - (double)(t_ind + 1 - info.t0) * s0 / k;
            double y_ind = index_y[t_ind * s0 + s_ind];
            int x_floor = (int)x_ind;
            x_floor = threshold(x_floor, info.Nx, 1);
            int x_ceil = x_floor < x_ind ? x_floor + 1 : x_ind;
            x_ceil = threshold(x_ceil, info.Nx, 1);
            int y_floor = (int)y_ind;
            y_floor = threshold(y_floor, info.Ny, 1);
            int y_ceil = y_floor < y_ind ? y_floor + 1 : y_ind;
            y_ceil = threshold(y_ceil, info.Ny, 1);
            double dis_x_c = x_ceil - x_ind;
            double dis_x_f = 1.0 - dis_x_c;
            double dis_y_c = y_ceil - y_ind;
            double dis_y_f = 1.0 - dis_y_c;
            mxArray* sai1 = mxGetCell(sai_m, (-t_ind + Nt - 1) + (-s_ind + s0 - 1 + 4) * Nt);
            double* sai_ = (double*)mxGetPr(sai1);
            double cost_ = 0;
            for (int i = 0; i < 3; i++) {
                temp1[i] = sai_[i * info.Nx * info.Ny + (x_floor - 1) * info.Ny + y_floor - 1] * dis_y_c * dis_x_c;
                temp2[i] = sai_[i * info.Nx * info.Ny + (x_ceil - 1) * info.Ny + y_floor - 1] * dis_y_c * dis_x_f;
                temp3[i] = sai_[i * info.Nx * info.Ny + (x_floor - 1) * info.Ny + y_ceil - 1] * dis_y_f * dis_x_c;
                temp4[i] = sai_[i * info.Nx * info.Ny + (x_ceil - 1) * info.Ny + y_ceil - 1] * dis_y_f * dis_x_f;
                temp[i] = temp1[i] + temp2[i] + temp3[i] + temp4[i];
                cost_ += abs(temp[i] - center[i]) / 3.0;
            }
            cost[t_ind * s0 + s_ind] = cost_ * cost_;
        }
    }
    double cost_sum = 0;
    for (int i = 0; i < s0 * Nt; i++) {
        cost_sum += cost[i];
    }
    delete[]index_x;
    delete[]index_y;
    delete[]index_x_temp;
    delete[]index_y_temp;
    delete[]cost;
    delete[]temp;
    delete[]temp1;
    delete[]temp2;
    delete[]temp3;
    delete[]temp4;
    return cost_sum;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    arguments info;
    setStructInfo(info, prhs);
    int Nt = mxGetScalar(Nt_m);
    int Ns = mxGetScalar(Ns_m);
    int kRes = mxGetScalar(kRes_m);
    double Rs = (Ns + 1) / 2;
    MyMatrix2d maskLeft((int)mxGetM(maskRight_m), (int)mxGetN(maskRight_m), mxGetPr(maskRight_m));
    vector<int> RightY, RightX;
    double** maskLeft_pr = maskLeft.getPr();
    MyMatrix2d dis((int)mxGetM(dis_m), (int)mxGetN(dis_m), mxGetPr(dis_m));
    double** dis_pr = dis.getPr();
    for (int x = 0; x < maskLeft.size(1); x++) {
        for (int y = 0; y < maskLeft.size(0); y++) {
            if (maskLeft_pr[y][x] == 1) {
                RightX.push_back(x);
                RightY.push_back(y);
            }
        }
    }
    disRefine_m = mxDuplicateArray(dis_m);
    MyMatrix2d disRefine((int)mxGetM(dis_m), (int)mxGetN(dis_m), mxGetPr(disRefine_m));
    double** disRefine_pr = disRefine.getPr();
    for (int i = 0; i < RightX.size(); i++) {
        int x = RightX[i];
        int y = RightY[i];
        double k0 = 1.0 / dis_pr[y][x];
        double b0 = (Nt * Rs * 1.0 + 1.0) / 2.0 - k0 * (x + 1); //保持delta与matlab一致
        double kmin = 0, kmax = 0;
        calRangek(k0, b0, Nt * Rs, kmin, kmax);
        vector<double> costDelta(201, 10000);
        int costInd = 0;
        for (double deltaTemp = -0.1; deltaTemp <= 0.1; deltaTemp += 0.001) {
            vector<double> kTemp(5, 0);
            vector<double> xTemp(5, 0);
            for (int i = 0; i < 5; i++) {
                kTemp[i] = k0 / (1 + i * deltaTemp * k0);
                xTemp[i] = x + i;
            }
            bool notAllZero = true, noNan = true;
            int count0 = 0;
            for (int i = 0; i < 5; i++) {
                if (kTemp[i] == 0) count0++;
                if (isnan(kTemp[i])) noNan = false;
            }
            if (count0 == 5) notAllZero = false;
            if (xTemp[0] < info.Nx && notAllZero && noNan) {
                double costTemp1 = indexCalDown(kTemp[0], xTemp[0], y, info, prhs);
                double costTemp2 = indexCalDown(kTemp[1], xTemp[1], y, info, prhs);
                double costTemp3 = indexCalDown(kTemp[2], xTemp[2], y, info, prhs);
                double costTemp4 = indexCalDown(kTemp[3], xTemp[3], y, info, prhs);
                double costTemp5 = indexCalDown(kTemp[4], xTemp[4], y, info, prhs);
                /*double costTemp1 = 1;
                double costTemp2 = 1;
                double costTemp3 = 1;
                double costTemp4 = 1;
                double costTemp5 = 1;*/
                costDelta[costInd] = costTemp1 + costTemp2 + costTemp3 + costTemp4 + costTemp5;
                costInd++;
            }
            else {
                costDelta[costInd] = 10000;
                costInd++;
            }
        }
        int deltaInd = min_element(costDelta.begin(), costDelta.end()) - costDelta.begin();
        double delta = -0.1 + 0.001 * deltaInd; //保持delta与matlab一致
        vector<double> cost(kRes + 1, 10000);
        costInd = 0;
        double kstep = (kmax - kmin) / (kRes - 1);
        for (int indK = 0; indK <= kRes; indK++) {
            double k = kmin + indK * kstep;
            vector<double> kTemp(5, 0);
            vector<double> xTemp(5, 0);
            for (int i = 0; i < 5; i++) {
                kTemp[i] = k / (1 + i * delta * k);
                xTemp[i] = x + i;
            }
            bool notAllZero = true, noNan = true;
            int count0 = 0;
            for (int i = 0; i < 5; i++) {
                if (kTemp[i] == 0) count0++;
                if (isnan(kTemp[i])) noNan = false;
            }
            if (count0 == 5) notAllZero = false;
            if (xTemp[0] < info.Nx && notAllZero && noNan) {
                double costTemp1 = indexCalDown(kTemp[0], xTemp[0], y, info, prhs);
                double costTemp2 = indexCalDown(kTemp[1], xTemp[1], y, info, prhs);
                double costTemp3 = indexCalDown(kTemp[2], xTemp[2], y, info, prhs);
                double costTemp4 = indexCalDown(kTemp[3], xTemp[3], y, info, prhs);
                double costTemp5 = indexCalDown(kTemp[4], xTemp[4], y, info, prhs);
                /*double costTemp1 = 1;
                double costTemp2 = 1;
                double costTemp3 = 1;
                double costTemp4 = 1;
                double costTemp5 = 1;*/
                cost[costInd] = costTemp1 * exp(-1) + costTemp2 * exp(-9.0 / 16.0) + costTemp3 * exp(-4.0 / 16.0) + costTemp4 * exp(-1.0 / 16.0) + costTemp5;
                costInd++;
            }
            else {
                cost[costInd] = 10000;
                costInd++;
            }
        }
        costInd = min_element(cost.begin(), cost.end()) - cost.begin();
        double k1 = kmin + kstep * costInd; //保持k1与matlab一致
        disRefine_pr[y][x] = threshold(1.0 / k1, info.dis_max, info.dis_min);
    }
    disRefine.Mat2MEXpr(mxGetPr(disRefine_m));
}
