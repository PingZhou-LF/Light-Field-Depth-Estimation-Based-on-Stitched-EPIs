#include "mex.h"
#include "math.h"
#include<iostream>
using namespace std;

#define sai      prhs[0]
#define csai     prhs[1]
#define flag_all prhs[2]
#define cost_half_shape prhs[3] //为了提供cost_half所需的Ny*Nx*dis_res的形状，全0
#define dis_step prhs[4]
#define info_m   prhs[5]
#define cost_half01 plhs[0]
#define dis_half01  plhs[1]
//#define args      plhs[2]
//#define indexX plhs[3]
//#define indexY plhs[4]
//#define pixel plhs[3]

struct arguments {
    int Nt, Ns, Nx, Ny, cc, t0, s0, res;
    double alpha1, alpha2, dpth, angle_min, angle_max, step, dis_max, dis_min;
    double* csai;
};
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
int count = 0;

template<typename T> T threshold(T input, T max, T min) {
    if (input > max) input = max;
    if (input < min) input = min;
    return input;
}

double indexCalUp(double k, int x0, int y0, const mxArray* prhs[], arguments info) {
    ::count++;
    int s0 = info.s0;
    int t0 = info.t0;
    int Nt = info.Nt;
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
    double* csai_d = (double*)mxGetPr(csai);
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
            mxArray* sai1 = mxGetCell(sai, (-t_ind + Nt - 1) + (-s_ind + s0 - 1) * Nt);
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

double indexCalDown(double k, int x0, int y0, const mxArray* prhs[], arguments info) {
    ::count++;
    int s0 = info.s0;
    int t0 = info.t0;
    int Nt = info.Nt;
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
    double* csai_d = (double*)mxGetPr(csai);
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
            mxArray* sai1 = mxGetCell(sai, (-t_ind + Nt - 1) + (-s_ind + s0 - 1 + 4) * Nt);
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
    //参数准备
    arguments info;
    setStructInfo(info, prhs);
    size_t num_t = mxGetM(sai);
    size_t num_s = mxGetN(sai);
    double dis_step_val = mxGetScalar(dis_step);
    int dis_res = int((info.dis_max - info.dis_min) / dis_step_val + 1);
    mexPrintf("dis_step=%f", dis_step_val);
    mexPrintf("dis_res=%d", dis_res);
    //数组准备
    double* flags = mxGetPr(flag_all);
    dis_half01 = mxCreateDoubleMatrix((mwSize)info.Ny, (mwSize)info.Nx, mxREAL);
    cost_half01 = mxDuplicateArray(cost_half_shape); //搞不出来3d数组，只能从matlab里传
    double* dis_half01_copy = (double*)mxGetPr(dis_half01);
    double* cost_half01_copy = (double*)mxGetPr(cost_half01);
    mxArray* cost = mxCreateDoubleMatrix((mwSize)dis_res, (mwSize)3, mxREAL);
    // auto cost_eigen = matlab2Eigen<double>(cost);
    double* cost_copy = (double*)mxGetPr(cost);
    double* dis_temp = new double[dis_res];
    for (int i = 0; i < dis_res; i++) {
        dis_temp[i] = info.dis_min + i * dis_step_val;
    }
    double* b = new double[dis_res];
    ::memset(b, 0, sizeof(b));

    //计算
    for (int y = 0; y < info.Ny; y++) {  //y为行，x为列，遍历先行再列，数组列优先存储
        for (int x = 0; x < info.Nx; x++) {
            double flag_xy = flags[x * info.Ny + y];
            if (flag_xy == 1.0) {
                for (int i = 0; i < dis_res; i++) {
                    double slope = 1.0 / dis_temp[i];
                    cost_copy[1 * dis_res + i] = slope;
                    cost_copy[2 * dis_res + i] = dis_temp[i];
                    //cost_copy[0 * dis_res + i] = 1.0;
                    cost_copy[0 * dis_res + i] = indexCalUp(slope, x, y, prhs, info);
                    b[i] = cost_copy[0 * dis_res + i];
                }
            }
            else {
                for (int i = 0; i < dis_res; i++) {
                    double slope = 1.0 / dis_temp[i];
                    cost_copy[1 * dis_res + i] = slope;
                    cost_copy[2 * dis_res + i] = dis_temp[i];
                    //cost_copy[0 * dis_res + i] = 1.0;
                    cost_copy[0 * dis_res + i] = indexCalDown(slope, x, y, prhs, info);
                    b[i] = cost_copy[0 * dis_res + i];
                }
            }
            int I = 0;
            double min_b = b[0];
            for (int i = 0; i < dis_res; i++) {
                if (b[i] < min_b) {
                    I = i;
                    min_b = b[i];
                }
            }
            dis_half01_copy[x * info.Ny + y] = cost_copy[2 * dis_res + I];
            for (int i = 0; i < dis_res; i++) {
                cost_half01_copy[i * info.Nx * info.Ny + x * info.Ny + y] = cost_copy[0 * dis_res + i];
            }

        }
    }
    mexPrintf("total call times = %d", ::count);
    int t_ind = 3;
    int s_ind = 4;
    //mxArray* temp = mxGetCell(sai, (-t_ind + info.Nt - 1) + (-s_ind + info.s0 - 1 + 4) * info.Nt);
    ////mxArray* temp = mxGetCell(sai, 6 + 1 * 9);
    //sai_xy = mxDuplicateArray(temp);
    delete[]dis_temp;
    delete[]b;

    //const int s0 = 5;
    //const int Nt = 9;
    //int radius = info.t0 - 1;
    //double* index_x = new double[s0 * Nt];
    //::memset(index_x, 0, sizeof(index_x));
    //double* index_y = new double[s0 * Nt];
    //::memset(index_y, 0, sizeof(index_y));
    //double* index_x_temp = new double[s0];
    //::memset(index_x_temp, 0, sizeof(index_x_temp));
    //double* index_y_temp = new double[s0];
    //::memset(index_y_temp, 0, sizeof(index_y_temp));
    //double* temp = new double[3];
    //memset(temp, 0, sizeof(temp));
    //double* temp1 = new double[3];
    //memset(temp1, 0, sizeof(temp1));
    //double* temp2 = new double[3];
    //memset(temp2, 0, sizeof(temp2));
    //double* temp3 = new double[3];
    //memset(temp3, 0, sizeof(temp3));
    //double* temp4 = new double[3];
    //memset(temp4, 0, sizeof(temp4));

    //int x0 = 107;
    //int y0 = 184;
    //double* center = new double[3];
    //double* csai_d = (double*)mxGetPr(csai);
    //for (int i = 0; i < 3; i++) {
    //    center[i] = csai_d[i * 512 * 512 + x0 * 512 + y0];
    //}
    //double k = 1 / dis_temp[1];
    //for (int delta_t = -radius; delta_t <= radius; delta_t++) {
    //    double index_x0 = x0 + (double)delta_t * s0 / k;
    //    double index_y_ = y0 + delta_t / k;
    //    for (int i = -s0 + 1; i <= 0; i++) {
    //        index_x_temp[i + s0 - 1] = index_x0 + i / k;
    //        index_y_temp[i + s0 - 1] = index_y_;
    //    }
    //    for (int i = 0; i < s0; i++) {
    //        index_x[i + s0 * (delta_t + radius)] = index_x_temp[i] + 1;
    //        index_y[i + s0 * (delta_t + radius)] = index_y_temp[i] + 1;
    //    }
    //}
    //double x_ind = index_x[t_ind * s0 + s_ind] - (t_ind + 1 - info.t0) * s0 / k;
    //double y_ind = index_y[t_ind * s0 + s_ind];
    //int x_floor = (int)x_ind;
    //x_floor = threshold(x_floor, info.Nx, 1);
    //int x_ceil = x_floor < x_ind ? x_floor + 1 : x_ind;
    //x_ceil = threshold(x_ceil, info.Nx, 1);
    //int y_floor = (int)y_ind;
    //y_floor = threshold(y_floor, info.Ny, 1);
    //int y_ceil = y_floor < y_ind ? y_floor + 1 : y_ind;
    //y_ceil = threshold(y_ceil, info.Ny, 1);
    //args = mxCreateDoubleMatrix((mwSize)12, (mwSize)1, mxREAL);
    ////indexX = mxCreateDoubleMatrix((mwSize)s0 * Nt, (mwSize)1, mxREAL);
    ////indexY = mxCreateDoubleMatrix((mwSize)s0 * Nt, (mwSize)1, mxREAL);
    ////pixel = mxCreateDoubleMatrix((mwSize)6, (mwSize)1, mxREAL);
    ////double* indexX_pr = (double*)mxGetPr(indexX);
    ////double* indexY_pr = (double*)mxGetPr(indexY);
    //double* args_pr = (double*)mxGetPr(args);
    ////double* pixel_pr = (double*)mxGetPr(pixel);
    //double dis_x_c = x_ceil - x_ind;
    //double dis_x_f = 1 - dis_x_c;
    //double dis_y_c = y_ceil - y_ind;
    //double dis_y_f = 1 - dis_y_c;
    //args_pr[0] = x_ceil; args_pr[1] = x_floor; args_pr[2] = y_ceil; args_pr[3] = y_floor;
    //args_pr[4] = k; args_pr[5] = x_ind; args_pr[6] = y_ind; args_pr[7] = dis_x_c;
    //args_pr[8] = dis_x_f; args_pr[9] = dis_y_c; args_pr[10] = dis_y_f; args_pr[11] = flags[y0 * info.Nx + x0];
    ///*for (int i = 0; i < s0 * Nt; i++){
    //    indexX_pr[i] = index_x[i];
    //    indexY_pr[i] = index_y[i];
    //}*/
    ///*mxArray* sai2 = mxGetCell(sai, (-t_ind + Nt - 1) + (-s_ind + s0 - 1 + 4) * Nt);
    //double* sai2_pr = (double*)mxGetPr(sai2);
    //for (int i = 0; i < 3; i++) {
    //    pixel_pr[i] = sai2_pr[i * 512 * 512 + (x_floor - 1) * 512 + y_floor - 1];
    //}
    //for (int i = 0; i < 3; i++) {
    //    pixel_pr[i + 3] = center[i];
    //}
    //double cost_ = 0;
    //for (int i = 0; i < 3; i++) {
    //    temp1[i] = sai2_pr[i * 512 * 512 + (x_floor - 1) * 512 + y_floor - 1] * dis_y_c * dis_x_c;
    //    temp2[i] = sai2_pr[i * 512 * 512 + (x_ceil - 1) * 512 + y_floor - 1] * dis_y_c * dis_x_f;
    //    temp3[i] = sai2_pr[i * 512 * 512 + (x_floor - 1) * 512 + y_ceil - 1] * dis_y_f * dis_x_c;
    //    temp4[i] = sai2_pr[i * 512 * 512 + (x_ceil - 1) * 512 + y_ceil - 1] * dis_y_f * dis_x_f;
    //    temp[i] = temp1[i] + temp2[i] + temp3[i] + temp4[i];
    //    cost_ += abs(temp[i] - center[i]) / 3;
    //}
    //args_pr[11] = cost_;*/
    //
    //delete[]index_x;
    //delete[]index_y;
    //delete[]index_x_temp;
    //delete[]index_y_temp;
    //delete[]temp;
    //delete[]temp1;
    //delete[]temp2;
    //delete[]temp3;
    //delete[]temp4;
    //delete[]center;
}

