#include "mex.h"
#include "math.h"
#include <algorithm>

using namespace std;


// 输入
#define	Im_in_remap         prhs[0]  // the remap image
#define	Nx_in      			prhs[1]  // 
#define	Ny_in               prhs[2]  // 
#define	Nt_in               prhs[3]  // 
#define	Ns_in               prhs[4]  // 
#define slope_in            prhs[5]  //
#define ChosenPtsMask       prhs[6]  // 1: 优化 0：不优化
#define Slope_Res           prhs[7]
#define Rlabel              prhs[8]
#define e_in                prhs[9]    


//输出
#define refined_slope       plhs[0]  // slope cost with sx-epis

//全局参数
int Nt, Ns, Nst, Nx, Ny, pixelNum, remap_Ny, remap_Nx, remap_pixelNum;
double* slope_cost;
double* re_slope;
double* slope;
double res;
double e;

int checkX(int x) {
    if (x < 0) return 0;
    else if (x >= Nx) return Nx - 1;
    else return x;
}
int checkY(int y) {
    if (y < 0) return 0;
    else if (y >= Ny) return Ny - 1;
    else return y;
}
double checkSign(double k0, double k1) {
    if ((k0 > 0 && k1 > 0) || (k0 < 0 && k1 < 0))
    {
        return k1;
    }
    if (k0 > 0 && k1 < 0)
    {
        return 0.001;
    }
    if (k0 < 0 && k1>0)
    {
        return -0.001;
    }


}
int checkST(int x) {
    if (x < 0) return 0;
    else if (x >= Nst) return Nst - 1;
    else return x;
}
double getNorm2(double array[], int size)
{
    double mean = 0.0;
    int num = 0;
    for (int i = 0; i < size; i++) {
        if (array[i] >= 0) {
            mean += array[i];
            num++;
        }
    }
    mean /= num;
    double var = 0;
    for (int i = 0; i < size; i++) {
        if (array[i] >= 0) {
            var += (array[i] - mean) * (array[i] - mean);
        }
    }
    if (num == 1) var = INT_MAX;
    else var /= (num - 1.0);
    return var;
}
double getNorm1(double array[], int size)
{
    double center = array[(size + 1) / 2];

    double norm1 = 0;
    for (int i = 0; i < size; i++) {
        norm1 += fabs(array[i] - center);
    }
    return norm1;
}
double ComputeDelta(double initial_slope[], int x, int y, double ChosenPts[], int ChosenPtsNum) {

    double delta;
    double mincost = 1000000;
    double newcost;
    double bestdelta = 0;
    double k = initial_slope[x * Ny + y];
    ;    for (delta = -0.1; delta <= 0.1; delta += 0.001) {
        newcost = 0.0;
        for (int i = 0; i < ChosenPtsNum; i++) {
            newcost += fabs((k / (1.0 + k * delta * (ChosenPts[i] - x))) - initial_slope[Ny * (int)ChosenPts[i] + y]);
        }
        if (newcost < mincost) {
            bestdelta = delta;
            mincost = newcost;
        }
    }
    return bestdelta;
}
void kRange(double k0, double* k1, double* k2, double beta) {

    int flag, p, q, l, s;
    int x = 0;
    int num_chain = Nst - 1;
    int* chain = new int[num_chain];

    for (int i = 0; i < num_chain; i++) {
        chain[i] = floor(k0 * (x + 1.0) + beta + 0.5) - floor(k0 * x + beta + 0.5);
        x++;
    }

    // q
    for (q = 1; q <= num_chain; q++) { // i 是待选周期
        flag = 1;
        for (int j = 1; j <= num_chain - q; j++) { // 遍历的坐标
            if (chain[j - 1] != chain[j - 1 + q]) {
                flag = 0;
                break;
            }
        }
        if (flag) {
            break;
        }
    }

    // p
    p = 0;
    for (int i = 1; i <= q; i++) {
        p += chain[i - 1];
    }

    // s
    int temp;
    for (s = 0; s <= q - 1; s++) {
        flag = 1;
        for (int i = 1; i <= q; i++) {
            temp = floor(p * (i - s) / q) - floor(p * (i - s - 1) / q);
            if (chain[i - 1] != temp) {
                flag = 0;
                break;
            }
        }

        if (flag)
            break;
    }

    // l
    for (l = 0; l <= q - 1; l++) {
        flag = 0;
        if (fmod(l * p + 1, q) == 0) {
            flag = 1;
            break;
        }
    }

    double Fs = s;
    double Ls = s + floor((num_chain - s) / q) * q;

    double Fsl = s + l - floor((s + l) / q) * q;
    double Lsl = s + l + floor((num_chain - s - l) / q) * q;

    *k1 = (double)p / q - 1.0 / (q * (Ls - Fsl));
    *k2 = (double)p / q + 1.0 / (q * (Lsl - Fs));

}
void sEPICost(double im_in_remap[], int x, int y, double k, double remap[], double* Ix, double* Iy, double* pur)
{
    double* mask = new double[Nst * Nx]();
    //
    for (int j = -(Nt - 1) / 2; j <= (Nt - 1) / 2; j++) { //j=t
        int t = -j;
        double y_ind = (double)t / k + y; //  每张单极图对应的y
        int   y_floor = floor(y_ind); // 以下是为了插出y_ind的所需数据
        int   y_1 = checkY(y_floor);
        int   y_2 = checkY(y_floor + 1);
        double y_1_w = 1 - (y_ind - y_floor); // 权重
        double y_2_w = 1 - y_1_w;

        for (int i = -(Ns - 1) / 2; i <= (Ns - 1) / 2; i++) { // i=s
            int s = -i;
            double x_ind = round(x + (t * Ns + s) / k) - Ns * t / k;
            int   x_floor = floor(x_ind);
            int   x_1 = checkX(x_floor);
            int   x_2 = checkX(x_floor + 1);
            double x_1_w = 1 - (x_ind - x_floor);
            double x_2_w = 1 - x_1_w;

            int x_1_index = i + (Ns - 1) / 2 + (x_1)*Ns;
            int y_1_index = j + (Nt - 1) / 2 + (y_1)*Nt;
            int x_2_index = i + (Ns - 1) / 2 + (x_2)*Ns;
            int y_2_index = j + (Nt - 1) / 2 + (y_2)*Nt;

            for (int c = 0; c < 3; c++) {
                remap[(i + (Ns - 1) / 2) * Ns + (j + (Nt - 1) / 2) + c * Nst] =
                    y_1_w * x_1_w * im_in_remap[y_1_index + x_1_index * remap_Ny + c * remap_pixelNum] +
                    y_2_w * x_1_w * im_in_remap[y_2_index + x_1_index * remap_Ny + c * remap_pixelNum] +
                    y_1_w * x_2_w * im_in_remap[y_1_index + x_2_index * remap_Ny + c * remap_pixelNum] +
                    y_2_w * x_2_w * im_in_remap[y_2_index + x_2_index * remap_Ny + c * remap_pixelNum];
            }

            int mask_y = -j * Ns - i + (Nst - 1) / 2;
            int mask_x = checkX(round(x + (t * Ns + s) / k));
            if (mask[mask_y + mask_x * Nst] == 0)
            {
                mask[mask_y + mask_x * Nst] = 1;
                (*pur)++;
            }
        }
    }
    delete[] mask;
}
void computeCost(double im_in_remap[], double kmin, double kmax, int y0, int x0, double xIdx[], int ChosenPtsNum, double DeltaDis, double maxdis)
{

    double* remap = new double[Nst * 3.0];
    double* Ix = new double[Nst * 3.0];
    double* Iy = new double[Nst * 3.0];
    double* w_dis = new double[1];
    double* pur = new double[1];

    double k0, ki;
    double kstep = (kmax - kmin) / ((double)res - 1.0);

    double norm1, norm2, norm3;
    double vals1, vals2, vals3;
    int xi;

    for (int j = 0; j < res; j++) {// j是斜率的索引
        k0 = kmin + j * kstep; //中心直线的斜率，和斜率索引有关
        vals1 = 0;
        //vals2 = 0;
        //vals3 = 0;
        *pur = 0;
        for (int i = 0; i < ChosenPtsNum; i++) {
            xi = xIdx[i];
            // 计算当前直线i 是否和（Curry，Currx）处于同一个超像素，并且Confi=1
            ki = k0 / (1.0 + k0 * DeltaDis * (xi - x0)); //计算当前直线的斜率
            sEPICost(im_in_remap, xi, y0, ki, remap, Ix, Iy, pur); //计算当前直线沿当前斜率方向的cost
            norm1 = getNorm1(remap, Nst) + getNorm1(&remap[Nst], Nst) + getNorm1(&remap[Nst * 2], Nst);
            vals1 += exp(-(xi - x0)* (xi - x0) / maxdis) * norm1 / 3.0;

        }

        //slope_cost[x0 * Ny + y0 + j * pixelNum] = ((*pur) / ((double)Nst * ChosenPtsNum)) * vals1;
        slope_cost[x0 * Ny + y0 + j * pixelNum] =  vals1;
    }
    delete[] remap;
    delete[] Ix;
    delete[] Iy;
    delete[] w_dis;
    delete[] pur;
}
int getMinIdx(int y, int x)
{
    int minIdx = 0;
    double minVal = INT_MAX;
    for (int i = 0; i < res; i++) {
        double sumVal = slope_cost[x * Ny + y + i * pixelNum];
        if (sumVal < minVal) {
            minVal = sumVal;
            minIdx = i;
        }
    }
    return minIdx;
}


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    double* im_in_remap_pt = mxGetPr(Im_in_remap);
    double* Nx_pt = mxGetPr(Nx_in);
    double* Ny_pt = mxGetPr(Ny_in);
    double* Nt_pt = mxGetPr(Nt_in);
    double* Ns_pt = mxGetPr(Ns_in);
    double* ChosenPtsMask_pt = mxGetPr(ChosenPtsMask);
    double* Slope_Res_pt = mxGetPr(Slope_Res);
    double* Rlabel_pt = mxGetPr(Rlabel);
    double* e_pt = mxGetPr(e_in);
    slope = mxGetPr(slope_in);

    /* 设置参数*/
    res = *Slope_Res_pt;
    e = *e_pt;
    Nt = *Nt_pt;
    Ns = *Ns_pt;
    Nst = Ns * Nt;

    Ny = *Ny_pt;
    Nx = *Nx_pt;
    pixelNum = Ny * Nx;

    remap_Ny = Ny * Nt;
    remap_Nx = Nx * Ns;
    remap_pixelNum = remap_Ny * remap_Nx;

    // Mine
    // 构建输出
    const mwSize dims2[] = { Ny,Nx };
    refined_slope = mxCreateNumericArray(2, dims2, mxDOUBLE_CLASS, mxREAL);
    re_slope = mxGetPr(refined_slope);

    // Mine
    double InitialSlope, SlopeMin,SlopeMax;
    double M, N;
    double* MMin = new double[1];
    double* MMax = new double[1];

    double DeltaDis;
    double label;
    int minIdx;
    int ChosenPtsNum;
    double maxdis;
    int ptIdx;

    slope_cost = new double[Ny * Nx * res];
    for (int j = 0; j < Ny; j++) {
       for (int i = 0; i < Nx; i++) { //对于像素点(j,x)  j-->y   i-->x    
            // 当前点要不要优化，不要优化直接跳到下一个
            if (ChosenPtsMask_pt[i * Ny + j] == 0) {
                re_slope[i * Ny + j] = slope[i * Ny + j];
                continue;
            }

            //计算要参与优化的点的坐标  如果没有要优化的点 直接跳到下一个
            label = Rlabel_pt[i * Ny + j];
            ChosenPtsNum = 0; //选中点的数量(不包含(j,i)本身)
            maxdis = 0; //距离

            for (int pti = 0; pti < Nx; pti++) {
                if ((Rlabel_pt[pti * Ny + j] == label) && (ChosenPtsMask_pt[pti * Ny + j] == 1) && (pti != i))
                {
                    ChosenPtsNum++;
                }
            }

            if (ChosenPtsNum == 0) { // 不优化
                re_slope[i * Ny + j] = slope[i * Ny + j];
                continue;
            }

            ChosenPtsNum++; //把(j,i)本身加上
            double* xIdx = new double[ChosenPtsNum];
            ptIdx = -1;
            for (int pti = 0; pti < Nx; pti++) {
                if ((Rlabel_pt[pti * Ny + j] == label) && (ChosenPtsMask_pt[pti * Ny + j] == 1))
                {
                    ptIdx++;
                    xIdx[ptIdx] = pti;
                    if ((pti - i)* (pti - i) > maxdis)
                    {
                        maxdis = (pti - i) * (pti - i);
                    }
                }
            }
            SlopeMin = 100000;
            SlopeMax = -100000;
            // 计算可变的斜率范围

                InitialSlope = slope[i* Ny + j]; // 点(j,i)初始斜率
                if (InitialSlope >= 1) {
                    M = 1.0 / InitialSlope; // 离散图中： 初始斜率
                    N = (double)i;
                    kRange(M, MMin, MMax, N); //离散图中：斜率范围
                    SlopeMin = 1.0 / *MMax;
                    SlopeMax = 1.0 / *MMin; // 极图中斜率范围（y=1/x是减函数）
                }
                else if (InitialSlope > 0) {
                    M = InitialSlope;
                    N = -(double)i * InitialSlope;
                    kRange(M, MMin, MMax, N); //离散图中：斜率范围
                    SlopeMin = *MMin;
                    SlopeMax = *MMax; // 极图中斜率范围（y=x是增函数）

                }
                else if (InitialSlope >= -1) {
                    M = -InitialSlope;
                    N = (double)i * InitialSlope;
                    kRange(M, MMin, MMax, N); //离散图中：斜率范围
                    SlopeMin = -*MMax;
                    SlopeMax = -*MMin; // 极图中斜率范围（y=-x是减函数）

                }
                else {
                    M = -1.0 / InitialSlope;
                    N = (double)i;
                    kRange(M, MMin, MMax, N); //离散图中：斜率范围
                    SlopeMin = -1.0 / *MMin;
                    SlopeMax = -1.0 / *MMax; // 极图中斜率范围（y=-1/x是增函数）

                }
                SlopeMin = checkSign(InitialSlope, SlopeMin - e);
                SlopeMax = checkSign(InitialSlope, SlopeMax + e);





            if (isinf(SlopeMin) || isinf(SlopeMax)) {
                re_slope[i * Ny + j] = slope[i * Ny + j];
                continue;
            }



            //计算当前邻域内的deltadis

            DeltaDis = ComputeDelta(slope, i, j, xIdx, ChosenPtsNum);

            //计算所有xIdx直线在[SlopeMin,SlopeMax]之间的Cost
            computeCost(im_in_remap_pt, SlopeMin, SlopeMax, j, i, xIdx, ChosenPtsNum, DeltaDis, maxdis);

            //计算最优cost
            minIdx = getMinIdx(j, i);

            re_slope[i * Ny + j] = SlopeMin + minIdx * (SlopeMax - SlopeMin) / (res - 1.0);

            delete[] xIdx;
        }
    }
    delete[] MMax;
    delete[] MMin;
    delete[] slope_cost;

}
