#include "mex.h"
#include "math.h"

using namespace std;

// 输入
#define	Im_in_remap         prhs[0]  // the remap image
#define	x_size   			prhs[1]  // spatial image width
#define	y_size              prhs[2]  // spatial image height
#define	AngleResolution     prhs[3]  // angular resolution
#define Angle_min           prhs[4]  // 最小角度
#define Angle_max           prhs[5]  // 最大角度
#define Angle_res           prhs[6]  // angle resolution

//输出
#define Slope_cost_sx       plhs[0]  // slope cost with sx-epis
#define Gx_cost_sx          plhs[1]
#define Gy_cost_sx          plhs[2]

//全局参数
int st_diameter;
int st_radius;
int st_size;
int height;
int width;
int pixelNum;
int remap_height;
int remap_width;
int remap_pixelNum;
double angle_min;
double angle_max;
int angle_res;
double angle_step;
int bsz = 5;               // spatial window size of bilateral filter for depth cost
double sigma_range = 0.01;
constexpr auto pi = 3.1415926;





int checkX(int x) {
    if (x < 0) return 0;
    else if (x >= width) return width - 1;
    else return x;
}

int checkY(int y) {
    if (y < 0) return 0;
    else if (y >= height) return height - 1;
    else return y;
}
int checkST(int x) {
    if (x < 0) return 0;
    else if (x >= st_diameter * st_diameter) return st_diameter * st_diameter - 1;
    else return x;
}


void stitching(double im_in_remap[], int x, int y, int angle_num, double remap[], double* Ix, double* Iy)
{
    double slope = tan((angle_min + angle_step * angle_num) * pi / 180);

    for (int j = -st_radius; j <= st_radius; j++) { //j=t
        float y_ind = (float)j / slope + y; //  每张单极图对应的y
        int   y_floor = floor(y_ind); // 以下是为了插出y_ind的所需数据
        int   y_1 = checkY(y_floor);
        int   y_2 = checkY(y_floor + 1);
        float y_1_w = 1 - (y_ind - y_floor); // 权重
        float y_2_w = 1 - y_1_w;

        for (int i = -st_radius; i <= st_radius; i++) { // i=s
            float x_ind = round(x + (j * st_diameter + i) / slope) - st_diameter * j / slope;
            int   x_floor = floor(x_ind);
            int   x_1 = checkX(x_floor);
            int   x_2 = checkX(x_floor + 1);
            float x_1_w = 1 - (x_ind - x_floor);
            float x_2_w = 1 - x_1_w;

            int x_1_index = i + st_radius + (x_1)*st_diameter;
            int y_1_index = j + st_radius + (y_1)*st_diameter;
            int x_2_index = i + st_radius + (x_2)*st_diameter;
            int y_2_index = j + st_radius + (y_2)*st_diameter;

            for (int c = 0; c < 3; c++) {
                remap[(i + st_radius) * st_diameter + (j + st_radius) + c * st_diameter * st_diameter] =
                    y_1_w * x_1_w * im_in_remap[y_1_index + x_1_index * remap_height + c * remap_pixelNum] +
                    y_2_w * x_1_w * im_in_remap[y_2_index + x_1_index * remap_height + c * remap_pixelNum] +
                    y_1_w * x_2_w * im_in_remap[y_1_index + x_2_index * remap_height + c * remap_pixelNum] +
                    y_2_w * x_2_w * im_in_remap[y_2_index + x_2_index * remap_height + c * remap_pixelNum];
            }


            // for Ix  【DONE】
            int x_1_left = checkX(x_1 - 1);
            int x_2_left = x_1;
            int x_1_index_left = i + st_radius + (x_1_left)*st_diameter;
            int x_2_index_left = i + st_radius + (x_2_left)*st_diameter;

            int x_1_right = x_1;
            int x_2_right = checkX(x_1 + 1);
            int x_1_index_right = i + st_radius + (x_1_right)*st_diameter;
            int x_2_index_right = i + st_radius + (x_2_right)*st_diameter;

            for (int c = 0; c < 3; c++) {
                Ix[(i + st_radius) * st_diameter + (j + st_radius) + c * st_diameter * st_diameter] =
                    ((y_1_w * x_1_w * im_in_remap[y_1_index + x_1_index_right * remap_height + c * remap_pixelNum] +
                        y_2_w * x_1_w * im_in_remap[y_2_index + x_1_index_right * remap_height + c * remap_pixelNum] +
                        y_1_w * x_2_w * im_in_remap[y_1_index + x_2_index_right * remap_height + c * remap_pixelNum] +
                        y_2_w * x_2_w * im_in_remap[y_2_index + x_2_index_right * remap_height + c * remap_pixelNum])
                        -
                        (y_1_w * x_1_w * im_in_remap[y_1_index + x_1_index_left * remap_height + c * remap_pixelNum] +
                            y_2_w * x_1_w * im_in_remap[y_2_index + x_1_index_left * remap_height + c * remap_pixelNum] +
                            y_1_w * x_2_w * im_in_remap[y_1_index + x_2_index_left * remap_height + c * remap_pixelNum] +
                            y_2_w * x_2_w * im_in_remap[y_2_index + x_2_index_left * remap_height + c * remap_pixelNum])) / ((x_2_left - x_1_left) + (x_2_right - x_1_right));
            }

            // for Iy

            // 先确定s和t

            int st_ind = (j + st_radius) * st_diameter + (i + st_radius);
            int st_up_ind = checkST(st_ind - 1);
            int st_down_ind = checkST(st_ind + 1);

            int t_up_ind = floor(st_up_ind / st_diameter) - st_radius;
            int s_up_ind = st_up_ind - (t_up_ind + st_radius) * st_diameter - st_radius;

            int t_down_ind = floor(st_down_ind / st_diameter) - st_radius;
            int s_down_ind = st_down_ind - (t_down_ind + st_radius) * st_diameter - st_radius;

            // up 点
            int x_up_1_index, x_up_2_index;
            if (t_up_ind != j)
            {
                float x_up_ind = x_ind - st_diameter / slope;
                int   x_up_floor = floor(x_up_ind);
                int   x_up_1 = checkX(x_up_floor);
                int   x_up_2 = checkX(x_up_floor + 1);
                x_up_1_index = s_up_ind + st_radius + (x_up_1)*st_diameter;
                x_up_2_index = s_up_ind + st_radius + (x_up_2)*st_diameter;
            }
            else
            {
                x_up_1_index = s_up_ind + st_radius + (x_1)*st_diameter;
                x_up_2_index = s_up_ind + st_radius + (x_2)*st_diameter;
            }

            float y_up_ind = (float)t_up_ind / slope + y;
            int   y_up_floor = floor(y_up_ind);
            int   y_up_1 = checkY(y_up_floor);
            int   y_up_2 = checkY(y_up_floor + 1);
            float y_up_1_w = 1 - (y_up_ind - y_up_floor);
            float y_up_2_w = 1 - y_up_1_w;
            int   y_up_1_index = t_up_ind + st_radius + (y_up_1)*st_diameter;
            int   y_up_2_index = t_up_ind + st_radius + (y_up_2)*st_diameter;


            // down 点
            int x_down_1_index, x_down_2_index;
            if (t_down_ind != j)
            {
                float x_down_ind = x_ind + st_diameter / slope;
                int   x_down_floor = floor(x_down_ind);
                int   x_down_1 = checkX(x_down_floor);
                int   x_down_2 = checkX(x_down_floor + 1);
                x_down_1_index = s_down_ind + st_radius + (x_down_1)*st_diameter;
                x_down_2_index = s_down_ind + st_radius + (x_down_2)*st_diameter;
            }
            else
            {
                x_down_1_index = s_down_ind + st_radius + (x_1)*st_diameter;
                x_down_2_index = s_down_ind + st_radius + (x_2)*st_diameter;
            }

            float y_down_ind = (float)t_down_ind / slope + y;
            int   y_down_floor = floor(y_down_ind);
            int   y_down_1 = checkY(y_down_floor);
            int   y_down_2 = checkY(y_down_floor + 1);
            float y_down_1_w = 1 - (y_down_ind - y_down_floor);
            float y_down_2_w = 1 - y_down_1_w;
            int y_down_1_index = t_down_ind + st_radius + (y_down_1)*st_diameter;
            int y_down_2_index = t_down_ind + st_radius + (y_down_2)*st_diameter;


            for (int c = 0; c < 3; c++) {
                Iy[(i + st_radius) * st_diameter + (j + st_radius) + c * st_diameter * st_diameter] =
                    ((y_up_1_w * x_1_w * im_in_remap[y_up_1_index + x_up_1_index * remap_height + c * remap_pixelNum] +
                        y_up_2_w * x_1_w * im_in_remap[y_up_2_index + x_up_1_index * remap_height + c * remap_pixelNum] +
                        y_up_1_w * x_2_w * im_in_remap[y_up_1_index + x_up_2_index * remap_height + c * remap_pixelNum] +
                        y_up_2_w * x_2_w * im_in_remap[y_up_2_index + x_up_2_index * remap_height + c * remap_pixelNum]) -

                        (y_down_1_w * x_1_w * im_in_remap[y_down_1_index + x_down_1_index * remap_height + c * remap_pixelNum] +
                            y_down_2_w * x_1_w * im_in_remap[y_down_2_index + x_down_1_index * remap_height + c * remap_pixelNum] +
                            y_down_1_w * x_2_w * im_in_remap[y_down_1_index + x_down_2_index * remap_height + c * remap_pixelNum] +
                            y_down_2_w * x_2_w * im_in_remap[y_down_2_index + x_down_2_index * remap_height + c * remap_pixelNum])) / (st_down_ind - st_up_ind);
            }

        }
    }
}

void stitching1(double im_in_remap[], int x, int y, int angle_num, double remap[], double* Ix, double* Iy)
{
    double slope = tan((angle_min + angle_step * angle_num) * pi / 180);

    for (int j = -st_radius; j <= st_radius; j++) { //j=t
        float y_ind = (float)j / slope + y; //  每张单极图对应的y
        int   y_floor = floor(y_ind); // 以下是为了插出y_ind的所需数据
        int   y_1 = checkY(y_floor);
        int   y_2 = checkY(y_floor + 1);
        float y_1_w = 1 - (y_ind - y_floor); // 权重
        float y_2_w = 1 - y_1_w;

        for (int i = -st_radius; i <= st_radius; i++) { // i=s
            float x_ind = round(x + (j * st_diameter + i) / slope) - st_diameter * j / slope;
            int   x_floor = floor(x_ind);
            int   x_1 = checkX(x_floor);
            int   x_2 = checkX(x_floor + 1);
            float x_1_w = 1 - (x_ind - x_floor);
            float x_2_w = 1 - x_1_w;

            int x_1_index = i + st_radius + (x_1)*st_diameter;
            int y_1_index = j + st_radius + (y_1)*st_diameter;
            int x_2_index = i + st_radius + (x_2)*st_diameter;
            int y_2_index = j + st_radius + (y_2)*st_diameter;

            for (int c = 0; c < 3; c++) {
                remap[(i + st_radius) * st_diameter + (j + st_radius) + c * st_diameter * st_diameter] =
                    y_1_w * x_1_w * im_in_remap[y_1_index + x_1_index * remap_height + c * remap_pixelNum] +
                    y_2_w * x_1_w * im_in_remap[y_2_index + x_1_index * remap_height + c * remap_pixelNum] +
                    y_1_w * x_2_w * im_in_remap[y_1_index + x_2_index * remap_height + c * remap_pixelNum] +
                    y_2_w * x_2_w * im_in_remap[y_2_index + x_2_index * remap_height + c * remap_pixelNum];
            }


            // for Ix  【DONE】
            int x_1_left = checkX(x_1 - 1);
            int x_2_left = x_1;
            int x_1_index_left = i + st_radius + (x_1_left)*st_diameter;
            int x_2_index_left = i + st_radius + (x_2_left)*st_diameter;

            int x_1_right = x_1;
            int x_2_right = checkX(x_1 + 1);
            int x_1_index_right = i + st_radius + (x_1_right)*st_diameter;
            int x_2_index_right = i + st_radius + (x_2_right)*st_diameter;

            for (int c = 0; c < 3; c++) {
                Ix[(i + st_radius) * st_diameter + (j + st_radius) + c * st_diameter * st_diameter] =
                    ((y_1_w * x_1_w * im_in_remap[y_1_index + x_1_index_right * remap_height + c * remap_pixelNum] +
                        y_2_w * x_1_w * im_in_remap[y_2_index + x_1_index_right * remap_height + c * remap_pixelNum] +
                        y_1_w * x_2_w * im_in_remap[y_1_index + x_2_index_right * remap_height + c * remap_pixelNum] +
                        y_2_w * x_2_w * im_in_remap[y_2_index + x_2_index_right * remap_height + c * remap_pixelNum])
                        -
                        (y_1_w * x_1_w * im_in_remap[y_1_index + x_1_index_left * remap_height + c * remap_pixelNum] +
                            y_2_w * x_1_w * im_in_remap[y_2_index + x_1_index_left * remap_height + c * remap_pixelNum] +
                            y_1_w * x_2_w * im_in_remap[y_1_index + x_2_index_left * remap_height + c * remap_pixelNum] +
                            y_2_w * x_2_w * im_in_remap[y_2_index + x_2_index_left * remap_height + c * remap_pixelNum])) / ((x_2_left - x_1_left) + (x_2_right - x_1_right));
            }

            // for Iy

            // 先确定s和t
            int t_up_ind = j;
            int t_down_ind = j;

            int s_up_ind = i == -st_radius ? i : i - 1;
            int s_down_ind = i == st_radius ? i : i + 1;

            // up 点
            int x_up_1_index = s_up_ind + st_radius + (x_1)*st_diameter;
            int x_up_2_index = s_up_ind + st_radius + (x_2)*st_diameter;

            float y_up_ind = (float)t_up_ind / slope + y;
            int   y_up_floor = floor(y_up_ind);
            int   y_up_1 = checkY(y_up_floor);
            int   y_up_2 = checkY(y_up_floor + 1);
            float y_up_1_w = 1 - (y_up_ind - y_up_floor);
            float y_up_2_w = 1 - y_up_1_w;
            int   y_up_1_index = t_up_ind + st_radius + (y_up_1)*st_diameter;
            int   y_up_2_index = t_up_ind + st_radius + (y_up_2)*st_diameter;

            // down 点
            int x_down_1_index = s_down_ind + st_radius + (x_1)*st_diameter;
            int x_down_2_index = s_down_ind + st_radius + (x_2)*st_diameter;

            float y_down_ind = (float)t_down_ind / slope + y;
            int   y_down_floor = floor(y_down_ind);
            int   y_down_1 = checkY(y_down_floor);
            int   y_down_2 = checkY(y_down_floor + 1);
            float y_down_1_w = 1 - (y_down_ind - y_down_floor);
            float y_down_2_w = 1 - y_down_1_w;
            int y_down_1_index = t_down_ind + st_radius + (y_down_1)*st_diameter;
            int y_down_2_index = t_down_ind + st_radius + (y_down_2)*st_diameter;


            for (int c = 0; c < 3; c++) {
                Iy[(i + st_radius) * st_diameter + (j + st_radius) + c * st_diameter * st_diameter] =
                    ((y_up_1_w * x_1_w * im_in_remap[y_up_1_index + x_up_1_index * remap_height + c * remap_pixelNum] +
                        y_up_2_w * x_1_w * im_in_remap[y_up_2_index + x_up_1_index * remap_height + c * remap_pixelNum] +
                        y_up_1_w * x_2_w * im_in_remap[y_up_1_index + x_up_2_index * remap_height + c * remap_pixelNum] +
                        y_up_2_w * x_2_w * im_in_remap[y_up_2_index + x_up_2_index * remap_height + c * remap_pixelNum]) -

                        (y_down_1_w * x_1_w * im_in_remap[y_down_1_index + x_down_1_index * remap_height + c * remap_pixelNum] +
                            y_down_2_w * x_1_w * im_in_remap[y_down_2_index + x_down_1_index * remap_height + c * remap_pixelNum] +
                            y_down_1_w * x_2_w * im_in_remap[y_down_1_index + x_down_2_index * remap_height + c * remap_pixelNum] +
                            y_down_2_w * x_2_w * im_in_remap[y_down_2_index + x_down_2_index * remap_height + c * remap_pixelNum])) / (s_down_ind - s_up_ind);
            }

        }
    }
}
void stitching2(double im_in_remap[], int x, int y, int angle_num, double remap[], double* Ix, double* Iy)
{
    double slope = tan((angle_min + angle_step * angle_num) * pi / 180);

    for (int j = -st_radius; j <= st_radius; j++) { //j=t
        int t = -j;
        float y_ind = (float)t / slope + y; //  每张单极图对应的y
        int   y_floor = floor(y_ind); // 以下是为了插出y_ind的所需数据
        int   y_1 = checkY(y_floor);
        int   y_2 = checkY(y_floor + 1);
        float y_1_w = 1 - (y_ind - y_floor); // 权重
        float y_2_w = 1 - y_1_w;

        for (int i = -st_radius; i <= st_radius; i++) { // i=s
            int s = -i;
            float x_ind = round(x + (t * st_diameter + s) / slope) - st_diameter * t / slope;
            int   x_floor = floor(x_ind);
            int   x_1 = checkX(x_floor);
            int   x_2 = checkX(x_floor + 1);
            float x_1_w = 1 - (x_ind - x_floor);
            float x_2_w = 1 - x_1_w;

            int x_1_index = i + st_radius + (x_1)*st_diameter;
            int y_1_index = j + st_radius + (y_1)*st_diameter;
            int x_2_index = i + st_radius + (x_2)*st_diameter;
            int y_2_index = j + st_radius + (y_2)*st_diameter;

            for (int c = 0; c < 3; c++) {
                remap[(i + st_radius) * st_diameter + (j + st_radius) + c * st_diameter * st_diameter] =
                    y_1_w * x_1_w * im_in_remap[y_1_index + x_1_index * remap_height + c * remap_pixelNum] +
                    y_2_w * x_1_w * im_in_remap[y_2_index + x_1_index * remap_height + c * remap_pixelNum] +
                    y_1_w * x_2_w * im_in_remap[y_1_index + x_2_index * remap_height + c * remap_pixelNum] +
                    y_2_w * x_2_w * im_in_remap[y_2_index + x_2_index * remap_height + c * remap_pixelNum];
            }


            // for Ix  【DONE】
            int x_1_left = checkX(x_1 - 1);
            int x_2_left = x_1;
            int x_1_index_left = i + st_radius + (x_1_left)*st_diameter;
            int x_2_index_left = i + st_radius + (x_2_left)*st_diameter;

            int x_1_right = x_1;
            int x_2_right = checkX(x_1 + 1);
            int x_1_index_right = i + st_radius + (x_1_right)*st_diameter;
            int x_2_index_right = i + st_radius + (x_2_right)*st_diameter;

            for (int c = 0; c < 3; c++) {
                Ix[(i + st_radius) * st_diameter + (j + st_radius) + c * st_diameter * st_diameter] =
                    ((y_1_w * x_1_w * im_in_remap[y_1_index + x_1_index_right * remap_height + c * remap_pixelNum] +
                        y_2_w * x_1_w * im_in_remap[y_2_index + x_1_index_right * remap_height + c * remap_pixelNum] +
                        y_1_w * x_2_w * im_in_remap[y_1_index + x_2_index_right * remap_height + c * remap_pixelNum] +
                        y_2_w * x_2_w * im_in_remap[y_2_index + x_2_index_right * remap_height + c * remap_pixelNum])
                        -
                        (y_1_w * x_1_w * im_in_remap[y_1_index + x_1_index_left * remap_height + c * remap_pixelNum] +
                            y_2_w * x_1_w * im_in_remap[y_2_index + x_1_index_left * remap_height + c * remap_pixelNum] +
                            y_1_w * x_2_w * im_in_remap[y_1_index + x_2_index_left * remap_height + c * remap_pixelNum] +
                            y_2_w * x_2_w * im_in_remap[y_2_index + x_2_index_left * remap_height + c * remap_pixelNum])) / ((x_2_left - x_1_left) + (x_2_right - x_1_right));
            }

            // for Iy

            // 先确定s和t
            int j_up_ind = j;
            int j_down_ind = j;

            int t_up_ind = t;
            int t_down_ind = t;

            int i_up_ind = i == -st_radius ? i : i - 1;
            int i_down_ind = i == st_radius ? i : i + 1;

            int s_up_ind = -i_up_ind;
            int s_down_ind = -i_down_ind;

            // up 点
            int x_up_1_index = i_up_ind + st_radius + (x_1)*st_diameter;
            int x_up_2_index = i_up_ind + st_radius + (x_2)*st_diameter;

            float y_up_ind = (float)t_up_ind / slope + y;
            int   y_up_floor = floor(y_up_ind);
            int   y_up_1 = checkY(y_up_floor);
            int   y_up_2 = checkY(y_up_floor + 1);
            float y_up_1_w = 1 - (y_up_ind - y_up_floor);
            float y_up_2_w = 1 - y_up_1_w;
            int   y_up_1_index = j_up_ind + st_radius + (y_up_1)*st_diameter;
            int   y_up_2_index = j_up_ind + st_radius + (y_up_2)*st_diameter;

            // down 点
            int x_down_1_index = i_down_ind + st_radius + (x_1)*st_diameter;
            int x_down_2_index = i_down_ind + st_radius + (x_2)*st_diameter;

            float y_down_ind = (float)t_down_ind / slope + y;
            int   y_down_floor = floor(y_down_ind);
            int   y_down_1 = checkY(y_down_floor);
            int   y_down_2 = checkY(y_down_floor + 1);
            float y_down_1_w = 1 - (y_down_ind - y_down_floor);
            float y_down_2_w = 1 - y_down_1_w;
            int y_down_1_index = j_down_ind + st_radius + (y_down_1)*st_diameter;
            int y_down_2_index = j_down_ind + st_radius + (y_down_2)*st_diameter;


            for (int c = 0; c < 3; c++) {
                Iy[(i + st_radius) * st_diameter + (j + st_radius) + c * st_diameter * st_diameter] =
                    ((y_up_1_w * x_1_w * im_in_remap[y_up_1_index + x_up_1_index * remap_height + c * remap_pixelNum] +
                        y_up_2_w * x_1_w * im_in_remap[y_up_2_index + x_up_1_index * remap_height + c * remap_pixelNum] +
                        y_up_1_w * x_2_w * im_in_remap[y_up_1_index + x_up_2_index * remap_height + c * remap_pixelNum] +
                        y_up_2_w * x_2_w * im_in_remap[y_up_2_index + x_up_2_index * remap_height + c * remap_pixelNum]) -

                        (y_down_1_w * x_1_w * im_in_remap[y_down_1_index + x_down_1_index * remap_height + c * remap_pixelNum] +
                            y_down_2_w * x_1_w * im_in_remap[y_down_2_index + x_down_1_index * remap_height + c * remap_pixelNum] +
                            y_down_1_w * x_2_w * im_in_remap[y_down_1_index + x_down_2_index * remap_height + c * remap_pixelNum] +
                            y_down_2_w * x_2_w * im_in_remap[y_down_2_index + x_down_2_index * remap_height + c * remap_pixelNum])) / (i_down_ind - i_up_ind);
            }

        }
    }
}
// compute mean and variance
double getVar(double array[], int size, double* mean)
{
    *mean = 0;
    int num = 0;
    for (int i = 0; i < size; i++) {
        if (array[i] >= 0) {
            *mean += array[i];
            num++;
        }
    }
    *mean /= num;
    double var = 0;
    for (int i = 0; i < size; i++) {
        if (array[i] >= 0) {
            var += (array[i] - *mean) * (array[i] - *mean);
        }
    }
    if (num == 1) var = INT_MAX;
    else var /= (num - 1.0);
    return var;
}

// 计算1-范数
double getNorm1(double array[], int size)
{
    double center = array[(size + 1) / 2];

    double norm1 = 0;
    for (int i = 0; i < size; i++) {
        norm1 += fabs(array[i] - center);
    }
    return norm1;
}

void computeCost(float slope_cost[], float gx_cost[], float gy_cost[], double im_in_remap[], double CSAI[], int x, int y)
{
    double* remap = new double[st_diameter * st_diameter * 3.0];
    double* Ix = new double[st_diameter * st_diameter * 3.0];
    double* Iy = new double[st_diameter * st_diameter * 3.0];

    double mean_p1[3];

    for (int angle_num = 0; angle_num < angle_res; angle_num++) {

        stitching2(im_in_remap, x, y, angle_num, remap, Ix, Iy);

        // Jake // if (orient[y+x*height] < -100*M_PI/180+0.1) { // no edge, apply traditional photo-consistency // Jake: variance of three color channels
        //double var1 = getVar(remap, st_size, mean_p1) + getVar(&remap[st_size], st_size, &mean_p1[1]) + getVar(&remap[st_size * 2], st_size, &mean_p1[2]);
        //slope_cost[y + x * height + angle_num * pixelNum] = sqrt(var1);

        double norm1 = getNorm1(remap, st_size) + getNorm1(&remap[st_size], st_size) + getNorm1(&remap[st_size * 2], st_size);
        slope_cost[y + x * height + angle_num * pixelNum] = norm1 / 3;

        double norm2 = getNorm1(Ix, st_size) + getNorm1(&Ix[st_size], st_size) + getNorm1(&Ix[st_size * 2], st_size);
        gx_cost[y + x * height + angle_num * pixelNum] = norm2 / 3;

        double norm3 = getNorm1(Iy, st_size) + getNorm1(&Iy[st_size], st_size) + getNorm1(&Iy[st_size * 2], st_size);
        gy_cost[y + x * height + angle_num * pixelNum] = norm3 / 3;




        //     for (int c = 0; c < 3; c++) {
  //          slope_cost[y + x * height + c * pixelNum + angle_num * 3 * pixelNum] = mean_p1[c];
  //          refocus[y + x * height + c * pixelNum + alpha_num * 3 * pixelNum + depth_res * pixelNum * 3] = mean_p1[c];
  //      }

 //       depth_cost[y + x * height + alpha_num * pixelNum + depth_res * pixelNum] = sqrt(var1);
        // Jake // continue;
    // Jake // }
    }
    delete[] remap;
    delete[] Ix;
    delete[] Iy;

}


void costFilter(float slope_cost[], float slope_cost_filtered[], float gx_cost[], float gx_cost_filtered[], float gy_cost[], float gy_cost_filtered[], double im_pinhole[])
{
    for (int i = 0; i < angle_res; i++) {
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                double sumVal_slope = 0;
                double weight_sum = 0;
                double sumVal_gx = 0;
                double sumVal_gy = 0;
                for (int u = -bsz / 2; u <= bsz / 2; u++) {
                    for (int v = -bsz / 2; v <= bsz / 2; v++) {
                        int xu = checkX(x + u);
                        int yv = checkY(y + v);

                        // compute bilateral weight
                        double color_dif = 0;
                        for (int c = 0; c < 3; c++)
                            color_dif += pow(im_pinhole[y + x * height + c * pixelNum] - im_pinhole[yv + xu * height + c * pixelNum], 2);
                        double weight = exp(-color_dif / (2 * pow(sigma_range, 2)));

                        weight_sum += weight;
                        // slope cost
                        sumVal_slope += weight * slope_cost[xu * height + yv + i * pixelNum];

                        // gx cost
                        sumVal_gx += weight * gx_cost[xu * height + yv + i * pixelNum];

                        // gy cost
                        sumVal_gy += weight * gy_cost[xu * height + yv + i * pixelNum];

                        // refocus cost
                        //double gra = 0;
                        //for (int c = 0; c < 3; c++) {
                        //    gra += fabs(refocus[yv + xu * height + c * pixelNum + i * 3 * pixelNum] - im_pinhole[yv + xu * height + c * pixelNum]);
                        //}
                        //sumVal_f += weight * gra;
                        //weight_sum_f += weight;
                    }
                }
                sumVal_slope = weight_sum != 0 ? sumVal_slope / weight_sum : 1;
                slope_cost_filtered[y + x * height + i * pixelNum] = sumVal_slope;

                sumVal_gx = weight_sum != 0 ? sumVal_gx / weight_sum : 1;
                gx_cost_filtered[y + x * height + i * pixelNum] = sumVal_gx;

                sumVal_gy = weight_sum != 0 ? sumVal_gy / weight_sum : 1;
                gy_cost_filtered[y + x * height + i * pixelNum] = sumVal_gy;
                //sumVal_c = weight_sum_c != 0 ? sumVal_c / weight_sum_c : 1;  // Jake: 1 for sharp forcus 
                //sumVal_f = weight_sum_f != 0 ? sumVal_f / weight_sum_f : 1;  // Jake: 1 for sharp forcus
                //depth_cost_sum[y+x*height+i*pixelNum] = sumVal_c + ratio*sumVal_f;

                //slope_cost_filtered[y + x * height + i * pixelNum] = sumVal_c;
               // slope_cost_filtered[y + x * height + i * pixelNum] = slope_cost[y + x * height + i * pixelNum];

                // Unlike the orignial defocus clue, which calculates the regional variance, w.r.t. to summed laplacian
                // here the defocus clue calculates the regional difference with the sharply focused pinhole image, 
                // therefore also the smaller, the more focused.
                //depth_var_defocus[y + x * height + i * pixelNum] = sumVal_f;
            }
        }
    }
}




void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    double* im_in_remap_pt = mxGetPr(Im_in_remap);
    double* width_pt = mxGetPr(x_size);
    double* height_pt = mxGetPr(y_size);
    double* st_diameter_pt = mxGetPr(AngleResolution);
    double* angle_min_pt = mxGetPr(Angle_min);
    double* angle_max_pt = mxGetPr(Angle_max);
    double* angle_res_pt = mxGetPr(Angle_res);


    /* 设置参数*/
    angle_min = *angle_min_pt;
    angle_max = *angle_max_pt;
    angle_res = *angle_res_pt;
    angle_step = (angle_max - angle_min) / (angle_res - 1.0);

    st_diameter = *st_diameter_pt;
    st_radius = (st_diameter - 1) / 2;
    st_size = st_diameter * st_diameter;

    height = *height_pt;
    width = *width_pt;
    pixelNum = height * width;

    remap_height = height * st_diameter;
    remap_width = width * st_diameter;
    remap_pixelNum = remap_height * remap_width;



    // Mine
    // 构建输出数组
    const mwSize dims2[] = { height, width, angle_res };
    Slope_cost_sx = mxCreateNumericArray(3, dims2, mxSINGLE_CLASS, mxREAL);
    Gx_cost_sx = mxCreateNumericArray(3, dims2, mxSINGLE_CLASS, mxREAL);
    Gy_cost_sx = mxCreateNumericArray(3, dims2, mxSINGLE_CLASS, mxREAL);

    //float* slope_cost_sx_pt = (float*)mxGetPr(Slope_cost_sx);
    float* slope_cost_sx_pt = new float[height * width * angle_res];
    float* slope_cost_sx_filtered_pt = (float*)mxGetPr(Slope_cost_sx);

    float* gx_cost_sx_pt = new float[height * width * angle_res];
    float* gx_cost_sx_filtered_pt = (float*)mxGetPr(Gx_cost_sx);

    float* gy_cost_sx_pt = new float[height * width * angle_res];
    float* gy_cost_sx_filtered_pt = (float*)mxGetPr(Gy_cost_sx);

    double* CSAI = new double[pixelNum * 3.0];
    // Mine



    // 中心子孔径图像
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            for (int c = 0; c < 3; c++) {
                CSAI[i * height + j + pixelNum * c] =
                    im_in_remap_pt[j * st_diameter + st_radius + (i * st_diameter + st_radius) * remap_height + remap_pixelNum * c];
            }
        }
    }


    for (int i = 0; i < width; i++) { //x
        for (int j = 0; j < height; j++) { //y    
            computeCost(slope_cost_sx_pt, gx_cost_sx_pt, gy_cost_sx_pt, im_in_remap_pt, CSAI, i, j);
        }
    }

    costFilter(slope_cost_sx_pt, slope_cost_sx_filtered_pt, gx_cost_sx_pt, gx_cost_sx_filtered_pt, gy_cost_sx_pt, gy_cost_sx_filtered_pt, CSAI);

    delete[] CSAI;
    delete[] slope_cost_sx_pt;
    delete[] gx_cost_sx_pt;
    delete[] gy_cost_sx_pt;

}