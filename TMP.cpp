#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include <windows.h>
#include <gdiplus.h>
#pragma comment(lib, "gdiplus.lib")

using namespace std;
using namespace Gdiplus;

#define M_PI 3.14159265358979323846
const int max_size = 1024;
complex<double> gra[max_size][max_size] = { 0 };
complex<double> foru_g[max_size][max_size] = { 0 };


void RaderReverse(complex<double>* arr, int N) {
    int j, k;
    for (int i = 1, j = N / 2; i < N - 1; ++i) {
        if (i < j) {
            complex<double> temp = arr[j];
            arr[j] = arr[i];
            arr[i] = temp;
        }
        k = N / 2;
        while (k <= j) {
            j = j - k;
            k = k / 2;
        }
        j = j + k;
    }
}

void resize(Bitmap* bmp, int height, int width)
{

    Color color;
    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            bmp->GetPixel(i, j, &color);
            gra[i][j] = complex<double>(color.GetRed(), 0.0);
        }
    }
    return;
}

void dctr(int row)
{
    int i = 0, j = 0, k = 0, l = 0;
    complex<double> up, down, product, x[max_size];
    complex<double>* W = new complex<double>[max_size];

    for (int i = 0; i < max_size; i++)
    {
        W[i] = exp(complex<double>(0, -2.0 * M_PI * i / max_size));
    }
    for (int i = 0; i < max_size; i++) {
        x[i] = foru_g[i][row];
    }
    RaderReverse(x, max_size);

    for (i = 0; i < log(max_size) / log(2); ++i)
    {
        l = 1 << i;
        int count = 0;
        complex<double>* A = new complex<double>[max_size];
        for (j = 0; j < max_size; j += 2 * l)
        {
            for (k = 0; k < l * 2; ++k)
            {

                if (k < l) {
                    A[count++] = x[j + k] + x[j + k + l] * W[max_size * k / 2 / l];
                }
                else {
                    A[count++] = x[j + k - l] - x[j + k] * W[max_size * (k - l) / 2 / l];
                }

            }
        }
        for (int t = 0; t < max_size; t++)x[t] = A[t];
    }
    for (i = 0; i < max_size; i++)foru_g[i][row] = x[i];
    for (i = 0; i < max_size; i++) {
        complex<double> tmp = complex<double>(x[i].real() * 2 / max_size, 0);
        foru_g[i][row] = tmp;

    }
    int ans = 0;
    for (i = 0; i < max_size; i++)
    {
        ans += foru_g[i][row].real();

    }
    foru_g[0][row] = complex<double>(ans / max_size, 0);
    delete[] W;
}

void dctc(int column)
{
    int i = 0, j = 0, k = 0, l = 0;
    complex<double> up, down, product, x[max_size];
    complex<double>* W = new complex<double>[max_size];

    for (int i = 0; i < max_size; i++)
    {
        W[i] = exp(complex<double>(0, -2.0 * M_PI * i / max_size));
        //if (column == 5)printf("%lf %lf\n", W[i].real(), W[i].imag());
    }

    for (int i = 0; i < max_size; i++) {
        x[i] = gra[column][i];
    }
    RaderReverse(x, max_size);

    for (i = 0; i < log(max_size) / log(2); ++i)
    {
        l = 1 << i;
        int count = 0;
        complex<double>* A = new complex<double>[max_size];
        for (j = 0; j < max_size; j += 2 * l)
        {
            for (k = 0; k < l * 2; ++k)
            {

                if (k < l) {
                    A[count++] = x[j + k] + x[j + k + l] * W[max_size * k / 2 / l];
                }
                else {
                    A[count++] = x[j + k - l] - x[j + k] * W[max_size * (k - l) / 2 / l];
                }
            }
        }
        for (int t = 0; t < max_size; t++)x[t] = A[t];
    }
    for (i = 0; i < max_size; i++) {
        complex<double> tmp = complex<double>(x[i].real(), 0);
        foru_g[column][i] = tmp;

    }
    int ans = 0;
    for (i = 0; i < max_size; i++)
    {
        ans += foru_g[column][i].real();

    }
    foru_g[column][0] = complex<double>(ans, 0);
    delete[] W;
}

void ndctr(int row)
{
    int i = 0, j = 0, k = 0, l = 0;
    complex<double> up, down, product, x[max_size];
    complex<double>* W = new complex<double>[max_size];

    for (int i = 0; i < max_size; i++)
    {
        W[i] = exp(complex<double>(0, 2.0 * M_PI * i / max_size));
    }
    for (int i = 0; i < max_size; i++) {
        x[i] = foru_g[i][row];
        if (i > max_size / 2)x[i] = 0;
    }
    RaderReverse(x, max_size);

    for (i = 0; i < log(max_size) / log(2); ++i)
    {
        l = 1 << i;
        int count = 0;
        complex<double>* A = new complex<double>[max_size];
        for (j = 0; j < max_size; j += 2 * l)
        {
            for (k = 0; k < l * 2; ++k)
            {

                if (k < l) {
                    A[count++] = x[j + k] + x[j + k + l] * W[max_size * k / 2 / l];
                }
                else {
                    A[count++] = x[j + k - l] - x[j + k] * W[max_size * (k - l) / 2 / l];
                }

            }
        }
        for (int t = 0; t < max_size; t++)x[t] = A[t];
    }
    complex<double> a(max_size, 0);
    for (i = 0; i < max_size; i++)foru_g[i][row] = x[i] / a / a;
    delete[] W;
}

void dct()
{
    for (int i = 0; i < max_size; i++)
    {
        dctc(i);
    }

    for (int i = 0; i < max_size; i++)
    {
        dctr(i);
    }
}

void ndctc(int column)
{
    int i = 0, j = 0, k = 0, l = 0;
    complex<double> up, down, product, x[max_size];
    complex<double>* W = new complex<double>[max_size];

    for (int i = 0; i < max_size; i++)
    {
        W[i] = exp(complex<double>(0, 2.0 * M_PI * i / max_size));
        //if (column == 5)printf("%lf %lf\n", W[i].real(), W[i].imag());
    }

    for (int i = 0; i < max_size; i++) {
        x[i] = foru_g[column][i];
        if (i > max_size / 2)x[i] = 0;
    }
    RaderReverse(x, max_size);

    for (i = 0; i < log(max_size) / log(2); ++i)
    {
        l = 1 << i;
        int count = 0;
        complex<double>* A = new complex<double>[max_size];
        for (j = 0; j < max_size; j += 2 * l)
        {
            for (k = 0; k < l * 2; ++k)
            {

                if (k < l) {
                    A[count++] = x[j + k] + x[j + k + l] * W[max_size * k / 2 / l];
                }
                else {
                    A[count++] = x[j + k - l] - x[j + k] * W[max_size * (k - l) / 2 / l];
                }
            }
        }
        for (int t = 0; t < max_size; t++)x[t] = A[t];
    }
    for (i = 0; i < max_size; i++)foru_g[column][i] = x[i];
    delete[] W;
}

void ndct()
{
    for (int i = 0; i < max_size; i++)
    {
        ndctc(i);
    }

    for (int i = 0; i < max_size; i++)
    {
        ndctr(i);
    }
}

double getmax(int a, int b)
{
    return a < b ? b : a;
}

double getmin(int a, int b)
{
    return a < b ? a : b;
}

int GetEncoderClsid(const WCHAR* format, CLSID* pClsid)
{
    UINT  num = 0;          // number of image encoders
    UINT  size = 0;         // size of the image encoder array in bytes

    ImageCodecInfo* pImageCodecInfo = NULL;

    GetImageEncodersSize(&num, &size);
    if (size == 0)
        return -1;  // Failure

    pImageCodecInfo = (ImageCodecInfo*)(malloc(size));
    if (pImageCodecInfo == NULL)
        return -1;  // Failure

    GetImageEncoders(num, size, pImageCodecInfo);

    for (UINT j = 0; j < num; ++j)
    {
        if (wcscmp(pImageCodecInfo[j].MimeType, format) == 0)
        {
            *pClsid = pImageCodecInfo[j].Clsid;
            free(pImageCodecInfo);
            return j;  // Success
        }
    }

    free(pImageCodecInfo);
    return -1;  // Failure
}


int main() {
    GdiplusStartupInput gdiplusstartupinput;
    ULONG_PTR gdiplustoken;
    GdiplusStartup(&gdiplustoken, &gdiplusstartupinput, NULL);
    int grey_vec[256] = { 0 };
    int num = -1;
    double maxn = -1, minn = 100000000;
    wstring infilename(L"lena512.bmp");
    Bitmap* bmp = new Bitmap(infilename.c_str());
    int height = bmp->GetHeight();
    int width = bmp->GetWidth();
    cout << "width " << width << ", height " << height << endl;
    Color color, tmp;
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            bmp->GetPixel(x, y, &color);
            int val = (int)(color.GetRed() * 0.299 + color.GetGreen() * 0.587 + color.GetBlue() * 0.114);
            num = val;
            grey_vec[num]++;
            tmp = Color(num, num, num);
            bmp->SetPixel(x, y, tmp);
        }
    }

    resize(bmp, height, width);
    dct();
    printf("%lf %lf\n", maxn, minn);
    CLSID pngClsid;
    GetEncoderClsid(L"image/bmp", &pngClsid);
    // bmp->Save(L"dct_origin.bmp", &pngClsid, NULL);
    Bitmap* ptr_bitmap = new Bitmap(max_size, max_size);
    for (int y = 3; y < max_size; y++)
        for (int x = 3; x < max_size; x++) {
            double real = foru_g[x][y].real(), im = foru_g[x][y].imag(), num1;
            num1 = (sqrt(real * real));
            maxn = getmax(maxn, num1);
            minn = getmin(minn, num1);
        }
    printf("%lf %lf\n", maxn, minn);
    for (int y = 3; y < max_size; y++) {
        for (int x = 3; x < max_size; x++) {
            double real = foru_g[x][y].real(), im = foru_g[x][y].imag(), num1;
            num1 = (sqrt(real * real));
            double val = (num1 - minn) / (maxn - minn);
            val *= 255;
            if (val < 0)val = 0;
            if (val > 255)val = 255;
            tmp = Color((int)val, (int)val, (int)val);
            ptr_bitmap->SetPixel(x, y, tmp);
        }
    }
    //ptr_bitmap->Save(L"dct.bmp", &pngClsid, NULL);
    minn = 1000000;
    maxn = -1111;
    for (int y = 0; y < max_size; y++)
        for (int x = 0; x < max_size; x++) {
            double real = foru_g[x][y].real(), im = foru_g[x][y].imag(), num1;
            num1 = log(sqrt(real * real) + 1);
            maxn = getmax(maxn, num1);
            minn = getmin(minn, num1);
        }
    printf("%lf %lf\n", maxn, minn);
    for (int y = 0; y < max_size; y++) {
        for (int x = 0; x < max_size; x++) {
            double real = foru_g[x][y].real(), im = foru_g[x][y].imag(), num1;
            num1 = log(sqrt(real * real) + 1);
            double val = (num1 - minn) / (maxn - minn);
            val *= 255;
            if (val < 0)val = 0;
            if (val > 255)val = 255;
            tmp = Color((int)val, (int)val, (int)val);
            ptr_bitmap->SetPixel(x, y, tmp);
        }
    }
    //ptr_bitmap->Save(L"dct_log.bmp", &pngClsid, NULL);

    ndct();
    minn = 1000000;
    maxn = -1;
    for (int y = 1; y < max_size; y++)
        for (int x = 1; x < max_size; x++) {
            double real = foru_g[x][y].real(), im = foru_g[x][y].imag(), num1;
            num1 = sqrt(real * real) * 100;
            //printf("%lf \n", num1);
            maxn = getmax(maxn, num1);
            minn = getmin(minn, num1);
        }
    printf("%lf %lf\n", minn, maxn);
    for (int y = 1; y < max_size; y++) {
        for (int x = 1; x < max_size; x++) {
            double real = foru_g[x][y].real(), im = foru_g[x][y].imag(), num1;
            num1 = (sqrt(real * real)) * 100;
            double val = (num1 - minn) / (maxn - minn);
            val *= 255;
            if (val < 0)val = 0;
            if (val > 255)val = 255;
            tmp = Color((int)val, (int)val, (int)val);
            ptr_bitmap->SetPixel(x, y, tmp);
        }
    }
    ptr_bitmap->Save(L"idct.bmp", &pngClsid, NULL);
    delete bmp;
    GdiplusShutdown(gdiplustoken);
    return 0;
}