#include "WaterMark.h"

WaterMark::WaterMark()
{
}

WaterMark::~WaterMark()
{
}
//原始水印位数
void WaterMark::setM(int m_m)
{
	m = m_m;
}

//码片速率
void WaterMark::setC(int m_c)
{
	c = m_c;
}


void WaterMark::setAlpha(double a)
{
	alpha = a;
}

void WaterMark::setKey(int k)
{
	key = k;
}

/*为VecA赋值*/
void WaterMark::createA()
{
	srand(time(NULL));//随机数种子值
	for (int i = 0; i < m; i++)
	{
		vecA.push_back(rand() % 2);//将A赋值为长度为m，值为0或1的随机序列
	}
}

/*为VecB赋值，并将写入到文件Wb.txt中*/
void WaterMark::createWB()
{
	//将二值图像转换得到的36*33的01矩阵读成长度为1188的原始的水印序列vecB
	const int SEQ_LENGTH = 1188;//36*33=1188
	double *B = new double[SEQ_LENGTH];
	ifstream Wfile;
	Wfile.open("D:\\firejq\\毕设\\Watermarking\\Txt\\AMesh\\flower.txt", ios_base::in);
	if (Wfile)
	{
		// read into memory
		Wfile.seekg(0, Wfile.end);
		int length = Wfile.tellg();
		Wfile.seekg(0, Wfile.beg);

		char *buffer = new char[length];
		Wfile.read(buffer, length);
		Wfile.close();

		// parse into array
		std::istringstream iss(buffer);
		int i = 0;
		while (i < SEQ_LENGTH)
		{
			iss >> B[i++];
		}
		delete[] buffer;
		// print or use it.
	}
	Wfile.close();

	for (int i = 0; i < SEQ_LENGTH/*vecA.size() * c*/; i++)
	{
		vecB.push_back(B[i]);
		if (vecB[i] == 0)
			vecB[i] = -1;
	}


	//将水印序列写入文件中
	ofstream Wbfile;
	Wbfile.open("D:\\firejq\\毕设\\Watermarking\\Txt\\Wb.txt", ios_base::out);
	if (Wbfile)
	{
		for (int i = 0; i < SEQ_LENGTH/*vecA.size() * c*/; i++)
		{
			Wbfile << vecB[i] << " ";
		}
	}
	Wbfile.close();
}

/*为P赋值，并写入到文件P.txt中*/
void WaterMark::createP()
{
	srand(key);//随机数种子值
	for (int i = 0; i < m*c; i++)
	{
		P.push_back(rand() % 2 == 0 ? -1 : 1);
	}

	//将P写入文件中
	ofstream Pfile;
	Pfile.open("D:\\firejq\\毕设\\Watermarking\\Txt\\P.txt", ios_base::out);
	if (Pfile)
	{
		for (int i = 0; i < m*c; i++)
		{
			Pfile << P[i] << " ";
		}
	}
	Pfile.close();
}