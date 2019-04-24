#include "WaterMark.h"

WaterMark::WaterMark()
{
}

WaterMark::~WaterMark()
{
}
//ԭʼˮӡλ��
void WaterMark::setM(int m_m)
{
	m = m_m;
}

//��Ƭ����
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

/*ΪVecA��ֵ*/
void WaterMark::createA()
{
	srand(time(NULL));//���������ֵ
	for (int i = 0; i < m; i++)
	{
		vecA.push_back(rand() % 2);//��A��ֵΪ����Ϊm��ֵΪ0��1���������
	}
}

/*ΪVecB��ֵ������д�뵽�ļ�Wb.txt��*/
void WaterMark::createWB()
{
	//����ֵͼ��ת���õ���36*33��01������ɳ���Ϊ1188��ԭʼ��ˮӡ����vecB
	const int SEQ_LENGTH = 1188;//36*33=1188
	double *B = new double[SEQ_LENGTH];
	ifstream Wfile;
	Wfile.open("D:\\firejq\\����\\Watermarking\\Txt\\AMesh\\flower.txt", ios_base::in);
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


	//��ˮӡ����д���ļ���
	ofstream Wbfile;
	Wbfile.open("D:\\firejq\\����\\Watermarking\\Txt\\Wb.txt", ios_base::out);
	if (Wbfile)
	{
		for (int i = 0; i < SEQ_LENGTH/*vecA.size() * c*/; i++)
		{
			Wbfile << vecB[i] << " ";
		}
	}
	Wbfile.close();
}

/*ΪP��ֵ����д�뵽�ļ�P.txt��*/
void WaterMark::createP()
{
	srand(key);//���������ֵ
	for (int i = 0; i < m*c; i++)
	{
		P.push_back(rand() % 2 == 0 ? -1 : 1);
	}

	//��Pд���ļ���
	ofstream Pfile;
	Pfile.open("D:\\firejq\\����\\Watermarking\\Txt\\P.txt", ios_base::out);
	if (Pfile)
	{
		for (int i = 0; i < m*c; i++)
		{
			Pfile << P[i] << " ";
		}
	}
	Pfile.close();
}