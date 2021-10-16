// Matrix.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include  <omp.h>

float* pfAlign16(float* pf)
{
	int32_t i = reinterpret_cast<int32_t>(pf);

	i += (16 - (i & 0x0F)) & 0x0F;     //остаток от деления на 16, исключаем 16
	return reinterpret_cast<float *>(i);
}

void Init1(float* pf)
{
	for (int i = 0; i < 16; i++)
		pf[i] = i;
}

void TraceMatrix4x4(float *pf)
{
	for (int iRow = 0; iRow < 4; iRow++)
	{
		for (int iCol = 0; iCol < 4; iCol++)
		{
			printf("%8.1f ", pf[(iRow << 2) + iCol]);//на всё число 8 позиций, на дробную часть 1 позиция (красивый вывод)
		}
		printf("\n");
	}
	printf("\n");
}

void Transpmatrix4x4(float* pf1, float* pf2)
{
	for (int iRow = 0; iRow < 4; iRow++)
	{
		for (int iCol = 0; iCol < 4; iCol++)
		{
			pf2[(iCol << 2) + iRow] = pf1[(iRow << 2) + iCol];//читаем по строкам, потому что работаем механизм упреждающего чтения
		}
	}
}

void Transpmatrix4x4_asm(float* pf1, float* pf2)// с использованием SSE
{
	_asm
	{
		mov eax, pf1                // в переменную eax записывается * pf1
		movaps xmm0, [eax]          // 3  2  1  0
		movaps xmm1, [eax + 16]     // 7  6  5  4 (смещение даём в байтах)
		movaps xmm2, [eax + 2 * 16] //11 10  9  8
		movaps xmm3, [eax + 3 * 16] //15 14 13 12

		mov eax, pf2
		movaps xmm4, xmm0           // 3  2  1  0
		shufps xmm0, xmm1, 044h     // 5  4  1  0  (01 00 01 00)
		shufps xmm4, xmm1, 0EEh     // 7  6  3  2  (11 10 11 10)

		//0,2,3,4 - те регистры, которые уже заняты

		movaps xmm1, xmm2           // 11  10   9   8
		shufps xmm2, xmm3, 044h     // 13  12   9   8  (01 00 01 00)
		shufps xmm1, xmm3, 0EEh     // 15  14  11  10  (11 10 11 10)

																//0,4, 2,1

		movaps xmm3, xmm0           //  5  4  1  0
		shufps xmm0, xmm2, 088h     // 12  8  4  0  (10 00 10 00)
		shufps xmm3, xmm2, 0DDh     // 13  9  5  1  (11 01 11 01)

		movaps[eax], xmm0
		movaps[eax + 16], xmm3

		movaps xmm3, xmm4           //  7   6  3  2
		shufps xmm4, xmm1, 088h     // 14  10  6  2  (10 00 10 00)
		shufps xmm3, xmm1, 0DDh     // 15  11  7  3  (11 01 11 01)

		movaps[eax + 2 * 16], xmm4
		movaps[eax + 3 * 16], xmm3
	}
}

void MultMatrix4x4(float* pf1, float* pf2, float* pfRes)
{
	int i, j, k;
	float f;
	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			f = 0.f;
			for (k = 0; k < 4; k++)
			{
				f += pf1[i * 4 + k] * pf2[k * 4 + j];
			}
			pfRes[i * 4 + j] = f;
		}
	}
}

void MultMatrix4x4_asm(float* pf1, float* pf2, float* pfRes)
{
	Transpmatrix4x4_asm(pf2, pf2);
	_asm
	{
		mov    eax,  pf2               //transposed matrix addres to EAX 
		movaps xmm4, [eax]          
		movaps xmm5, [eax + 16]     
		movaps xmm6, [eax + 2 * 16] 
		movaps xmm7, [eax + 3 * 16] 

		mov    ecx, 4
		loop1:
			mov    eax, pf1
			add    pf1, 16
			movaps xmm0, [eax]       


			movaps xmm1, xmm0
			dpps   xmm1, xmm4, 11110001b
			movaps xmm2, xmm0
			dpps   xmm2, xmm5, 11110001b
			movaps xmm3, xmm0
			dpps   xmm3, xmm6, 11110001b
			dpps   xmm0, xmm7, 11110001b

			shufps xmm2, xmm1, 00000000b
			shufps xmm2, xmm2, 00000010b

			shufps xmm0, xmm3, 00000000b
			shufps xmm0, xmm0, 00000010b

 			shufps xmm2, xmm0, 11001100b

			mov    eax, pfRes
			movaps [eax], xmm2
			add    pfRes, 16
		loop loop1
	}

}

int main()
{
	printf("omp_get_wtick() = %lg%lg\n\n", omp_get_wtick());
	//for (int i = 0; i < 10; i++)
	//{
	//    float* pf = new float[10];
	//    float* pfA = pfAlign16(pf); //храним исходный указатель для удаления памяти
	//    printf("pf = %p; pfA = %p\n", pf, pfA);
	//    //delete[] pf;
	//}
	float* pf = new float[3 * 4 * 4 + 16];//3 матрицы 4*4 с выравниванием на 64 байта
	float *pf1 = pfAlign16(pf); //выровненный указатель
	float* pf2 = pf1 + 16; //&(pf1[16])
	float* pf3 = pf1 + 2 * 16;

	Init1(pf1);
	TraceMatrix4x4(pf1);

	Transpmatrix4x4(pf1, pf2);
	TraceMatrix4x4(pf2);

	Transpmatrix4x4_asm(pf1, pf3);
	TraceMatrix4x4(pf3);

	int iCount = 10000000;
	double dStart = omp_get_wtime(), dEnd;//Для замерки времени работы, возвращет результат в секундах
	for (int i = 0; i < iCount; i++)
		Transpmatrix4x4(pf1, pf2);
	dEnd = omp_get_wtime();
	printf("%d = Transpmatrix4x4 %lg\n", iCount, dEnd - dStart);

	dStart = omp_get_wtime();
	for (int i = 0; i < iCount; i++)
		Transpmatrix4x4_asm(pf1, pf2);
	dEnd = omp_get_wtime();
	printf("%d = Transpmatrix4x4_asm %lg\n", iCount, dEnd - dStart);

	Init1(pf1);
	Init1(pf2);
	MultMatrix4x4(pf1, pf2, pf3);
	TraceMatrix4x4(pf3);

	std::cout << std::endl << "------------------------------------------------" << std::endl;
	Init1(pf1);
	Init1(pf2);
	MultMatrix4x4_asm(pf1, pf2, pf3);
	TraceMatrix4x4(pf3);

	dStart = omp_get_wtime();
	for (int i = 0; i < iCount; i++)
		MultMatrix4x4_asm(pf1, pf2, pf3);
	dEnd = omp_get_wtime();
	printf("%d = MultMatrix4x4 %lg\n", iCount, dEnd - dStart);

	dStart = omp_get_wtime();
	for (int i = 0; i < iCount; i++)
		Transpmatrix4x4_asm(pf1, pf2);
	dEnd = omp_get_wtime();
	printf("%d = MultMatrix4x4_asm %lg\n", iCount, dEnd - dStart);

	delete[] pf;
	//std::cout << "Hello World!" << std::endl;
	return 0;
}