#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mat.h"

int main()
{
        srand(time(NULL));
        
        Mat M(4,4);

        // Filling the matrix with random elements between 0 and 1
        for (uint8_t i = 0; i < M.rows(); ++i) {
                for (uint8_t j = 0; j < M.cols(); ++j) {
                        M.at(i,j) = (rand() % 100) / 100.0;
                }
        }

        // Assigning the first element to 0 to check that pivot element is
        // correctly chosen
        M.at(0,0) = 0;

        Mat M_i = M.inv();
        Mat M_tr = M.transpose();

        printf("A Matrix : \n");
        M.print();
        printf("\nIt's Inverse : \n");
        M_i.print();
        printf("\nIt's Transpose : \n");
        M_tr.print();

        Mat result = M * M_i;
        printf("\nProduct of M and M_i : \n");
        result.print();
        printf("\n");

        float a1(1), b1(1), c1(1), d1(1);
        float a2(-1), b2(1), c2(2), d2(0);
        float a3(1), b3(-1), c3(-1), d3(1);

        Mat A(3,3);
        A.at(0,0) = a1; A.at(0,1) = b1; A.at(0,2) = c1;
        A.at(1,0) = a2; A.at(1,1) = b2; A.at(1,2) = c2;
        A.at(2,0) = a3; A.at(2,1) = b3; A.at(2,2) = c3;

        Vec3<float> b(d1, d2, d3);

        Mat A_i = A.inv();
        Mat x = A_i * b;
        printf("Solution : \n");
        x.print();
        printf("\n");

        return 0;
}

void _Error_Handler(const char *file, size_t line)
{
        printf("File : %s, Line : %d\n", file, line);
}
