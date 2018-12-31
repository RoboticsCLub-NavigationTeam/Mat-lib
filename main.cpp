#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mat.h"

int main()
{
        srand(time(NULL));
        
        Mat M(5,5);

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

        return 0;
}

void _Error_Handler(const char *file, size_t line)
{
        printf("File : %s, Line : %d\n", file, line);
}
