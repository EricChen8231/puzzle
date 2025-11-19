#include <iostream>
#include <vector>
#include <string>
#include "sudoku.h"

using namespace std;

int main()
{
    vector<vector<int>> test = {
        {1, 3, 4, 6, 7, 8, 2, 5, 0},
        {1, 3, 4, 6, 7, 9, 2, 5, 0},
        {1, 3, 4, 6, 8, 9, 2, 5, 0},
        {1, 3, 4, 7, 8, 9, 2, 5, 0},
        {1, 3, 6, 7, 8, 9, 2, 5, 0},
        {1, 4, 6, 7, 8, 9, 2, 5, 0},
        {3, 4, 6, 7, 8, 9, 2, 5, 0}};

    vector<vector<int>> grid = {
        {-1, -1, -1, -1, -1, -1, -1, 2, -1},
        {-1, -1, -1, -1, 2, -1, -1, -1, 5},
        {-1, 2, -1, -1, -1, -1, -1, -1, -1},
        {-1, -1, 0, -1, -1, -1, -1, -1, -1},
        {-1, -1, -1, -1, -1, -1, -1, -1, -1},
        {-1, -1, -1, 2, -1, -1, -1, -1, -1},
        {-1, -1, -1, -1, 0, -1, -1, -1, -1},
        {-1, -1, -1, -1, -1, 2, -1, -1, -1},
        {-1, -1, -1, -1, -1, -1, 5, -1, -1}};

    for (int i = 0; i < 7; i++)
    {
        SudokuSolver solver(grid, test[i]);
        string name = "file" + to_string(0) + ".txt";
        solver.solve(name);
        solver.printSolution();
        }

        return 0;
}   