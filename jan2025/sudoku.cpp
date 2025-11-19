    #include "sudoku.h"
    #include <iostream>
    #include <cmath>
    #include <numeric>
    #include <algorithm>
    #include <fstream>
    #include <unordered_map>
    #include <bitset>
    using namespace std;
    #define GRID_SIZE 9
    std::bitset<GRID_SIZE> rowUsed[GRID_SIZE], colUsed[GRID_SIZE], boxUsed[GRID_SIZE];
    long long counter = 0;


    SudokuSolver::SudokuSolver(const vector<vector<int>> &initialGrid, const vector<int> digit)
    {
        for (vector<int>::size_type i = 0; i < digit.size(); ++i){
            digitToIndex[digit[i]] = i;

        }
        grid = initialGrid;

        // Initialize the bitsets using the custom digits
        for (int i = 0; i < GRID_SIZE; i++)
        {
            for (int j = 0; j < GRID_SIZE; j++)
            {
                if (grid[i][j] != -1)
                {
                    int num = grid[i][j];
                    int index = digitToIndex[num];
                    rowUsed[i].set(index);                    
                    colUsed[j].set(index);                    
                    boxUsed[(i / 3) * 3 + (j / 3)].set(index);
                }
            }
        }
        digits = digit;

        maxGCD = 12345679;
        bestMiddleRow = vector<int>(GRID_SIZE, -1);
        solutions = vector<vector<int>>(GRID_SIZE, vector<int>(GRID_SIZE, 0));
    }

    int SudokuSolver::gcd(int a, int b)
    {
        return std::gcd(a, b);
    }
    // Calculate the GCD of the numbers formed by the rows
    int SudokuSolver::calculateGCDOfRows(int num)
    {
        int rowGCD = 0;
        for (int i = 0; i < num; ++i)
        {
            int rowNumber = 0;
            for (int j = 0; j < GRID_SIZE; ++j)
            {
                rowNumber = rowNumber * 10 + grid[i][j];
            }
            rowGCD = (i == 0) ? rowNumber : gcd(rowGCD, rowNumber);
            if  (maxGCD >= rowGCD){
                return 0 ;
            }
        }
        return rowGCD;
    }

    bool SudokuSolver::isValid(int row, int col, int num1)
    {
        int num = digitToIndex[num1];
        int boxIndex = (row / 3) * 3 + col / 3; 
        return !rowUsed[row][num] && !colUsed[col][num] && !boxUsed[boxIndex][num];
    }

    void SudokuSolver::placeNumber(int row, int col, int num1)
    {
        int num = digitToIndex[num1];
        int boxIndex = (row / 3) * 3 + col / 3;
        rowUsed[row].set(num);
        colUsed[col].set(num);
        boxUsed[boxIndex].set(num);
    }

    void SudokuSolver::removeNumber(int row, int col, int num1)
    {
        int num = digitToIndex[num1];
        int boxIndex = (row / 3) * 3 + col / 3;
        rowUsed[row].reset(num);
        colUsed[col].reset(num);
        boxUsed[boxIndex].reset(num);
    }

    void SudokuSolver::backtrack(int row, int col)
    {
        counter++;
        if (counter % 1000000000 == 0)
        {
            cout << "Recursion count: " << counter << endl;
            for (const auto &row : grid)
            {
                for (int num : row)
                {
                    cout << num << " ";
                }
                cout << endl;
            }
        }
        // if (grid[0][0] == 1){
        //     return;
        // }
        // if (grid[0][0] == 4)
        // {
        //     return;
        // }
        // if (grid[0][0] == 3)
        // {
        //     return;
        // }
        // if (grid[0][0] == 6)
        // {
        //     return;
        // }
        if (row == GRID_SIZE)
        {
            int temp = calculateGCDOfRows(9);
            if (temp > maxGCD)
            {
                maxGCD = temp;
                solutions = grid;
                for (const auto &row : solutions)
                {
                    for (int num : row)
                    {
                        cout << num << " ";
                    }
                    cout << endl;
                }
                cout <<"Current Max GCD: "<< maxGCD;
                cout << endl;
                saveSol(file1);
            }
            return;
        }

        if (col == GRID_SIZE)
        {
            if (row > 1 ){
                int temp = calculateGCDOfRows(row);
                if (maxGCD >= temp)
                {
                    return;
                }
            }
            backtrack(row + 1, 0);
            return;
        }

        if (grid[row][col] != -1)
        {
            backtrack(row, col + 1);
            return;
        }
        for (int num : digits)
        {

            if (isValid(row, col, num))
            {
                grid[row][col] = num;
                placeNumber(row, col, num);
                backtrack(row, col + 1);
                grid[row][col] = -1;
                removeNumber(row, col, num);
            }
        }
        return;
    }

    void SudokuSolver::saveSol(const std::string &filename)
    {
        ofstream file(filename);
        if (file.is_open())
        {
            for (const auto &row : solutions)
            {
                for (int num : row)
                {
                    if (num == -1)
                        file << ". ";
                    else
                        file << num << " ";
                }
                file << endl;
            }
            file << endl;
            file.close();
        }
        else
        {
            cerr << "Unable to open file: " << filename << endl;
        }
    }

    void SudokuSolver::printSolution()
    {
        cout << "Solution found: " << endl;
        for (const auto &row : solutions)
        {
            for (int num : row)
            {
                cout << num << " ";
            }
            cout << endl;
        }
    }

    void SudokuSolver::solve(const string &file)
    {
        file1 = file; 
        backtrack(0, 0);
        printSolution();
        saveSol(file1);
    }