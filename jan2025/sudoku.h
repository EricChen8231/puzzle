#ifndef SUDOKU_H
#define SUDOKU_H

#include <vector>
#include <fstream>
#include <unordered_map>


using namespace std;

class SudokuSolver
{
private:
    vector<vector<int>> grid;
    vector<vector<int>> empty;
    vector<vector<int>> solboard;
    vector<int> digits;
    int maxGCD;
    vector<int> bestMiddleRow;
    vector<vector<int>> solutions;
    unordered_map<int, int> digitToIndex;
    string file1;

public:
    SudokuSolver(const vector<vector<int>> &initialGrid, const vector<int> digits);
    void printSolution();
    bool isValid(int row, int col, int num);
    void placeNumber(int row, int col, int num1);
    void removeNumber(int row, int col, int num1);
    int calculateGCDOfRows(int row);
    int gcd(int a, int b);
    void backtrack(int row, int col);
    void solve(const string &filename);
    void saveSol(const string &file);
};

#endif // SUDOKU_H
