#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <unordered_set>
#include <algorithm>

// Directions for word search: right, down-right, down, down-left, left, up-left, up, up-right
const int dx[8] = {1, 1, 0, -1, -1, -1, 0, 1};
const int dy[8] = {0, 1, 1, 1, 0, -1, -1, -1};

// Structure to represent a word placement
struct WordPlacement {
    std::string word;
    int row, col, dir;
};

class InverseWordSearch {
private:
    int width, height;
    std::vector<std::string> requiredWords;
    std::vector<std::string> forbiddenWords;
    std::vector<std::vector<char>> grid;
    std::vector<std::vector<std::vector<char>>> solutions;
    bool findAllSolutions;
    // Set of strings to represent solutions for faster duplicate checking
    std::set<std::string> solutionHashes;

public:
    InverseWordSearch(int w, int h, 
                      const std::vector<std::string>& required, 
                      const std::vector<std::string>& forbidden,
                      bool findAll) 
        : width(w), height(h), requiredWords(required), forbiddenWords(forbidden), 
          findAllSolutions(findAll) {
        // Initialize grid with placeholder characters
        grid.resize(height, std::vector<char>(width, '.'));
    }

    bool solve() {
        // Clear any existing solutions
        solutions.clear();
        solutionHashes.clear();
        
        // Special case: If no required words, just fill with valid letters
        if (requiredWords.empty()) {
            return fillEmptyCells(0, 0);
        }
        
        // Pre-process forbidden words for faster checking
        // Extract single-letter forbidden words for quick checking
        std::unordered_set<char> forbiddenLetters;
        std::vector<std::string> multiLetterForbidden;
        
        for (const auto& word : forbiddenWords) {
            if (word.length() == 1) {
                forbiddenLetters.insert(word[0]);
            } else {
                multiLetterForbidden.push_back(word);
            }
        }
        
        // Try to place all the required words
        std::vector<WordPlacement> placements;
        bool result = placeWords(0, placements, forbiddenLetters, multiLetterForbidden);
        
        return result || !solutions.empty();
    }

    const std::vector<std::vector<std::vector<char>>>& getSolutions() const {
        return solutions;
    }

private:
    // Generate a string hash of the current grid for duplicate checking
    std::string gridHash() const {
        std::string hash;
        hash.reserve(width * height);
        for (const auto& row : grid) {
            hash.append(row.begin(), row.end());
        }
        return hash;
    }

    // Try to place a word at a specific position and direction
    bool canPlaceWord(const std::string& word, int row, int col, int dir) const {
        // Check if the word fits within the grid boundaries
        int endRow = row + (word.length() - 1) * dy[dir];
        int endCol = col + (word.length() - 1) * dx[dir];
        
        if (endRow < 0 || endRow >= height || endCol < 0 || endCol >= width) {
            return false;
        }
        
        // Check if the word can be placed at this position
        for (size_t i = 0; i < word.length(); i++) {
            int r = row + i * dy[dir];
            int c = col + i * dx[dir];
            
            // If the cell is already filled and doesn't match the word's letter, can't place
            if (grid[r][c] != '.' && grid[r][c] != word[i]) {
                return false;
            }
        }
        
        return true;
    }
    
    // Place a word in the grid
    void placeWord(const std::string& word, int row, int col, int dir) {
        for (size_t i = 0; i < word.length(); i++) {
            int r = row + i * dy[dir];
            int c = col + i * dx[dir];
            grid[r][c] = word[i];
        }
    }
    
    // Remove a word from the grid (reset to '.' where no other word overlaps)
    void removeWord(const std::string& word, int row, int col, int dir, 
                   const std::vector<WordPlacement>& placements) {
        // First place '.' in every cell the word occupies
        for (size_t i = 0; i < word.length(); i++) {
            int r = row + i * dy[dir];
            int c = col + i * dx[dir];
            grid[r][c] = '.';
        }
        
        // Then restore other words' letters
        for (const auto& placement : placements) {
            // Skip the current word being removed
            if (placement.word == word && placement.row == row && 
                placement.col == col && placement.dir == dir) {
                continue;
            }
            
            // Restore letters from other placements
            for (size_t i = 0; i < placement.word.length(); i++) {
                int r = placement.row + i * dy[placement.dir];
                int c = placement.col + i * dx[placement.dir];
                grid[r][c] = placement.word[i];
            }
        }
    }
    
    // Check if a word exists in the grid
    bool wordExists(const std::string& word) const {
        // For single-letter words, need to check if they appear in the grid
        if (word.length() == 1) {
            // Check if this letter appears anywhere in the grid
            for (int row = 0; row < height; row++) {
                for (int col = 0; col < width; col++) {
                    if (grid[row][col] == word[0]) {
                        return true;
                    }
                }
            }
            return false;
        }
        
        // First letter lookup optimization - create a list of positions where the first letter appears
        std::vector<std::pair<int, int>> firstLetterPositions;
        for (int row = 0; row < height; row++) {
            for (int col = 0; col < width; col++) {
                if (grid[row][col] == word[0]) {
                    firstLetterPositions.emplace_back(row, col);
                }
            }
        }
        
        // Now check only those positions
        for (const auto& pos : firstLetterPositions) {
            int row = pos.first;
            int col = pos.second;
            
            for (int dir = 0; dir < 8; dir++) {
                // Check word at this position and direction
                int endRow = row + (word.length() - 1) * dy[dir];
                int endCol = col + (word.length() - 1) * dx[dir];
                
                // Skip if out of bounds
                if (endRow < 0 || endRow >= height || endCol < 0 || endCol >= width) {
                    continue;
                }
                
                bool match = true;
                for (size_t i = 1; i < word.length(); i++) { // Start from 1 since we already matched the first letter
                    int r = row + i * dy[dir];
                    int c = col + i * dx[dir];
                    
                    if (grid[r][c] != word[i]) {
                        match = false;
                        break;
                    }
                }
                
                if (match) return true;
            }
        }
        
        return false;
    }
    
    // Check if the grid is valid (all required words present, no forbidden words)
    bool isValidGrid() const {
        // Check required words
        for (const auto& word : requiredWords) {
            if (!wordExists(word)) {
                return false;
            }
        }
        
        // Check forbidden words
        for (const auto& word : forbiddenWords) {
            if (wordExists(word)) {
                return false;
            }
        }
        
        return true;
    }
    
    // Fill empty cells with valid letters
    bool fillEmptyCells(int row, int col, 
                       const std::unordered_set<char>& forbiddenLetters = std::unordered_set<char>(),
                       const std::vector<std::string>& multiLetterForbidden = std::vector<std::string>()) {
        // If we've filled the entire grid, check if it's valid
        if (row == height) {
            if (isValidGrid()) {
                // Use hash for faster duplicate detection
                std::string hash = gridHash();
                
                if (solutionHashes.find(hash) == solutionHashes.end()) {
                    solutionHashes.insert(hash);
                    solutions.push_back(grid);
                    return !findAllSolutions;  // Stop if we only need one solution
                }
            }
            return false;
        }
        
        // Calculate next position
        int nextRow = row;
        int nextCol = col + 1;
        if (nextCol == width) {
            nextRow++;
            nextCol = 0;
        }
        
        // If the cell is already filled, move to the next one
        if (grid[row][col] != '.') {
            return fillEmptyCells(nextRow, nextCol, forbiddenLetters, multiLetterForbidden);
        }
        
        // Try each possible letter
        for (char c = 'a'; c <= 'z'; c++) {
            // Quick check for forbidden single letters
            if (forbiddenLetters.find(c) != forbiddenLetters.end()) {
                continue; // Skip this letter
            }
            
            grid[row][col] = c;
            
            // Check if this letter would create a forbidden word
            bool createsForbiddenWord = false;
            
            // Only check multi-letter forbidden words when we're near the end of filling the grid
            // This optimization avoids unnecessary checks when the grid is still mostly empty
            if ((row >= height-1) || (row == height-1 && col >= width-2)) {
                for (const auto& word : multiLetterForbidden) {
                    if (wordExists(word)) {
                        createsForbiddenWord = true;
                        break;
                    }
                }
            }
            
            if (!createsForbiddenWord) {
                bool result = fillEmptyCells(nextRow, nextCol, forbiddenLetters, multiLetterForbidden);
                if (result && !findAllSolutions) {
                    return true;  // If we only need one solution and found it, return immediately
                }
                // If we need all solutions, or no solution was found, continue with the next letter
            }
        }
        
        // Backtrack
        grid[row][col] = '.';
        return false;
    }
    
    // Recursive function to place required words in the grid
    bool placeWords(int wordIndex, std::vector<WordPlacement>& placements,
                   const std::unordered_set<char>& forbiddenLetters,
                   const std::vector<std::string>& multiLetterForbidden) {
        // If all words are placed, fill the remaining cells
        if (wordIndex == requiredWords.size()) {
            return fillEmptyCells(0, 0, forbiddenLetters, multiLetterForbidden);
        }
        
        const std::string& word = requiredWords[wordIndex];
        bool foundSolution = false;
        
        // Try all possible positions and directions for this word
        for (int row = 0; row < height && (!foundSolution || findAllSolutions); row++) {
            for (int col = 0; col < width && (!foundSolution || findAllSolutions); col++) {
                for (int dir = 0; dir < 8 && (!foundSolution || findAllSolutions); dir++) {
                    if (canPlaceWord(word, row, col, dir)) {
                        // Place the word
                        placeWord(word, row, col, dir);
                        placements.push_back({word, row, col, dir});
                        
                        // Try to place the next word
                        bool result = placeWords(wordIndex + 1, placements, forbiddenLetters, multiLetterForbidden);
                        foundSolution = foundSolution || result;
                        
                        // If we're only looking for one solution and we found it, return early
                        if (result && !findAllSolutions) {
                            return true;
                        }
                        
                        // Backtrack - remove the word and try another position
                        placements.pop_back();
                        removeWord(word, row, col, dir, placements);
                    }
                }
            }
        }
        
        return foundSolution;
    }
};

// Function to print a grid
void printGrid(std::ostream& out, const std::vector<std::vector<char>>& grid) {
    out << "Board: " << std::endl;
    for (const auto& row : grid) {
        out << "  ";
        for (char c : row) {
            out << c;
        }
        out << std::endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " input_file output_file mode" << std::endl;
        return 1;
    }
    
    std::string inputFile = argv[1];
    std::string outputFile = argv[2];
    std::string mode = argv[3];
    
    // Read input
    std::ifstream inFile(inputFile);
    if (!inFile) {
        std::cerr << "Error opening input file: " << inputFile << std::endl;
        return 1;
    }
    
    int width, height;
    inFile >> width >> height;
    inFile.ignore(); // Skip the newline
    
    std::vector<std::string> requiredWords;
    std::vector<std::string> forbiddenWords;
    
    std::string line;
    while (std::getline(inFile, line)) {
        if (line.empty()) {
            continue;
        }
        
        char type = line[0];
        std::string word = line.substr(2);
        
        // Trim trailing whitespace
        word.erase(word.find_last_not_of(" \n\r\t") + 1);
        
        if (type == '+') {
            requiredWords.push_back(word);
        } else if (type == '-') {
            forbiddenWords.push_back(word);
        }
    }
    
    inFile.close();
    
    // Sort required words by length (descending) for better pruning
    std::sort(requiredWords.begin(), requiredWords.end(), 
              [](const std::string& a, const std::string& b) {
                  return a.length() > b.length();
              });
    
    // Check if the puzzle is solvable
    bool impossiblePuzzle = false;
    for (const auto& word : requiredWords) {
        if (word.length() > std::max(width, height)) {
            impossiblePuzzle = true;
            break;
        }
    }
    
    // Solve the puzzle
    std::vector<std::vector<std::vector<char>>> solutions;
    
    if (!impossiblePuzzle) {
        bool findAllSols = (mode == "all_solutions");
        InverseWordSearch solver(width, height, requiredWords, forbiddenWords, findAllSols);
        solver.solve();
        solutions = solver.getSolutions();
    }
    
    // Write output
    std::ofstream outFile(outputFile);
    if (!outFile) {
        std::cerr << "Error opening output file: " << outputFile << std::endl;
        return 1;
    }
    
    if (solutions.empty()) {
        outFile << "No solutions found" << std::endl;
    } else if (mode == "all_solutions") {
        outFile << solutions.size() << " solution(s)" << std::endl;
        for (const auto& sol : solutions) {
            printGrid(outFile, sol);
        }
    } else {
        printGrid(outFile, solutions[0]);
    }
    
    outFile.close();
    
    return 0;
}