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
    size_t length; // Cache the word length for efficiency
};

class InverseWordSearch {
private:
    int width, height;
    std::vector<std::string> requiredWords;
    std::vector<std::string> forbiddenWords;
    std::vector<char> grid; // Flat grid for better cache locality
    std::vector<std::vector<std::vector<char>>> solutions; // Changed from 2D to 3D vector
    bool findAllSolutions;
    // Set of strings to represent solutions for faster duplicate checking
    std::set<std::string> solutionHashes;

    // Helper function to convert 2D coordinates to 1D index
    inline int flatIndex(int row, int col) const {
        return row * width + col;
    }

public:
    InverseWordSearch(int w, int h,
                      const std::vector<std::string>& required,
                      const std::vector<std::string>& forbidden,
                      bool findAll)
        : width(w), height(h), requiredWords(required), forbiddenWords(forbidden),
          findAllSolutions(findAll) {
        // Initialize flat grid with placeholder characters
        grid.resize(width * height, '.');
    }

    bool solve() {
        // Clear any existing solutions
        solutions.clear();
        solutionHashes.clear();
        
        // Try to place all the required words
        std::vector<WordPlacement> placements;
        return placeWords(0, placements);
    }

    const std::vector<std::vector<std::vector<char>>>& getSolutions() const {
        return solutions;
    }

private:
    // Generate a string hash of the current grid for duplicate checking
    std::string gridHash() const {
        // Create hash directly from flat grid
        return std::string(grid.begin(), grid.end());
    }

    // Try to place a word at a specific position and direction
    inline bool canPlaceWord(const std::string& word, int row, int col, int dir) const {
        const size_t wordLen = word.length();
        
        // Pre-compute end coordinates
        const int endRow = row + (wordLen - 1) * dy[dir];
        const int endCol = col + (wordLen - 1) * dx[dir];
        
        // Combined boundary check
        if (endRow < 0 || endRow >= height || endCol < 0 || endCol >= width) {
            return false;
        }
        
        // Calculate the starting flat index and the direction offset
        const int startIdx = flatIndex(row, col);
        const int dirOffset = dy[dir] * width + dx[dir];
        
        // Check if the word can be placed at this position
        int idx = startIdx;
        for (size_t i = 0; i < wordLen; i++) {
            // If the cell is already filled and doesn't match the word's letter, can't place
            if (grid[idx] != '.' && grid[idx] != word[i]) {
                return false;
            }
            idx += dirOffset;
        }
        
        return true;
    }

    // Place a word in the grid
    inline void placeWord(const std::string& word, int row, int col, int dir) {
        const size_t wordLen = word.length(); // Cache word length
        
        for (size_t i = 0; i < wordLen; i++) {
            const int r = row + i * dy[dir];
            const int c = col + i * dx[dir];
            grid[flatIndex(r, c)] = word[i];
        }
    }

    // Remove a word from the grid (reset to '.' where no other word overlaps)
    void removeWord(const std::string& word, int row, int col, int dir, 
                    const std::vector<WordPlacement>& placements) {
        const size_t wordLen = word.length(); // Cache word length
        
        // First place '.' in every cell the word occupies
        for (size_t i = 0; i < wordLen; i++) {
            const int r = row + i * dy[dir];
            const int c = col + i * dx[dir];
            grid[flatIndex(r, c)] = '.';
        }
        
        // Then restore other words' letters
        for (const auto& placement : placements) {
            // Skip the current word being removed
            if (placement.word == word && placement.row == row && 
                placement.col == col && placement.dir == dir) {
                continue;
            }
            
            // Restore letters from other placements
            for (size_t i = 0; i < placement.length; i++) {
                const int r = placement.row + i * dy[placement.dir];
                const int c = placement.col + i * dx[placement.dir];
                grid[flatIndex(r, c)] = placement.word[i];
            }
        }
    }

    // Check if a word exists in the grid
    bool wordExists(const std::string& word) const {
        const size_t wordLen = word.length(); // Cache word length
        const char firstChar = word[0]; // Cache first character
        
        // For single-letter words, need to check if they appear in the grid
        if (wordLen == 1) {
            // Check if this letter appears anywhere in the grid
            for (int idx = 0; idx < width * height; idx++) {
                if (grid[idx] == firstChar) {
                    return true;
                }
            }
            return false;
        }
        
        // First letter lookup optimization - create a list of positions where the first letter appears
        std::vector<std::pair<int, int>> firstLetterPositions;
        firstLetterPositions.reserve(width * height / 10); // Rough estimate
        
        for (int row = 0; row < height; row++) {
            for (int col = 0; col < width; col++) {
                if (grid[flatIndex(row, col)] == firstChar) {
                    firstLetterPositions.emplace_back(row, col);
                }
            }
        }
        
        // Now check only those positions
        for (const auto& pos : firstLetterPositions) {
            const int row = pos.first;
            const int col = pos.second;
            
            for (int dir = 0; dir < 8; dir++) {
                // Check word at this position and direction
                const int endRow = row + (wordLen - 1) * dy[dir];
                const int endCol = col + (wordLen - 1) * dx[dir];
                
                // Skip if out of bounds
                if (endRow < 0 || endRow >= height || endCol < 0 || endCol >= width) {
                    continue;
                }
                
                bool match = true;
                for (size_t i = 1; i < wordLen; i++) { // Start from 1 since we already matched the first letter
                    const int r = row + i * dy[dir];
                    const int c = col + i * dx[dir];
                    
                    if (grid[flatIndex(r, c)] != word[i]) {
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
                    
                    // Convert flat grid to 2D for output
                    std::vector<std::vector<char>> gridSolution(height, std::vector<char>(width));
                    for (int r = 0; r < height; r++) {
                        for (int c = 0; c < width; c++) {
                            gridSolution[r][c] = grid[flatIndex(r, c)];
                        }
                    }
                    
                    solutions.push_back(gridSolution);
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
        if (grid[flatIndex(row, col)] != '.') {
            return fillEmptyCells(nextRow, nextCol, forbiddenLetters, multiLetterForbidden);
        }
        
        // Try each possible letter
        for (char c = 'a'; c <= 'z'; c++) {
            // Quick check for forbidden single letters
            if (forbiddenLetters.find(c) != forbiddenLetters.end()) {
                continue; // Skip this letter
            }
            
            grid[flatIndex(row, col)] = c;
            
            // Check if placing this letter would create a forbidden word
            // This optimization was causing issues - check all forbidden words always
            bool createsForbiddenWord = false;
            for (const auto& word : forbiddenWords) {
                if (wordExists(word)) {
                    createsForbiddenWord = true;
                    break;
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
        grid[flatIndex(row, col)] = '.';
        return false;
    }

    // Recursive function to place required words in the grid
    bool placeWords(int wordIndex, std::vector<WordPlacement>& placements) {
        // If all words are placed, fill the remaining cells
        if (wordIndex == requiredWords.size()) {
            return fillEmptyCells(0, 0);
        }
        
        const std::string& word = requiredWords[wordIndex];
        const size_t wordLen = word.length();
        bool foundSolution = false;
        
        // Fixed-size array to track changes - no heap allocations
        const int MAX_CHANGES = 20;
        std::pair<int, char> changes[MAX_CHANGES];
        int numChanges;
        
        // Pre-compute direction offsets in flat grid coordinates
        // This eliminates repeated calculations in the inner loops
        const int dirOffsets[8] = {
            1,              // right
            width + 1,      // down-right
            width,          // down
            width - 1,      // down-left
            -1,             // left
            -width - 1,     // up-left
            -width,         // up
            -width + 1      // up-right
        };
        
        // Try all possible positions and directions
        for (int dir = 0; dir < 8; dir++) {
            // Get pre-computed direction offset
            const int dirOffset = dirOffsets[dir];
            
            // Boundary calculations
            const int dirMaxOffset = static_cast<int>(wordLen) - 1;
            int minRow = (dy[dir] < 0) ? dirMaxOffset : 0;
            int minCol = (dx[dir] < 0) ? dirMaxOffset : 0;
            int maxRow = height - ((dy[dir] > 0) ? dirMaxOffset : 0);
            int maxCol = width - ((dx[dir] > 0) ? dirMaxOffset : 0);
            
            for (int row = minRow; row < maxRow; row++) {
                for (int col = minCol; col < maxCol; col++) {
                    // Fast in-place canPlaceWord check
                    int idx = flatIndex(row, col);
                    bool canPlace = true;
                    
                    // Check if all cells are available
                    for (size_t i = 0; i < wordLen; i++) {
                        if (grid[idx] != '.' && grid[idx] != word[i]) {
                            canPlace = false;
                            break;
                        }
                        idx += dirOffset;
                    }
                    
                    if (canPlace) {
                        // Place word and track changes
                        numChanges = 0;
                        idx = flatIndex(row, col);
                        
                        for (size_t i = 0; i < wordLen; i++) {
                            if (grid[idx] != word[i]) {
                                changes[numChanges++] = {idx, grid[idx]};
                                grid[idx] = word[i];
                            }
                            idx += dirOffset;
                        }
                        
                        // Add to placements
                        placements.push_back({word, row, col, dir, wordLen});
                        
                        // Try to place the next word
                        bool result = placeWords(wordIndex + 1, placements);
                        foundSolution = foundSolution || result;
                        
                        if (result && !findAllSolutions) {
                            return true;
                        }
                        
                        // Backtrack - restore only changed cells
                        placements.pop_back();
                        for (int i = 0; i < numChanges; i++) {
                            grid[changes[i].first] = changes[i].second;
                        }
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
    // Reserve space for vectors based on typical values
    requiredWords.reserve(20);
    forbiddenWords.reserve(20);
    
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