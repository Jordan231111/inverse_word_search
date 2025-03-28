#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <memory>
#include <chrono>

// --- Trie Data Structure for efficient word checking ---
class TrieNode {
public:
    bool isEndOfWord;
    TrieNode* children[26];
    
    TrieNode() : isEndOfWord(false) {
        for (int i = 0; i < 26; i++) {
            children[i] = nullptr;
        }
    }
    
    ~TrieNode() {
        for (int i = 0; i < 26; i++) {
            delete children[i];
        }
    }
};

class Trie {
private:
    TrieNode* root;
    
public:
    Trie() : root(new TrieNode()) {}
    
    ~Trie() {
        delete root;
    }
    
    void insert(const std::string& word) {
        TrieNode* current = root;
        
        for (char c : word) {
            int index = c - 'a';
            if (!current->children[index]) {
                current->children[index] = new TrieNode();
            }
            current = current->children[index];
        }
        
        current->isEndOfWord = true;
    }
    
    bool contains(const std::string& word) const {
        TrieNode* current = root;
        
        for (char c : word) {
            int index = c - 'a';
            if (!current->children[index]) {
                return false;
            }
            current = current->children[index];
        }
        
        return current->isEndOfWord;
    }
    
    // Check if there's any word with the given prefix
    bool hasPrefix(const std::string& prefix) const {
        TrieNode* current = root;
        
        for (char c : prefix) {
            int index = c - 'a';
            if (!current->children[index]) {
                return false;
            }
            current = current->children[index];
        }
        
        return true;
    }
};

// --- Direction Arrays (Standard 8 Directions) ---
// N, NE, E, SE, S, SW, W, NW
const int dr[] = {-1, -1, 0, 1, 1, 1, 0, -1};
const int dc[] = {0, 1, 1, 1, 0, -1, -1, -1};

// --- Global Variables ---
int grid_width;
int grid_height;
std::vector<std::string> required_words;
std::unordered_set<std::string> forbidden_words; // Using unordered_set for faster lookups
std::vector<std::vector<std::string>> solutions; // Store found grids
bool find_one_solution = false;
size_t max_word_len = 0; // Maximum length of any word (required or forbidden)
Trie required_trie;   // Trie for required words
Trie forbidden_trie;  // Trie for forbidden words

// --- Function Declarations ---
bool readInputFile(const std::string& filename);
void solve(std::vector<std::string>& current_grid, int row, int col);
bool isGridValid(const std::vector<std::string>& grid);
bool checkCompletedWords(const std::vector<std::string>& grid);
bool isForbiddenWordFormed(const std::vector<std::string>& grid, int r, int c);
bool isRequiredWordMissing(const std::vector<std::string>& grid);
void printGrid(const std::vector<std::string>& grid, std::ostream& os);

// --- Main Function ---
int main(int argc, char* argv[]) {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_file> <one_solution|all_solutions>" << std::endl;
        return 1;
    }

    std::string input_filename = argv[1];
    std::string output_filename = argv[2];
    std::string mode = argv[3];

    if (mode == "one_solution") {
        find_one_solution = true;
    } else if (mode != "all_solutions") {
        std::cerr << "Error: Mode must be 'one_solution' or 'all_solutions'." << std::endl;
        return 1;
    }

    if (!readInputFile(input_filename)) {
        return 1;
    }

    // Initialize grid with '.' placeholder
    std::vector<std::string> grid(grid_height, std::string(grid_width, '.'));
    
    // Start the recursive backtracking
    solve(grid, 0, 0);

    std::ofstream output_file(output_filename);
    if (!output_file) {
        std::cerr << "Error: Could not open output file '" << output_filename << "'." << std::endl;
        return 1;
    }

    if (solutions.empty()) {
        output_file << "No solutions found" << std::endl;
    } else {
        if (!find_one_solution) {
            output_file << solutions.size() << " solution(s)" << std::endl;
            for (const auto& sol : solutions) {
                printGrid(sol, output_file);
            }
        } else {
            // Print just the first solution for one_solution mode
            printGrid(solutions[0], output_file);
        }
    }

    output_file.close();
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    // Print execution time to console
    std::cerr << "Execution time: " << duration.count() << " ms" << std::endl;
    std::cerr << "Found " << solutions.size() << " solution(s)" << std::endl;
    
    return 0;
}

// --- Function Implementations ---

bool readInputFile(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error: Could not open input file '" << filename << "'." << std::endl;
        return false;
    }

    if (!(infile >> grid_width >> grid_height)) {
        std::cerr << "Error: Could not read grid dimensions from file '" << filename << "'." << std::endl;
        infile.close();
        return false;
    }
    
    if (grid_width <= 0 || grid_height <= 0) {
        std::cerr << "Error: Invalid grid dimensions (" << grid_width << "x" << grid_height << ")." << std::endl;
        infile.close();
        return false;
    }

    char type;
    std::string word;
    max_word_len = 0;

    while (infile >> type >> word) {
        if (word.empty()) continue;

        // Update max word length
        max_word_len = std::max(max_word_len, word.length());

        if (type == '+') {
            required_words.push_back(word);
            required_trie.insert(word);
        } else if (type == '-') {
            forbidden_words.insert(word);
            forbidden_trie.insert(word);
        } else {
            std::cerr << "Warning: Invalid type character '" << type << "' in file '" << filename << "'. Skipping line." << std::endl;
        }
    }

    if (infile.bad()) {
        std::cerr << "Error reading from input file '" << filename << "'." << std::endl;
        infile.close();
        return false;
    }

    infile.close();
    return true;
}

void solve(std::vector<std::string>& current_grid, int row, int col) {
    // Early exit if we only need one solution and found it
    if (find_one_solution && !solutions.empty()) {
        return;
    }

    // Base Case: Grid is completely filled
    if (row == grid_height) {
        // Final validation
        if (checkCompletedWords(current_grid)) {
            solutions.push_back(current_grid);
        }
        return;
    }

    // Calculate next cell coordinates
    int next_row = row;
    int next_col = col + 1;
    if (next_col == grid_width) {
        next_row = row + 1;
        next_col = 0;
    }

    // Recursive Step: Try filling current_grid[row][col] with each letter
    for (char ch = 'a'; ch <= 'z'; ++ch) {
        current_grid[row][col] = ch;

        // Pruning: Check if placing this character causes any forbidden word to form
        if (!isForbiddenWordFormed(current_grid, row, col)) {
            // If placing 'ch' is currently valid, recurse to the next cell
            solve(current_grid, next_row, next_col);

            // If we found a solution in the recursive call and only need one, return
            if (find_one_solution && !solutions.empty()) {
                return;
            }
        }
    }

    // Backtrack: Reset the current cell before returning up the call stack
    current_grid[row][col] = '.';
}

// Check if any forbidden word is formed by placing a character at (r, c)
bool isForbiddenWordFormed(const std::vector<std::string>& grid, int r, int c) {
    if (forbidden_words.empty()) {
        return false; // No forbidden words to check
    }

    // Check all 8 directions radiating from the current cell (r, c)
    for (int dir = 0; dir < 8; ++dir) {
        // Check up to max_word_len characters in this direction
        for (size_t len = 1; len <= max_word_len; ++len) {
            std::string word;
            bool valid_word = true;
            
            // Extract word in this direction starting at position (r, c)
            for (size_t i = 0; i < len; ++i) {
                int curr_r = r - (len - 1 - i) * dr[dir];
                int curr_c = c - (len - 1 - i) * dc[dir];
                
                // Check bounds
                if (curr_r < 0 || curr_r >= grid_height || curr_c < 0 || curr_c >= grid_width) {
                    valid_word = false;
                    break;
                }
                
                // If we encounter an empty cell, word is not complete
                if (grid[curr_r][curr_c] == '.') {
                    valid_word = false;
                    break;
                }
                
                word += grid[curr_r][curr_c];
            }
            
            // If we have a valid word of correct length and it's in our forbidden list, return true
            if (valid_word && forbidden_words.count(word)) {
                return true;
            }
        }
    }
    
    return false;
}

// Check if all required words are in the grid
bool isRequiredWordMissing(const std::vector<std::string>& grid) {
    for (const std::string& word : required_words) {
        bool found = false;
        
        // Try each possible starting cell
        for (int r = 0; r < grid_height && !found; ++r) {
            for (int c = 0; c < grid_width && !found; ++c) {
                // Check all 8 directions
                for (int dir = 0; dir < 8 && !found; ++dir) {
                    int end_r = r + (word.length() - 1) * dr[dir];
                    int end_c = c + (word.length() - 1) * dc[dir];
                    
                    // Check if word would fit in bounds
                    if (end_r < 0 || end_r >= grid_height || end_c < 0 || end_c >= grid_width) {
                        continue;
                    }
                    
                    // Check the word
                    bool match = true;
                    for (size_t i = 0; i < word.length(); ++i) {
                        int curr_r = r + i * dr[dir];
                        int curr_c = c + i * dc[dir];
                        
                        if (grid[curr_r][curr_c] != word[i]) {
                            match = false;
                            break;
                        }
                    }
                    
                    if (match) {
                        found = true;
                    }
                }
            }
        }
        
        if (!found) {
            return true; // Missing required word
        }
    }
    
    return false; // All required words found
}

// Final validation of a completed grid
bool checkCompletedWords(const std::vector<std::string>& grid) {
    // 1. Check all required words are present
    if (isRequiredWordMissing(grid)) {
        return false;
    }
    
    // 2. Double check for forbidden words (redundant with pruning but good sanity check)
    if (!forbidden_words.empty()) {
        for (const std::string& word : forbidden_words) {
            for (int r = 0; r < grid_height; ++r) {
                for (int c = 0; c < grid_width; ++c) {
                    for (int dir = 0; dir < 8; ++dir) {
                        int end_r = r + (word.length() - 1) * dr[dir];
                        int end_c = c + (word.length() - 1) * dc[dir];
                        
                        // Check if word would fit in bounds
                        if (end_r < 0 || end_r >= grid_height || end_c < 0 || end_c >= grid_width) {
                            continue;
                        }
                        
                        // Check the word
                        bool match = true;
                        for (size_t i = 0; i < word.length(); ++i) {
                            int curr_r = r + i * dr[dir];
                            int curr_c = c + i * dc[dir];
                            
                            if (grid[curr_r][curr_c] != word[i]) {
                                match = false;
                                break;
                            }
                        }
                        
                        if (match) {
                            return false; // Found forbidden word
                        }
                    }
                }
            }
        }
    }
    
    return true; // All checks passed
}

// Helper function to print the grid in the required format
void printGrid(const std::vector<std::string>& grid, std::ostream& os) {
    os << "Board:\n";
    for (const auto& row : grid) {
        os << "  " << row << "\n";
    }
}
