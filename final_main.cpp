#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_set>
#include <algorithm>
#include <chrono>
#include <array>
#include <bitset>
#include <memory> // For std::unique_ptr
#include <unordered_map> // For Trie children

// --- Direction Arrays (Standard 8 Directions) ---
// N, NE, E, SE, S, SW, W, NW
const int dr[] = {-1, -1, 0, 1, 1, 1, 0, -1};
const int dc[] = {0, 1, 1, 1, 0, -1, -1, -1};

// --- Trie Implementation for Forbidden Words ---
struct TrieNode {
    std::unordered_map<char, std::unique_ptr<TrieNode>> children;
    bool isEndOfWord = false;

    TrieNode() = default; // Use default constructor
};

class Trie {
public:
    Trie() : root(std::make_unique<TrieNode>()) {}

    void insert(const std::string& word) {
        TrieNode* current = root.get();
        for (char ch : word) {
            // If the character node doesn't exist, create it
            if (current->children.find(ch) == current->children.end()) {
                current->children[ch] = std::make_unique<TrieNode>();
            }
            current = current->children[ch].get();
        }
        current->isEndOfWord = true;
    }

    // Checks if any substring of the given segment is a forbidden word in the Trie
    bool containsForbiddenWord(const std::string& segment) const {
        if (!root || segment.empty()) {
            return false;
        }

        for (size_t i = 0; i < segment.length(); ++i) {
            TrieNode* current = root.get();
            for (size_t j = i; j < segment.length(); ++j) {
                char ch = segment[j];
                if (current->children.find(ch) == current->children.end()) {
                    // No path for this character, break inner loop
                    break;
                }
                current = current->children[ch].get();
                if (current->isEndOfWord) {
                    // Found a forbidden word ending at index j, starting at or before i
                    return true;
                }
            }
        }
        return false; // No forbidden word found as a substring
    }

    // TODO: Add a method to check if a segment contains a forbidden word

private:
    std::unique_ptr<TrieNode> root;
};

// --- Global Variables ---
int grid_width;
int grid_height;
std::vector<std::string> required_words;
Trie forbidden_trie; // Global Trie for forbidden words
std::vector<std::vector<std::string>> solutions;
bool find_one_solution = false;
size_t max_forbidden_len = 0;
size_t min_forbidden_len = 1000;
size_t max_required_len = 0;

// --- Function Declarations ---
bool readInputFile(const std::string& filename);
void solve(std::vector<std::string>& current_grid, int row, int col);
bool isGridValid(const std::vector<std::string>& grid);
bool findWordInGrid(const std::vector<std::string>& grid, const std::string& word);
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
    max_forbidden_len = 0;
    min_forbidden_len = 1000;
    max_required_len = 0;
    bool has_forbidden = false;

    while (infile >> type >> word) {
        if (word.empty()) continue;

        if (type == '+') {
            required_words.push_back(word);
            max_required_len = std::max(max_required_len, word.length());
        } else if (type == '-') {
            if (!word.empty()) { // Avoid inserting empty strings
                forbidden_trie.insert(word);
                // Also insert the reverse for backward checks
                std::string reversed_word = word;
                std::reverse(reversed_word.begin(), reversed_word.end());
                forbidden_trie.insert(reversed_word);

                max_forbidden_len = std::max(max_forbidden_len, word.length());
                min_forbidden_len = std::min(min_forbidden_len, word.length());
                has_forbidden = true;
            }
        } else {
            std::cerr << "Warning: Invalid type character '" << type << "' in file '" << filename << "'. Skipping line." << std::endl;
        }
    }
    
    if (!has_forbidden) {
        min_forbidden_len = 0; // If no forbidden words, set min length to 0
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
        // Final validation (redundant if pruning is perfect, but good failsafe)
        if (isGridValid(current_grid)) {
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

        // --- Pruning Step: Check if placing 'ch' creates a forbidden word --- 
        bool creates_forbidden = false;
        if (min_forbidden_len > 0) { // Only check if there are forbidden words
            std::string segment; // Moved segment declaration outside dir loop
            for (int dir = 0; dir < 8; ++dir) {
                segment.clear(); // Clear segment for each direction
                // Build segment backwards from (row, col) up to max_forbidden_len
                for (size_t k = 0; k < max_forbidden_len; ++k) {
                    int current_r = row - k * dr[dir];
                    int current_c = col - k * dc[dir];

                    // Check bounds
                    if (current_r < 0 || current_r >= grid_height || current_c < 0 || current_c >= grid_width) {
                        break; // Went off grid
                    }

                    char grid_char = current_grid[current_r][current_c];
                    if (grid_char == '.') {
                        // If we hit an empty cell before reaching min_forbidden_len, 
                        // this segment cannot form a forbidden word yet in this direction.
                        // However, a shorter forbidden word might have already been formed.
                        // Let the containsForbiddenWord handle prefixes correctly.
                        break; 
                    }
                    segment += grid_char; // Append character (builds backward)
                }
                
                if (segment.empty()) continue; // Skip if segment is empty

                std::reverse(segment.begin(), segment.end()); // Reverse to get correct order
                
                // Check segment and its relevant substrings using the Trie
                if (segment.length() >= min_forbidden_len && forbidden_trie.containsForbiddenWord(segment)) {
                    creates_forbidden = true;
                    break; // Found forbidden word in this direction, no need to check others
                }
            }
        }
        // --- End Pruning Step ---

        // If placing 'ch' is currently valid (did not create forbidden word), recurse
        if (!creates_forbidden) {
            solve(current_grid, next_row, next_col);

            // If we found a solution and only need one, stop backtracking
            if (find_one_solution && !solutions.empty()) {
                return;
            }
        }
    }

    // Backtrack: Reset the current cell before returning
    current_grid[row][col] = '.';
}

// Check if a specific word exists in the grid
bool findWordInGrid(const std::vector<std::string>& grid, const std::string& word) {
    if (word.empty()) return true;
    if (word.length() > std::max(grid_width, grid_height)) return false; // Quick size check
    
    for (int r = 0; r < grid_height; ++r) {
        for (int c = 0; c < grid_width; ++c) {
            // Quick first character check
            if (grid[r][c] != word[0]) continue;
            
            // Try all 8 directions
            for (int dir = 0; dir < 8; ++dir) {
                size_t k;
                for (k = 0; k < word.length(); ++k) {
                    int current_r = r + k * dr[dir];
                    int current_c = c + k * dc[dir];
                    
                    // Check bounds
                    if (current_r < 0 || current_r >= grid_height || current_c < 0 || current_c >= grid_width) {
                        break;
                    }
                    
                    // Check character match
                    if (grid[current_r][current_c] != word[k]) {
                        break;
                    }
                }
                
                // If we matched the whole word
                if (k == word.length()) {
                    return true;
                }
            }
        }
    }
    
    return false; // Word not found
}

// Validate the completed grid
bool isGridValid(const std::vector<std::string>& grid) {
    // Check all required words
    for (const std::string& req_word : required_words) {
        if (!findWordInGrid(grid, req_word)) {
            return false; // Missing a required word
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
