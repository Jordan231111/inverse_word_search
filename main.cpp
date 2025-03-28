#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <set>
#include <stdexcept> // For exception handling (though not actively used here)
#include <algorithm> // Required for std::max
#include <vector>    // Make sure vector is included if not implicitly

// --- Global Variables ---
int grid_width;
int grid_height;
std::vector<std::string> required_words;
std::set<std::string> forbidden_words; // Using set for faster lookups
std::vector<std::vector<std::string>> solutions; // Store found grids
bool find_one_solution = false;
size_t max_forbidden_len = 0; // Optimization: Max length of a forbidden word
size_t min_forbidden_len = 1000; // Optimization: Min length (init high)

// --- Direction Arrays (Standard 8 Directions) ---
// N, NE, E, SE, S, SW, W, NW
const int dr[] = {-1, -1, 0, 1, 1, 1, 0, -1};
const int dc[] = {0, 1, 1, 1, 0, -1, -1, -1};

// --- Function Declarations ---
bool readInputFile(const std::string& filename);
void solve(std::vector<std::string>& current_grid, int row, int col);
bool isGridValid(const std::vector<std::string>& grid, int height, int width);
bool findWordInGrid(const std::vector<std::string>& grid, const std::string& word, int height, int width);
// *** NEW PRUNING FUNCTION ***
bool checkForbiddenPlacement(const std::vector<std::string>& grid, int r, int c, int height, int width);
void printGrid(const std::vector<std::string>& grid, std::ostream& os);

// --- Main Function ---
int main(int argc, char* argv[]) {
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
                printGrid(sol, output_file); // Call helper to print board
            }
        } else {
            // Print just the first solution for one_solution mode
             printGrid(solutions[0], output_file); // Call helper
        }
    }

    output_file.close();
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
    min_forbidden_len = 1000; // Re-initialize min length
    bool has_forbidden = false;

    while (infile >> type >> word) {
        if (word.empty()) continue;

        if (type == '+') {
            required_words.push_back(word);
        } else if (type == '-') {
            forbidden_words.insert(word);
            max_forbidden_len = std::max(max_forbidden_len, word.length());
            min_forbidden_len = std::min(min_forbidden_len, word.length());
            has_forbidden = true;
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

// *** THE NEW, MORE CORRECT PRUNING FUNCTION ***
// Checks if placing the character at grid[r][c] (which is already done)
// has COMPLETED any forbidden word passing through or ending/starting at (r,c).
bool checkForbiddenPlacement(const std::vector<std::string>& grid, int r, int c, int height, int width) {
    if (forbidden_words.empty() || min_forbidden_len == 0) { // Quick exit if no forbidden words
        return false;
    }

    std::string current_segment; // Buffer to build word segments

    // Check all 8 directions radiating from the current cell (r, c)
    for (int dir = 0; dir < 8; ++dir) {
        // For each direction, check potential forbidden words centered at (r, c)
        // Iterate through possible lengths (from min to max forbidden length)
        for (size_t len = min_forbidden_len; len <= max_forbidden_len; ++len) {

             // Iterate through all possible starting positions OF THE FORBIDDEN WORD
             // such that the character at grid[r][c] corresponds to a character within that forbidden word.
             for(size_t k = 0; k < len; ++k) { // k is the index within the potential forbidden word that corresponds to (r,c)
                int start_r = r - k * dr[dir];
                int start_c = c - k * dc[dir];

                // Optimization: If the calculated start is way off, continue
                // We only need to check if the start allows the END of the word to be within bounds.
                int end_r = start_r + (len - 1) * dr[dir];
                int end_c = start_c + (len - 1) * dc[dir];
                if (start_r < 0 || start_r >= height || start_c < 0 || start_c >= width ||
                    end_r < 0 || end_r >= height || end_c < 0 || end_c >= width) {
                    continue; // This alignment isn't possible within bounds
                }

                // Build the segment of length 'len' starting from (start_r, start_c)
                current_segment.clear();
                current_segment.reserve(len);
                bool possible = true;
                for (size_t i = 0; i < len; ++i) {
                    int current_r = start_r + i * dr[dir];
                    int current_c = start_c + i * dc[dir];

                     // We already checked start and end bounds, but intermediate could technically
                     // wrap around in a weird grid shape (not applicable here, but safe)
                     // Bounds check primarily ensures we access valid grid indices.
                     // if (current_r < 0 || current_r >= height || current_c < 0 || current_c >= width) {
                     //    possible = false; break; // Should not happen due to start/end check
                     // }

                    char grid_char = grid[current_r][current_c];
                    if (grid_char == '.') { // If any part of the segment is still empty '.', it's not a completed word
                        possible = false;
                        break;
                    }
                    current_segment += grid_char;
                }

                // If we successfully built a full segment of the correct length
                if (possible && current_segment.length() == len) {
                    // Check if this completed segment is a forbidden word
                    if (forbidden_words.count(current_segment)) {
                        //std::cout << "Pruning: Found forbidden '" << current_segment << "' by placing char at (" << r << "," << c << ")\n";
                        return true; // Yes, this placement creates a forbidden word
                    }
                }
             } // End loop k (index within forbidden word)
        } // End loop len (length of forbidden word)
    } // End loop dir (direction)

    return false; // No forbidden words were completed by this placement
}


void solve(std::vector<std::string>& current_grid, int row, int col) {
    // Early exit if we only need one solution and found it
    if (find_one_solution && !solutions.empty()) {
        return;
    }

    // Base Case: Grid is completely filled
    if (row == grid_height) {
        // Final validation (redundant if pruning is perfect, but good failsafe)
        if (isGridValid(current_grid, grid_height, grid_width)) {
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

        // *** Pruning Step: Use the new, more robust check ***
        if (!checkForbiddenPlacement(current_grid, row, col, grid_height, grid_width)) {
            // If placing 'ch' is currently valid, recurse to the next cell
            solve(current_grid, next_row, next_col);

            // If we found a solution in the recursive call and only need one, stop backtracking here
            if (find_one_solution && !solutions.empty()) {
                 current_grid[row][col] = '.'; // Backtrack before returning
                 return;
             }
        }
        // else: Pruning failed, the loop will try the next character after backtracking.
    }

    // Backtrack: Reset the current cell before returning up the call stack
    current_grid[row][col] = '.';
}

// Checks if a *completed* grid contains all required and no forbidden words
bool isGridValid(const std::vector<std::string>& grid, int height, int width) {
    // 1. Check for all required words
    for (const std::string& req_word : required_words) {
        if (!findWordInGrid(grid, req_word, height, width)) {
            return false; // Missing a required word
        }
    }

    // 2. Double-check for forbidden words (mostly redundant if pruning works)
    if (!forbidden_words.empty()) { // Only check if there are forbidden words
        for (const std::string& f_word : forbidden_words) {
             // Optimization: skip check if word is longer than possible in grid
             if(f_word.length() > std::max(width, height) && f_word.length() > 1) { // diagonal check needed? len > sqrt(w^2+h^2)
                 int max_dim = std::max(width, height);
                 if (f_word.length() > max_dim) continue; // Simplified check
             }

             if (findWordInGrid(grid, f_word, height, width)) {
                return false; // Found a forbidden word
            }
        }
    }

    return true; // All checks passed
}

// Finds a specific word in a grid (8 directions)
bool findWordInGrid(const std::vector<std::string>& grid, const std::string& word, int height, int width) {
    if (word.empty()) return true; // Define behavior for empty word search
    if (height == 0 || width == 0) return false;

    // Check if word is even possible given grid dimensions
     size_t max_possible_len = 0;
     if (height > 0 && width > 0) {
        max_possible_len = std::max({(size_t)width, (size_t)height, (size_t)std::min(width,height)}); // Approx max diagonal
     }
     if (word.length() > max_possible_len && word.length() > 1) { // Be careful with length 1 words
         // More accurate check: word length cannot exceed grid dimensions
         if (word.length() > (size_t)width && word.length() > (size_t)height) return false;
     }


    for (int r = 0; r < height; ++r) {
        for (int c = 0; c < width; ++c) {
            // Optimization: If first char doesn't match, skip direction checks
            if (grid[r][c] != word[0]) continue;

            for (int dir = 0; dir < 8; ++dir) {
                size_t k;
                for (k = 0; k < word.length(); ++k) {
                    int current_r = r + k * dr[dir];
                    int current_c = c + k * dc[dir];

                    // Check bounds
                    if (current_r < 0 || current_r >= height || current_c < 0 || current_c >= width) {
                        break; // Out of bounds
                    }
                    // Check character match (ensure grid cell isn't '.')
                    if (grid[current_r][current_c] == '.' || grid[current_r][current_c] != word[k]) {
                        break; // Mismatch or empty cell
                    }
                }
                // If we matched the whole word (k reached word.length())
                if (k == word.length()) {
                    return true;
                }
            }
        }
    }

    return false; // Word not found
}

// Helper function to print the grid in the required format
void printGrid(const std::vector<std::string>& grid, std::ostream& os) {
    os << "Board:\n"; // Add the "Board:" label
    for (const auto& row : grid) {
        os << "  " << row << "\n"; // Add indentation
    }
}