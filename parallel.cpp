#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <atomic>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <functional>
#include <future>
#include <memory>
#include <deque>
#include <type_traits>
#include <chrono>

// Directions for word search
const int dx[8] = {1, 1, 0, -1, -1, -1, 0, 1};
const int dy[8] = {0, 1, 1, 1, 0, -1, -1, -1};

// Structure to represent a word placement
struct WordPlacement {
    std::string word;
    int row, col, dir;
    size_t length;
};

// Print grid to output stream
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

// C++11 compatible make_unique implementation
template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

// Simple thread-safe work queue
class WorkQueue {
private:
    std::queue<std::function<void()>> tasks;
    std::mutex mutex;
    std::condition_variable condition;
    std::atomic<bool> stop;

public:
    WorkQueue() : stop(false) {}

    template<class F>
    void push(F&& f) {
        {
            std::unique_lock<std::mutex> lock(mutex);
            tasks.push(std::forward<F>(f));
        }
        condition.notify_one();
    }

    void waitForTask(std::function<void()>& task) {
        std::unique_lock<std::mutex> lock(mutex);
        condition.wait(lock, [this] { return stop || !tasks.empty(); });
        
        if (stop && tasks.empty()) {
            return;
        }
        
        if (!tasks.empty()) {
            task = std::move(tasks.front());
            tasks.pop();
        }
    }

    bool try_pop(std::function<void()>& task) {
        std::lock_guard<std::mutex> lock(mutex);
        if (tasks.empty()) return false;
        task = std::move(tasks.front());
        tasks.pop();
        return true;
    }

    void shutdown() {
        {
            std::unique_lock<std::mutex> lock(mutex);
            stop = true;
        }
        condition.notify_all();
    }

    bool empty() {
        std::lock_guard<std::mutex> lock(mutex);
        return tasks.empty();
    }
};

// Simple thread pool with work stealing
class ThreadPool {
private:
    std::vector<std::thread> workers;
    WorkQueue queue;
    std::atomic<bool> stop;
    
public:
    ThreadPool(size_t threads = 0) : stop(false) {
        if (threads == 0) {
            threads = std::thread::hardware_concurrency();
            if (threads == 0) threads = 4;
        }
        
        for (size_t i = 0; i < threads; ++i) {
            workers.emplace_back([this] {
                std::function<void()> task;
                
                while (!stop) {
                    queue.waitForTask(task);
                    
                    if (task) {
                        task();
                        task = nullptr;
                    } else if (stop) {
                        break;
                    }
                }
            });
        }
    }
    
    ~ThreadPool() {
        stop = true;
        queue.shutdown();
        
        for (auto& worker : workers) {
            if (worker.joinable()) {
                worker.join();
            }
        }
    }
    
    template<class F, class... Args>
    std::future<typename std::result_of<F(Args...)>::type> 
    enqueue(F&& f, Args&&... args) {
        using return_type = typename std::result_of<F(Args...)>::type;
        
        auto task = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );
        
        std::future<return_type> result = task->get_future();
        
        queue.push([task] { (*task)(); });
        
        return result;
    }
    
    size_t num_threads() const {
        return workers.size();
    }
};

// Optimized inverse word search solver
class OptimizedInverseWordSearch {
private:
    int width, height;
    std::vector<std::string> requiredWords;
    std::vector<std::string> forbiddenWords;
    bool findAllSolutions;
    std::atomic<bool> solutionFound;
    std::mutex solutionsMutex;
    std::vector<std::vector<std::vector<char>>> solutions;
    std::set<std::string> solutionHashes;
    int dirOffsets[8]; // Pre-computed direction offsets

    // Helper function to convert 2D coordinates to 1D index
    inline int flatIndex(int row, int col) const {
        return row * width + col;
    }
    
public:
    OptimizedInverseWordSearch(int w, int h,
                      const std::vector<std::string>& required,
                      const std::vector<std::string>& forbidden,
                      bool findAll)
        : width(w), height(h), requiredWords(required), forbiddenWords(forbidden),
          findAllSolutions(findAll), solutionFound(false) {
        
        // Pre-compute direction offsets for faster grid traversal
        for (int dir = 0; dir < 8; dir++) {
            dirOffsets[dir] = dy[dir] * width + dx[dir];
        }
    }

    // Optimize word order for more efficient search
    void optimizeWordOrder() {
        const size_t wordCount = requiredWords.size();
        if (wordCount <= 1) return;
        
        // Score words based on constraint potential and rarity
        std::vector<std::pair<int, std::string>> scoredWords;
        scoredWords.reserve(wordCount);
        
        for (const auto& word : requiredWords) {
            int score = 0;
            
            // Longer words are placed first
            score += word.length() * 5;
            
            // Count unique letters (higher uniqueness = higher score)
            bool seen[26] = {false};
            int uniqueCount = 0;
            
            for (char c : word) {
                int idx = c - 'a';
                if (idx >= 0 && idx < 26 && !seen[idx]) {
                    seen[idx] = true;
                    uniqueCount++;
                    
                    // Rarer letters get higher scores
                    switch (c) {
                        case 'q': case 'x': case 'z': case 'j': score += 10; break;
                        case 'k': case 'v': case 'w': case 'y': score += 8; break;
                        case 'b': case 'f': case 'g': case 'h': case 'p': score += 5; break;
                        default: score += 1; break;
                    }
                }
            }
            
            score += uniqueCount * 3;
            scoredWords.push_back(std::make_pair(score, word));
        }
        
        // Sort words by score (descending)
        std::sort(scoredWords.begin(), scoredWords.end(), 
                 [](const std::pair<int, std::string>& a, const std::pair<int, std::string>& b) {
                     // If scores differ, sort by score
                     if (a.first != b.first) {
                         return a.first > b.first;
                     }
                     // If scores are equal, sort by length
                     return a.second.length() > b.second.length();
                 });
        
        // Update requiredWords with sorted order
        requiredWords.clear();
        for (const auto& pair : scoredWords) {
            requiredWords.push_back(pair.second);
        }
    }

    bool solve() {
        // Clear any existing solutions
        solutions.clear();
        solutionHashes.clear();
        solutionFound.store(false);
        
        // Create a thread pool with hardware concurrency
        ThreadPool pool;
        int numThreads = pool.num_threads();
        
        // First calculate all possible placements for the first word
        struct FirstPlacement {
            int row, col, dir;
        };
        
        std::vector<FirstPlacement> firstPlacements;
        const std::string& firstWord = requiredWords[0];
        const size_t firstWordLen = firstWord.length();
        
        for (int dir = 0; dir < 8; dir++) {
            // Calculate valid bounds for this direction
            const int dirMaxOffset = static_cast<int>(firstWordLen) - 1;
            int minRow = (dy[dir] < 0) ? dirMaxOffset : 0;
            int minCol = (dx[dir] < 0) ? dirMaxOffset : 0;
            int maxRow = height - ((dy[dir] > 0) ? dirMaxOffset : 0);
            int maxCol = width - ((dx[dir] > 0) ? dirMaxOffset : 0);
            
            for (int row = minRow; row < maxRow; row++) {
                for (int col = minCol; col < maxCol; col++) {
                    firstPlacements.push_back({row, col, dir});
                }
            }
        }
        
        // If too few placements, just solve sequentially
        if (firstPlacements.size() < 20) {
            return solveSequential();
        }

        // Divide placements into batches for each thread
        int placementsPerThread = (firstPlacements.size() + numThreads - 1) / numThreads;
        std::vector<std::future<bool>> futures;
        
        for (size_t i = 0; i < firstPlacements.size(); i += placementsPerThread) {
            size_t end = std::min(i + placementsPerThread, firstPlacements.size());
            
            futures.push_back(pool.enqueue(
                [this, &firstPlacements, i, end]() {
                    for (size_t j = i; j < end; j++) {
                        // Skip if solution already found
                        if (!findAllSolutions && solutionFound.load(std::memory_order_acquire)) {
                            return true;
                        }
                        
                        const auto& placement = firstPlacements[j];
                        
                        // Initialize grid
                        std::vector<char> grid(width * height, '.');
                        std::vector<WordPlacement> placements;
                        
                        // Place first word
                        const std::string& word = requiredWords[0];
                        const size_t wordLen = word.length();
                        
                        int idx = flatIndex(placement.row, placement.col);
                        const int dirOffset = dirOffsets[placement.dir];
                        
                        for (size_t k = 0; k < wordLen; k++) {
                            grid[idx] = word[k];
                            idx += dirOffset;
                        }
                        
                        // Add to placements
                        placements.push_back({word, placement.row, placement.col, placement.dir, wordLen});
                        
                        // Try to place remaining words
                        if (placeWords(1, placements, grid)) {
                            return true;
                        }
                    }
                    
                    return false;
                }
            ));
        }
        
        // Wait for results
        bool success = false;
        for (auto& future : futures) {
            success = future.get() || success;
            
            if (!findAllSolutions && success) {
                break;
            }
        }
        
        return success;
    }
    
    bool solveSequential() {
        std::vector<char> grid(width * height, '.');
        std::vector<WordPlacement> placements;
        
        const std::string& firstWord = requiredWords[0];
        const size_t firstWordLen = firstWord.length();
        bool found = false;
        
        for (int dir = 0; dir < 8 && !found; dir++) {
            const int dirMaxOffset = static_cast<int>(firstWordLen) - 1;
            int minRow = (dy[dir] < 0) ? dirMaxOffset : 0;
            int minCol = (dx[dir] < 0) ? dirMaxOffset : 0;
            int maxRow = height - ((dy[dir] > 0) ? dirMaxOffset : 0);
            int maxCol = width - ((dx[dir] > 0) ? dirMaxOffset : 0);
            
            for (int row = minRow; row < maxRow && !found; row++) {
                for (int col = minCol; col < maxCol && !found; col++) {
                    // Reset grid
                    std::fill(grid.begin(), grid.end(), '.');
                    placements.clear();
                    
                    // Place first word
                    int idx = flatIndex(row, col);
                    const int dirOffset = dirOffsets[dir];
                    
                    for (size_t i = 0; i < firstWordLen; i++) {
                        grid[idx] = firstWord[i];
                        idx += dirOffset;
                    }
                    
                    placements.push_back({firstWord, row, col, dir, firstWordLen});
                    
                    if (placeWords(1, placements, grid)) {
                        found = true;
                        if (!findAllSolutions) {
                            break;
                        }
                    }
                }
            }
        }
        
        return found;
    }

    const std::vector<std::vector<std::vector<char>>>& getSolutions() const {
        return solutions;
    }

private:
    // Check if word exists in the grid
    bool wordExists(const std::string& word, const std::vector<char>& grid) const {
        const size_t wordLen = word.length();
        
        // Optimize for single letter words
        if (wordLen == 1) {
            char c = word[0];
            for (int i = 0; i < width * height; i++) {
                if (grid[i] == c) {
                    return true;
                }
            }
            return false;
        }

        char firstChar = word[0];
        
        // Find occurrences of first letter
        for (int row = 0; row < height; row++) {
            for (int col = 0; col < width; col++) {
                if (grid[flatIndex(row, col)] == firstChar) {
                    // Check all directions
                    for (int dir = 0; dir < 8; dir++) {
                        int endRow = row + dy[dir] * (wordLen - 1);
                        int endCol = col + dx[dir] * (wordLen - 1);
                        
                        // Skip out of bounds
                        if (endRow < 0 || endRow >= height || endCol < 0 || endCol >= width) {
                            continue;
                        }
                        
                        // Check if word matches
                        bool match = true;
                        for (size_t i = 1; i < wordLen; i++) {
                            int r = row + dy[dir] * i;
                            int c = col + dx[dir] * i;
                            if (grid[flatIndex(r, c)] != word[i]) {
                                match = false;
                                break;
                            }
                        }
                        
                        if (match) {
                            return true;
                        }
                    }
                }
            }
        }
        
        return false;
    }

    // Try to place a word at a specific position
    bool canPlaceWord(const std::string& word, int row, int col, int dir, const std::vector<char>& grid) const {
        const size_t wordLen = word.length();
        
        // Check bounds
        int endRow = row + dy[dir] * (wordLen - 1);
        int endCol = col + dx[dir] * (wordLen - 1);
        
        if (endRow < 0 || endRow >= height || endCol < 0 || endCol >= width) {
            return false;
        }
        
        // Check if word can be placed
        int idx = flatIndex(row, col);
        int dirOffset = dirOffsets[dir];
        
        for (size_t i = 0; i < wordLen; i++) {
            if (grid[idx] != '.' && grid[idx] != word[i]) {
                return false;
            }
            idx += dirOffset;
        }
        
        return true;
    }

    // Check if grid is valid
    bool isValidGrid(const std::vector<char>& grid) const {
        // Check required words
        for (const auto& word : requiredWords) {
            if (!wordExists(word, grid)) {
                return false;
            }
        }
        
        // Check forbidden words
        for (const auto& word : forbiddenWords) {
            if (wordExists(word, grid)) {
                return false;
            }
        }
        
        return true;
    }

    // Fill empty cells
    bool fillEmptyCells(int row, int col, std::vector<char>& grid) {
        // Check for early termination
        if (!findAllSolutions && solutionFound.load(std::memory_order_acquire)) {
            return false;
        }
        
        // If we've filled the entire grid, check if it's valid
        if (row >= height) {
            if (isValidGrid(grid)) {
                // Convert to string for hash
                std::string hash(grid.begin(), grid.end());
                
                // Lock for updating solutions
                std::lock_guard<std::mutex> lock(solutionsMutex);
                
                if (solutionHashes.find(hash) == solutionHashes.end()) {
                    solutionHashes.insert(hash);
                    
                    // Convert flat grid to 2D for output
                    std::vector<std::vector<char>> solution(height, std::vector<char>(width));
                    for (int r = 0; r < height; r++) {
                        for (int c = 0; c < width; c++) {
                            solution[r][c] = grid[flatIndex(r, c)];
                        }
                    }
                    
                    solutions.push_back(solution);
                    
                    if (!findAllSolutions) {
                        solutionFound.store(true, std::memory_order_release);
                    }
                    
                    return true;
                }
            }
            
            return false;
        }
        
        // Calculate next position
        int nextRow = row;
        int nextCol = col + 1;
        if (nextCol >= width) {
            nextRow++;
            nextCol = 0;
        }
        
        // If cell is already filled, move to next
        if (grid[flatIndex(row, col)] != '.') {
            return fillEmptyCells(nextRow, nextCol, grid);
        }
        
        // Try each letter
        for (char c = 'a'; c <= 'z'; c++) {
            grid[flatIndex(row, col)] = c;
            
            // Check if this creates any forbidden words
            bool validPlacement = true;
            for (const auto& word : forbiddenWords) {
                if (word.length() == 1 && word[0] == c) {
                    validPlacement = false;
                    break;
                }
                
                // Only check words that contain this letter (optimization)
                if (word.find(c) != std::string::npos) {
                    if (wordExists(word, grid)) {
                        validPlacement = false;
                        break;
                    }
                }
            }
            
            if (validPlacement && fillEmptyCells(nextRow, nextCol, grid)) {
                return true;
            }
        }
        
        // Backtrack
        grid[flatIndex(row, col)] = '.';
        return false;
    }

    // Place words recursively
    bool placeWords(int wordIndex, std::vector<WordPlacement>& placements, std::vector<char>& grid) {
        // Check for early termination
        if (!findAllSolutions && solutionFound.load(std::memory_order_acquire)) {
            return false;
        }
        
        // If all words placed, fill remaining cells
        if (wordIndex >= static_cast<int>(requiredWords.size())) {
            return fillEmptyCells(0, 0, grid);
        }
        
        const std::string& word = requiredWords[wordIndex];
        const size_t wordLen = word.length();
        bool foundSolution = false;
        
        // Track changes for efficient backtracking
        struct Change {
            int pos;
            char prevChar;
        };
        
        std::vector<Change> changes;
        changes.reserve(wordLen);
        
        // Try all possible positions
        for (int dir = 0; dir < 8; dir++) {
            // Boundary calculations
            const int dirMaxOffset = static_cast<int>(wordLen) - 1;
            int minRow = (dy[dir] < 0) ? dirMaxOffset : 0;
            int minCol = (dx[dir] < 0) ? dirMaxOffset : 0;
            int maxRow = height - ((dy[dir] > 0) ? dirMaxOffset : 0);
            int maxCol = width - ((dx[dir] > 0) ? dirMaxOffset : 0);
            
            for (int row = minRow; row < maxRow; row++) {
                for (int col = minCol; col < maxCol; col++) {
                    // Check if word can be placed
                    if (canPlaceWord(word, row, col, dir, grid)) {
                        // Place word and track changes
                        changes.clear();
                        int idx = flatIndex(row, col);
                        const int dirOffset = dirOffsets[dir];
                        
                        for (size_t i = 0; i < wordLen; i++) {
                            if (grid[idx] != word[i]) {
                                changes.push_back({idx, grid[idx]});
                                grid[idx] = word[i];
                            }
                            idx += dirOffset;
                        }
                        
                        // Add to placements
                        placements.push_back({word, row, col, dir, wordLen});
                        
                        // Try to place next word
                        bool result = placeWords(wordIndex + 1, placements, grid);
                        
                        if (result) {
                            foundSolution = true;
                            if (!findAllSolutions) {
                                return true;
                            }
                        }
                        
                        // Backtrack - restore changed cells
                        placements.pop_back();
                        for (const auto& change : changes) {
                            grid[change.pos] = change.prevChar;
                        }
                    }
                }
            }
        }
        
        return foundSolution;
    }
};

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
    requiredWords.reserve(20);
    forbiddenWords.reserve(20);
    
    while (std::getline(inFile, line)) {
        if (line.empty()) {
            continue;
        }
        
        char type = line[0];
        std::string word = line.substr(2);
        
        // Trim trailing whitespace
        size_t end = word.find_last_not_of(" \n\r\t");
        if (end != std::string::npos) {
            word = word.substr(0, end + 1);
        }
        
        if (type == '+') {
            requiredWords.push_back(word);
        } else if (type == '-') {
            forbiddenWords.push_back(word);
        }
    }

    inFile.close();

    // Check if puzzle is solvable
    bool impossiblePuzzle = false;
    for (const auto& word : requiredWords) {
        if (word.length() > std::max(width, height)) {
            impossiblePuzzle = true;
            break;
        }
    }

    // Solve puzzle
    std::vector<std::vector<std::vector<char>>> solutions;

    if (!impossiblePuzzle) {
        bool findAllSols = (mode == "all_solutions");
        
        auto startTime = std::chrono::high_resolution_clock::now();
        
        OptimizedInverseWordSearch solver(width, height, requiredWords, forbiddenWords, findAllSols);
        solver.optimizeWordOrder();
        bool success = solver.solve();
        
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
        
        std::cout << "Solved in " << duration.count() << " ms" << std::endl;
        
        if (success) {
            solutions = solver.getSolutions();
        }
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