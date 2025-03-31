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
#include <string_view> // Modern C++ for non-owning string references

// Directions for word search (8 directions: right, right-down, down, down-left, left, left-up, up, up-right)
const int directionOffsetX[8] = {1, 1, 0, -1, -1, -1, 0, 1};
const int directionOffsetY[8] = {0, 1, 1, 1, 0, -1, -1, -1};

// Structure to represent a word placement
struct WordPlacement {
    std::string word;
    int startRow{0}, startCol{0}, direction{0}; // Default member initialization
    size_t length{0};
};

// Print grid to output stream
void printGrid(std::ostream& outputStream, const std::vector<std::vector<char>>& grid) noexcept {
    outputStream << "Board: " << std::endl;
    for (const auto& row : grid) {
        outputStream << "  ";
        for (const char cellValue : row) {
            outputStream << cellValue;
        }
        outputStream << std::endl;
    }
}

// Simple thread-safe work queue
class WorkQueue {
private:
    std::queue<std::function<void()>> taskQueue;
    std::mutex queueMutex;
    std::condition_variable queueCondition;
    std::atomic<bool> isShuttingDown{false};

public:
    WorkQueue() = default; // Use defaulted constructor instead of empty one

    template<class F>
    void push(F&& task) {
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            taskQueue.push(std::forward<F>(task));
        }
        queueCondition.notify_one();
    }

    void waitForTask(std::function<void()>& task) {
        std::unique_lock<std::mutex> lock(queueMutex);
        queueCondition.wait(lock, [this] { return isShuttingDown || !taskQueue.empty(); });
        
        if (isShuttingDown && taskQueue.empty()) {
            return;
        }
        
        if (!taskQueue.empty()) {
            task = std::move(taskQueue.front());
            taskQueue.pop();
        }
    }

    [[nodiscard]] bool try_pop(std::function<void()>& task) noexcept {
        std::lock_guard<std::mutex> lock(queueMutex);
        if (taskQueue.empty()) return false;
        task = std::move(taskQueue.front());
        taskQueue.pop();
        return true;
    }

    void shutdown() noexcept {
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            isShuttingDown = true;
        }
        queueCondition.notify_all();
    }

    [[nodiscard]] bool empty() noexcept {
        std::lock_guard<std::mutex> lock(queueMutex);
        return taskQueue.empty();
    }
};

// Simple thread pool with work stealing
class ThreadPool {
private:
    std::vector<std::thread> workerThreads;
    WorkQueue workQueue;
    std::atomic<bool> isShuttingDown{false};
    
public:
    explicit ThreadPool(size_t numThreads = 0) {
        if (numThreads == 0) {
            numThreads = std::thread::hardware_concurrency();
            if (numThreads == 0) numThreads = 4;
        }
        
        for (size_t threadIdx = 0; threadIdx < numThreads; ++threadIdx) {
            workerThreads.emplace_back([this] {
                std::function<void()> task;
                
                while (!isShuttingDown) {
                    workQueue.waitForTask(task);
                    
                    if (task) {
                        task();
                        task = nullptr;
                    } else if (isShuttingDown) {
                        break;
                    }
                }
            });
        }
    }
    
    // Rule of five: declare destructor and delete other special members
    ~ThreadPool() {
        isShuttingDown = true;
        workQueue.shutdown();
        
        for (auto& workerThread : workerThreads) {
            if (workerThread.joinable()) {
                workerThread.join();
            }
        }
    }
    
    // Delete copy/move operations for thread safety
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;
    ThreadPool(ThreadPool&&) = delete;
    ThreadPool& operator=(ThreadPool&&) = delete;
    
    template<class F, class... Args>
    [[nodiscard]] std::future<std::invoke_result_t<F, Args...>> // Using invoke_result instead of result_of
    enqueue(F&& func, Args&&... args) {
        using ReturnType = std::invoke_result_t<F, Args...>;
        
        auto taskPackage = std::make_shared<std::packaged_task<ReturnType()>>(
            std::bind(std::forward<F>(func), std::forward<Args>(args)...)
        );
        
        std::future<ReturnType> resultFuture = taskPackage->get_future();
        
        workQueue.push([taskPackage] { (*taskPackage)(); });
        
        return resultFuture;
    }
    
    [[nodiscard]] size_t num_threads() const noexcept {
        return workerThreads.size();
    }
};

// Optimized inverse word search solver
class OptimizedInverseWordSearch {
private:
    int gridWidth, gridHeight;
    std::vector<std::string> requiredWords;
    std::vector<std::string> forbiddenWords;
    bool findAllSolutions;
    std::atomic<bool> solutionFound{false};
    std::mutex solutionsMutex;
    std::vector<std::vector<std::vector<char>>> solutions;
    std::set<std::string> solutionHashes;
    int directionOffsets[8] = {}; // Pre-computed direction offsets with zero initialization

    // Helper function to convert 2D coordinates to 1D index
    [[nodiscard]] inline int flatIndex(int row, int col) const noexcept {
        return row * gridWidth + col;
    }
    
public:
    OptimizedInverseWordSearch(int width, int height,
                      const std::vector<std::string>& required,
                      const std::vector<std::string>& forbidden,
                      bool findAll)
        : gridWidth(width), gridHeight(height), 
          requiredWords(required), forbiddenWords(forbidden),
          findAllSolutions(findAll) {
        
        // Pre-compute direction offsets for faster grid traversal
        for (int direction = 0; direction < 8; direction++) {
            directionOffsets[direction] = directionOffsetY[direction] * gridWidth + directionOffsetX[direction];
        }
    }

    // Optimize word order for more efficient search
    void optimizeWordOrder() {
        const size_t wordCount = requiredWords.size();
        if (wordCount <= 1) return;
        
        // Score words based on constraint potential and rarity
        std::vector<std::pair<int, std::string>> scoredWords;
        scoredWords.reserve(wordCount);
        
        for (const auto& currentWord : requiredWords) {
            int scoreValue = 0;
            
            // Longer words are placed first
            scoreValue += currentWord.length() * 5;
            
            // Count unique letters (higher uniqueness = higher score)
            bool seenLetters[26] = {false};
            int uniqueLetterCount = 0;
            
            for (const char letter : currentWord) {
                int letterIdx = letter - 'a';
                if (letterIdx >= 0 && letterIdx < 26 && !seenLetters[letterIdx]) {
                    seenLetters[letterIdx] = true;
                    uniqueLetterCount++;
                    
                    // Rarer letters get higher scores
                    switch (letter) {
                        case 'q': case 'x': case 'z': case 'j': scoreValue += 10; break;
                        case 'k': case 'v': case 'w': case 'y': scoreValue += 8; break;
                        case 'b': case 'f': case 'g': case 'h': case 'p': scoreValue += 5; break;
                        default: scoreValue += 1; break;
                    }
                }
            }
            
            scoreValue += uniqueLetterCount * 3;
            scoredWords.push_back(std::make_pair(scoreValue, currentWord));
        }
        
        // Sort words by score (descending)
        std::sort(scoredWords.begin(), scoredWords.end(), 
                 [](const auto& a, const auto& b) {
                     // If scores differ, sort by score
                     if (a.first != b.first) {
                         return a.first > b.first;
                     }
                     // If scores are equal, sort by length
                     return a.second.length() > b.second.length();
                 });
        
        // Update requiredWords with sorted order
        requiredWords.clear();
        for (const auto& [score, word] : scoredWords) { // Using structured binding
            requiredWords.push_back(word);
        }
    }

    [[nodiscard]] bool solve() {
        // Clear any existing solutions
        solutions.clear();
        solutionHashes.clear();
        solutionFound.store(false);
        
        // Create a thread pool with hardware concurrency
        ThreadPool threadPool;
        int numThreads = threadPool.num_threads();
        
        // First calculate all possible placements for the first word
        struct FirstPlacement {
            int startRow, startCol, direction;
        };
        
        std::vector<FirstPlacement> firstPlacements;
        const std::string& firstWord = requiredWords[0];
        const size_t firstWordLen = firstWord.length();
        
        for (int direction = 0; direction < 8; direction++) {
            // Calculate valid bounds for this direction
            const int dirMaxOffset = static_cast<int>(firstWordLen) - 1;
            int minRow = (directionOffsetY[direction] < 0) ? dirMaxOffset : 0;
            int minCol = (directionOffsetX[direction] < 0) ? dirMaxOffset : 0;
            int maxRow = gridHeight - ((directionOffsetY[direction] > 0) ? dirMaxOffset : 0);
            int maxCol = gridWidth - ((directionOffsetX[direction] > 0) ? dirMaxOffset : 0);
            
            for (int row = minRow; row < maxRow; row++) {
                for (int col = minCol; col < maxCol; col++) {
                    firstPlacements.push_back({row, col, direction});
                }
            }
        }
        
        // If too few placements, just solve sequentially
        if (firstPlacements.size() < 20) [[unlikely]] {
            return solveSequential();
        }

        // Divide placements into batches for each thread
        int placementsPerThread = (firstPlacements.size() + numThreads - 1) / numThreads;
        std::vector<std::future<bool>> futures;
        
        for (size_t batchStart = 0; batchStart < firstPlacements.size(); batchStart += placementsPerThread) {
            size_t batchEnd = std::min(batchStart + placementsPerThread, firstPlacements.size());
            
            futures.push_back(threadPool.enqueue(
                [this, &firstPlacements, batchStart, batchEnd]() {
                    for (size_t placementIdx = batchStart; placementIdx < batchEnd; placementIdx++) {
                        // Skip if solution already found
                        if (!findAllSolutions && solutionFound.load(std::memory_order_acquire)) [[unlikely]] {
                            return true;
                        }
                        
                        const auto& placement = firstPlacements[placementIdx];
                        
                        // Initialize grid
                        std::vector<char> grid(gridWidth * gridHeight, '.');
                        std::vector<WordPlacement> placements;
                        
                        // Place first word
                        const std::string& word = requiredWords[0];
                        const size_t wordLen = word.length();
                        
                        int gridPos = flatIndex(placement.startRow, placement.startCol);
                        const int dirOffset = directionOffsets[placement.direction];
                        
                        for (size_t charIdx = 0; charIdx < wordLen; charIdx++) {
                            grid[gridPos] = word[charIdx];
                            gridPos += dirOffset;
                        }
                        
                        // Add to placements
                        placements.push_back({word, placement.startRow, placement.startCol, placement.direction, wordLen});
                        
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
            
            if (!findAllSolutions && success) [[unlikely]] {
                break;
            }
        }
        
        return success;
    }
    
    [[nodiscard]] bool solveSequential() {
        std::vector<char> grid(gridWidth * gridHeight, '.');
        std::vector<WordPlacement> placements;
        
        const std::string& firstWord = requiredWords[0];
        const size_t firstWordLen = firstWord.length();
        bool found = false;
        
        for (int direction = 0; direction < 8 && !found; direction++) {
            const int dirMaxOffset = static_cast<int>(firstWordLen) - 1;
            int minRow = (directionOffsetY[direction] < 0) ? dirMaxOffset : 0;
            int minCol = (directionOffsetX[direction] < 0) ? dirMaxOffset : 0;
            int maxRow = gridHeight - ((directionOffsetY[direction] > 0) ? dirMaxOffset : 0);
            int maxCol = gridWidth - ((directionOffsetX[direction] > 0) ? dirMaxOffset : 0);
            
            for (int row = minRow; row < maxRow && !found; row++) {
                for (int col = minCol; col < maxCol && !found; col++) {
                    // Reset grid
                    std::fill(grid.begin(), grid.end(), '.');
                    placements.clear();
                    
                    // Place first word
                    int gridPos = flatIndex(row, col);
                    const int dirOffset = directionOffsets[direction];
                    
                    for (size_t charIdx = 0; charIdx < firstWordLen; charIdx++) {
                        grid[gridPos] = firstWord[charIdx];
                        gridPos += dirOffset;
                    }
                    
                    placements.push_back({firstWord, row, col, direction, firstWordLen});
                    
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

    [[nodiscard]] const std::vector<std::vector<std::vector<char>>>& getSolutions() const noexcept {
        return solutions;
    }

private:
    // Check if word exists in the grid
    [[nodiscard]] bool wordExists(std::string_view word, const std::vector<char>& grid) const {
        const size_t wordLen = word.length();
        
        // Optimize for single letter words
        if (wordLen == 1) [[unlikely]] {
            char letterToFind = word[0];
            for (int gridPos = 0; gridPos < gridWidth * gridHeight; gridPos++) {
                if (grid[gridPos] == letterToFind) {
                    return true;
                }
            }
            return false;
        }

        char firstChar = word[0];
        
        // Find occurrences of first letter
        for (int row = 0; row < gridHeight; row++) {
            for (int col = 0; col < gridWidth; col++) {
                if (grid[flatIndex(row, col)] == firstChar) {
                    // Check all directions
                    for (int direction = 0; direction < 8; direction++) {
                        int endRow = row + directionOffsetY[direction] * (wordLen - 1);
                        int endCol = col + directionOffsetX[direction] * (wordLen - 1);
                        
                        // Skip out of bounds
                        if (endRow < 0 || endRow >= gridHeight || endCol < 0 || endCol >= gridWidth) {
                            continue;
                        }
                        
                        // Check if word matches
                        bool match = true;
                        for (size_t charIdx = 1; charIdx < wordLen; charIdx++) {
                            int checkRow = row + directionOffsetY[direction] * charIdx;
                            int checkCol = col + directionOffsetX[direction] * charIdx;
                            if (grid[flatIndex(checkRow, checkCol)] != word[charIdx]) {
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
    [[nodiscard]] bool canPlaceWord(std::string_view word, int row, int col, int direction, const std::vector<char>& grid) const {
        const size_t wordLen = word.length();
        
        // Check bounds
        int endRow = row + directionOffsetY[direction] * (wordLen - 1);
        int endCol = col + directionOffsetX[direction] * (wordLen - 1);
        
        if (endRow < 0 || endRow >= gridHeight || endCol < 0 || endCol >= gridWidth) {
            return false;
        }
        
        // Check if word can be placed
        int gridPos = flatIndex(row, col);
        int dirOffset = directionOffsets[direction];
        
        for (size_t charIdx = 0; charIdx < wordLen; charIdx++) {
            if (grid[gridPos] != '.' && grid[gridPos] != word[charIdx]) {
                return false;
            }
            gridPos += dirOffset;
        }
        
        return true;
    }

    // Check if grid is valid
    [[nodiscard]] bool isValidGrid(const std::vector<char>& grid) const {
        // Check required words
        for (const auto& requiredWord : requiredWords) {
            if (!wordExists(requiredWord, grid)) {
                return false;
            }
        }
        
        // Check forbidden words
        for (const auto& forbiddenWord : forbiddenWords) {
            if (wordExists(forbiddenWord, grid)) {
                return false;
            }
        }
        
        return true;
    }

    // Fill empty cells
    bool fillEmptyCells(int row, int col, std::vector<char>& grid) {
        // Check for early termination
        if (!findAllSolutions && solutionFound.load(std::memory_order_acquire)) [[unlikely]] {
            return false;
        }
        
        // If we've filled the entire grid, check if it's valid
        if (row >= gridHeight) [[unlikely]] {
            if (isValidGrid(grid)) {
                // Convert to string for hash
                std::string gridHash(grid.begin(), grid.end());
                
                // Lock for updating solutions
                std::lock_guard<std::mutex> lock(solutionsMutex);
                
                if (solutionHashes.find(gridHash) == solutionHashes.end()) {
                    solutionHashes.insert(gridHash);
                    
                    // Convert flat grid to 2D for output
                    std::vector<std::vector<char>> solution(gridHeight, std::vector<char>(gridWidth));
                    for (int r = 0; r < gridHeight; r++) {
                        for (int c = 0; c < gridWidth; c++) {
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
        if (nextCol >= gridWidth) {
            nextRow++;
            nextCol = 0;
        }
        
        // If cell is already filled, move to next
        if (grid[flatIndex(row, col)] != '.') {
            return fillEmptyCells(nextRow, nextCol, grid);
        }
        
        // Try each letter
        for (char letter = 'a'; letter <= 'z'; letter++) {
            grid[flatIndex(row, col)] = letter;
            
            // Check if this creates any forbidden words
            bool validPlacement = true;
            for (const auto& forbiddenWord : forbiddenWords) {
                if (forbiddenWord.length() == 1 && forbiddenWord[0] == letter) {
                    validPlacement = false;
                    break;
                }
                
                // Only check words that contain this letter (optimization)
                if (forbiddenWord.find(letter) != std::string::npos) {
                    if (wordExists(forbiddenWord, grid)) {
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
        if (!findAllSolutions && solutionFound.load(std::memory_order_acquire)) [[unlikely]] {
            return false;
        }
        
        // If all words placed, fill remaining cells
        if (wordIndex >= static_cast<int>(requiredWords.size())) [[unlikely]] {
            return fillEmptyCells(0, 0, grid);
        }
        
        const std::string& currentWord = requiredWords[wordIndex];
        const size_t wordLen = currentWord.length();
        bool foundSolution = false;
        
        // Track changes for efficient backtracking
        struct CellChange {
            int position;
            char previousChar;
        };
        
        // Use thread_local storage to avoid repeated allocation/deallocation during recursion
        static thread_local std::vector<CellChange> changes;
        changes.clear();
        changes.reserve(wordLen);
        
        // Try all possible positions
        for (int direction = 0; direction < 8; direction++) {
            // Boundary calculations
            const int dirMaxOffset = static_cast<int>(wordLen) - 1;
            int minRow = (directionOffsetY[direction] < 0) ? dirMaxOffset : 0;
            int minCol = (directionOffsetX[direction] < 0) ? dirMaxOffset : 0;
            int maxRow = gridHeight - ((directionOffsetY[direction] > 0) ? dirMaxOffset : 0);
            int maxCol = gridWidth - ((directionOffsetX[direction] > 0) ? dirMaxOffset : 0);
            
            for (int row = minRow; row < maxRow; row++) {
                for (int col = minCol; col < maxCol; col++) {
                    // Check if word can be placed
                    if (canPlaceWord(currentWord, row, col, direction, grid)) {
                        // Place word and track changes
                        changes.clear();
                        int gridPos = flatIndex(row, col);
                        const int dirOffset = directionOffsets[direction];
                        
                        for (size_t charIdx = 0; charIdx < wordLen; charIdx++) {
                            if (grid[gridPos] != currentWord[charIdx]) {
                                changes.push_back({gridPos, grid[gridPos]});
                                grid[gridPos] = currentWord[charIdx];
                            }
                            gridPos += dirOffset;
                        }
                        
                        // Add to placements
                        placements.push_back({currentWord, row, col, direction, wordLen});
                        
                        // Try to place next word
                        bool result = placeWords(wordIndex + 1, placements, grid);
                        
                        if (result) {
                            foundSolution = true;
                            if (!findAllSolutions) [[unlikely]] {
                                return true;
                            }
                        }
                        
                        // Backtrack - restore changed cells
                        placements.pop_back();
                        for (const auto& change : changes) {
                            grid[change.position] = change.previousChar;
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

    std::string_view inputFilePath{argv[1]};
    std::string_view outputFilePath{argv[2]};
    std::string_view solveMode{argv[3]};

    // Read input
    std::ifstream inputFile{std::string(inputFilePath)};
    if (!inputFile) {
        std::cerr << "Error opening input file: " << inputFilePath << std::endl;
        return 1;
    }

    int gridWidth, gridHeight;
    inputFile >> gridWidth >> gridHeight;
    inputFile.ignore(); // Skip the newline

    std::vector<std::string> requiredWords;
    std::vector<std::string> forbiddenWords;

    std::string line;
    requiredWords.reserve(20);
    forbiddenWords.reserve(20);
    
    while (std::getline(inputFile, line)) {
        if (line.empty()) {
            continue;
        }
        
        char wordType = line[0];
        std::string word = line.substr(2);
        
        // Trim trailing whitespace
        size_t endPos = word.find_last_not_of(" \n\r\t");
        if (endPos != std::string::npos) {
            word = word.substr(0, endPos + 1);
        }
        
        if (wordType == '+') {
            requiredWords.push_back(word);
        } else if (wordType == '-') {
            forbiddenWords.push_back(word);
        }
    }

    inputFile.close();

    // Check if puzzle is solvable
    bool impossiblePuzzle = false;
    for (const auto& word : requiredWords) {
        if (word.length() > std::max(gridWidth, gridHeight)) {
            impossiblePuzzle = true;
            break;
        }
    }

    // Solve puzzle
    std::vector<std::vector<std::vector<char>>> solutions;

    if (!impossiblePuzzle) {
        bool findAllSolutions = (solveMode == "all_solutions");
        
        auto startTime = std::chrono::high_resolution_clock::now();
        
        OptimizedInverseWordSearch solver(gridWidth, gridHeight, requiredWords, forbiddenWords, findAllSolutions);
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
    std::ofstream outputFile{std::string(outputFilePath)};
    if (!outputFile) {
        std::cerr << "Error opening output file: " << outputFilePath << std::endl;
        return 1;
    }

    if (solutions.empty()) {
        outputFile << "No solutions found" << std::endl;
    } else if (solveMode == "all_solutions") {
        outputFile << solutions.size() << " solution(s)" << std::endl;
        for (const auto& solution : solutions) {
            printGrid(outputFile, solution);
        }
    } else {
        printGrid(outputFile, solutions[0]);
    }

    outputFile.close();

    return 0;
}