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

// Directions arrays must be defined
const int dx[8] = {1, 1, 0, -1, -1, -1, 0, 1};
const int dy[8] = {0, 1, 1, 1, 0, -1, -1, -1};

// Define WordPlacement structure
struct WordPlacement {
    std::string word;
    int row, col, dir;
    size_t length;
};

// Define printGrid function
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

// Work-stealing thread pool implementation
class ThreadPool {
private:
    // Structure for a task
    struct Task {
        std::function<void()> func;
        Task(std::function<void()> f) : func(std::move(f)) {}
    };

    // Work-stealing queue for each thread
    class WorkStealingQueue {
    private:
        std::deque<std::shared_ptr<Task>> queue;
        std::mutex mutex;

    public:
        void push(std::shared_ptr<Task> task) {
            std::lock_guard<std::mutex> lock(mutex);
            queue.push_front(std::move(task));
        }

        bool pop(std::shared_ptr<Task>& task) {
            std::lock_guard<std::mutex> lock(mutex);
            if (queue.empty()) return false;
            task = std::move(queue.front());
            queue.pop_front();
            return true;
        }

        bool steal(std::shared_ptr<Task>& task) {
            std::lock_guard<std::mutex> lock(mutex);
            if (queue.empty()) return false;
            task = std::move(queue.back());
            queue.pop_back();
            return true;
        }

        bool empty() {
            std::lock_guard<std::mutex> lock(mutex);
            return queue.empty();
        }
    };

    std::vector<std::unique_ptr<WorkStealingQueue>> queues;
    std::vector<std::thread> threads;
    std::atomic<bool> stop{false};
    std::atomic<size_t> index{0};

    void workerThread(size_t id) {
        while (!stop) {
            std::shared_ptr<Task> task;
            if (popTask(task, id)) {
                // Execute the task
                task->func();
            } else {
                // Try to steal tasks from other queues
                if (!stealTask(task)) {
                    // If no tasks to steal, yield to other threads
                    std::this_thread::yield();
                } else {
                    task->func();
                }
            }
        }
    }

    bool popTask(std::shared_ptr<Task>& task, size_t id) {
        return queues[id]->pop(task);
    }

    bool stealTask(std::shared_ptr<Task>& task) {
        const size_t count = queues.size();
        for (size_t i = 0; i < count; ++i) {
            const size_t idx = (index.fetch_add(1, std::memory_order_relaxed) + i) % count;
            if (queues[idx]->steal(task)) {
                return true;
            }
        }
        return false;
    }

public:
    ThreadPool(size_t threadCount = std::thread::hardware_concurrency()) {
        queues.resize(threadCount);
        for (size_t i = 0; i < threadCount; ++i) {
            queues[i] = std::make_unique<WorkStealingQueue>();
        }

        for (size_t i = 0; i < threadCount; ++i) {
            threads.emplace_back([this, i] { workerThread(i); });
        }
    }

    ~ThreadPool() {
        stop = true;
        for (auto& thread : threads) {
            if (thread.joinable()) thread.join();
        }
    }

    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args) -> std::future<std::invoke_result_t<F, Args...>> {
        using return_type = std::invoke_result_t<F, Args...>;
        
        auto task = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );
        
        std::future<return_type> result = task->get_future();
        
        auto wrapper = std::make_shared<Task>([task]() { (*task)(); });
        
        size_t idx = index.fetch_add(1, std::memory_order_relaxed) % queues.size();
        queues[idx]->push(std::move(wrapper));
        
        return result;
    }
};

// Modify the InverseWordSearch class to support parallelization
class ParallelInverseWordSearch {
private:
    int width, height;
    std::vector<std::string> requiredWords;
    std::vector<std::string> forbiddenWords;
    bool findAllSolutions;
    std::atomic<bool> solutionFound{false};
    std::mutex solutionsMutex;
    std::vector<std::vector<std::vector<char>>> solutions;
    std::set<std::string> solutionHashes;

    // Helper function to convert 2D coordinates to 1D index
    inline int flatIndex(int row, int col) const {
        return row * width + col;
    }

public:
    ParallelInverseWordSearch(int w, int h,
                      const std::vector<std::string>& required,
                      const std::vector<std::string>& forbidden,
                      bool findAll)
        : width(w), height(h), requiredWords(required), forbiddenWords(forbidden),
          findAllSolutions(findAll) {
    }

    void optimizeWordOrder() {
        // Same implementation as the original
        const size_t wordCount = requiredWords.size();
        if (wordCount <= 1) return;
        
        // Static letter rarity scores (higher = rarer)
        // Using uint8_t for better cache usage
        static const uint8_t letterScore[26] = {
            18, 98, 79, 67, 13, 91, 90, 49, 49, 99,
            97, 66, 79, 50, 48, 90, 99, 54, 51, 37,
            80, 96, 84, 99, 91, 99
        };
        
        // For each word length, find where groups start and end
        size_t currentLength = requiredWords[0].length();
        size_t groupStart = 0;
        
        // Temporary storage for a single group - allocated once and reused
        std::vector<std::pair<uint16_t, uint16_t>> scoreIndices; // (score, index)
        
        for (size_t i = 1; i <= wordCount; i++) {
            if (i == wordCount || requiredWords[i].length() != currentLength) {
                const size_t groupSize = i - groupStart;
                
                // Only need to process groups with multiple words
                if (groupSize > 1) {
                    // Resize our temporary vector for this group
                    scoreIndices.resize(groupSize);
                    
                    // Calculate scores for this group
                    for (size_t j = 0; j < groupSize; j++) {
                        const std::string& word = requiredWords[groupStart + j];
                        
                        // Fixed-size array is much faster than unordered_set for this case
                        bool seen[26] = {false};
                        uint8_t uniqueCount = 0;
                        uint16_t wordScore = 0;
                        
                        // Process each character once
                        for (char c : word) {
                            uint8_t idx = c - 'a';
                            if (!seen[idx]) {
                                seen[idx] = true;
                                uniqueCount++;
                                wordScore += letterScore[idx];
                            }
                        }
                        
                        // Final score with uniqueness bonus
                        scoreIndices[j] = {
                            static_cast<uint16_t>(wordScore + (uniqueCount * 10)), 
                            static_cast<uint16_t>(j)
                        };
                    }
                    
                    // Sort by score (higher first)
                    std::sort(scoreIndices.begin(), scoreIndices.end(),
                        [](const auto& a, const auto& b) {
                            return a.first > b.first;
                        }
                    );
                    
                    // Create a temporary copy of the reordered group
                    std::vector<std::string> reordered(groupSize);
                    for (size_t j = 0; j < groupSize; j++) {
                        reordered[j] = requiredWords[groupStart + scoreIndices[j].second];
                    }
                    
                    // Copy back to original array
                    std::copy(reordered.begin(), reordered.end(), 
                              requiredWords.begin() + groupStart);
                }
                
                // Move to next group
                if (i < wordCount) {
                    currentLength = requiredWords[i].length();
                    groupStart = i;
                }
            }
        }
    }

    bool solve() {
        // Clear any existing solutions
        solutions.clear();
        solutionHashes.clear();
        solutionFound.store(false);
        
        // Determine number of threads to use
        unsigned int numThreads = std::thread::hardware_concurrency();
        if (numThreads == 0) numThreads = 4; // Default if hardware_concurrency fails
        
        // Create thread pool
        ThreadPool pool(numThreads);
        std::vector<std::future<bool>> futures;
        
        // Create tasks for parallel execution - start with different positions and directions for the first word
        for (int dir = 0; dir < 8; dir++) {
            for (int row = 0; row < height; row++) {
                for (int col = 0; col < width; col++) {
                    // Skip invalid placements for the first word
                    const std::string& word = requiredWords[0];
                    const size_t wordLen = word.length();
                    
                    // Pre-compute end coordinates
                    const int endRow = row + (wordLen - 1) * dy[dir];
                    const int endCol = col + (wordLen - 1) * dx[dir];
                    
                    // Skip if out of bounds
                    if (endRow < 0 || endRow >= height || endCol < 0 || endCol >= width) {
                        continue;
                    }
                    
                    // Create a new task
                    futures.push_back(pool.enqueue(
                        [this, row, col, dir]() {
                            // Initialize grid for this thread
                            std::vector<char> grid(width * height, '.');
                            std::vector<WordPlacement> placements;
                            
                            // Try to place the first word
                            const std::string& word = requiredWords[0];
                            const size_t wordLen = word.length();
                            
                            // Place the first word
                            for (size_t i = 0; i < wordLen; i++) {
                                const int r = row + i * dy[dir];
                                const int c = col + i * dx[dir];
                                grid[flatIndex(r, c)] = word[i];
                            }
                            
                            // Add to placements
                            placements.push_back({word, row, col, dir, wordLen});
                            
                            // Try to place the remaining words
                            bool result = placeWordsParallel(1, placements, grid);
                            return result;
                        }
                    ));
                }
            }
        }
        
        // Wait for results
        bool anySuccess = false;
        for (auto& future : futures) {
            anySuccess = future.get() || anySuccess;
            
            // If we only need one solution and found it, break early
            if (!findAllSolutions && anySuccess) {
                break;
            }
        }
        
        return anySuccess;
    }

    const std::vector<std::vector<std::vector<char>>>& getSolutions() const {
        return solutions;
    }

private:
    // Generate a string hash of the grid for duplicate checking
    std::string gridHash(const std::vector<char>& grid) const {
        return std::string(grid.begin(), grid.end());
    }

    // Check if a word can be placed
    inline bool canPlaceWord(const std::string& word, int row, int col, int dir, 
                            const std::vector<char>& grid) const {
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
    inline void placeWord(const std::string& word, int row, int col, int dir, 
                         std::vector<char>& grid) {
        const size_t wordLen = word.length();
        
        for (size_t i = 0; i < wordLen; i++) {
            const int r = row + i * dy[dir];
            const int c = col + i * dx[dir];
            grid[flatIndex(r, c)] = word[i];
        }
    }

    // Remove a word from the grid
    void removeWord(const std::string& word, int row, int col, int dir, 
                   const std::vector<WordPlacement>& placements,
                   std::vector<char>& grid) {
        const size_t wordLen = word.length();
        
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
    bool wordExists(const std::string& word, const std::vector<char>& grid) const {
        const size_t wordLen = word.length();
        const char firstChar = word[0];
        
        // For single-letter words, need to check if they appear in the grid
        if (wordLen == 1) {
            for (int idx = 0; idx < width * height; idx++) {
                if (grid[idx] == firstChar) {
                    return true;
                }
            }
            return false;
        }
        
        // Find positions of the first letter
        std::vector<std::pair<int, int>> firstLetterPositions;
        firstLetterPositions.reserve(width * height / 10);
        
        for (int row = 0; row < height; row++) {
            for (int col = 0; col < width; col++) {
                if (grid[flatIndex(row, col)] == firstChar) {
                    firstLetterPositions.emplace_back(row, col);
                }
            }
        }
        
        // Check from those positions
        for (const auto& pos : firstLetterPositions) {
            const int row = pos.first;
            const int col = pos.second;
            
            for (int dir = 0; dir < 8; dir++) {
                const int endRow = row + (wordLen - 1) * dy[dir];
                const int endCol = col + (wordLen - 1) * dx[dir];
                
                if (endRow < 0 || endRow >= height || endCol < 0 || endCol >= width) {
                    continue;
                }
                
                bool match = true;
                for (size_t i = 1; i < wordLen; i++) {
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

    // Check if the grid is valid
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
    bool fillEmptyCells(int row, int col, std::vector<char>& grid,
                       const std::unordered_set<char>& forbiddenLetters = std::unordered_set<char>(),
                       const std::vector<std::string>& multiLetterForbidden = std::vector<std::string>()) {
        // If we already found a solution and don't need all, return early
        if (!findAllSolutions && solutionFound.load(std::memory_order_acquire)) {
            return false;
        }

        // If we've filled the entire grid, check if it's valid
        if (row == height) {
            if (isValidGrid(grid)) {
                // Create hash for duplicate detection
                std::string hash = gridHash(grid);
                
                // Lock for updating shared solution data
                std::lock_guard<std::mutex> lock(solutionsMutex);
                
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
                    
                    // Set solution found flag if we only need one
                    if (!findAllSolutions) {
                        solutionFound.store(true, std::memory_order_release);
                        return true;
                    }
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
            return fillEmptyCells(nextRow, nextCol, grid, forbiddenLetters, multiLetterForbidden);
        }
        
        // Try each possible letter
        for (char c = 'a'; c <= 'z'; c++) {
            // Skip forbidden letters
            if (forbiddenLetters.find(c) != forbiddenLetters.end()) {
                continue;
            }
            
            grid[flatIndex(row, col)] = c;
            
            // Check if placing this letter would create a forbidden word
            bool createsForbiddenWord = false;
            for (const auto& word : forbiddenWords) {
                if (wordExists(word, grid)) {
                    createsForbiddenWord = true;
                    break;
                }
            }
            
            if (!createsForbiddenWord) {
                bool result = fillEmptyCells(nextRow, nextCol, grid, forbiddenLetters, multiLetterForbidden);
                if (result) {
                    return true;
                }
            }
        }
        
        // Backtrack
        grid[flatIndex(row, col)] = '.';
        return false;
    }

    // Parallel version of placeWords
    bool placeWordsParallel(int wordIndex, std::vector<WordPlacement>& placements, std::vector<char>& grid) {
        // If we already found a solution and don't need all, return early
        if (!findAllSolutions && solutionFound.load(std::memory_order_acquire)) {
            return false;
        }

        // If all words are placed, fill the remaining cells
        if (wordIndex == requiredWords.size()) {
            return fillEmptyCells(0, 0, grid);
        }
        
        const std::string& word = requiredWords[wordIndex];
        const size_t wordLen = word.length();
        bool foundSolution = false;
        
        // Fixed-size array to track changes
        const int MAX_CHANGES = 20;
        std::pair<int, char> changes[MAX_CHANGES];
        int numChanges;
        
        // Pre-compute direction offsets
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
                        bool result = placeWordsParallel(wordIndex + 1, placements, grid);
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

// Update main to use the parallel solver
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

    // Solve the puzzle using the parallel solver
    std::vector<std::vector<std::vector<char>>> solutions;

    if (!impossiblePuzzle) {
        bool findAllSols = (mode == "all_solutions");
        ParallelInverseWordSearch solver(width, height, requiredWords, forbiddenWords, findAllSols);
        solver.optimizeWordOrder(); // Call optimizeWordOrder before solving
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