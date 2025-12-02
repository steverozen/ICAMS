/*
 * Traced version of segment.cpp with debug output
 *
 * Compile:
 *   g++ -std=c++11 -o segment_traced segment_traced.cpp
 *
 * Example usage:
 *   ./segment_traced ATATAT ATATGG
 *   ./segment_traced ATATAT ATATGG --trace
 *   ./segment_traced CGCGCG CGCGAA --trace
 *   ./segment_traced TTTTT TTTGGG --trace
 *   ./segment_traced ACGTACGT ACGTACGTNN --trace
 *
 * Example inputs explained:
 *   1. Indel="ATATAT", Context="ATATGG"
 *      Perfect AT repeat, continues 2x in context
 *      Expected: unit=AT, internal_reps=4, spacer=0, prime3_reps=4
 *
 *   2. Indel="CGCGCG", Context="CGCGAA"
 *      CG repeat, continues 2x in context
 *      Expected: unit=CG, internal_reps=4, spacer=0, prime3_reps=4
 *
 *   3. Indel="TTTTT", Context="TTTGGG"
 *      T repeat (5 copies), continues 3x in context
 *      Expected: unit=T, internal_reps=4, spacer=0, prime3_reps=3
 *
 *   4. Indel="ACGTACGT", Context="ACGTACGTNN"
 *      ACGT repeat (2 copies), continues 2x in context
 *      Expected: unit=ACGT, internal_reps=4, spacer=0, prime3_reps=8
 */

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstring>

using namespace std;

// Global trace flag - default to true for detailed output
bool TRACE = true;

/*
 * Comparison function for sorting segmentations by priority
 *
 * Segmentation vector format: [unit_length, internal_reps, spacer_length, prime3_reps]
 *
 * Priority order (best segmentation has):
 * 1. HIGHEST 3' context repeats [3] - most important
 * 2. HIGHEST internal repeats [1]
 * 3. LOWEST spacer length [2] - prefer clean repeats
 * 4. LOWEST unit length [0] - prefer simpler repeat units
 */
bool compareSegmentations(std::vector<int> v1, std::vector<int> v2, bool trace = false){
    if(trace){
        cout << "  [COMPARE] v1=[" << v1[0] << "," << v1[1] << "," << v1[2] << "," << v1[3] << "]"
             << " vs v2=[" << v2[0] << "," << v2[1] << "," << v2[2] << "," << v2[3] << "]";
    }

    // Priority 1: HIGHEST 3' context repeats (most important)
    // Indels are more likely in regions that continue repeating into flanking sequence
    if(v1[3] > v2[3]){
        if(trace) cout << " -> v1 wins (higher prime3_reps)" << endl;
        return true;   // v1 is better
    }else if(v1[3] < v2[3]){
        if(trace) cout << " -> v2 wins (higher prime3_reps)" << endl;
        return false;  // v2 is better
    } else{
        // Priority 2: HIGHEST internal repeats (when 3' repeats are tied)
        // More internal repetition suggests stronger repeat-mediated mechanism
        if(v1[1] > v2[1]){
            if(trace) cout << " -> v1 wins (higher internal_reps)" << endl;
            return true;
        }else if(v1[1] < v2[1]){
            if(trace) cout << " -> v2 wins (higher internal_reps)" << endl;
            return false;
        }else{
            // Priority 3: LOWEST spacer length (when internal repeats are tied)
            // Cleaner repeat patterns (less leftover bases) are preferred
            if(v1[2] < v2[2]){
                if(trace) cout << " -> v1 wins (lower spacer)" << endl;
                return true;
            }else if(v1[2] > v2[2]){
                if(trace) cout << " -> v2 wins (lower spacer)" << endl;
                return false;
            }else{
                // Priority 4: LOWEST unit length (when spacer is tied)
                // Simpler repeat units (e.g., "AT" better than "ATAT")
                if(v1[0] < v2[0]){
                    if(trace) cout << " -> v1 wins (lower unit_length)" << endl;
                    return true;
                }else if(v1[0] > v2[0]){
                    if(trace) cout << " -> v2 wins (lower unit_length)" << endl;
                    return false;
                }else{
                    if(trace) cout << " -> tie (equal)" << endl;
                    return false;  // Equal segmentations
                }
            }
        }
    }
}


/*
 * LTRS: Longest Tandem Repeat Search
 *
 * Counts how many consecutive times a repeat unit appears in a target sequence
 *
 * Parameters:
 *   str1_start to str1_end: The repeat unit to search for
 *   str2_start to str2_end: The target sequence to search in
 *   trace: If true, print debug information
 *
 * Returns: Total length of tandem repeats found (not the count, but total bases)
 *
 * Example:
 *   unit = "AT", target = "ATATGC"
 *   Returns: 4 (two copies of "AT" = 4 bases)
 */
int ltrs(std::string::iterator str1_start, std::string::iterator str1_end,
         std::string::iterator str2_start, std::string::iterator str2_end, bool trace = false){

    // Extract strings for debugging
    string unit(str1_start, str1_end);
    string target(str2_start, str2_end);

    if(trace){
        cout << "  [LTRS] Searching for unit=\"" << unit << "\" in target=\"" << target << "\"" << endl;
    }

    int increment_unit = std::distance(str1_start, str1_end);  // Length of repeat unit
    int coverage = 0;  // Total bases covered by tandem repeats
    bool equal_unit = true; // assume all substrings are identical

    int repeat_count = 0; // For trace output

    // Iterate through target sequence, checking if repeat unit matches
    for(auto i = str2_start; i != str2_end;){
        equal_unit = true;  // Reset for each potential repeat
        for(auto j = str1_start; j != str1_end; j++, i++){
            if(i == str2_end){
                equal_unit = false;
                break;
            }
            equal_unit = equal_unit & (*i == *j);

            if(!equal_unit){
                break;  // Unit doesn't match, stop checking
            }
        }

        if(!equal_unit){ // if is_equal_unit is false, break out.
            break;  // No more tandem repeats found
        }
        coverage += increment_unit;  // Add this repeat unit to coverage
        repeat_count++;

        if(trace){
            cout << "    Found repeat #" << repeat_count << ", coverage now: " << coverage << endl;
        }
    }

    if(trace){
        cout << "  [LTRS] Result: " << coverage << " bases (" << repeat_count << " repeats)" << endl;
    }

    return coverage;
}


/*
 * segmentSingle: Find optimal repeat unit segmentation for an indel
 *
 * This function tries every possible repeat unit size (1 to string length)
 * and scores each segmentation based on:
 * - How well the unit repeats internally
 * - How well it repeats in the flanking context
 * - How little "spacer" (leftover bases) remains
 *
 * Parameters:
 *   string: The indel sequence to segment
 *   context: The flanking sequence (3' context)
 *   trace: If true, print debug information
 *
 * Returns: Best segmentation as [unit_length, internal_reps, spacer_length, prime3_reps]
 */
std::vector<int> segmentSingle(std::string &string, std::string &context, bool trace = false){

    if(trace){
        cout << "\n[SEGMENT_SINGLE] Starting segmentation" << endl;
        cout << "  Input string: \"" << string << "\" (length=" << string.size() << ")" << endl;
        cout << "  Input context: \"" << context << "\" (length=" << context.size() << ")" << endl;
        cout << endl;
    }

    // Get the total length of the indel sequence
    int string_size = string.size();

    // Create scoring matrix: one row per possible unit size, 4 columns per row
    // Column 0: unit length
    // Column 1: internal repeats (bases)
    // Column 2: spacer length (leftover bases)
    // Column 3: 3' context repeats (bases)
    std::vector<std::vector<int> > scores(string_size, std::vector<int>(4));

    int i = 0;

    if(trace){
        cout << "Trying all possible unit sizes:" << endl;
    }

    // Try every possible repeat unit size (from 1 to string length)
    for(auto unit_iter = string.begin() + 1; unit_iter <= string.end(); ++unit_iter, ++i){

        std::string unit(string.begin(), unit_iter);

        if(trace){
            cout << "\n--- Unit size " << (i+1) << ": \"" << unit << "\" ---" << endl;
        }

        // 1. Set unit length (1, 2, 3, ... string_size)
        scores[i][0] = i + 1;

        // 2. Count how many times this unit repeats INSIDE the indel
        //    Uses LTRS to find tandem repeats after the first occurrence of the unit
        scores[i][1] = ltrs(string.begin(), unit_iter, unit_iter, string.end(), trace);

        // 3. Calculate spacer length (leftover bases that don't fit the repeat pattern)
        //    spacer = total_length - unit_length - internal_repeat_length
        scores[i][2] = string_size - scores[i][0] - scores[i][1];

        // 4. Count how many times this unit repeats in the 3' FLANKING CONTEXT
        //    This is crucial: indels often occur where repeats continue into flanking region
        scores[i][3] = ltrs(string.begin(), unit_iter, context.begin(), context.end(), trace);

        if(trace){
            cout << "  Score: [unit_len=" << scores[i][0]
                 << ", internal_reps=" << scores[i][1]
                 << ", spacer=" << scores[i][2]
                 << ", prime3_reps=" << scores[i][3] << "]" << endl;
        }
    }

    if(trace){
        cout << "\n[SORTING] Sorting segmentations by priority..." << endl;
        cout << "  Format: [unit_length, internal_reps, spacer_length, prime3_reps]" << endl;
        cout << "  Priority: 1) Highest prime3_reps  2) Highest internal_reps  3) Lowest spacer  4) Lowest unit_length" << endl;
    }

    // Sort all segmentations using priority comparison function
    // Best segmentation will be first (index 0)
    std::sort(scores.begin(), scores.end(),
              [trace](const vector<int>& a, const vector<int>& b) {
                  return compareSegmentations(a, b, trace);
              });

    if(trace){
        cout << "\n[RESULT] Best segmentation:" << endl;
        cout << "  Unit length: " << scores[0][0] << endl;
        cout << "  Internal repeats (bases): " << scores[0][1] << endl;
        cout << "  Spacer length: " << scores[0][2] << endl;
        cout << "  3' context repeats (bases): " << scores[0][3] << endl;

        // Show the actual unit
        std::string best_unit(string.begin(), string.begin() + scores[0][0]);
        cout << "  Unit sequence: \"" << best_unit << "\"" << endl;
    }

    // Return the best segmentation
    return scores[0];
}


/*
 * Main function for command-line testing
 */
int main(int argc, char* argv[]){
    // Parse command line arguments
    if(argc < 3){
        cerr << "Usage: " << argv[0] << " <indel_sequence> <context_sequence> [--trace]" << endl;
        cerr << "\nExamples:" << endl;
        cerr << "  " << argv[0] << " ATATAT ATATGG" << endl;
        cerr << "  " << argv[0] << " ATATAT ATATGG --trace" << endl;
        cerr << "  " << argv[0] << " CGCGCG CGCGAA --trace" << endl;
        cerr << "  " << argv[0] << " TTTTT TTTGGG --trace" << endl;
        return 1;
    }

    string indel_seq = argv[1];
    string context_seq = argv[2];

    // Check for --trace flag
    for(int i = 3; i < argc; i++){
        if(strcmp(argv[i], "--trace") == 0){
            TRACE = true;
        }
    }

    cout << "===========================================" << endl;
    cout << "Indel Segmentation Analysis" << endl;
    cout << "===========================================" << endl;

    // Run segmentation
    vector<int> result = segmentSingle(indel_seq, context_seq, TRACE);

    // Print summary
    cout << "\n===========================================" << endl;
    cout << "SUMMARY" << endl;
    cout << "===========================================" << endl;
    cout << "Indel sequence: " << indel_seq << endl;
    cout << "Context sequence: " << context_seq << endl;

    string unit(indel_seq.begin(), indel_seq.begin() + result[0]);
    cout << "\nBest segmentation:" << endl;
    cout << "  Unit: \"" << unit << "\"" << endl;
    cout << "  Unit length: " << result[0] << endl;
    cout << "  Internal repeats: " << result[1] << " bases" << endl;
    cout << "  Spacer: " << result[2] << " bases" << endl;
    cout << "  3' context repeats: " << result[3] << " bases" << endl;

    // Calculate repeat counts
    int internal_count = result[1] / result[0];
    int context_count = result[3] / result[0];
    int total_count = 1 + internal_count + context_count;

    cout << "\nRepeat counts:" << endl;
    cout << "  Initial: 1 copy (the unit itself)" << endl;
    cout << "  Internal: " << internal_count << " additional copies" << endl;
    cout << "  Context: " << context_count << " copies in flanking region" << endl;
    cout << "  Total: " << total_count << " copies" << endl;

    return 0;
}
