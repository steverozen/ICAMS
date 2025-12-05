/*
 * Rcpp interface for indel segmentation
 *
 * This file provides Rcpp wrappers for the segmentSingle function,
 * allowing direct C++ calls from R without using system2.
 */

/* This code was adapated from indelsig.tools.lib, 
https://github.com/Nik-Zainal-Group/indelsig.tools.lib,
which was released under this license:

# MIT License

Copyright (c) 2021 indelsig.tools.lib authors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <Rcpp.h>
#include <vector>
#include <string>
#include <algorithm>

using namespace Rcpp;

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
bool compareSegmentations(std::vector<int> v1, std::vector<int> v2) {
    // Priority 1: HIGHEST 3' context repeats (most important)
    if (v1[3] > v2[3]) {
        return true;   // v1 is better
    } else if (v1[3] < v2[3]) {
        return false;  // v2 is better
    } else {
        // Priority 2: HIGHEST internal repeats (when 3' repeats are tied)
        if (v1[1] > v2[1]) {
            return true;
        } else if (v1[1] < v2[1]) {
            return false;
        } else {
            // Priority 3: LOWEST spacer length (when internal repeats are tied)
            if (v1[2] < v2[2]) {
                return true;
            } else if (v1[2] > v2[2]) {
                return false;
            } else {
                // Priority 4: LOWEST unit length (when spacer is tied)
                if (v1[0] < v2[0]) {
                    return true;
                } else if (v1[0] > v2[0]) {
                    return false;
                } else {
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
 *
 * Returns: Total length of tandem repeats found (not the count, but total bases)
 */
int ltrs(std::string::iterator str1_start, std::string::iterator str1_end,
         std::string::iterator str2_start, std::string::iterator str2_end) {

    int increment_unit = std::distance(str1_start, str1_end);  // Length of repeat unit
    int coverage = 0;  // Total bases covered by tandem repeats
    bool equal_unit = true;

    // Iterate through target sequence, checking if repeat unit matches
    for (auto i = str2_start; i != str2_end;) {
        equal_unit = true;  // Reset for each potential repeat
        for (auto j = str1_start; j != str1_end; j++, i++) {
            if (i == str2_end) {
                equal_unit = false;
                break;
            }
            equal_unit = equal_unit & (*i == *j);

            if (!equal_unit) {
                break;  // Unit doesn't match, stop checking
            }
        }

        if (!equal_unit) {
            break;  // No more tandem repeats found
        }
        coverage += increment_unit;  // Add this repeat unit to coverage
    }

    return coverage;
}

/*
 * segmentSimple: Find optimal repeat unit segmentation for an indel
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
 *
 * Returns: Best segmentation as [unit_length, internal_reps, spacer_length, prime3_reps]
 */
std::vector<int> segmentSimple(std::string &string, std::string &context) {

    // Get the total length of the indel sequence
    int string_size = string.size();

    // Create scoring matrix: one row per possible unit size, 4 columns per row
    // Column 0: unit length
    // Column 1: internal repeats (bases)
    // Column 2: spacer length (leftover bases)
    // Column 3: 3' context repeats (bases)
    std::vector<std::vector<int> > scores(string_size, std::vector<int>(4));

    int i = 0;

    // Try every possible repeat unit size (from 1 to string length)
    for (auto unit_iter = string.begin() + 1; unit_iter <= string.end(); ++unit_iter, ++i) {

        // 1. Set unit length (1, 2, 3, ... string_size)
        scores[i][0] = i + 1;

        // 2. Count how many times this unit repeats INSIDE the indel
        scores[i][1] = ltrs(string.begin(), unit_iter, unit_iter, string.end());

        // 3. Calculate spacer length (leftover bases that don't fit the repeat pattern)
        scores[i][2] = string_size - scores[i][0] - scores[i][1];

        // 4. Count how many times this unit repeats in the 3' FLANKING CONTEXT
        scores[i][3] = ltrs(string.begin(), unit_iter, context.begin(), context.end());
    }

    // Sort all segmentations using priority comparison function
    // Best segmentation will be first (index 0)
    std::sort(scores.begin(), scores.end(),
              [](const std::vector<int>& a, const std::vector<int>& b) {
                  return compareSegmentations(a, b);
              });

    // Return the best segmentation
    return scores[0];
}

//' Segment a single indel using Rcpp interface
//'
//' This function provides direct C++ interface for indel segmentation without
//' calling an external binary via system2.
//'
//' @param ins_or_del Character indicating insertion ("i") or deletion ("d")
//' @param string A single indel sequence to segment (character scalar)
//' @param context A single flanking context sequence (character scalar)
//'
//' @return A list with the following elements:
//'   \item{unit}{Repeat unit sequence (character)}
//'   \item{unit_length}{Unit length (integer)}
//'   \item{internal_rep}{Internal repeat region (character)}
//'   \item{internal_reps}{Internal repeat count (integer)}
//'   \item{spacer}{Spacer sequence (character)}
//'   \item{spacer_length}{Spacer length (integer)}
//'   \item{prime3_rep}{3' flanking repeat region (character)}
//'   \item{prime3_reps}{3' flanking repeat count (integer)}
//'   \item{original_reps}{Original repeat count (integer)}
//'
//' @export
// [[Rcpp::export]]
Rcpp::List segment_simple_cpp(std::string ins_or_del, std::string string, std::string context) {

    // Call the segmentation function
    std::vector<int> result = segmentSimple(string, context);

    int unit_len = result[0];
    int internal_bases = result[1];
    int spacer_len = result[2];
    int prime3_bases = result[3];

    // Extract string components
    std::string unit = string.substr(0, unit_len);

    std::string internal_rep = "";
    if (internal_bases > 0) {
        internal_rep = string.substr(unit_len, internal_bases);
    }

    std::string spacer = "";
    if (spacer_len > 0) {
        int spacer_start = unit_len + internal_bases;
        spacer = string.substr(spacer_start, spacer_len);
    }

    std::string prime3_rep = "";
    if (prime3_bases > 0) {
        prime3_rep = context.substr(0, prime3_bases);
    }

    // Convert bases to repeat counts
    int internal_reps = 0;
    int prime3_reps = 0;
    if (unit_len != 0) {
        internal_reps = internal_bases / unit_len;
        prime3_reps = prime3_bases / unit_len;
    }

    // Calculate original_reps based on ins_or_del
    int original_reps;
    if (ins_or_del == "i") {
        original_reps = prime3_reps;
    } else {
        if (spacer_len == 0) {
            original_reps = internal_reps + 1 + prime3_reps;
        } else {
            original_reps = prime3_reps;
        }
    }

    // Return as named list
    return List::create(
        Named("unit") = unit,
        Named("unit_length") = unit_len,
        Named("internal_rep") = internal_rep,
        Named("internal_reps") = internal_reps,
        Named("spacer") = spacer,
        Named("spacer_length") = spacer_len,
        Named("prime3_rep") = prime3_rep,
        Named("prime3_reps") = prime3_reps,
        Named("original_reps") = original_reps
    );
}
