#ifndef included_Utilities_hpp
#define included_Utilities_hpp

#include "Utilities.h"

#include <vector>

namespace Utilities {

/************************************************************************
 * templated quicksort routines                                          *
 ************************************************************************/
template <class T> void quicksort(std::vector<T> &x) {
    if (x.size() <= 1u)
        return;
    T *arr = &x[0];
    bool test;
    long int i, ir, j, jstack, k, l, istack[100];
    T a, tmp_a;
    jstack = 0;
    l = 0;
    ir = x.size() - 1;
    while (1) {
        if (ir - l < 7) { // Insertion sort when subarray small enough.
            for (j = l + 1; j <= ir; j++) {
                a = arr[j];
                test = true;
                for (i = j - 1; i >= 0; i--) {
                    if (arr[i] < a) {
                        arr[i + 1] = a;
                        test = false;
                        break;
                    }
                    arr[i + 1] = arr[i];
                }
                if (test) {
                    i = l - 1;
                    arr[i + 1] = a;
                }
            }
            if (jstack == 0)
                return;
            ir = istack
                [jstack]; // Pop stack and begin a new round of partitioning.
            l = istack[jstack - 1];
            jstack -= 2;
        } else {
            k = (l + ir) /
                2; // Choose median of left, center and right elements as partitioning
                // element a. Also rearrange so that a(l) < a(l+1) < a(ir).
            tmp_a = arr[k];
            arr[k] = arr[l + 1];
            arr[l + 1] = tmp_a;
            if (arr[l] > arr[ir]) {
                tmp_a = arr[l];
                arr[l] = arr[ir];
                arr[ir] = tmp_a;
            }
            if (arr[l + 1] > arr[ir]) {
                tmp_a = arr[l + 1];
                arr[l + 1] = arr[ir];
                arr[ir] = tmp_a;
            }
            if (arr[l] > arr[l + 1]) {
                tmp_a = arr[l];
                arr[l] = arr[l + 1];
                arr[l + 1] = tmp_a;
            }
            // Scan up to find element > a
            j = ir;
            a = arr[l + 1]; // Partitioning element.
            for (i = l + 2; i <= ir; i++) {
                if (arr[i] < a)
                    continue;
                while (arr[j] > a) // Scan down to find element < a.
                    j--;
                if (j < i)
                    break; // Pointers crossed. Exit with partitioning complete.
                tmp_a = arr[i]; // Exchange elements of both arrays.
                arr[i] = arr[j];
                arr[j] = tmp_a;
            }
            arr[l + 1] = arr[j]; // Insert partitioning element in both arrays.
            arr[j] = a;
            jstack += 2;
            // Push pointers to larger subarray on stack, process smaller subarray immediately.
            if (ir - i + 1 >= j - l) {
                istack[jstack] = ir;
                istack[jstack - 1] = i;
                ir = j - 1;
            } else {
                istack[jstack] = j - 1;
                istack[jstack - 1] = l;
                l = i;
            }
        }
    }
}
template <class T1, class T2>
void quicksort(std::vector<T1> &x, std::vector<T2> &y) {
    if (x.size() <= 1u)
        return;
    T1 *arr = &x[0];
    T2 *brr = &y[0];
    bool test;
    long int i, ir, j, jstack, k, l, istack[100];
    T1 a, tmp_a;
    T2 b, tmp_b;
    jstack = 0;
    l = 0;
    ir = x.size() - 1;
    while (1) {
        if (ir - l < 7) { // Insertion sort when subarray small enough.
            for (j = l + 1; j <= ir; j++) {
                a = arr[j];
                b = brr[j];
                test = true;
                for (i = j - 1; i >= 0; i--) {
                    if (arr[i] < a) {
                        arr[i + 1] = a;
                        brr[i + 1] = b;
                        test = false;
                        break;
                    }
                    arr[i + 1] = arr[i];
                    brr[i + 1] = brr[i];
                }
                if (test) {
                    i = l - 1;
                    arr[i + 1] = a;
                    brr[i + 1] = b;
                }
            }
            if (jstack == 0)
                return;
            ir = istack
                [jstack]; // Pop stack and begin a new round of partitioning.
            l = istack[jstack - 1];
            jstack -= 2;
        } else {
            k = (l + ir) /
                2; // Choose median of left, center and right elements as partitioning
                // element a. Also rearrange so that a(l) ? a(l+1) ? a(ir).
            tmp_a = arr[k];
            arr[k] = arr[l + 1];
            arr[l + 1] = tmp_a;
            tmp_b = brr[k];
            brr[k] = brr[l + 1];
            brr[l + 1] = tmp_b;
            if (arr[l] > arr[ir]) {
                tmp_a = arr[l];
                arr[l] = arr[ir];
                arr[ir] = tmp_a;
                tmp_b = brr[l];
                brr[l] = brr[ir];
                brr[ir] = tmp_b;
            }
            if (arr[l + 1] > arr[ir]) {
                tmp_a = arr[l + 1];
                arr[l + 1] = arr[ir];
                arr[ir] = tmp_a;
                tmp_b = brr[l + 1];
                brr[l + 1] = brr[ir];
                brr[ir] = tmp_b;
            }
            if (arr[l] > arr[l + 1]) {
                tmp_a = arr[l];
                arr[l] = arr[l + 1];
                arr[l + 1] = tmp_a;
                tmp_b = brr[l];
                brr[l] = brr[l + 1];
                brr[l + 1] = tmp_b;
            }
            // Scan up to find element > a
            j = ir;
            a = arr[l + 1]; // Partitioning element.
            b = brr[l + 1];
            for (i = l + 2; i <= ir; i++) {
                if (arr[i] < a)
                    continue;
                while (arr[j] > a) // Scan down to find element < a.
                    j--;
                if (j < i)
                    break; // Pointers crossed. Exit with partitioning complete.
                tmp_a = arr[i]; // Exchange elements of both arrays.
                arr[i] = arr[j];
                arr[j] = tmp_a;
                tmp_b = brr[i];
                brr[i] = brr[j];
                brr[j] = tmp_b;
            }
            arr[l + 1] = arr[j]; // Insert partitioning element in both arrays.
            arr[j] = a;
            brr[l + 1] = brr[j];
            brr[j] = b;
            jstack += 2;
            // Push pointers to larger subarray on stack, process smaller subarray immediately.
            if (ir - i + 1 >= j - l) {
                istack[jstack] = ir;
                istack[jstack - 1] = i;
                ir = j - 1;
            } else {
                istack[jstack] = j - 1;
                istack[jstack - 1] = l;
                l = i;
            }
        }
    }
}
template <class T> void unique(std::vector<T> &x) {
    if (x.size() <= 1)
        return;
    // First perform a quicksort
    quicksort(x);
    // Next remove duplicate entries
    size_t pos = 1;
    for (size_t i = 1; i < x.size(); i++) {
        if (x[i] != x[pos - 1]) {
            x[pos] = x[i];
            pos++;
        }
    }
    if (pos < x.size())
        x.resize(pos);
}

} // namespace Utilities

#endif
