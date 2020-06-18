#include <sort.h>

void merge(void **arr, int (*comparator)(const void *, const void *), int l, int m, int r) {
  int i = 0, j = 0, k = l;
  int n1 = m - l + 1;
  int n2 = r - m;

  void *L[n1], *R[n2];

  // Copy to subarrays
  for (i = 0; i < n1; i++)
    L[i] = arr[l + i];
  for (j = 0; j < n2; j++)
    R[j] = arr[m + 1 + j];

  while (i < n1 && j < n2) {
    if ((*comparator)(L[i], R[i]) <= 0) {
      arr[k] = L[i];
      i++;
    } else {
      arr[k] = R[j];
      j++;
    }
    k++;
  }

  while (i < n1) {
    arr[k] = L[i];
    i++;
    k++;
  }

  while (j < n2) {
    arr[k] = R[j];
    j++;
    k++;
  }
}

void merge_sort(void **arr, int (*comparator)(const void *, const void *),
                int l, int r) {
  if (l < r) {
    int m = l + (r - l) / 2;
    // Sort first and second halves
    merge_sort(arr, comparator, l, m);
    merge_sort(arr, comparator, m + 1, r);
    merge(arr, comparator, l, m, r);
  }
}