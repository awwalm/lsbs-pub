import time
import bisect
import tracemalloc



def lcs_hs(s1, s2):
    """
    Hunt-Szymanski algorithm implemented by Ray Gardner (rdg, raygard)

    See Also: 
        `Source code on GitHub <https://github.com/raygard/lcs_diff_demo>`_
    """
    m, n = len(s1), len(s2)

    # Step 1: build linked lists
    matchlist = [[] for k in range(m + 1)]
    # Note line numbers in reverse order
    aa = sorted(zip(s1, range(1, m+1)), key=lambda t: (t[0], -t[1]))
    bb = sorted(zip(s2, range(1, n+1)), key=lambda t: (t[0], -t[1]))
    ai = bi = 0
    while ai < m and bi < n:
        av, bv = aa[ai][0], bb[bi][0]
        if av < bv:
            ai += 1
        elif av > bv:
            bi += 1
        else:
            k = aa[ai][1]
            while bi < n and bb[bi][0] == bv:
                matchlist[k] += [bb[bi][1]]
                bi += 1
            ai += 1
            while ai < m and aa[ai][0] == av:
                matchlist[aa[ai][1]] = matchlist[k]
                ai += 1

    # Step 2: initialize the THRESH array
    thresh = [n+1] * (m + 1)
    thresh[0] = 0

    # Step 3: compute successive THRESH values
    link = [None] * (m+1)
    for i in range(1, m+1):
        for j in matchlist[i]:
            #find k such that thresh[k-1] < j <= thresh[k]
            k = bisect.bisect_left(thresh, j)
            #assert thresh[k-1] < j <= thresh[k]
            if j < thresh[k]:
                thresh[k] = j
                link[k] = (i, j, link[k-1])
                #print(f'dmatch({i}, {j})')
                #assert A[i-1] == B[j-1]

    # Step 4: recover longest common subsequence pairs in reverse order
    k = 0
    while k < m and thresh[k+1] != n + 1:
        k += 1
    p = link[k]
    # v will hold (i,j) pairs
    v = []
    while p != None:
        v.append(p[:2])
        p = p[2]
    v.reverse()

    #print(f'lcslen: {len(v)=}')
    return v, s1, thresh

def get_lcs_hs(a: str, b: str):
    """Helper for obtaining the LCS via Hunt-Szymanski algorithm."""
    
    result = lcs_hs(a, b)
    r0, r1, r2 = result
    for k in range(len(r0)):
        r0[k] = a[r0[k][0]-1]
    return r0
    # return {"lcs": r0, "time": t2, "memory": mem, "meta": [r2]}
