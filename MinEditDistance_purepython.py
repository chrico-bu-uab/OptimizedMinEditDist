import re


def amean(*Z):
    return sum(Z) / len(Z)


def gmean(*Z):
    result = 1
    for z in Z:
        result *= z
    return result**(1.0 / len(Z))


def hmean(*Z):
    if min(Z) > 0:
        return len(Z) / sum(1 / z for z in Z)
    return 0


def MinEditDistance(x, y, C, editmatrices=None, normalize=False, isDamerau=False):

    if editmatrices is not None:
        for edit, matrix in editmatrices.items():
            if edit in ('sub', 'rev'):
                assert matrix.shape == (26, 26), matrix.shape
                for i in range(26):
                    c = chr(i + 97)
                    assert matrix[c][c] == 0, (c, matrix[c][c])
            else:
                assert matrix.shape == (27, 26), matrix.shape

    x = '#' + re.sub('[^a-z]', '', x.lower())
    y = '#' + re.sub('[^a-z]', '', y.lower())
    m = len(x)
    n = len(y)

    # initialize edit distance matrix
    D = [[0] * n for _ in range(m)]
    for i in range(m):
        D[i][0] = i
    for i in range(n):
        D[0][i] = i

    # initialize backtrace matrix
    backtrace = [[[0, 0, 0, 0] for _ in range(n - 1)] for _ in range(m - 1)]

    # compute recurrence relations
    for i in range(1, m):
        for j in range(1, n):
            swap = None
            if editmatrices is None:
                down = D[i - 1][j] + C[0]
                left = D[i][j - 1] + C[1]
                diag = D[i - 1][j - 1] + C[2] * (x[i] != y[j])
                if isDamerau and (i > 1) and (j > 1) and (
                        x[i - 1] == y[j]) and (x[i] == y[j - 1]):
                    swap = D[i - 2][j - 2] + C[3]
            else:
                # compute distances using lookup tables
                down = D[i - 1][j] + C[0] * editmatrices['del'][x[i]][y[j]]
                left = D[i][j - 1] + C[1] * editmatrices['add'][x[i]][y[j]]
                diag = D[i - 1][j - 1] + C[2] * editmatrices['sub'][x[i]][y[j]]
                if isDamerau and (i > 1) and (j > 1) and (
                        x[i - 1] == y[j]) and (x[i] == y[j - 1]):
                    swap = D[i - 2][j -
                                    2] + C[3] * editmatrices['rev'][x[i]][y[j]]
            # assign minimum distance to current cell
            if swap is None:
                D[i][j] = min(left, down, diag)
            else:
                D[i][j] = min(left, down, diag, swap)
            # backtrace from current cell to previous cells
            backtrace[i - 1][j - 1][0] = D[i][j] == left
            backtrace[i - 1][j - 1][1] = D[i][j] == down
            backtrace[i - 1][j - 1][2] = D[i][j] == diag
            backtrace[i - 1][j - 1][3] = D[i][j] == swap

    if normalize:
        # initialize harmonic mean matrix
        H = [[0] * n for _ in range(m)]
        for i in range(m):
            H[i][0] = i
        for i in range(n):
            H[0][i] = i
        for i in range(1, m):
            H[i][1] = hmean(C[0] + H[i - 1][1], C[1] + H[i][0],
                            C[2] + H[i - 1][0])
        for i in range(1, n):
            H[1][i] = hmean(C[0] + H[1][i - 1], C[1] + H[0][i],
                            C[2] + H[0][i - 1])
        for i in range(2, m):
            for j in range(2, n):
                H[i][j] = hmean(C[0] + H[i - 1][j], C[1] + H[i][j - 1],
                                C[2] + H[i - 1][j - 1], C[3] + H[i - 2][j - 2])
        for i in range(m):
            for j in range(n):
                D[i][j] /= H[m - 1][n - 1] * hmean(*C)

    backtrace = [[('-' if x[0] else ' ') + ('|' if x[1] else ' ') +
                  ('/' if x[2] else ' ') + ('\\' if x[3] else ' ')
                  for x in y] for y in backtrace]

    # return value is the final minimum edit distnce D(m,n)
    out = D[m - 1][n - 1]

    return out, D, backtrace

