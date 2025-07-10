import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def scaleArray(arr, s):
    for i in range(len(arr)):
        arr[i] *=s

def findStart(arr, cod):
    re = []
    for i in range(1, len(arr)):
        if arr[i - 1] < cod and arr[i] >= cod:
            re.append(i)
    return re


def findEnd(arr, cod):
    re = []
    for i in range(1, len(arr)):
        if arr[i - 1] >= cod and arr[i] < cod:
            re.append(i - 1)
    return re


def find_ranges(arr, cod):
    # cod = tpcod[type]
    leng = len(arr)
    if cod <= min(arr):
        return [[0, leng - 1]]  # all
    elif cod > max(arr):
        return []  # none

    _s = findStart(arr, cod)  # start position ids
    _e = findEnd(arr, cod)  # end position ids
    # print(_s)
    # print(_e)
    if not _s:  # no start, must has one end
        assert 1 == len(_e)
        return [[0, _e[0]]]  # add 0
    if not _e:  # no end, must has one start
        assert 1 == len(_s)
        return [[_s[0], leng - 1]]  # add leng-1

    if _e[0] < _s[0]:
        _s.insert(0, 0)
    if _e[-1] < _s[-1]:
        _e.append(leng - 1)

    assert len(_s) == len(_e)
    intervals = []
    for i, si in enumerate(_s):
        ts = _s[i]
        te = _e[i]
        assert ts <= te
        intervals.append([ts, te])
    return intervals

tpcod9 = (None, 3, 7, 4, 9, 1, 6, 5, 8, 2)  # grade 1-9 for x1 x2
tpcod7 = (None, 3, 7, 4, 1, 6, 5, 2)  # grade 1-7 for x3 x4

if __name__ == "__main__":
    R = 30
    C = 120
    assert C > R
    data = pd.read_excel("/Users/huahao/Documents/labAAA research/research/digital fabrication/bamboo-sawdust/sol5.xlsx", header=None, usecols="A,B")
    dat = data.to_numpy()
    x1 = np.reshape(dat[:, 0], (R, C), 'F')
    x2 = np.reshape(dat[:, 1], (C, R))

    shifY = 40
    lnwei = 1  #1 for show(); 0.2 for PDF
    # draw grid
    for row in range(R+1):
        plt.plot([0, C], [row, row], color='k', lw=lnwei, alpha=0.07)
        plt.plot([0, C], [row + shifY, row + shifY], color='k', lw=lnwei, alpha=0.07)
        plt.plot([0, C], [row + 2 * shifY, row + 2 * shifY], color='k', lw=lnwei, alpha=0.07)
    for col in range(C+1):
        plt.plot([col, col], [0, R], color='k', lw=lnwei, alpha=0.07)
        plt.plot([col, col], [0 + shifY, R + shifY], color='k', lw=lnwei, alpha=0.07)
        plt.plot([col, col], [0 + 2 * shifY, R + 2 * shifY], color='k', lw=lnwei, alpha=0.07)

    for row in range(R):  # draw x1
        _arr = [round(v * 9) for v in x1[row, :]]  # float to 0-9
        cy = (R - 1 - row)  # tricky Y
        for k in range(1, 10):  # type
            ints = find_ranges(_arr, tpcod9[k])
            y = cy + (1 - (k - 0.5) / 9)  # trick
            for lin in ints:
                plt.plot([lin[0], lin[1] + 1], [y + 2 * shifY, y + 2 * shifY], color=plt.cm.tab20c(k), lw=lnwei, solid_capstyle="butt", alpha=0.7)
                plt.plot([lin[0], lin[1] + 1], [y + shifY, y + shifY], color=plt.cm.tab20c(k), solid_capstyle="butt", lw=lnwei)

    for col in range(C):  # draw x2
        _arr = [round(v * 9) for v in x2[col, :]]  # float to 0-9
        for k in range(1, 10):  # type
            ints = find_ranges(_arr, tpcod9[k])
            x = col + (k - 0.5) / 9
            for lin in ints:
                plt.plot([x, x], [R - lin[0], R - (lin[1] + 1)], color=plt.cm.tab20c(k), solid_capstyle="butt", lw=lnwei)  # tricky Y
                plt.plot([x, x], [2 * shifY + R - lin[0], 2 * shifY + 30 - (lin[1] + 1)], color=plt.cm.tab20c(k), solid_capstyle="butt", lw=lnwei, alpha=0.7)  # tricky Y

    plt.axis('equal')
    plt.axis('off')
    plt.show()
    # plt.savefig("strips5.pdf", format="pdf", bbox_inches="tight", pad_inches=0)
