import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from drawgraph.findStrip5 import find_ranges, tpcod7

R = 30
C = 120
assert C > R
data = pd.read_excel("/Users/huahao/Documents/labAAA research/research/digital fabrication/bamboo-sawdust/sol5c.xlsx", header=None, usecols="A,B")
dat = data.to_numpy()
x3 = np.reshape(dat[:, 0], (R, C), 'F')
x4 = np.reshape(dat[:, 1], (R, C), 'F')
lnwei = 1 #  1 for show(); 0.2 for PDF

# draw grid
for row in range(R + 1):
    plt.plot([0, C], [row, row], color='k', lw=lnwei, alpha=0.07)
for col in range(C + 1):
    plt.plot([col, col], [0, R], color='k', lw=lnwei, alpha=0.07)

for k in range(1, 8):  # type  1-7
    sk = (k+3)/7  #(k + 4) / 9
    for i in range(1 - R, C):
        # draw x3
        if 0 >= i:  # head
            _arr = [round(x3[R - 1 + i - j, j] * 9) for j in range(R + i)]
        elif C - R <= i:  # tail
            _arr = [round(x3[R - 1 - j, i + j] * 9) for j in range(C - i)]
        else:  # middle
            _arr = [round(x3[R - 1 - j, i + j] * 9) for j in range(R)]
        ints = find_ranges(_arr, tpcod7[k])
        for ln in ints:
            x1, x2, y1, y2 = max(0, i) + ln[0] - 1 + sk, max(0, i) + ln[1] + 1, max(0, -i) + ln[0], max(0, -i) + ln[1] + 2 - sk
            a1, a2, b1, b2 = max(0, i) + ln[0], max(0, i) + ln[1] + sk, max(0, -i) + ln[0] + 1 - sk, max(0, -i) + ln[1] + 1
            if k <= 3:
                st = [x1, y1] if 0.5 < x1 else [a1, b1]     #corner 1
                ed = [x2, y2] if y2 < R - 0.5 else [a2, b2] #corner 2
            else:
                st = [a1, b1] if y1 > 0.5 else [x1, y1]      #corner 3
                ed = [a2, b2] if x2 < C - 0.5 else [x2, y2]  #corner 4
            plt.plot([st[0], ed[0]], [st[1], ed[1]], color=plt.cm.tab20c(k+1), lw=lnwei, solid_capstyle="butt", alpha=0.7)

        # draw x4
        if 0 >= i:  # head
            _arr = [round(x4[j - i, j] * 9) for j in range(R + i)]
        elif C - R <= i:  # tail
            _arr = [round(x4[j, i + j] * 9) for j in range(C - i)]
        else:  # middle
            _arr = [round(x4[j, i + j] * 9) for j in range(R)]
        ints = find_ranges(_arr, tpcod7[k])
        for ln in ints:
            x1, x2, y1, y2 = max(0, i) + ln[0] - 1 + sk, max(0, i) + ln[1] + 1, min(0, i) + R - ln[0], min(0, i) + R - ln[1] - 2 + sk
            a1, a2, b1, b2 = max(0, i) + ln[0], max(0, i) + ln[1] + sk, min(0, i) + R - ln[0] - 1 + sk, min(0, i) + R - ln[1] - 1
            if k <= 3:
                st = [x1, y1] if 0.5 < x1 else [a1, b1] #corner 1
                ed = [x2, y2] if y2 > 0.5 else [a2, b2] #corner 2
            else:
                st = [a1, b1] if y1 < R - 0.5 else [x1, y1]  #corner 3
                ed = [a2, b2] if x2 < C - 0.5 else [x2, y2]  #corner 4
            plt.plot([st[0], ed[0]], [st[1], ed[1]], color=plt.cm.tab20c(k+1), lw=lnwei, solid_capstyle="butt", alpha=0.7)

plt.axis('equal')
plt.axis('off')
plt.show()
# plt.savefig("strips5c.pdf", format="pdf", bbox_inches="tight", pad_inches=0)
