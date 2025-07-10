import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from drawgraph.findStrip5 import find_ranges, tpcod9, tpcod7, scaleArray

R = 30
C = 120
assert C > R
data = pd.read_excel("/Users/huahao/Documents/labAAA research/research/digital fabrication/bamboo-sawdust/sol6.xlsx", header=None, usecols="A:D")
dat = data.to_numpy()
x1 = np.reshape(dat[:, 0], (R, C), 'F')
x2 = np.reshape(dat[:, 1], (C, R))
x3 = np.reshape(dat[:, 2], (R, C), 'F')
x4 = np.reshape(dat[:, 3], (R, C), 'F')
lnwei = 0.2   #1 for show(); 0.2 for PDF
num_fila = [0, 0, 0, 0]
x1lns = []
x2lns = []
x3lns = []
x4lns = []

# draw grid
for row in range(R + 1):
    plt.plot([0, C], [row, row], color='k', lw=lnwei, alpha=0.07)
for col in range(C + 1):
    plt.plot([col, col], [0, R], color='k', lw=lnwei, alpha=0.07)

for row in range(R):  # draw x1
    _arr = [round(v * 9) for v in x1[row, :]]  # float to 0-9
    cy = (R - 1 - row)  # tricky Y
    for k in range(1, 10):  # type
        ints = find_ranges(_arr, tpcod9[k])
        y = cy + (1 - (k - 0.5) / 9)  # trick
        num_fila[0] += len(ints)
        for lin in ints:
            x1lns.append([lin[0], y, lin[1] + 1, y])
            plt.plot([lin[0], lin[1] + 1], [y, y], color=plt.cm.tab20c(k), solid_capstyle="butt", lw=lnwei, alpha=0.7)

for col in range(C):  # draw x2
    _arr = [round(v * 9) for v in x2[col, :]]  # float to 0-9
    for k in range(1, 10):  # type
        ints = find_ranges(_arr, k)
        x = col + (k - 0.5) / 9
        num_fila[1] += len(ints)
        for lin in ints:
            x2lns.append([x, R - lin[0], x, R - (lin[1] + 1)])
            plt.plot([x, x], [R - lin[0], R - (lin[1] + 1)], color=plt.cm.tab20c(k), solid_capstyle="butt", lw=lnwei, alpha=0.7)  # tricky Y

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
        num_fila[2] += len(ints)
        for ln in ints:
            x1, x2, y1, y2 = max(0, i) + ln[0] - 1 + sk, max(0, i) + ln[1] + 1, max(0, -i) + ln[0], max(0, -i) + ln[1] + 2 - sk
            a1, a2, b1, b2 = max(0, i) + ln[0], max(0, i) + ln[1] + sk, max(0, -i) + ln[0] + 1 - sk, max(0, -i) + ln[1] + 1
            if k <= 3:
                st = [x1, y1] if 0.5 < x1 else [a1, b1]     #corner 1
                ed = [x2, y2] if y2 < R - 0.5 else [a2, b2] #corner 2
            else:
                st = [a1, b1] if y1 > 0.5 else [x1, y1]      #corner 3
                ed = [a2, b2] if x2 < C - 0.5 else [x2, y2]  #corner 4
            x3lns.append([*st, *ed])
            plt.plot([st[0], ed[0]], [st[1], ed[1]], color=plt.cm.tab20c(k+1), lw=lnwei, solid_capstyle="butt", alpha=0.7)

        # draw x4
        if 0 >= i:  # head
            _arr = [round(x4[j - i, j] * 9) for j in range(R + i)]
        elif C - R <= i:  # tail
            _arr = [round(x4[j, i + j] * 9) for j in range(C - i)]
        else:  # middle
            _arr = [round(x4[j, i + j] * 9) for j in range(R)]
        ints = find_ranges(_arr, tpcod7[k])
        num_fila[3] += len(ints)
        for ln in ints:
            x1, x2, y1, y2 = max(0, i) + ln[0] - 1 + sk, max(0, i) + ln[1] + 1, min(0, i) + R - ln[0], min(0, i) + R - ln[1] - 2 + sk
            a1, a2, b1, b2 = max(0, i) + ln[0], max(0, i) + ln[1] + sk, min(0, i) + R - ln[0] - 1 + sk, min(0, i) + R - ln[1] - 1
            if k <= 3:
                st = [x1, y1] if 0.5 < x1 else [a1, b1] #corner 1
                ed = [x2, y2] if y2 > 0.5 else [a2, b2] #corner 2
            else:
                st = [a1, b1] if y1 < R - 0.5 else [x1, y1]  #corner 3
                ed = [a2, b2] if x2 < C - 0.5 else [x2, y2]  #corner 4
            x4lns.append([*st, *ed])
            plt.plot([st[0], ed[0]], [st[1], ed[1]], color=plt.cm.tab20c(k+1), lw=lnwei, solid_capstyle="butt", alpha=0.7)

# for ln in x1lns:
#     scaleArray(ln, 40)  # 40*40*14mm
# for ln in x2lns:
#     scaleArray(ln, 40)  # 40*40*14mm
# for ln in x3lns:
#     scaleArray(ln, 40)  # 40*40*14mm
# for ln in x4lns:
#     scaleArray(ln, 40)  # 40*40*14mm
# with pd.ExcelWriter("xe.xlsx") as writer:
#     pd.DataFrame(x1lns).to_excel(writer, sheet_name="x1", header=False, index=False)
#     pd.DataFrame(x2lns).to_excel(writer, sheet_name="x2", header=False, index=False)
#     pd.DataFrame(x3lns).to_excel(writer, sheet_name="x3", header=False, index=False)
#     pd.DataFrame(x4lns).to_excel(writer, sheet_name="x4", header=False, index=False)

print(num_fila)
print(sum(num_fila))
plt.axis('equal')
plt.axis('off')
plt.show()
# plt.savefig("strips6_4.pdf", format="pdf", bbox_inches="tight", pad_inches=0)
