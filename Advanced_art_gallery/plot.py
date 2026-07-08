import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
import sys

def visualize():
    try:
        with open("input.txt", "r") as f_in, open("output.txt", "r") as f_out:
            lines_in = f_in.read().split()
            lines_out = f_out.read().split()
    except FileNotFoundError:
        print("Ensure input.txt and output.txt exist in the directory.")
        return

    ptr_in = 0
    ptr_out = 0
    
    if len(lines_in) == 0: return
    t = int(lines_in[ptr_in])
    ptr_in += 1

    for tc in range(t):
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_aspect('equal')
        ax.set_title(f"Test Case {tc+1} - Triangulation & Guards")

        # Read Outer Boundary
        v0 = int(lines_in[ptr_in])
        ptr_in += 1
        outer_x, outer_y = [], []
        for _ in range(v0):
            outer_x.append(float(lines_in[ptr_in]))
            outer_y.append(float(lines_in[ptr_in+1]))
            ptr_in += 2
        outer_x.append(outer_x[0]) 
        outer_y.append(outer_y[0])
        ax.plot(outer_x, outer_y, 'k-', linewidth=2, label="Outer Boundary")

        # Read Holes
        h = int(lines_in[ptr_in])
        ptr_in += 1
        for _ in range(h):
            vh = int(lines_in[ptr_in])
            ptr_in += 1
            hx, hy = [], []
            for _ in range(vh):
                hx.append(float(lines_in[ptr_in]))
                hy.append(float(lines_in[ptr_in+1]))
                ptr_in += 2
            hx.append(hx[0])
            hy.append(hy[0])
            ax.plot(hx, hy, 'r-', linewidth=2, label="Hole Boundary" if _ == 0 else "")

        # Read Triangles from Output
        while ptr_out < len(lines_out) and lines_out[ptr_out] != "GUARDS":
            tx = [float(lines_out[ptr_out]), float(lines_out[ptr_out+2]), float(lines_out[ptr_out+4]), float(lines_out[ptr_out])]
            ty = [float(lines_out[ptr_out+1]), float(lines_out[ptr_out+3]), float(lines_out[ptr_out+5]), float(lines_out[ptr_out+1])]
            ax.plot(tx, ty, 'b-', alpha=0.5, linewidth=0.8)
            ptr_out += 6
            
        ptr_out += 1 # skip "GUARDS"
        
        # Read Point Guards from Output
        k = int(lines_out[ptr_out])
        ptr_out += 1
        gx, gy = [], []
        for _ in range(k):
            gx.append(float(lines_out[ptr_out]))
            gy.append(float(lines_out[ptr_out+1]))
            ptr_out += 2
            
        ax.scatter(gx, gy, c='green', marker='o', s=80, zorder=5, label=f"Guards (k={k})")
        
        if lines_out[ptr_out] == "END_TC":
            ptr_out += 1

        ax.legend(loc='upper right')
        plt.savefig(f"plot_{tc+1}.png")
        plt.close()

if __name__ == "__main__":
    visualize()