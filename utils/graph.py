import matplotlib.pyplot as plt
import sys

parallel = [0, 1]
y1=[2,1]
y2=[2,1]
y3=[2,1]
y4=[2,1]

fig = plt.figure()
ax1=fig.add_subplot(2, 2, 1)
ax2=fig.add_subplot(2, 2, 2)
ax3=fig.add_subplot(2, 2, 3)
ax4=fig.add_subplot(2, 2, 4)
ax1.plot(parallel, y1)
ax2.plot(parallel, y2)
ax3.plot(parallel, y3)
ax4.plot(parallel, y4)
ax1.set_xlabel('Parallel')
ax2.set_xlabel('Parallel')
ax3.set_xlabel('Parallel')
ax4.set_xlabel('Parallel')
ax1.set_ylabel('Gap Score')
ax2.set_ylabel('Gap Score')
ax3.set_ylabel('Gap Score')
ax4.set_ylabel('Gap Score')
plt.tight_layout()
plt.show()
