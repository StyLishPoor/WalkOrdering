import matplotlib.pyplot as plt
a=[2.061, 1.633, 1.535, 1.4856, 1.48301, 1.49243, 1.51853, 1.538, 1.559, 1.595, 1.76489, 1.96141, 2.1626, 2.3562] 
b=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
c=[0.35, 0.356, 0.360, 0.358, 0.355, 0.350, 0.353, 0.352, 0.351, 0.355, 0.348, 0.345, 0.348, 0.348]
d=c*10
e=[]
f=[]

base=0.407+1.23
for i in range(len(b)):
  tmp = a[i] + c[i]
  e.append(base/tmp)

base=a[0]+c[0]
for i in range(len(b)):
  tmp = a[i] + c[i]
  f.append(base/tmp)

fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(b,e,label='RW + Reordering + PageRank')
#ax.plot(x_ratio,rw,label='MHRW')
ax.set_xlabel('Fractoin Rate (%)')
ax.set_ylabel('Speed Up Compared to Original')
plt.legend(loc='best')
plt.xlim(1, b[-1])
plt.show()

fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(b,f,label='RW + Reordering + PageRank')
#ax.plot(x_ratio,rw,label='MHRW')
ax.set_xlabel('Fractoin Rate (%)')
ax.set_ylabel('Speed Up Compared to Normal')
plt.legend(loc='best')
plt.xlim(1, b[-1])
plt.show()
