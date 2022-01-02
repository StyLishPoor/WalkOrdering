import matplotlib.pyplot as plt
sequential_gorder_time=[2.061, 1.619, 1.511, 1.448, 1.45491, 1.45442,  1.44389, 1.44332, 1.45356, 1.44026, 1.4342, 1.46291, 1.48103, 1.49167, 1.48389, 1.49562, 1.49573, 1.512, 1.53208]
sequential_gorder_gap=[3.767, 4.0397, 4.2115, 4.32238, 4.42812, 4.48082, 4.58728, 4.60394, 4.62421, 4.68753, 4.72806, 4.75331, 4.76817, 4.8206, 4.837103, 4.83794, 4.8461, 4.92727, 4.91087]
para_gorder_time=[2.061, 1.17432, 0.88217, 0.749043, 0.644664, 0.596755, 0.605193, 0.56183, 0.550571, 0.537648, 0.54139, 0.523439, 0.52322, 0.527824, 0.513943, 0.507879, 0.498702, 0.507543, 0.50051]
para_gorder_gap=[3.767, 4.011, 4.1866, 4.3691, 4.40702, 4.48972, 4.5758, 4.63483, 4.70423, 4.70387, 4.70938, 4.76297, 4.79415, 4.8324, 4.86376, 4.87744, 4.87352, 4.93414, 4.96465]
b=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
test_para_gorder_gap=[3.767, 4.011, 4.1866, 4.3691, 4.40702, 4.48972, 4.5758, 4.63483, 4.70423, 4.70387, 4.70938, 4.76297, 4.79415, 4.8324, 4.86376, 4.87744, 4.87352, 4.93414, 4.96465, 7.03576, 6.44966, 5.63594]
b_test=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
min_sequential_gorder_time=[2.061, 1.619, 1.511, 1.448, 1.45491, 1.45442,  1.44389, 1.44332, 1.45356, 1.44026,  1.48389, 1.53208]
min_para_gorder_time=[2.061, 1.17432, 0.88217, 0.749043, 0.644664, 0.605193, 0.596755, 0.56183, 0.550571, 0.537648, 0.513943, 0.50051]
bb=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 19]
c=[0.35, 0.356, 0.360, 0.358, 0.355, 0.350, 0.353, 0.352, 0.351, 0.355, 0.348, 0.345, 0.348, 0.348]
d=c*10
e=[]
f=[]

base=0.407+1.23
for i in range(len(bb)):
  tmp1 = min_sequential_gorder_time[i] + c[i]
  tmp2 = min_para_gorder_time[i] + c[i]
  e.append(base/tmp1)
  f.append(base/tmp2)


fig=plt.figure(figsize=(100, 100))
ax=fig.add_subplot(111)
#ax.bar(b_test,sequential_gorder_gap,label='Naive')
ax.bar(b_test,test_para_gorder_gap,label='Better')
#ax.plot(x_ratio,rw,label='MHRW')
ax.set_xlabel('Thread Num')
ax.set_ylabel('Gap Score')
plt.legend(loc='best')
plt.xlim(1, b_test[-1])
plt.show()

#fig=plt.figure()
#ax=fig.add_subplot(111)
#ax.plot(b,sequential_gorder_time,label='Naive')
#ax.plot(b,para_gorder_time,label='Better')
#ax.set_xlabel('Thread Num')
#ax.set_ylabel('RW + Reorder Time (sec)')
#plt.legend(loc='best')
#plt.xlim(1, b[-1])
#plt.show()

fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(bb,e,label='Naive')
ax.plot(bb,f,label='Better')
ax.set_xlabel('Thread Num')
ax.set_ylabel('Total Speed Up')
plt.legend(loc='best')
plt.xlim(1, b[-1])
plt.show()
