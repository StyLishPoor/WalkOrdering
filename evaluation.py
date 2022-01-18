import matplotlib.pyplot as plt
import sys
parallel=list(range(1, int(sys.argv[1])+1))
total_time=[]
base_gapscore=[]
sqrt025_time=[]
sqrt025_gapscore=[]
sqrt050_time=[]
sqrt050_gapscore=[]
sqrt075_time=[]
sqrt075_gapscore=[]

# gapscoreとreordering timeを格納
for i in range(1, int(sys.argv[1])+1):
  filename = str(i) + "-.ans"
  with open(filename) as f:
    row = [s.strip() for s in f.readlines()]
    ans=row[0].split()
    base_gapscore.append(float(ans[0]))
    total_time.append(float(ans[1]))

# pagerank演算の時間を格納
exec_list = []
for i in range(1, int(sys.argv[1])+1):
  filename = str(i) + ".exectime"
  with open(filename) as f:
    tmp = f.readlines()
    exectime = 0
    for ptime in tmp:
      exectime += float(str(ptime))
    exectime = exectime / len(tmp)
    exec_list.append(exectime)
    total_time[i-1] += exectime

filename="rw-time.txt"
rwtime_ave = 0
with open(filename) as f:
  tmp = f.readlines()
  for rwtime in tmp:
    rwtime_ave += float(str(rwtime))
  rwtime_ave = rwtime_ave / len(tmp)

filename="original.exectime"
original_exec_time = 0
with open(filename) as f:
  tmp = f.readlines()
  for rwtime in tmp:
    original_exec_time += float(str(rwtime))
  original_exec_time = original_exec_time / len(tmp)

filename="random.exectime"
random_exec_time = 0
with open(filename) as f:
  tmp = f.readlines()
  for rwtime in tmp:
    random_exec_time += float(str(rwtime))
  random_exec_time = random_exec_time / len(tmp)


original_time = rwtime_ave + original_exec_time
random_time = rwtime_ave + random_exec_time

#with open(filename) as f:
#  tmp = f.readline()
#  original_time += float(str(tmp))
#filename = "original.exectime"
#with open(filename) as f:
  #tmp = f.readline()
  #original_time += float(str(tmp))
#for i in range(1, int(sys.argv[1])+1):
#  filename = str(i) + "-025.ans"
#  with open(filename) as f:
#    row = [s.strip() for s in f.readlines()]
#    ans=row[0].split()
#    sqrt025_gapscore.append(float(ans[0]))
#    sqrt025_time.append(float(ans[1]))
#
#for i in range(1, int(sys.argv[1])+1):
#  filename = str(i) + "-050.ans"
#  with open(filename) as f:
#    row = [s.strip() for s in f.readlines()]
#    ans=row[0].split()
#    sqrt050_gapscore.append(float(ans[0]))
#    sqrt050_time.append(float(ans[1]))
#
#for i in range(1, int(sys.argv[1])+1):
#  filename = str(i) + "-075.ans"
#  with open(filename) as f:
#    row = [s.strip() for s in f.readlines()]
#    ans=row[0].split()
#    sqrt075_gapscore.append(float(ans[0]))
#    sqrt075_time.append(float(ans[1]))

speed_up = []
speed_up_random = []
for time in total_time:
  speed_up.append(original_time/time)
  speed_up_random.append(random_time/time)
print("-----")
print(base_gapscore)
print(total_time)
print("-----")

fig=plt.figure(figsize=(100, 100))
ax=fig.add_subplot(111)
ax.plot(parallel, base_gapscore, label="gorder")
#ax.plot(parallel, sqrt025_gapscore, label="0.25")
#ax.plot(parallel, sqrt050_gapscore, label="0.50")
#ax.plot(parallel, sqrt075_gapscore, label="0.75")
ax.set_xlabel('Parallel')
ax.set_ylabel('Gap Score')
plt.legend(loc='best')
plt.xlim(1, parallel[-1])
plt.show()

fig=plt.figure(figsize=(100, 100))
ax=fig.add_subplot(111)
ax.plot(parallel, total_time, label="gorder")
#ax.plot(parallel, sqrt025_time, label="0.25")
#ax.plot(parallel, sqrt050_time, label="0.50")
#ax.plot(parallel, sqrt075_time, label="0.75")
ax.set_xlabel('Parallel')
ax.set_ylabel('Total Elapsed Time (sec)')
plt.legend(loc='best')
plt.xlim(1, parallel[-1])
plt.show()

fig=plt.figure(figsize=(100, 100))
ax=fig.add_subplot(111)
ax.plot(parallel, exec_list, label="gorder")
#ax.plot(parallel, sqrt025_time, label="0.25")
#ax.plot(parallel, sqrt050_time, label="0.50")
#ax.plot(parallel, sqrt075_time, label="0.75")
ax.set_xlabel('Parallel')
ax.set_ylabel('PageRank Time (sec)')
plt.legend(loc='best')
plt.xlim(1, parallel[-1])
plt.show()

fig=plt.figure(figsize=(100, 100))
ax=fig.add_subplot(111)
ax.plot(parallel, speed_up, label="original id")
ax.plot(parallel, speed_up_random, label="random id")
#ax.plot(parallel, sqrt025_time, label="0.25")
#ax.plot(parallel, sqrt050_time, label="0.50")
#ax.plot(parallel, sqrt075_time, label="0.75")
ax.set_xlabel('Parallel')
ax.set_ylabel('Speed Up')
plt.legend(loc='best')
plt.xlim(1, parallel[-1])
plt.show()
