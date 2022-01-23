import matplotlib.pyplot as plt
import sys
total_time=[]
reordering_time = []
base_gapscore=[]
graph_name=["ca", "roadNet", "dblp", "youtube"]
parameter_list = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

# gapscoreとreordering timeを格納
for name in graph_name:
  tmp_time = []
  tmp_gapscore = []
  #for i in range(1, int(sys.argv[1])+1):
  for i in parameter_list:
    filename = str(name) + "-" + str(i) + ".ans"
    with open(filename) as f:
      row = [s.strip() for s in f.readlines()]
      ans=row[0].split()
      tmp_gapscore.append(float(ans[0]))
      tmp_time.append(float(ans[1]))

  total_time.append(tmp_time)
  reordering_time.append(tmp_time)
  base_gapscore.append(tmp_gapscore)


# pagerank演算の時間を格納
exec_list = []
for index, name in enumerate(graph_name):
  tmp_exec_list = []
  for index2, i in enumerate(parameter_list):
  #for i in range(1, int(sys.argv[1])+1):
    filename = str(name) + "-" + str(i) + ".exectime"
    with open(filename) as f:
      tmp = f.readlines()
      exectime = []
      for ptime in tmp:
        exectime.append(float(str(ptime)))
      exectime.sort()
      mean_exectime = exectime[int(len(exectime)*0.2):int(len(exectime)*0.8)]
      tmp_exec_list.append(float(sum(mean_exectime))/len(mean_exectime))
      total_time[index][index2] += float(sum(mean_exectime))/len(mean_exectime)
  exec_list.append(tmp_exec_list)

rw_list = []
for name in graph_name:
  filename=str(name)+".rwtime"
  rwtime_ave = 0
  with open(filename) as f:
    tmp = f.readlines()
    for rwtime in tmp:
      rwtime_ave += float(str(rwtime))
    rwtime_ave = rwtime_ave / len(tmp)
    rw_list.append(rwtime_ave)

gap_time = []
for i,gtime in enumerate(reordering_time):
  tmp_gap_time = []
  for t in gtime:
    tmp_gap_time.append(t-rw_list[i])
  gap_time.append(tmp_gap_time)

gap_time_by_normal = []
for i, gtime in enumerate(gap_time):
  tmp = []
  for t in gtime:
    tmp.append(100*t/gtime[0])
  gap_time_by_normal.append(tmp)

original_exec_list = []
for name in graph_name:
  filename=str(name)+"-original.exectime"
  original_exec_time = 0
  with open(filename) as f:
    tmp = f.readlines()
    for rwtime in tmp:
      original_exec_time += float(str(rwtime))
    original_exec_time = original_exec_time / len(tmp)
    original_exec_list.append(original_exec_time)

random_exec_list = []
for name in graph_name:
  filename=str(name)+"-random.exectime"
  random_exec_time = 0
  with open(filename) as f:
    tmp = f.readlines()
    for rwtime in tmp:
      random_exec_time += float(str(rwtime))
    random_exec_time = random_exec_time / len(tmp)
    random_exec_list.append(random_exec_time)

original_time_list = []
random_time_list = []
for i in range(len(original_exec_list)):
  original_time_list.append(rw_list[i] + original_exec_list[i])
  random_time_list.append(rw_list[i] + random_exec_list[i])

speed_up = []
speed_up_random = []
for i, graph_time in enumerate(total_time):
  speed_up_tmp = []
  speed_up_random_tmp = []
  for time in graph_time:
    speed_up_tmp.append(original_time_list[i]/time)
    speed_up_random_tmp.append(random_time_list[i]/time)
  speed_up.append(speed_up_tmp)
  speed_up_random.append(speed_up_random_tmp)

fig=plt.figure(figsize=(100, 100))
ax1=fig.add_subplot(2, 2, 1)
ax2=fig.add_subplot(2, 2, 2)
ax3=fig.add_subplot(2, 2, 3)
ax4=fig.add_subplot(2, 2, 4)
ax1.plot(parameter_list, base_gapscore[0], label=graph_name[0])
ax2.plot(parameter_list, base_gapscore[1], label=graph_name[1])
ax3.plot(parameter_list, base_gapscore[2], label=graph_name[2])
ax4.plot(parameter_list, base_gapscore[3], label=graph_name[3])
#ax1.plot(parameter_list, base_gapscore[0], label=graph_name[0])
#ax2.plot(parameter_list, base_gapscore[1], label=graph_name[1])
#ax3.plot(parameter_list, base_gapscore[2], label=graph_name[2])
#ax4.plot(parameter_list, base_gapscore[3], label=graph_name[3])
#ax1.plot(parameter_list, base_gapscore[0])
#ax2.plot(parameter_list, base_gapscore[1])
#ax3.plot(parameter_list, base_gapscore[2])
#ax4.plot(parameter_list, base_gapscore[3])
#ax1.set_title(graph_name[0])
#ax2.set_title(graph_name[1])
#ax3.set_title(graph_name[2])
#ax4.set_title(graph_name[3])
ax1.set_xlabel('parameter')
ax2.set_xlabel('parameter')
ax3.set_xlabel('parameter')
ax4.set_xlabel('parameter')
ax1.set_ylabel('Gap Score')
ax2.set_ylabel('Gap Score')
ax3.set_ylabel('Gap Score')
ax4.set_ylabel('Gap Score')
ax1.legend(loc='best')
ax2.legend(loc='best')
ax3.legend(loc='best')
ax4.legend(loc='best')
#ax1.get_legend().remove()
#ax2.get_legend().remove()
#ax3.get_legend().remove()
#ax4.get_legend().remove()
#plt.xlim(1, parameter_list[-1])
#plt.tight_layout()
plt.show()

fig=plt.figure(figsize=(100, 100))
ax1=fig.add_subplot(2, 2, 1)
ax2=fig.add_subplot(2, 2, 2)
ax3=fig.add_subplot(2, 2, 3)
ax4=fig.add_subplot(2, 2, 4)
ax1.plot(parameter_list, reordering_time[0], label=graph_name[0])
ax2.plot(parameter_list, reordering_time[1], label=graph_name[1])
ax3.plot(parameter_list, reordering_time[2], label=graph_name[2])
ax4.plot(parameter_list, reordering_time[3], label=graph_name[3])
ax1.set_xlabel('Parameter')
ax2.set_xlabel('Parameter')
ax3.set_xlabel('Parameter')
ax4.set_xlabel('Parameter')
ax1.set_ylabel('Reordering Time(sec)')
ax2.set_ylabel('Reordering Time(sec)')
ax3.set_ylabel('Reordering Time(sec)')
ax4.set_ylabel('Reordering Time(sec)')
ax1.legend(loc='best')
ax2.legend(loc='best')
ax3.legend(loc='best')
ax4.legend(loc='best')
#ax1.get_legend().remove()
#ax2.get_legend().remove()
#ax3.get_legend().remove()
#ax4.get_legend().remove()
#plt.xlim(1, parameter_list[-1])
#plt.tight_layout()
plt.show()

fig=plt.figure(figsize=(100, 100))
ax1=fig.add_subplot(2, 2, 1)
ax2=fig.add_subplot(2, 2, 2)
ax3=fig.add_subplot(2, 2, 3)
ax4=fig.add_subplot(2, 2, 4)
ax1.plot(parameter_list, gap_time[0], label=graph_name[0])
ax2.plot(parameter_list, gap_time[1], label=graph_name[1])
ax3.plot(parameter_list, gap_time[2], label=graph_name[2])
ax4.plot(parameter_list, gap_time[3], label=graph_name[3])
ax1.set_xlabel('Parameter')
ax2.set_xlabel('Parameter')
ax3.set_xlabel('Parameter')
ax4.set_xlabel('Parameter')
ax1.set_ylabel('Wait Time(sec)')
ax2.set_ylabel('Wait Time(sec)')
ax3.set_ylabel('Wait Time(sec)')
ax4.set_ylabel('Wait Time(sec)')
ax1.legend(loc='best')
ax2.legend(loc='best')
ax3.legend(loc='best')
ax4.legend(loc='best')
#ax1.get_legend().remove()
#ax2.get_legend().remove()
#ax3.get_legend().remove()
#ax4.get_legend().remove()
#plt.xlim(1, parameter_list[-1])
#plt.tight_layout()
plt.show()

fig=plt.figure(figsize=(100, 100))
ax1=fig.add_subplot(2, 2, 1)
ax2=fig.add_subplot(2, 2, 2)
ax3=fig.add_subplot(2, 2, 3)
ax4=fig.add_subplot(2, 2, 4)
ax1.plot(parameter_list, gap_time_by_normal[0], label=graph_name[0])
ax2.plot(parameter_list, gap_time_by_normal[1], label=graph_name[1])
ax3.plot(parameter_list, gap_time_by_normal[2], label=graph_name[2])
ax4.plot(parameter_list, gap_time_by_normal[3], label=graph_name[3])
ax1.set_xlabel('Parameter')
ax2.set_xlabel('Parameter')
ax3.set_xlabel('Parameter')
ax4.set_xlabel('Parameter')
ax1.set_ylabel('Wait Compared To Normal Gorder')
ax2.set_ylabel('Wait Compared To Normal Gorder')
ax3.set_ylabel('Wait Compared To Normal Gorder')
ax4.set_ylabel('Wait Compared To Normal Gorder')
ax1.legend(loc='best')
ax2.legend(loc='best')
ax3.legend(loc='best')
ax4.legend(loc='best')
#ax1.get_legend().remove()
#ax2.get_legend().remove()
#ax3.get_legend().remove()
#ax4.get_legend().remove()
#plt.xlim(1, parameter_list[-1])
#plt.tight_layout()
plt.show()

fig=plt.figure(figsize=(100, 100))
ax1=fig.add_subplot(2, 2, 1)
ax2=fig.add_subplot(2, 2, 2)
ax3=fig.add_subplot(2, 2, 3)
ax4=fig.add_subplot(2, 2, 4)
ax1.plot(parameter_list, total_time[0], label=graph_name[0])
ax2.plot(parameter_list, total_time[1], label=graph_name[1])
ax3.plot(parameter_list, total_time[2], label=graph_name[2])
ax4.plot(parameter_list, total_time[3], label=graph_name[3])
#ax1.plot(parameter_list, total_time[0])
#ax2.plot(parameter_list, total_time[1])
#ax3.plot(parameter_list, total_time[2])
#ax4.plot(parameter_list, total_time[3])
#ax1.set_title(graph_name[0])
#ax2.set_title(graph_name[1])
#ax3.set_title(graph_name[2])
#ax4.set_title(graph_name[3])
ax1.set_xlabel('Parameter')
ax2.set_xlabel('Parameter')
ax3.set_xlabel('Parameter')
ax4.set_xlabel('Parameter')
ax1.set_ylabel('Total Elapsed Time(sec)')
ax2.set_ylabel('Total Elapsed Time(sec)')
ax3.set_ylabel('Total Elapsed Time(sec)')
ax4.set_ylabel('Total Elapsed Time(sec)')
#ax1.get_legend().remove()
#ax2.get_legend().remove()
#ax3.get_legend().remove()
#ax4.get_legend().remove()
ax1.legend(loc='best')
ax2.legend(loc='best')
ax3.legend(loc='best')
ax4.legend(loc='best')
#plt.xlim(1, parameter_list[-1])
#plt.tight_layout()
plt.show()

fig=plt.figure(figsize=(100, 100))
ax1=fig.add_subplot(2, 2, 1)
ax2=fig.add_subplot(2, 2, 2)
ax3=fig.add_subplot(2, 2, 3)
ax4=fig.add_subplot(2, 2, 4)
ax1.plot(parameter_list, exec_list[0], label=graph_name[0])
ax2.plot(parameter_list, exec_list[1], label=graph_name[1])
ax3.plot(parameter_list, exec_list[2], label=graph_name[2])
ax4.plot(parameter_list, exec_list[3], label=graph_name[3])
#ax1.plot(parameter_list, exec_list[0])
#ax2.plot(parameter_list, exec_list[1])
#ax3.plot(parameter_list, exec_list[2])
#ax4.plot(parameter_list, exec_list[3])
#ax1.set_title(graph_name[0])
#ax2.set_title(graph_name[1])
#ax3.set_title(graph_name[2])
#ax4.set_title(graph_name[3])
ax1.set_xlabel('parameter')
ax2.set_xlabel('parameter')
ax3.set_xlabel('parameter')
ax4.set_xlabel('parameter')
ax1.set_ylabel('PageRank Time(sec)')
ax2.set_ylabel('PageRank Time(sec)')
ax3.set_ylabel('PageRank Time(sec)')
ax4.set_ylabel('PageRank Time(sec)')
#ax1.get_legend().remove()
#ax2.get_legend().remove()
#ax3.get_legend().remove()
#ax4.get_legend().remove()
ax1.legend(loc='best')
ax2.legend(loc='best')
ax3.legend(loc='best')
ax4.legend(loc='best')
#plt.xlim(1, parameter_list[-1])
#plt.tight_layout()
plt.show()

fig=plt.figure(figsize=(100, 100))
ax1=fig.add_subplot(2, 2, 1)
ax2=fig.add_subplot(2, 2, 2)
ax3=fig.add_subplot(2, 2, 3)
ax4=fig.add_subplot(2, 2, 4)
ax1.plot(parameter_list, speed_up[0], label=graph_name[0]+" original id")
ax2.plot(parameter_list, speed_up[1], label=graph_name[1]+" original id")
ax3.plot(parameter_list, speed_up[2], label=graph_name[2]+" original id")
ax4.plot(parameter_list, speed_up[3], label=graph_name[3]+" original id")
ax1.plot(parameter_list, speed_up_random[0], label=graph_name[0]+" random id")
ax2.plot(parameter_list, speed_up_random[1], label=graph_name[1]+" random id")
ax3.plot(parameter_list, speed_up_random[2], label=graph_name[2]+" random id")
ax4.plot(parameter_list, speed_up_random[3], label=graph_name[3]+" random id")
#ax1.plot(parameter_list, speed_up[0], label="original id")
#ax2.plot(parameter_list, speed_up[1], label="original id")
#ax3.plot(parameter_list, speed_up[2], label="original id")
#ax4.plot(parameter_list, speed_up[3], label="original id")
#ax1.plot(parameter_list, speed_up_random[0], label="random id")
#ax2.plot(parameter_list, speed_up_random[1], label="random id")
#ax3.plot(parameter_list, speed_up_random[2], label="random id")
#ax4.plot(parameter_list, speed_up_random[3], label="random id")
#ax1.set_title(graph_name[0])
#ax2.set_title(graph_name[1])
#ax3.set_title(graph_name[2])
#ax4.set_title(graph_name[3])
ax1.set_xlabel('parameter')
ax2.set_xlabel('parameter')
ax3.set_xlabel('parameter')
ax4.set_xlabel('parameter')
ax1.set_ylabel('Speed Up')
ax2.set_ylabel('Speed Up')
ax3.set_ylabel('Speed Up')
ax4.set_ylabel('Speed Up')
ax1.legend(loc='best')
ax2.legend(loc='best')
ax3.legend(loc='best')
ax4.legend(loc='best')
#plt.xlim(1, parameter_list[-1])
#plt.tight_layout()
plt.show()

speed_up_by_base = []
for data in speed_up:
  tmp = []
  for data_speed in data:
    tmp.append(data_speed/data[0])
  speed_up_by_base.append(tmp)

speed_up_by_base_random = []
for data in speed_up_random:
  tmp = []
  for data_speed in data:
    tmp.append(data_speed/data[0])
  speed_up_by_base_random.append(tmp)

fig=plt.figure(figsize=(100, 100))
ax1=fig.add_subplot(2, 2, 1)
ax2=fig.add_subplot(2, 2, 2)
ax3=fig.add_subplot(2, 2, 3)
ax4=fig.add_subplot(2, 2, 4)
ax1.plot(parameter_list, speed_up_by_base[0], label=graph_name[0])
ax2.plot(parameter_list, speed_up_by_base[1], label=graph_name[1])
ax3.plot(parameter_list, speed_up_by_base[2], label=graph_name[2])
ax4.plot(parameter_list, speed_up_by_base[3], label=graph_name[3])
#ax1.plot(parameter_list, speed_up_by_base_random[0], label=graph_name[0]+" random id")
#ax2.plot(parameter_list, speed_up_by_base_random[1], label=graph_name[1]+" random id")
#ax3.plot(parameter_list, speed_up_by_base_random[2], label=graph_name[2]+" random id")
#ax4.plot(parameter_list, speed_up_by_base_random[3], label=graph_name[3]+" random id")
#ax1.plot(parameter_list, speed_up[0], label="original id")
#ax2.plot(parameter_list, speed_up[1], label="original id")
#ax3.plot(parameter_list, speed_up[2], label="original id")
#ax4.plot(parameter_list, speed_up[3], label="original id")
#ax1.plot(parameter_list, speed_up_random[0], label="random id")
#ax2.plot(parameter_list, speed_up_random[1], label="random id")
#ax3.plot(parameter_list, speed_up_random[2], label="random id")
#ax4.plot(parameter_list, speed_up_random[3], label="random id")
#ax1.set_title(graph_name[0])
#ax2.set_title(graph_name[1])
#ax3.set_title(graph_name[2])
#ax4.set_title(graph_name[3])
ax1.set_xlabel('Parameter')
ax2.set_xlabel('Parameter')
ax3.set_xlabel('Parameter')
ax4.set_xlabel('Parameter')
ax1.set_ylabel('Speed Up Compared to Normal Gorder')
ax2.set_ylabel('Speed Up Compared to Normal Gorder')
ax3.set_ylabel('Speed Up Compared to Normal Gorder')
ax4.set_ylabel('Speed Up Compared to Normal Gorder')
ax1.legend(loc='best')
ax2.legend(loc='best')
ax3.legend(loc='best')
ax4.legend(loc='best')
#plt.xlim(1, parameter_list[-1])
#plt.tight_layout()
plt.show()
##fig=plt.figure(figsize=(100, 100))
##ax=fig.add_subplot(111)
##ax.plot(parameter_list, total_time, label="gorder")
###ax.plot(parameter_list, sqrt025_time, label="0.25")
###ax.plot(parameter_list, sqrt050_time, label="0.50")
###ax.plot(parameter_list, sqrt075_time, label="0.75")
##ax.set_xlabel('parameter_list')
##ax.set_ylabel('Total Elapsed Time (sec)')
##plt.legend(loc='best')
##plt.xlim(1, parameter_list[-1])
##plt.show()
##
##fig=plt.figure(figsize=(100, 100))
##ax=fig.add_subplot(111)
##ax.plot(parameter_list, exec_list, label="gorder")
###ax.plot(parameter_list, sqrt025_time, label="0.25")
###ax.plot(parameter_list, sqrt050_time, label="0.50")
###ax.plot(parameter_list, sqrt075_time, label="0.75")
##ax.set_xlabel('parameter_list')
##ax.set_ylabel('PageRank Time (sec)')
##plt.legend(loc='best')
##plt.xlim(1, parameter_list[-1])
##plt.show()
##
##fig=plt.figure(figsize=(100, 100))
##ax=fig.add_subplot(111)
##ax.plot(parameter_list, speed_up, label="original id")
##ax.plot(parameter_list, speed_up_random, label="random id")
###ax.plot(parameter_list, sqrt025_time, label="0.25")
###ax.plot(parameter_list, sqrt050_time, label="0.50")
###ax.plot(parameter_list, sqrt075_time, label="0.75")
##ax.set_xlabel('parameter_list')
##ax.set_ylabel('Speed Up')
##plt.legend(loc='best')
##plt.xlim(1, parameter_list[-1])
##plt.show()
