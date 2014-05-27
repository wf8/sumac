
import numpy as np
import matplotlib.pyplot as plt

# some test data
supermatrix = np.array([
                       [100, 90, 60, 0, 0],
                       [95, 0, 50, 0, 40],
                       [90, 40, 0, 20, 30],
                       [100, 90, 60, 0, 0],
                       [95, 0, 50, 0, 40],
                       [95, 0, 50, 0, 40],
                       [100, 90, 60, 0, 0],
                       [95, 0, 50, 0, 40],
                       [90, 40, 0, 20, 30],
                       [100, 90, 60, 0, 0],
                       [100, 90, 60, 0, 0],
                       [95, 0, 50, 0, 40],
                       [100, 90, 60, 0, 0],
                       [95, 0, 50, 0, 40],
                       [90, 40, 0, 20, 30],
                       [100, 90, 60, 0, 0],
                       [95, 0, 50, 0, 40],
                       [90, 40, 0, 20, 30],
                       [90, 40, 0, 20, 30],
                       [100, 90, 60, 0, 0],
                       [95, 0, 50, 0, 40],
                       [90, 40, 0, 20, 30],
                       [95, 0, 50, 0, 40],
                       [100, 90, 60, 0, 0],
                       [95, 0, 50, 0, 40],
                       [90, 40, 0, 20, 30],
                       [100, 90, 60, 0, 0],
                       [95, 0, 50, 0, 40],
                       [90, 40, 0, 20, 30],
                       [90, 40, 0, 20, 30],
                       [100, 90, 60, 0, 0],
                       [95, 0, 50, 0, 40],
                       [90, 40, 0, 20, 30],
                       [100, 90, 60, 0, 0],
                       [95, 0, 50, 0, 40],
                       [90, 40, 0, 20, 30]])

# set up figure
fig, ax = plt.subplots()
ax.set_position([0.4, 0.1, .3, .7])
# fig.set_size_inches(8.5, 11)

# add data
heatmap = ax.pcolor(supermatrix, cmap=plt.cm.Blues, alpha=0.8)


# put the labels in the middle of each cell
ax.set_yticks(np.arange(supermatrix.shape[0]) + 0.5, minor=False)
ax.set_xticks(np.arange(supermatrix.shape[1]) + 0.5, minor=False)

# move x axis label stuff to top
ax.invert_yaxis()
ax.xaxis.tick_top()

# add the labels
genes = ['1','2','3','4','5'] 
otus = ['R semonivii','Chylismia brevipes', 'Hesperolinon somethiun',
'R semonivii','Chylismia brevipes', 'Hesperolinon somethiun',
'R semonivii','Chylismia brevipes', 'Hesperolinon somethiun',
'R semonivii','Chylismia brevipes', 'Hesperolinon somethiun',
'R semonivii','Chylismia brevipes', 'Hesperolinon somethiun',
'R semonivii','Chylismia brevipes', 'Hesperolinon somethiun',
'R semonivii','Chylismia brevipes', 'Hesperolinon somethiun',
'R semonivii','Chylismia brevipes', 'Hesperolinon somethiun',
'R semonivii','Chylismia brevipes', 'Hesperolinon somethiun',
'R semonivii','Chylismia brevipes', 'Hesperolinon somethiun',
'R semonivii','Chylismia brevipes', 'Hesperolinon somethiun',
'R semonivii','Chylismia brevipes', 'Hesperolinon somethiun',
'R semonivii','Chylismia brevipes', 'Hesperolinon somethiun',
'R semonivii','Chylismia brevipes', 'Hesperolinon somethiun',
'R semonivii','Chylismia brevipes', 'Hesperolinon somethiun',
'R semonivii','Chylismia brevipes', 'Hesperolinon somethiun'
]

ax.set_xticklabels(genes, minor=False, family="Arial", size=10)
ax.set_yticklabels(otus, minor=False, family="Arial", size=8)

# rotate the gene names
#plt.xticks(rotation=90)

# remove junk from axes
ax.grid(False)
ax.set_frame_on(False)

# turn off all the ticks
for t in ax.xaxis.get_major_ticks():
    t.tick1On = False
    t.tick2On = False
for t in ax.yaxis.get_major_ticks():
    t.tick1On = False
    t.tick2On = False

# add color legend
axcolor = fig.add_axes([0.425, .90, 0.25, 0.015])
cbar = fig.colorbar(heatmap, cax=axcolor, ticks=[0,50, 100], orientation='horizontal')
cbar.ax.set_xticklabels(['0%', '50%', '100%'], family="Arial", size=10)


plt.savefig("plot.pdf")
