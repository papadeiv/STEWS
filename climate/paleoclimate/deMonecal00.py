from pandas import read_csv
from matplotlib import pyplot as plt

# Render LaTeX to Matplotlib (https://matplotlib.org/stable/gallery/text_labels_and_annotations/tex_demo.html)
# needs extra packages (https://stackoverflow.com/questions/11354149/python-unable-to-render-tex-in-matplotlib)
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    })

timeseries = read_csv('SSTODP658C.csv', sep = ',', engine = 'python', header=0, index_col=0)
print(type(timeseries))
print(timeseries.head)

fig = plt.figure(figsize=[12.8,9.6], dpi=200, layout='tight')
plt.gca().invert_xaxis()
plt.plot(timeseries, linewidth = 4, color='green')
plt.xlabel(r'yr', fontsize = 40)
plt.ylabel(r'SST', fontsize = 40)
plt.xticks(fontsize = 35, weight = 'bold')
plt.yticks(fontsize = 35, weight = 'bold')
plt.title('desertification of North Africa', fontsize = 40)
plt.savefig("../../results/climate/paleoclimate/deMonecal00.png", transparent=True, dpi=600)
