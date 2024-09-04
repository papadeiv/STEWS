from pandas import read_csv
from matplotlib import pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    })

timeseries = read_csv('YeastPopulation.csv', sep = ',', engine = 'python', header=0, index_col=0)
print(type(timeseries))
print(timeseries.head)

fig = plt.figure(figsize=[12.8,9.6], dpi=200, layout='tight')
plt.plot(timeseries, linewidth = 4, color = 'orange')
plt.yscale("log")
plt.xlabel(r'days', fontsize = 40)
plt.ylabel(r'population density (cell/$\mu$l)', fontsize = 40)
plt.xticks(fontsize = 35, weight = 'bold')
plt.yticks(fontsize = 35, weight = 'bold')
plt.title('yeast population collapse', fontsize = 40)
plt.savefig("../results/ecosystems/Dai12.png", transparent=True, dpi=600)
