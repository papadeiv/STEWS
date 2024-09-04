from pandas import read_csv
from matplotlib import pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    })

timeseries = read_csv('CarbonateSite1218.csv', sep = ',', engine = 'python', header=0, index_col=0)
print(type(timeseries))
print(timeseries.head)

fig = plt.figure(figsize=[12.8,9.6], dpi=200, layout='tight')
plt.gca().invert_xaxis()
plt.plot(timeseries, linewidth = 4)
plt.xlabel(r'Myr', fontsize = 40)
plt.ylabel(r'$CaCO_{3}$ concentration', fontsize = 40)
plt.xticks(fontsize = 35, weight = 'bold')
plt.yticks(fontsize = 35, weight = 'bold')
plt.title('end of Greenhouse Earth', fontsize = 40)
plt.savefig("../../results/climate/paleoclimate/Tripati05.png", transparent=True, dpi=600)
