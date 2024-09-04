from pandas import read_csv
from matplotlib import pyplot as plt

# Render LaTeX to Matplotlib (https://matplotlib.org/stable/gallery/text_labels_and_annotations/tex_demo.html)
# needs extra packages (https://stackoverflow.com/questions/11354149/python-unable-to-render-tex-in-matplotlib)
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    })

timeseries = read_csv('S&P500FinancialCrisis.csv', sep = ',', engine = 'python', header=0, parse_dates = [0], index_col=0)
print(type(timeseries))
print(timeseries.head)

fig = plt.figure(figsize=[12.8,9.6], dpi=200, layout='tight')
plt.plot(timeseries, linewidth = 4, color='purple')
plt.xlabel(r'yr', fontsize = 20)
plt.ylabel(r'S\&P 500 index', fontsize = 40)
plt.xticks(['2008-01', '2009-01', '2010-01'], ['2008', '2009', '2010'], fontsize = 35)
plt.yticks(fontsize = 35)
plt.title('2009 global financial crisis', fontsize = 40)
plt.savefig("../results/financial/Diks18.png", transparent=True, dpi=600)
