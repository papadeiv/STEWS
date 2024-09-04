from pandas import read_csv
from matplotlib import pyplot as plt

# Function for custom axis tick format
def format_func(value, tick_number):
    # find number of multiples of pi/2
    N = int(np.round(2 * value / np.pi))
    if N == 0:
        return "0"
    elif N == 1:
        return r"$\pi/2$"
    elif N == 2:
        return r"$\pi$"
    elif N % 2 > 0:
        return r"${0}\pi/2$".format(N)
    else:
        return r"${0}\pi$".format(N // 2)

ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
fig

# Render LaTeX to Matplotlib (https://matplotlib.org/stable/gallery/text_labels_and_annotations/tex_demo.html)
# needs extra packages (https://stackoverflow.com/questions/11354149/python-unable-to-render-tex-in-matplotlib)
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    })

timeseries = read_csv('HFcompoundPhaseTransition.csv', sep = ',', engine = 'python', header=0, index_col=0)
print(type(timeseries))
print(timeseries.head)

fig = plt.figure(figsize=[12.8,9.6], dpi=200, layout='tight')
plt.plot(timeseries, linewidth = 4, color='gray')
plt.xscale("log")
plt.xlabel(r'temperature (K)', fontsize = 40)
plt.ylabel(r'resonance (THz)', fontsize = 40)
plt.xticks(fontsize = 35, weight = 'bold')
plt.yticks(fontsize = 35, weight = 'bold')
plt.title('heavy Fermi compound phase transition', fontsize = 40)
plt.savefig("../results/physics/Yang23.png", transparent=True, dpi=600)
