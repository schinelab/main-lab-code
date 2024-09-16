import matplotlib.pyplot as plt
import seaborn as sns


def plot_shit(df, xcol, ycol, title, ylabel, xlabel, filename):
    sns.set(style="whitegrid")

    sns.lineplot(data=df, x=xcol, y=ycol, marker='o')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    # Set additional x-ticks
    xticks = pd.date_range(start=df[xcol].min(), end=df[xcol].max(), freq='12H')
    plt.xticks(ticks=xticks, rotation=90)

    plt.savefig(f'{filename}.png', bbox_inches='tight')
    plt.show()
    sns.set(style="whitegrid")