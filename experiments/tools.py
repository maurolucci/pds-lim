import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import functools
from math import ceil, isnan
import numpy as np
import perfprof as pp


def show_solved_instances(data, solver):
    data2 = data[data.solver == solver]
    labels = range(0, 6, 1)
    palette = sns.color_palette("Spectral", n_colors=len(labels))
    sns.scatterplot(
        x=data2.instance, y=data2.k, hue=data2.result, hue_order=labels, palette=palette
    )
    plt.xticks(rotation=90)
    plt.legend(title=solver, facecolor="white", loc=2)
    plt.title("Number of solved instances by " + solver)
    plt.show()


def merge2_dataframes(data, solvers):
    datas = [data[data.solver == solver] for solver in solvers]
    data2 = datas[0].copy()
    for i in range(1, len(solvers) - 1):
        data2 = data2.merge(
            datas[i],
            on=["instance", "vertices", "edges", "propagating_vertices", "k"],
            suffixes=["", "_" + solvers[i]],
        )
    data2 = data2.merge(
        datas[len(solvers) - 1],
        on=["instance", "vertices", "edges", "propagating_vertices", "k"],
        suffixes=["_" + solvers[0], "_" + solvers[len(solvers) - 1]],
    )
    return data2


def safe_greater(x, y, rtol):
    if np.isclose(x, y, rtol, atol=0.0):
        return False
    return x > y


def safe_lower(x, y, rtol):
    if np.isclose(x, y, rtol, atol=0.0):
        return False
    return x < y


# First criterior: most solved instances
# Second criterion: lowest execution time
# Third criterion: lowest upper bound
# Fourth criterion: greatest lower bound
def compare(t1, t2):
    if safe_greater(t1[1], t2[1], rtol=0.0):
        return -1
    elif safe_lower(t1[1], t2[1], rtol=0.0):
        return 1
    if safe_lower(t1[2], t2[2], rtol=1e-2):
        return -1
    elif safe_greater(t1[2], t2[2], rtol=1e-2):
        return 1
    if safe_lower(t1[3], t2[3], rtol=1e-5):
        return -1
    elif safe_greater(t1[3], t2[3], rtol=1e-5):
        return 1
    if safe_greater(t1[4], t2[4], rtol=1e-5):
        return -1
    elif safe_lower(t1[4], t2[4], rtol=1e-5):
        return 1
    return 0


def get_tuple(serie, solver):
    return (
        serie["solver_" + solver],
        serie["result_" + solver],
        serie["time_opt_" + solver],
        serie["upper_bound_" + solver],
        serie["lower_bound_" + solver],
    )


def get_winner(serie, solvers):
    ls = [get_tuple(serie, solver) for solver in solvers]
    ls = sorted(ls, key=functools.cmp_to_key(compare))
    return ls[0][0]


def get_winner_or_tie(serie, solvers):
    if compare(get_tuple(serie, solvers[0]), get_tuple(serie, solvers[1])) == 0:
        return "Tie"
    else:
        return get_winner(serie, solvers)


def show_best_solver(data, solvers, ties=False):

    # Get winner solver
    data = merge2_dataframes(data, solvers)
    if ties and len(solvers) == 2:
        data["winner"] = data.apply(get_winner_or_tie, solvers=solvers, axis=1)
    else:
        data["winner"] = data.apply(get_winner, solvers=solvers, axis=1)
    data2 = data.groupby("winner", as_index=False).agg(number=("instance", "count"))
    data2.sort_values(by="winner")

    # Pie chart
    plt.pie(data2.number, labels=data2.winner, autopct="%1.1f%%")
    plt.title("Percentage of instances won by solver")
    plt.show()

    # Scatterplot
    plt.figure(figsize=(6, 6))
    ax = plt.gca()
    ax.set_aspect(0.9)
    ax.set_axisbelow(True)
    sns.scatterplot(
        x=data.instance,
        y=data.k,
        hue=data.winner,
        hue_order=np.sort(data.winner.unique()),
        style=data.winner,
    )
    plt.grid(alpha=0.5, zorder=-1, axis="y")
    plt.title("Winning solver by network and number of channels")
    plt.xticks(rotation=90)
    plt.legend()
    plt.show()


def show_cumulative_time(data, solvers, log_scale: bool = True, interval=None):
    if interval is not None:
        data1 = data[(data["L"] / data["L_star"]).apply(lambda x: x in interval)]
    else:
        data1 = data
    ax = sns.ecdfplot(
        data=data1[data1.solver.isin(solvers)],
        x="t_solver",
        hue=data1[data1.solver.isin(solvers)].solver,
        log_scale=log_scale,
    )
    ax.set(
        title="Cummulative percentage of solved runs over execution time by solver",
        xlabel="Time (seconds)",
        ylabel="Runs achieving optimality (%)",
    )
    ax.legend_.set_title(None)
    plt.show()

    # sns.boxplot(
    #     data=data1[data1.solver.isin(solvers)],
    #     x="t_solver",
    #     y="solver",
    #     hue=data1[data1.solver.isin(solvers)].solver,
    #     showfliers=False,
    #     legend=False,
    # ).set(
    #     title="Boxplot of execution time by solver", xlabel="Time (seconds)", ylabel=""
    # )
    # plt.show()


def show_execution_time(data, solvers, grid: bool = False, log_scale: bool = False):
    data = data[data.solver.isin(solvers)]
    if not grid:
        # sns.set_theme(rc={"figure.figsize": (10, 5)})
        sns.pointplot(
            data=data,
            x="instance",
            y="t_solver",
            hue="solver",
            markers=["o", "x", "s", "+"],
            linestyles=["-", "--", "-.", ":"],
            errorbar=None,
            linewidth=1.5,
        )
        plt.title("Average execution time")
        if log_scale:
            plt.yscale("log")
        plt.ylabel("time (s)")
        plt.xticks(rotation=90)
        plt.legend()
        plt.show()

        # sns.set_theme(rc={"figure.figsize": (10, 5)})
        # sns.pointplot(
        #     data=data,
        #     x="instance",
        #     y="t_solver",
        #     hue="solver",
        #     markers=["o", "x", "s", "+"],
        #     linestyle="none",
        #     linewidth=1.5,
        #     errorbar=("pi", 50),
        #     dodge=True,
        # )
        # plt.style.use("default")
        # plt.title("Interquartile range of execution time")
        # if log_scale:
        #     plt.yscale("log")
        # plt.ylabel("time (s)")
        # plt.xticks(rotation=90)
        # plt.show()

        # sns.set_theme(rc={"figure.figsize": (10, 5)})
        # sns.pointplot(
        #     data=data,
        #     x="instance",
        #     y="t_solver",
        #     hue="solver",
        #     markers=["o", "x", "s", "+"],
        #     linestyles=["-", "--", "-.", ":"],
        #     linewidth=1.25,
        #     errorbar=("pi", 50),
        #     dodge=True,
        # )
        # plt.style.use("default")
        # plt.title("Execution time")
        # if log_scale:
        #     plt.yscale("log")
        # plt.ylabel("time (s)")
        # plt.xticks(rotation=90)
        # plt.show()
    else:
        g = sns.FacetGrid(
            data=data,
            col="k_class",
            hue="solver",
            hue_kws={
                "markers": ["o", "x", "s", "+"],
                "linestyles": ["-", "--", "-.", ":"],
            },
            margin_titles=True,
        )
        g.map(
            sns.pointplot,
            "instance",
            "t_solver",
            markers="^",
            linewidth=1,
            order=data.instance.unique(),
            errorbar=None,
            # errorbar=("pi", 50),
            # dodge=True,
        )
        g.set_xticklabels(rotation=90)
        g.fig.suptitle("Execution time", y=1.02)
        g.set_axis_labels("instance", "time (s)")
        if log_scale:
            g.set(yscale="log")
        g.set_titles(col_template="k/k* {col_name}")
        g.add_legend()


def show_gap_1(data, solvers, grid: bool = False):
    data = data[data.solver.isin(solvers)]
    if not grid:
        sns.set_theme(rc={"figure.figsize": (10, 5)})
        sns.pointplot(
            data=data,
            x="instance",
            y="gap",
            hue="solver",
            markers=["o", "x", "s", "+"],
            linestyles=["-", "--", "-.", ":"],
            linewidth=1.5,
        )
        plt.title("Gap")
        plt.ylabel("gap (%)")
        plt.xticks(rotation=90)
        plt.show()
    else:
        g = sns.FacetGrid(
            data=data,
            col="L_class",
            hue="solver",
            hue_kws={
                "markers": ["o", "x", "s", "+"],
                "linestyles": ["-", "--", "-.", ":"],
            },
            margin_titles=True,
        )
        g.map(
            sns.pointplot,
            "instance",
            "gap",
            markers="^",
            linewidth=1,
            order=data.instance.unique(),
        )
        g.set_xticklabels(rotation=90)
        g.fig.suptitle("Gap", y=1.02)
        g.set_axis_labels("instance", "gap (%)")
        g.set_titles(col_template="L/L* {col_name}")
        g.add_legend()


def show_performance_profile(data, solvers, log_scale: bool = False):
    data = data[data.solver.isin(solvers)]
    data2 = data[["instance", "solver", "time_all"]].pivot_table(
        index="instance", columns="solver", values="time_all"
    )
    palette = ["o-C0", "x--C1", "s-.C2", "+:C3", "o-C4", "o--C5", "o-.C6"]
    pp.perfprof(data2, palette, markersize=4, markevery=[0])
    plt.legend(data2.columns)
    if log_scale:
        plt.xscale("log")
    plt.show()
