#!/usr/bin/env python3

"""Combine the six retained Stairway Plot demographic scenarios."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
DATA_DIRECTORY = (
    REPOSITORY_ROOT / "data" / "1.5.demography" / "stairway"
)
RESULT_DIRECTORY = REPOSITORY_ROOT / "results" / "demography"

SCENARIOS = [
    {
        "mutation": "1.01e8",
        "generation": "50y",
        "filename": (
            "pop1_1.01e8_50y_demography_analysis.final.summary"
        ),
        "color": "firebrick",
    },
    {
        "mutation": "1.01e8",
        "generation": "100y",
        "filename": (
            "pop1_1.01e8_100y_demography_analysis.final.summary"
        ),
        "color": "darkred",
    },
    {
        "mutation": "4.2e8",
        "generation": "50y",
        "filename": (
            "pop1_4.2e8_50y_demography_analysis.final.summary"
        ),
        "color": "gold",
    },
    {
        "mutation": "4.2e8",
        "generation": "100y",
        "filename": (
            "pop1_4.2e8_100y_demography_analysis.final.summary"
        ),
        "color": "goldenrod",
    },
    {
        "mutation": "5.2e8",
        "generation": "50y",
        "filename": (
            "pop1_5.2e8_50y_demography_analysis.final.summary"
        ),
        "color": "royalblue",
    },
    {
        "mutation": "5.2e8",
        "generation": "100y",
        "filename": (
            "pop1_5.2e8_100y_demography_analysis.final.summary"
        ),
        "color": "navy",
    },
]

MUTATION_LABELS = {
    "1.01e8": r"1.01$\times 10^{-8}$",
    "4.2e8": r"4.2$\times 10^{-8}$",
    "5.2e8": r"5.2$\times 10^{-8}$",
}


def read_summary(filename: str) -> pd.DataFrame:
    """Read one retained Stairway Plot final summary."""

    summary_file = DATA_DIRECTORY / filename

    if not summary_file.is_file():
        raise FileNotFoundError(
            f"Required Stairway Plot summary not found: {summary_file}"
        )

    return pd.read_csv(summary_file, sep="\t", comment="#")


def export_minimum_summary() -> None:
    """Export minimum year and minimum median Ne for each scenario."""

    records = []

    for scenario in SCENARIOS:
        dataframe = read_summary(scenario["filename"])

        records.append(
            {
                "file": (
                    f"{scenario['mutation']}_"
                    f"{scenario['generation']}"
                ),
                "min_year": dataframe["year"].min(),
                "min_Ne_median": dataframe["Ne_median"].min(),
            }
        )

    output = pd.DataFrame(records)
    output.to_csv(
        RESULT_DIRECTORY / "min_Ne_summary.csv",
        index=False,
    )

    print(output)
    print(
        "Summary exported to "
        f"{RESULT_DIRECTORY / 'min_Ne_summary.csv'}"
    )


def plot_scenarios() -> None:
    """Generate the combined six-scenario demographic-history figure."""

    plt.figure(figsize=(13, 8))

    for scenario in SCENARIOS:
        dataframe = read_summary(scenario["filename"])

        years_ka = dataframe["year"] / 1000
        median_ne = dataframe["Ne_median"] / 1000
        lower_ci = dataframe["Ne_2.5%"] / 1000
        upper_ci = dataframe["Ne_97.5%"] / 1000

        label = (
            f"{MUTATION_LABELS[scenario['mutation']]}_"
            f"{scenario['generation']}"
        )

        plt.fill_between(
            years_ka,
            lower_ci,
            upper_ci,
            color=scenario["color"],
            alpha=0.15,
        )

        plt.plot(
            years_ka,
            median_ne,
            linewidth=2,
            color=scenario["color"],
            alpha=0.8,
            label=label,
        )

    plt.axvspan(
        19,
        26,
        color="grey",
        alpha=0.2,
        label="LGM",
    )

    plt.axvline(
        11.7,
        color="green",
        linestyle="--",
        linewidth=2,
        alpha=0.7,
        label="Holocene start (11.7 ka)",
    )

    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(0.05, 20000)
    plt.ylim(0.1, 300000)

    plt.xlabel("ka Before Present", fontsize=15)
    plt.ylabel(r"$N_e$ ($\times 10^{3}$)", fontsize=15)

    plt.xticks(
        [0.1, 1, 10, 100, 1000, 10000],
        ["0.1", "1", "10", "100", "1000", "10000"],
    )

    plt.yticks(
        [0.1, 1, 10, 100, 1000, 10000, 100000],
        ["0.1", "1", "10", "100", "1000", "10000", "100000"],
    )

    plt.legend(
        fontsize=12,
        loc="upper left",
        ncol=2,
        framealpha=0.95,
    )

    plt.tight_layout()

    output_file = RESULT_DIRECTORY / "demography_paleo_events.png"

    plt.savefig(output_file, dpi=300)
    print(f"Figure exported to {output_file}")

    plt.show()


def main() -> None:
    """Run summary export and combined plotting."""

    RESULT_DIRECTORY.mkdir(parents=True, exist_ok=True)

    export_minimum_summary()
    plot_scenarios()


if __name__ == "__main__":
    main()
