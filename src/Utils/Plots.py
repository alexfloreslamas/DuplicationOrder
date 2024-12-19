import os
import pandas as pd
import matplotlib.pyplot as plt


def validate_columns(df: pd.DataFrame, required_columns: list[str]) -> None:
    """Validates that the required columns exist in the DataFrame."""
    for col in required_columns:
        if col not in df.columns:
            raise ValueError(f"Missing required column in DataFrame: {col}")


def plot(df: pd.DataFrame, save_path: str | None = None) -> None:
    """
    Plot comparisons from a given DataFrame.

    Args:
        df (pd.DataFrame): The input DataFrame containing metrics data.
        save_path (str | None): Directory to save plots. If None, plots are displayed.
    """

    df['og'] = df['og'].astype(str)     # Ensure 'og' is treated as strings for even distribution on the X-axis

    config = {
        "fig_size": (10, 6),
        "x_label": 'OG',
        "y_labels": ["Precision", "Recall", "Contradiction"],
        "tree_labels": [
            r"$\left(\tilde{T}, \hat{T}\right)$",
            r"$\left(T^*, \hat{T}\right)$"
        ]
    }

    # Define metrics and titles dynamically
    heads: list[tuple[str, str]] = [
        ("precision1", "precision2"), ("recall1", "recall2"), ("contradiction1", "contradiction2")
    ]
    titles = [
        f"{metric} Comparison: {metric}{config['tree_labels'][0]} vs {metric}{config['tree_labels'][1]}"
        for metric in config["y_labels"]
    ]

    # Validate DataFrame columns
    required_columns = ['og'] + [col for head in heads for col in head]
    validate_columns(df, required_columns)

    # Ensure the save path exists
    if save_path:
        os.makedirs(save_path, exist_ok=True)

    # Plot each comparison
    for idx, (head, ylabel) in enumerate(zip(heads, config["y_labels"])):
        plt.figure(figsize=config["fig_size"])
        plt.plot(df['og'], df[head[0]], label=f"{ylabel}{config['tree_labels'][0]}", marker='o')
        plt.plot(df['og'], df[head[1]], label=f"{ylabel}{config['tree_labels'][1]}", marker='x')
        plt.title(titles[idx])
        plt.xlabel(config["x_label"])
        plt.ylabel(ylabel)
        plt.legend()
        plt.grid(True)
        plt.xticks(rotation=90)
        plt.tight_layout()

        if save_path:
            file: str = f"{save_path}plot_{ylabel.replace(' ', '_')}.png"
            plt.savefig(file)
            print(f"Plot saved: {file}")
        else:
            plt.show()


# def main():
#     input_data_path: str = "../../output/results.tsv"
#     output_path: str = "../../output/plots/"
#     df = pd.read_csv(input_data_path, sep='\t')
#     plot(df, output_path)
#
#
# if __name__ == "__main__":
#     main()
