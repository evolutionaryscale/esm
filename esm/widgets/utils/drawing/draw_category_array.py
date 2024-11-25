import random
from typing import Sequence

import ipywidgets as widgets
import matplotlib.colors as mcolors
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from IPython.display import clear_output
from matplotlib.colors import Normalize


def draw_data_array(
    output: widgets.Output,
    data_array: Sequence[int | float],
    categories: list[str] = [],
    category_color_mapping: dict = {},
    pixel_height: int = 100,
    cmap="tab20",
    randomize_cmap=False,
    normalize_cmap=False,
    use_legend=True,
    highlighted_ranges: list[tuple[int, int, int | str | tuple]] = [],
):
    cmap_ = plt.get_cmap(cmap)

    def generate_color_palette(categories, category_color_mapping):
        if len(categories) == 0:
            # Continuous data, use colorbar
            return [], {}, []

        category_to_index = {category: idx for idx, category in enumerate(categories)}

        if normalize_cmap:
            cmap_colors = [cmap_(i / len(categories)) for i in range(len(categories))]
        else:
            cmap_colors = [cmap_(i) for i in range(cmap_.N)]

        cmap_colors_rgb = [mcolors.to_rgb(color) for color in cmap_colors]
        if randomize_cmap:
            rng = random.Random(42)
            rng.shuffle(cmap_colors_rgb)
        rgb_colors = cmap_colors_rgb[: len(categories)]

        for category, color in category_color_mapping.items():
            rgb_colors[category_to_index[category]] = mcolors.hex2color(color)

        if len(rgb_colors) < len(categories):
            raise ValueError("Not enough colors to match the number of categories.")

        return categories, category_to_index, rgb_colors

    categories, category_to_index, rgb_colors = generate_color_palette(
        categories, category_color_mapping
    )
    data_array_ = np.array(data_array)

    with output:
        clear_output(wait=True)

        fig, ax = plt.subplots(figsize=(12, pixel_height / 100))

        for idx, value in enumerate(data_array_):
            if len(categories) > 0:
                category_name = categories[value]
                color = rgb_colors[category_to_index[category_name]]
            else:
                max_value = max(data_array_)
                color = cmap_(value / max_value)

            rect = patches.Rectangle(
                (idx, 0),
                1,
                1,
                linewidth=1,
                edgecolor=mcolors.to_rgba("gray", alpha=0.1),
                facecolor=color,
            )
            ax.add_patch(rect)

        # Add highlighted ranges with bounding boxes
        if highlighted_ranges:
            for start, end, color in highlighted_ranges:
                rect = patches.Rectangle(
                    (start, -0.1),
                    end - start + 1,
                    1.2,
                    linewidth=2,
                    edgecolor=color,
                    facecolor="none",
                )
                ax.add_patch(rect)

        ax.set_xlim(0, len(data_array_))
        ax.set_ylim(-0.1, 1.1)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis("off")
        ax.grid(which="minor", color="gray", linestyle="-", linewidth=2)
        ax.tick_params(bottom=False, left=False, labelbottom=False, labelleft=False)

        if use_legend:
            if len(categories) == 0:
                sm = plt.cm.ScalarMappable(
                    cmap=cmap,
                    norm=Normalize(vmin=min(data_array_), vmax=max(data_array_)),
                )
                sm.set_array([])
                plt.colorbar(sm, ax=ax, orientation="horizontal", label="Value")
            else:
                legend_patches = [
                    patches.Patch(
                        color=rgb_colors[category_to_index[category]], label=category
                    )
                    for category in categories
                ]
                ax.legend(
                    handles=legend_patches,
                    bbox_to_anchor=(1.05, 1),
                    loc="upper left",
                    borderaxespad=0.0,
                )

        plt.show()
