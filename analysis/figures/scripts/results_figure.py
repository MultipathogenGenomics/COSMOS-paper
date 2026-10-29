from analysis.common.paths import RESULTS_DIR
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import seaborn as sns
from importlib import import_module
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.transforms import blended_transform_factory
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from analysis.figures.scripts.process import combine_data_for_results_figure
from analysis.common.mosquito_ids import shorted_genus
import geodatasets
import geopandas as gpd
from matplotlib.patches import Circle
from shapely.geometry import LineString
from matplotlib.transforms import Affine2D
from matplotlib.patches import ConnectionPatch

from analysis.common.annu_plot_groups import (
    ANNU_ABUNDANCE_THRESHOLD_FOR_COLOUR, AGG_SITE_FOR_AVERAGE_BUFFER,
    ANNU_MIN_COUNT_THRESHOLD, add_annu_site_clusters, annu_abundance_filter,
    average_annu_sites, write_annu_site_averages,
)

A4_WIDTH_INCHES = 8.27

COASTAL_PALETTE = {
    1: "#0072B2",
    0: "#F28E2B",
    None: "#D9D9D9",
}

GDR_COLOUR = "#8B3A3A"
MVEV_COLOUR = "#D62728"

# MVEV_COLOUR = "#DC267F"

def draw_points_as_averaged_gradients(ax: plt.Axes, averaged_sites: pl.DataFrame, colours: list[str], radius: float) -> None:
    # for each site, at the central lat/lon, draw a circle with colour based on the mean ratio
    cmap = sns.blend_palette(colours, as_cmap=True)
    for mean_lon, mean_lat, mean_ratio in averaged_sites.select(
        "site_longitude", "site_latitude", "mean_saus_ratio"
    ).iter_rows():
        colour = cmap(mean_ratio)
        circle = plt.Circle(
            (mean_lon, mean_lat),
            radius=radius,
            facecolor=colour,
            edgecolor="white",
            fill=True,
            linewidth=0.5,
            zorder=20,
        )
        ax.add_patch(circle)

def add_small_bar_chart_to_points(ax: plt.Axes, traps_to_add_bars_to: list, df: pl.DataFrame, bar_width: float = 0.1, bar_height: float = 0.7, outline_colour: str = "black") -> None:
    df = df.filter(pl.col("sample").is_in(traps_to_add_bars_to)) # REMOVE if want to show all bars
    TRAPS_TO_ADD_BARS_TO = {
        "23-013-29692Cos_S14": (-0.5, 0.25),  # 143.564511
        "23-013-29587Cos_S8": (-0.5, 0.8),   # 144.758824
        "23-013-29478Cos_S17": (-0.5, 0.3),  # 146.035165; shared site
        '23-013-29549Cos_S6': (-0.5, 0.3),
        "23-013-29548Cos_S5": (0.5, 0.6),   # 146.065377
        "23-013-30091Cos_S3": (-0.5, 0.25),   # 146.396361
        "23-013-29644Cos_S12": (-0.3, 0.3),  # 147.531248
        "23-013-29488Cos_S2": (1, 1),
        "23-013-29588Cos_S10": (0.01, 0.25),
    }
    between_bar_buffer = 0.05
    for trap in traps_to_add_bars_to:
        bar_offset_x, bar_offset_y = TRAPS_TO_ADD_BARS_TO[trap]
        # find all traps at that site
        central_trap_lon = df.filter(pl.col("sample") == trap).select("longitude").to_numpy()[0][0]
        central_trap_lat = df.filter(pl.col("sample") == trap).select("latitude").to_numpy()[0][0]
        site_traps = df.filter(
            pl.col("longitude").is_between(central_trap_lon - AGG_SITE_FOR_AVERAGE_BUFFER, central_trap_lon + AGG_SITE_FOR_AVERAGE_BUFFER)
        ).filter(
            pl.col("latitude").is_between(central_trap_lat - AGG_SITE_FOR_AVERAGE_BUFFER, central_trap_lat + AGG_SITE_FOR_AVERAGE_BUFFER)
        ).filter(
            pl.col("Culex annulirostris proportion") > 0
        ).filter(
            pl.col("avg_subtype_ratio") > 0
        )
        # order by season then week
        site_traps = site_traps.sort(["season", "weeknum"])
        # for each trap, draw a small bar
        for i, (trap, annu_prop, saus_proportion) in enumerate(site_traps.select("sample", "Culex annulirostris proportion", "avg_subtype_ratio").iter_rows()):
            if trap in traps_to_add_bars_to:
                outline = True
            else:
                outline = False
            bar_x = central_trap_lon + bar_offset_x + i * bar_width + i * between_bar_buffer
            bar_y = central_trap_lat + bar_offset_y
            annu_height = annu_prop * bar_height
            non_annu_height = bar_height - annu_height
            saus_height = saus_proportion * annu_height
            aus_height = annu_height - saus_height
            # bar for SAUS in orange, bar for AUS in blue, bar for non-annu in grey, stacked
            saus_patch = plt.Rectangle(
                (bar_x, bar_y),
                bar_width,
                saus_height,
                facecolor=COASTAL_PALETTE[0],
                edgecolor=COASTAL_PALETTE[0],
                linewidth=0.5,
                zorder=30,
            )
            aus_patch = plt.Rectangle(
                (bar_x, bar_y + saus_height),
                bar_width,
                aus_height,
                facecolor=COASTAL_PALETTE[1],
                edgecolor=COASTAL_PALETTE[1],
                linewidth=0.5,
                zorder=30,
            )
            non_annu_patch = plt.Rectangle(
                (bar_x, bar_y + annu_height),
                bar_width,
                non_annu_height,
                facecolor="lightgrey",
                edgecolor="lightgrey",
                linewidth=0.5,
                zorder=30,
            )
            ax.add_patch(saus_patch)
            ax.add_patch(aus_patch)
            ax.add_patch(non_annu_patch)
            # if outline is True, draw a black rectangle around the bar
            if outline:
                outline_patch = plt.Rectangle(
                    (bar_x, bar_y),
                    bar_width,
                    bar_height,
                    facecolor="none",
                    edgecolor=outline_colour,
                    linewidth=1.2,
                    zorder=40,
                )
                ax.add_patch(outline_patch)
        # draw a line connecting the centre of the trap to the bottom centre of the bars
        bar_attachment_y = central_trap_lat + bar_offset_y
        if bar_offset_y < 0:
            bar_attachment_y += bar_height
        bar_attachment_x = (
            central_trap_lon
            + bar_offset_x
            + len(site_traps) * bar_width / 2
            + (len(site_traps) - 1) * between_bar_buffer / 2
        )
        line = ConnectionPatch(
            xyA=(central_trap_lon, central_trap_lat),
            xyB=(bar_attachment_x, bar_attachment_y),
            coordsA="data",
            coordsB="data",
            axesA=ax,
            axesB=ax,
            color="black",
            linewidth=0.5,
            zorder=25,
        )
        ax.add_artist(line)

def draw_circle_around_traps(ax, traps_to_circle, df, outline_colour="black"):
    trap_df = df.filter(
        pl.col("sample").is_in(traps_to_circle)
    ).select(
        ["longitude", "latitude", "avg_subtype_ratio"]
    )
    cmap = sns.blend_palette(
        [COASTAL_PALETTE[1], COASTAL_PALETTE[0]],
        as_cmap=True,
    )
    for longitude, latitude, saus_prop in trap_df.iter_rows():
        colour = cmap(saus_prop)
        circle = plt.Circle(
            (longitude, latitude),
            radius=0.2,
            facecolor='none', # not colouring as now the sites are averages
            edgecolor=outline_colour,
            fill=True,
            linewidth=1,
            zorder=20,
        )
        ax.add_patch(circle)

def plot_cx_annulirostris_map(fig, subplot_spec, map_file: str, metadata_file: str, df: pl.DataFrame, traps_to_circle: list = None, traps_to_bars: list = None, outline_colour: str = "black", release_lineage_points=None, sydney_inset_bounds=None, site_averages_file=None) -> tuple:
    cryptic_process = import_module("analysis.cryptic_mosquito_population_structure.Cx_annulirostris.scripts.plotting.process")
    cryptic_plot = import_module("analysis.cryptic_mosquito_population_structure.Cx_annulirostris.scripts.plotting.plot")
    subtype_colours = {
        "S-AUS": COASTAL_PALETTE[0],
        "AUS": COASTAL_PALETTE[1],
    }
    panel_config = {
        "label": "Culex annulirostris",
        "species_abundance_column": "Culex annulirostris",
        "file_template": "data/cryptic_subpops/{dataset}.Cx_annulirostris_collated_results.tsv",
        "ratio_columns": ["avg_subtype_ratio", "non_subtype_ratio"],
        "legend_labels": ["S-AUS", "AUS"],
        "colours": [subtype_colours["S-AUS"], subtype_colours["AUS"]],
        "gradient_colours": [subtype_colours["AUS"], subtype_colours["S-AUS"]],
        "min_count_threshold": ANNU_MIN_COUNT_THRESHOLD,
        "complement_ratio_column": "avg_subtype_ratio",
        "complement_output_column": "non_subtype_ratio",
        "draw_type": "gradient"
    }
    if release_lineage_points is None:
        subpopulation_data = cryptic_process.prepare_cryptic_subpopulation_data(
            map_file,
            metadata_file,
            [panel_config],
            ["COSMOS1", "COSMOS2", "COSMOS3", "COSMOS4"],
        )
        panel = subpopulation_data["panels"][0]
        map_boundary = subpopulation_data["map_boundary"]
    else:
        map_boundary = cryptic_process.load_nsw_act_boundary(map_file).to_crs("EPSG:4326")
        point_data = release_lineage_points.to_pandas()
        points = gpd.GeoDataFrame(point_data, geometry=gpd.points_from_xy(
            point_data["longitude"], point_data["latitude"]), crs="EPSG:4326")
        panel = {"points": points,
                 "gradient_colours": panel_config["gradient_colours"]}

    samples_to_plot = set(df.filter(annu_abundance_filter())["sample"].to_list())
    points_to_plot = panel["points"][
        panel["points"]["sample"].isin(samples_to_plot)
    ]
    print(
        f"Cx. annulirostris samples plotted: {len(points_to_plot)} "
        f"(abundance threshold={ANNU_ABUNDANCE_THRESHOLD_FOR_COLOUR})"
    )

    # Draw and export the very same groups, so publication tables cannot drift.
    plotted_samples = pl.from_pandas(points_to_plot[
        ["sample", "longitude", "latitude", "avg_subtype_ratio"]
    ].reset_index(drop=True)).join(
        df.select("sample", "Culex annulirostris proportion", "location", "location-site", "MVEV"),
        on="sample", how="left", validate="1:1",
    )
    averaged_sites = average_annu_sites(plotted_samples)
    if site_averages_file is not None:
        write_annu_site_averages(averaged_sites, site_averages_file)

    map_ax = fig.add_subplot(subplot_spec)
    map_min_x, _, map_max_x, _ = map_boundary.total_bounds
    map_xlim = (map_min_x, map_max_x - 0.2)
    map_boundary.boundary.plot(ax=map_ax, color="black", linewidth=0.5)
    map_ylim = map_ax.get_ylim()
    # cryptic_plot.draw_points_as_gradients(
    #     map_ax,
    #     points_to_plot,
    #     panel["ratio_columns"],
    #     panel["gradient_colours"],
    #     radius=0.2,
    # )
    draw_points_as_averaged_gradients(
        map_ax,
        averaged_sites,
        panel["gradient_colours"],
        radius=0.2,
    )
    mainland_xlim = (140.9, 153.8)
    mainland_ylim = (-37.6, -28.0)

    map_ax.set_xlim(mainland_xlim)
    map_ax.set_ylim(mainland_ylim)

    map_ax.set_aspect("equal")
    map_ax.set_axis_off()

    # draw on the GDR coords from AK
    gdr_coords = np.array([
        [152.0310561857285, -29.05884438306106],
        [151.8419028204735, -29.74480814986759],
        [150.9466485608887, -32.0576776944991],
        [150.3174422786903, -33.71480730074187],
        [149.7800386106203, -34.75207284189212],
        [149.2961935232750, -35.36848469567541],
        [148.8017118226620, -36.21266565329890],
        [148.5659720513401, -36.82561305143643],
    ])
    map_ax.plot(
        gdr_coords[:, 0],
        gdr_coords[:, 1],
        color=GDR_COLOUR,
        linewidth=3,
        alpha=0.3,
        zorder=10,
    )
    gdr_labels = [
        ("Great", 2, 11, 0),
        ("Dividing", 1, -6, 3),
        ("Range", 0, -10, 0),
    ]

    for label, segment_index, along_offset, rotation_offset in gdr_labels:

        p1 = gdr_coords[segment_index + 1]
        p2 = gdr_coords[segment_index]
        midpoint = (p1 + p2) / 2

        p1_display = map_ax.transData.transform(p1)
        p2_display = map_ax.transData.transform(p2)

        dx = p2_display[0] - p1_display[0]
        dy = p2_display[1] - p1_display[1]
        length = np.hypot(dx, dy)

        tangent_x = dx / length
        tangent_y = dy / length
        normal_x = -dy / length
        normal_y = dx / length

        angle = np.degrees(np.arctan2(dy, dx)) + rotation_offset
        distance_from_line = 8

        map_ax.annotate(
            label,
            xy=midpoint,
            xytext=(
                normal_x * distance_from_line + tangent_x * along_offset,
                normal_y * distance_from_line + tangent_y * along_offset,
            ),
            textcoords="offset points",
            ha="center",
            va="center",
            rotation=angle,
            rotation_mode="anchor",
            color=GDR_COLOUR,
            fontsize=8,
            zorder=11,
        )

    # circle traps with MVEV detections
    if traps_to_circle:
        traps_to_circle = [
            sample for sample in traps_to_circle
            if sample in samples_to_plot
        ]
        draw_circle_around_traps(map_ax, traps_to_circle, df, outline_colour=outline_colour)

    if traps_to_bars:
        add_small_bar_chart_to_points(map_ax, traps_to_bars, df, outline_colour='none')

    map_position = map_ax.get_position()
    inset_left = 0.8
    inset_bottom = -0.05
    inset_width = 0.28
    inset_height = 0.38

    inset_ax = fig.add_axes([
        map_position.x0 + map_position.width * inset_left,
        map_position.y0 + map_position.height * inset_bottom,
        map_position.width * inset_width,
        map_position.height * inset_height,
    ])
    inset_ax.set_box_aspect(1)
    map_boundary.boundary.plot(ax=inset_ax, color="black", linewidth=0.5)
    # cryptic_plot.draw_points_as_gradients(
    #     inset_ax,
    #     points_to_plot,
    #     panel["ratio_columns"],
    #     panel["gradient_colours"],
    #     radius=0.04,
    # )
    draw_points_as_averaged_gradients(
        inset_ax,
        averaged_sites,
        panel["gradient_colours"],
        radius=0.04,
    )
    if sydney_inset_bounds is None:
        xmin, xmax, ymin, ymax = cryptic_plot.get_sydney_inset_limits(
            panel["points"], panel["points"].crs)
    else:
        xmin, xmax, ymin, ymax = sydney_inset_bounds
    inset_ax.set_xlim(xmin, xmax)
    inset_ax.set_ylim(ymin, ymax)
    inset_ax.set_aspect("equal")
    inset_ax.set_xticks([])
    inset_ax.set_yticks([])
    inset_ax.set_title("", fontsize=8)
    mark_inset(map_ax, inset_ax, loc1=3, loc2=1, fc="none", ec="black", lw=0.8)

    # # FLAT GLOBE INSET
    # inset_left = -0.2
    # inset_bottom = 0.7
    # inset_width = 0.64
    # inset_height = 0.32

    # world_inset_ax = fig.add_axes([
    #     map_position.x0 + map_position.width * inset_left,
    #     map_position.y0 + map_position.height * inset_bottom,
    #     map_position.width * inset_width,
    #     map_position.height * inset_height,
    # ])

    # world_inset_ax.set_box_aspect(0.5)
    # world = gpd.read_file(
    #     geodatasets.get_path("naturalearth.land")
    # )
    # world.boundary.plot(
    #     ax=world_inset_ax,
    #     color="black",
    #     linewidth=0.3,
    # )

    # # nsw_world_crs = map_boundary.to_crs(world.crs)
    # # nsw_world_crs.plot(
    # #     ax=world_inset_ax,
    # #     color="none",
    # #     edgecolor="#FF073A",
    # #     linewidth=0.4,
    # #     zorder=5,
    # # )
    # # draw a an arrow pointing to centre of NSW with an angle
    # nsw_centre_lon = 147.0
    # nsw_centre_lat = -32.5

    # world_inset_ax.annotate(
    #     "",
    #     xy=(nsw_centre_lon, nsw_centre_lat),
    #     xytext=(nsw_centre_lon + 20, nsw_centre_lat + 30),
    #     arrowprops={
    #         "arrowstyle": "-|>",
    #         "color": "#FF073A",
    #         "linewidth": 3,
    #         "shrinkA": 0,
    #         "shrinkB": 0,
    #         "mutation_scale": 10,
    #         "joinstyle": "miter",
    #         "capstyle": "butt",
    #     },
    #     zorder=10,
    # )

    # world_inset_ax.set_xlim(-180, 180)
    # world_inset_ax.set_ylim(-90, 90)
    # world_inset_ax.set_aspect("equal")
    # world_inset_ax.set_xticks([])
    # world_inset_ax.set_yticks([])

    # ROUND GLOBE INSET
    inset_left = -0.15
    inset_bottom = -0.04
    inset_width = 0.32
    inset_height = 0.32

    tilt_degrees = 12

    world_inset_ax = fig.add_axes([
        map_position.x0 + map_position.width * inset_left,
        map_position.y0 + map_position.height * inset_bottom,
        map_position.width * inset_width,
        map_position.height * inset_height,
    ])
    world = gpd.read_file(geodatasets.get_path("naturalearth.land"))
    # Centre approximately over the Indo-Pacific
    globe_crs = (
        "+proj=ortho "
        "+lat_0=-25 "
        "+lon_0=134 "
        "+datum=WGS84 "
        "+units=m "
        "+no_defs"
    )
    globe_zoom = 1.5 # 1 would be whole hemisphere
    world_globe = world.to_crs(globe_crs)
    world_globe.plot(
        ax=world_inset_ax,
        facecolor="white",
        edgecolor="black",
        linewidth=0.3,
    )
    world_inset_ax.set_box_aspect(1)
    world_inset_ax.set_aspect("equal")
    world_inset_ax.set_axis_off()

    globe_rotation = (
        Affine2D()
        .rotate_deg_around(0, 0, tilt_degrees)
        + world_inset_ax.transData
    )

    # draw the outline of NSW in neon red
    neon_red = "#FF073A"
    grey = "dimgrey"
    nsw_globe = map_boundary.to_crs(globe_crs)
    nsw_globe.boundary.plot(
        ax=world_inset_ax,
        color=neon_red,
        facecolor=neon_red,
        linewidth=0,
        zorder=10,
    )
    # draw on some longitudes
    longitudes_to_draw = [140, 155]
    longitude_transform = blended_transform_factory(
        map_ax.transData,  # x position uses longitude
        map_ax.transAxes,  # y position uses 0–1 axes coordinates
    )
    for i, long in enumerate(longitudes_to_draw):
        horizontal_alignment = "right" if i == 0 else "left"
        meridian = gpd.GeoSeries(
            [
                LineString(
                    [(long, latitude) for latitude in np.linspace(-90, 90, 361)]
                )
            ],
            crs="EPSG:4326",
        ).to_crs(globe_crs)

        meridian.plot(
            ax=world_inset_ax,
            color=grey,
            linewidth=0.5,
            linestyle="--",
            alpha=0.7,
            zorder=2,
        )
        # at lat = -50, annotate
        label_point = gpd.GeoSeries(
            gpd.points_from_xy([long], [-50]),
            crs="EPSG:4326",
        ).to_crs(globe_crs).iloc[0]

        # world_inset_ax.annotate(
        #     f"{long}°E",
        #     xy=(label_point.x, label_point.y),
        #     color=neon_red,
        #     fontsize=10,
        #     ha=horizontal_alignment,
        #     va="top",
        #     xytext=(0, -2),
        #     textcoords="offset points",
        #     clip_on=False,
        #     zorder=10,
        # )
        map_ax.plot(
            [long, long],
            [0.7, 0.675],
            transform=longitude_transform,
            color=grey,
            linewidth=1,
            clip_on=False,
            zorder=20,
        )

        map_ax.annotate(
            f"{long}°E",
            xy=(long, 0.7),
            xycoords=longitude_transform,
            xytext=(0, 2),
            textcoords="offset points",
            ha="center",
            va="bottom",
            color=grey,
            fontsize=9,
            clip_on=False,
            zorder=20,
        )

    for collection in world_inset_ax.collections:
        collection.set_transform(globe_rotation)

    globe_radius = 6_371_000  # approximately Earth’s radius in metres
    globe_outline = Circle(
        (0.5, 0.5),
        0.495,
        facecolor="none",
        edgecolor="black",
        linewidth=0.5,
        transform=world_inset_ax.transAxes,
        clip_on=False,
        zorder=100,
    )
    world_inset_ax.add_patch(globe_outline)
    for collection in world_inset_ax.collections:
        collection.set_clip_path(globe_outline)

    visible_radius = globe_radius / globe_zoom
    world_inset_ax.set_xlim(-visible_radius, visible_radius)
    world_inset_ax.set_ylim(-visible_radius, visible_radius)

    # annotate on the main map, near the top a tick with 140E and 155E

        
    return map_ax

def plot_SDI_boxplot(ax, df: pl.DataFrame):
    # colour boxes by Coastal and Urban/Rural status
    coastal_urban_rural_palette = {
        "1 - Rural": COASTAL_PALETTE[1],
        "1 - Met_Fringe": COASTAL_PALETTE[1],
        "1 - Urban": COASTAL_PALETTE[1],
        "0 - Rural": COASTAL_PALETTE[0],
        "0 - Met_Fringe": COASTAL_PALETTE[0],
        "0 - Urban": COASTAL_PALETTE[0],
        "1 - Met": COASTAL_PALETTE[1],
        "0 - Met": COASTAL_PALETTE[0],
    }
    category_col = "Coastal_urban_rural"
    plot_source_df = (
        df.with_row_index("_x")
        .filter(
            pl.col("Simpsons_Diversity_Index").is_not_null()
            & ~pl.col(category_col).str.starts_with("Unknown")
        )
    )
    # going to replace the above colours with the shade of the proportion of S-AUS
    cmap = sns.blend_palette(
        [COASTAL_PALETTE[1], COASTAL_PALETTE[0]],
        as_cmap=True,
    )
    category_saus_prop = (
        plot_source_df.filter(
            pl.col("avg_subtype_ratio").is_not_null()
            & pl.col("avg_subtype_ratio").is_finite()
        )
        .group_by(category_col, maintain_order=True)
        .agg(
            pl.col("avg_subtype_ratio")
            .mean()
            .alias("mean_S_AUS_proportion")
        )
    )
    coastal_urban_rural_palette = {
        category: cmap(proportion)
        for category, proportion in zip(
            category_saus_prop[category_col].to_list(),
            category_saus_prop["mean_S_AUS_proportion"].to_list(),
        )
        if proportion is not None
    }
    # compute box positions along the x axis as the midpoint of the category
    category_positions = (
        plot_source_df.group_by(category_col, maintain_order=True)
        .agg(
            (
                (pl.col("_x").min() + pl.col("_x").max() + 1) / 2
            ).alias("midpoint")
        )
    )
    box_positions = dict(
        zip(
            category_positions[category_col].to_list(),
            category_positions["midpoint"].to_list(),
        )
    )
    plot_df = (
        plot_source_df.with_columns(
            pl.col(category_col)
            .replace_strict(box_positions)
            .cast(pl.Float64)
            .alias("_box_position")
        )
        .to_pandas()
    )
    sns.boxplot(
        data=plot_df,
        x="_box_position",
        y="Simpsons_Diversity_Index",
        hue=category_col,
        palette=coastal_urban_rural_palette,
        native_scale=True,
        dodge=False,
        legend=False,
        width=0.5,
        fliersize=0,
        ax=ax,
    )
    sns.stripplot(
        data=plot_df,
        x="_box_position",
        y="Simpsons_Diversity_Index",
        native_scale=True,
        color="black",
        size=2,
        jitter=0.15,
        alpha=0.5,
        ax=ax,
    )
    # The first and last boxes can extend slightly beyond the heatmap bounds.
    for artist in [*ax.lines, *ax.patches, *ax.collections]:
        artist.set_clip_on(False)

    # Match the full width of the heatmaps beneath it.
    ax.set_xlim(0, df.height)
    ordered_categories = category_positions[category_col].to_list()
    ax.set_xticks([box_positions[category] for category in ordered_categories])
    ax.set_xticklabels(
        [
            category.split(" - ", maxsplit=1)[-1].replace("Met_Fringe", "Fringe")
            for category in ordered_categories
        ],
        fontsize=8,
    )

    # label x tick labels along top
    ax.set_xlabel("")
    ax.set_ylabel("Diversity\n(1 - D)", fontsize=9, rotation=0, va="center", ha="right")
    ax.tick_params(
        axis="x",
        top=True,
        labeltop=True,
        bottom=False,
        labelbottom=False,
    )
    ax.grid(axis="y", linestyle="-", alpha=0.5)
    add_vertical_lines(ax, df)
    ax.set_ylim(0, 1)
    # set ticks at .25, .5, .75
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticklabels(["0", "0.25", "0.5", "0.75", "1"], fontsize=8)
    # remove all spines
    for spine in ax.spines.values():
        spine.set_visible(False)
    # add a bracket at the top of the plot to indicate Coastal vs Inland
    bracket_info = {
        (0, 3): "Inland",
        (3, 6): "Coastal",
    }
    bottom_left_y = ax.get_ylim()[1] + 0.55
    bracket_ax = ax.inset_axes([0, bottom_left_y, 1, 0.2])
    add_brackets_horizontal(bracket_ax, bracket_info, brakcet_colours = ["black", "black"])

def plot_mosquito_heatmap(ax, df: pl.DataFrame):
    # first take the top 10 species by mean abundance across all samples
    mosquito_genera = [
        "Aedes",
        "Anopheles",
        "Culex",
        "Toxorhynchites",
        "Ochlerotatus",
        "Haemagogus",
        "Mansonia",
        "Psorophora",
        "Coquillettidia",
        "Uranotaenia",
        "Tripteroides",
        "Verrallina",
        "Sabethes",
    ]
    species_cols = [col for col in df.columns if any(genus in col for genus in mosquito_genera)]
    # get the mean abundance for each species across all samples
    mean_abundance = df.select(species_cols).mean().to_dict(as_series=False)
    # sort species by mean abundance
    sorted_species = sorted(mean_abundance, key=mean_abundance.get, reverse=True)
    top_species = sorted_species[:10]
    # now doing custom order
    top_species = ["Culex annulirostris proportion", "Aedes notoscriptus", "Culex australicus", "Anopheles annulipes", "Culex pipiens", "Mansonia uniformis", "Aedes vigilax", "Culex sitiens", "Ochlerotatus procax", "Coquillettidia linealis"]
    # draw heatmap of top species
    heatmap_data = df.select(top_species).to_numpy().T
    sns.heatmap(
        heatmap_data,
        ax=ax,
        cmap="Greys",
        cbar=False,
        yticklabels=False,
    )
    # Colour Cx. annulirostris by subtype ratio, with colour strength set by
    # its relative abundance in the sample.
    row = top_species.index("Culex annulirostris proportion")
    subtype_ratio = df["avg_subtype_ratio"].to_numpy()
    abundance = heatmap_data[row, :]
    s_aus_colour = COASTAL_PALETTE[0]
    aus_colour = COASTAL_PALETTE[1]
    cmap = sns.blend_palette(
        [aus_colour, s_aus_colour],
        as_cmap=True,
    )
    if ANNU_ABUNDANCE_THRESHOLD_FOR_COLOUR == 0:
        passes_abundance_threshold = abundance >= 0
    else:
        passes_abundance_threshold = (
            abundance > ANNU_ABUNDANCE_THRESHOLD_FOR_COLOUR
        )
    valid = (
        np.isfinite(subtype_ratio)
        & np.isfinite(abundance)
        & passes_abundance_threshold
    )
    subtype_colours = cmap(np.clip(np.nan_to_num(subtype_ratio), 0, 1))
    colour_strength = np.clip(np.nan_to_num(abundance), 0, 1)[:, None]
    cell_colours = np.ones((df.height, 4))
    cell_colours[:, :3] = (
        1 - colour_strength * (1 - subtype_colours[:, :3])
    )
    cell_colours[:, 3] = 1

    # Recolour the existing seaborn cells so this row has exactly the same
    # dimensions as every other heatmap row.
    heatmap_mesh = ax.collections[0]
    heatmap_colours = heatmap_mesh.get_facecolors().reshape(
        heatmap_data.shape[0],
        heatmap_data.shape[1],
        4,
    )
    heatmap_colours[row, valid, :] = cell_colours[valid, :]
    heatmap_mesh.set_array(None)
    heatmap_mesh.set_facecolors(heatmap_colours.reshape(-1, 4))
    # no ticks on the x-axis
    ax.set_xticks([])
    # label y-axis with species names
    ax.set_yticks(np.arange(len(top_species)) + 0.5)
    ytick_labels = [
        "Cx. annulirostris" if species == "Culex annulirostris proportion"
        else shorted_genus(species)
        for species in top_species
    ]
    # ytick_labels = [label.replace("annulirostris", "annu.") for label in ytick_labels]
    ax.set_yticklabels(ytick_labels, rotation=0, fontsize=9, fontstyle="italic")
    # ax.set_ylabel("Vectors", fontsize=11)

    add_outline_to_heatmap(ax)
    add_vertical_lines(ax, df)

def add_outline_to_heatmap(ax):
    # add very thin line outline
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.5)
        spine.set_color("black")
        spine.set_alpha(0.3)

def plot_virus_nc2_heatmap(ax, df: pl.DataFrame):
    virus_acronyms = {
        "stratfordvirus": "STRV",
        "rossrivervirus": "RRV",
        "murrayvalleyencephalitisvirus": "MVEV",
        "kunjinvirus": "KUN",
        "kokoberavirus": "KOK",
        "edgehillvirus": "EHV",
        "japaneseencephalitisvirus": "JEV",
        "barmahforestvirus": "BFV",
        "sindbisvirus": "SINDV",
        "umbrevirus": "LSV",
    }
    virus_cols = [virus_acronyms[v] for v in virus_acronyms.keys()]
    virus_cols = [col for col in virus_cols if col in df.columns]
    # custom order
    virus_cols = ["MVEV", "JEV", "EHV", "STRV", "KOK", "KUN", "RRV", "BFV", "SINDV", "LSV"]
    heatmap_data = (
        df.select(
            pl.col(virus_cols).cast(pl.Float64, strict=False)
        )
        .fill_null(0)
        .to_numpy()
        .T
    )
    # filter to where nc2 > 0.05
    heatmap_data = np.where(heatmap_data < 0.05, 0, heatmap_data)
    sns.heatmap(
        heatmap_data,
        ax=ax,
        cmap="Greys",
        cbar=False,
        yticklabels=False,
    )

    # going to make the top MVEV row colouring go from white to MVEV colour, instead of greyscale
    mvev_row = virus_cols.index("MVEV")
    mvev_overlay = np.full_like(
        heatmap_data,
        np.nan,
        dtype=float,
    )
    mvev_overlay[mvev_row, :] = heatmap_data[mvev_row, :]
    mvev_cmap = sns.blend_palette(
        ["white", MVEV_COLOUR],
        as_cmap=True,
    )
    sns.heatmap(
        mvev_overlay,
        mask=np.isnan(mvev_overlay),
        ax=ax,
        cmap=mvev_cmap,
        vmin=0,
        vmax=1,
        cbar=False,
        xticklabels=False,
        yticklabels=False,
    )
    # add outline to the heatmap
    add_outline_to_heatmap(ax)

    # no ticks on the x-axis
    ax.set_xticks([])
    # virus acronyms as y-axis labels
    ax.set_yticks(np.arange(len(virus_cols)) + 0.5)
    ax.set_yticklabels(virus_cols, rotation=0, fontsize=9)
    ax.set_ylabel("Viruses", fontsize=11)
    # draw vertical lines at the boundaries between coastal_urban_rural groups
    add_vertical_lines(ax, df)

    # to the left of the heatmap, between y tick labels and y axis label, I'm going to add a symbol to indicate anthrigenic or encephalitis
    virus_symptom_info = {
        "MVEV": "E",
        "JEV": "E",
        "EHV": None,
        "STRV": None,
        "KOK": "E",
        "KUN": "E",
        "RRV": "A",
        "BFV": "A",
        "SINDV": "A",
        "LSV": "E"
    }
    marker_info = {
        "A": {"label": "Arthrogenic", "xpos": -0.16},
        "E": {"label": "Encephalitic", "xpos": -0.20},
    }
    marker_transform = blended_transform_factory(ax.transAxes, ax.transData)
    for i, virus in enumerate(virus_cols):
        fc = MVEV_COLOUR if virus == "MVEV" else "black"
        symptom = virus_symptom_info[virus]
        for symptom_code, info in marker_info.items():
            ax.scatter(
                x=info["xpos"],
                y=i + 0.5,
                marker="o",
                facecolor=fc if symptom == symptom_code else "none",
                edgecolor="black",
                linewidth=0.5,
                s=30,
                zorder=10,
                clip_on=False,
                transform=marker_transform,
            )
    # add label at bottom of the markers, 45 rotated, with the text "Arthrogenic" and "Encephalitic" below the markers
    for symptom_code, info in marker_info.items():
        ax.text(
            x=info["xpos"],
            y=len(virus_cols) + 0.05,
            s=info["label"],
            ha="right",
            va="top",
            fontsize=8,
            rotation=45,
            rotation_mode="anchor",
            transform=marker_transform,
            clip_on=False,
        )

def add_vertical_lines(ax, df: pl.DataFrame):
    coastal_urban_rural_boundaries = df.select("Coastal_urban_rural").to_numpy().flatten()
    boundary_indices = np.where(coastal_urban_rural_boundaries[:-1] != coastal_urban_rural_boundaries[1:])[0] + 1
    for boundary_index in boundary_indices:
        ax.axvline(boundary_index, color="grey", linewidth=0.5, alpha=1)

def add_brackets(ax, bracket_info):
    n_rows = max(end for start, end in bracket_info)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, n_rows)
    ax.set_axis_off()

    bracket_x = 0.45
    bracket_arm_end = 0.15
    label_x = 0.60
    bracket_buffer = 0.12
    for (start, end), label in bracket_info.items():
        y_bottom = start + bracket_buffer
        y_top = end - bracket_buffer
        y_mid = (y_bottom + y_top) / 2
        ax.plot(
            [bracket_x, bracket_x],
            [y_top, y_bottom],
            color="black",
            linewidth=0.8,
            clip_on=False,
        )
        ax.plot(
            [bracket_x, bracket_arm_end],
            [y_top, y_top],
            color="black",
            linewidth=0.8,
            clip_on=False,
        )
        ax.plot(
            [bracket_x, bracket_arm_end],
            [y_bottom, y_bottom],
            color="black",
            linewidth=0.8,
            clip_on=False,
        )
        ax.text(
            label_x,
            y_mid,
            label,
            ha="left",
            va="center",
            fontsize=8,
            clip_on=False,
        )

def add_brackets_horizontal(ax, bracket_info, brakcet_colours = None):
    n_cols = max(end for start, end in bracket_info)
    ax.set_xlim(0, n_cols)
    ax.set_ylim(0, 1)
    ax.set_axis_off()

    bracket_y = 0.45
    bracket_arm_end = 0
    label_y = 0.60
    bracket_buffer = 0.12
    for (start, end), label in bracket_info.items():
        if brakcet_colours:
            colour = brakcet_colours.pop(0)
        else:
            colour = "black"
        x_left = start + bracket_buffer
        x_right = end - bracket_buffer
        x_mid = (x_left + x_right) / 2
        ax.plot(
            [x_left, x_right],
            [bracket_y, bracket_y],
            color=colour,
            linewidth=0.8,
            clip_on=False,
        )
        ax.plot(
            [x_left, x_left],
            [bracket_y, bracket_arm_end],
            color=colour,
            linewidth=0.8,
            clip_on=False,
        )
        ax.plot(
            [x_right, x_right],
            [bracket_y, bracket_arm_end],
            color=colour,
            linewidth=0.8,
            clip_on=False,
        )
        ax.text(
            x_mid,
            label_y,
            label,
            ha="center",
            va="bottom",
            fontsize=11,
            clip_on=False,
        )

def determine_heatmap_start_end_col_gdr(
    df: pl.DataFrame,
) -> tuple[int, int]:
    values = df["E_or_W_of_GDR"].drop_nulls().to_numpy()
    initial_state = values[0]  # W
    final_state = values[-1]   # E
    transition_start = np.where(values == final_state)[0][0]
    transition_end = np.where(values == initial_state)[0][-1] + 1
    return transition_start + 1, transition_end

def connect_gdr_to_virus_heatmap(
    fig,
    map_ax,
    virus_ax,
    vector_ax,
    gdr_coords,
    start_col=71,
    end_col=77,
    colour="lightgrey",
):
    linewidth = 1
    linestyle = '--'
    # Northern endpoint of the GDR line.
    top_gdr_index = np.argmax(gdr_coords[:, 1])
    top_gdr_lon, top_gdr_lat = gdr_coords[top_gdr_index]
    # Draw a short vertical tick floating above the GDR endpoint.
    tick_gap = 0.8
    tick_height = 0.3
    tick_bottom = top_gdr_lat + tick_gap
    tick_top = tick_bottom + tick_height
    map_ax.plot(
        [top_gdr_lon, top_gdr_lon],
        [tick_bottom, tick_top],
        color=colour,
        linewidth=linewidth,
        clip_on=False,
        zorder=20,
    )
    # Rows 71–77 are one-based table positions. Heatmap cell 71 begins at
    # x=70, and cell 77 ends at x=77.
    bracket_left = start_col - 1
    bracket_right = end_col
    bracket_midpoint = (bracket_left + bracket_right) / 2
    # Seaborn heatmaps have an inverted y-axis. The first ylim value is the
    # lower edge of the displayed heatmap.
    virus_bottom = virus_ax.get_ylim()[0]
    bracket_y = virus_bottom + 0.25
    # Horizontal line covering columns 71–77.
    virus_ax.plot(
        [bracket_left, bracket_right],
        [bracket_y, bracket_y],
        color=colour,
        linewidth=linewidth + 1,
        clip_on=False,
        zorder=20,
    )
    # Short central tick extending from the bracket towards the map.
    bracket_tick_end = bracket_y + 0.35
    virus_ax.plot(
        [bracket_midpoint, bracket_midpoint],
        [bracket_y, bracket_tick_end],
        color=colour,
        linewidth=linewidth,
        clip_on=False,
        zorder=20,
    )
    # Connect the outward ends of the two vertical ticks.
    connector = ConnectionPatch(
        xyA=(top_gdr_lon, tick_top),
        coordsA=map_ax.transData,
        axesA=map_ax,
        xyB=(bracket_midpoint, bracket_tick_end),
        coordsB=virus_ax.transData,
        axesB=virus_ax,
        color=colour,
        linewidth=linewidth,
        linestyle=linestyle,
        clip_on=False,
        zorder=15,
    )
    fig.add_artist(connector)

    # gonna also shade the area of the GDR on each heatmap
    for ax in [virus_ax, vector_ax]:
        ax.axvspan(
            start_col - 1,
            end_col,
            color=colour,
            alpha=0.2,
            zorder=50,
        )

def plot_results_figure(
        metadata_file: str,
        SDI_results_file: str,
        mosquito_abundance_file: str,
        virus_nc2_file: str,
        NSW_map_file: str,
        out_file: str,
        prepared_data=None,
        release_lineage_points=None,
        sydney_inset_bounds=None,
        site_averages_file=None,
):
    fig = plt.figure(
        figsize=(7.5, 7),
        constrained_layout=False,
    )
    heatmap_prop = 0.4
    vector_len = (heatmap_prop / 20) * 10
    virus_len = (heatmap_prop / 20) * 10

    ## USE BELOW INSTEAD IF WANT WIDER MAP
    # outer_grid = fig.add_gridspec(
    #     nrows=4,
    #     ncols=4,
    #     width_ratios=[0.1, 0.53, 0.03, 0.34],
    #     height_ratios=[0.08, vector_len, virus_len, 0.57],
    #     left=0.08,
    #     right=0.95,
    #     top=0.92,
    #     bottom=0.03,
    #     wspace=0.02,
    #     hspace=0.1,
    #     )
    # SDI_boxplot_spec = outer_grid[0, 1]
    # mosquito_heatmap_spec = outer_grid[1, 1]
    # virus_heatmap_spec = outer_grid[2, 1]
    # NSW_map_spec = outer_grid[3, 0:2]
    # MVEV_tree_spec = outer_grid[1:3, 3]
    # vector_bracket_spec = outer_grid[1, 2]
    # virus_bracket_spec = outer_grid[2, 2]

    outer_grid = fig.add_gridspec(
        nrows=4,
        ncols=3,
        width_ratios=[0.95, 0.05, 0],#.37],
        height_ratios=[0.08, vector_len, virus_len, 0.52],
        left=0.2,
        right=0.9,
        top=0.92,
        bottom=0.03,
        wspace=0.02,
        hspace=0.1,
        )
    SDI_boxplot_spec = outer_grid[0, 0]
    mosquito_heatmap_spec = outer_grid[1, 0]
    virus_heatmap_spec = outer_grid[2, 0]
    NSW_map_spec = outer_grid[3, 0]
    MVEV_tree_spec = outer_grid[1:3, 2]
    vector_bracket_spec = outer_grid[1, 1]
    virus_bracket_spec = outer_grid[2, 1]

    # parse data - this is ordered
    if prepared_data is None:
        df = combine_data_for_results_figure(
            metadata_file=metadata_file,
            SDI_results_file=SDI_results_file,
            mosquito_abundance_file=mosquito_abundance_file,
            virus_nc2_file=virus_nc2_file,
            map_file=NSW_map_file
        )
    else:
        df = prepared_data
    assert df.height == 219
    assert df["sample"].n_unique() == 219

    # plot SDI boxplot
    ax_sdi = fig.add_subplot(SDI_boxplot_spec)
    plot_SDI_boxplot(ax_sdi, df)

    # plot mosquito heatmap
    ax_mosquito = fig.add_subplot(mosquito_heatmap_spec)
    plot_mosquito_heatmap(ax_mosquito, df)

    # add mosquito brackets
    mosquito_bracket_info = {
        (0, 5): "Coastal",
        (5, 10): "Widely\ndistributed",
    }
    mosquito_bracket_ax = fig.add_subplot(vector_bracket_spec)
    add_brackets(mosquito_bracket_ax, mosquito_bracket_info)

    # add the Vectors ylabel at the midpoint of the mosquito heatmap and boxpplots
    box_pos = ax_sdi.get_position()
    mos_pos = ax_mosquito.get_position()
    combined_y_mid = (box_pos.y1 + mos_pos.y0) / 2
    ax_mosquito.set_ylabel("Vectors", fontsize=11)
    label_transform = blended_transform_factory(
        ax_mosquito.transAxes,
        fig.transFigure,
    )
    ax_mosquito.yaxis.set_label_coords(
        -0.32,
        combined_y_mid,
        transform=label_transform,
    )

    # plot virus heatmap
    ax_virus = fig.add_subplot(virus_heatmap_spec)
    plot_virus_nc2_heatmap(ax_virus, df)
    # add virus brackets
    virus_bracket_info = {
        (0, 1): "Orthobunya-\nvirus",
        (1, 5): "Alphavirus",
        (5, 11): "Orthoflavi-\nvirus",
    }
    virus_bracket_ax = fig.add_subplot(virus_bracket_spec)
    add_brackets(virus_bracket_ax, virus_bracket_info)

    # plot NSW map
    traps_to_circle = df.filter(pl.col("MVEV") > 0.05).select("sample").to_series().to_list()
    traps_to_bars = df.filter(pl.col("ct_MVEV").is_not_null()).select("sample").to_series().to_list()
    ax_map = plot_cx_annulirostris_map(
        fig,
        NSW_map_spec,
        NSW_map_file,
        metadata_file,
        df=df,
        traps_to_circle=traps_to_circle,
        # traps_to_bars=traps_to_circle, # ask TG about those COSMOS did not detect
        outline_colour=MVEV_COLOUR, release_lineage_points=release_lineage_points,
        sydney_inset_bounds=sydney_inset_bounds, site_averages_file=site_averages_file)

    # add the missing MVEV trap as a dot on the map
    add_the_missing_MVEV_trap_as_dot(ax_map, df, trap_id="23-020623-30091_S21", colour=MVEV_COLOUR)

    # join GDR to virus heatmap
    gdr_coords = np.array([
        [152.0310561857285, -29.05884438306106],
        [151.8419028204735, -29.74480814986759],
        [150.9466485608887, -32.0576776944991],
        [150.3174422786903, -33.71480730074187],
        [149.7800386106203, -34.75207284189212],
        [149.2961935232750, -35.36848469567541],
        [148.8017118226620, -36.21266565329890],
        [148.5659720513401, -36.82561305143643],
    ])
    start_col, end_col = determine_heatmap_start_end_col_gdr(df)
    connect_gdr_to_virus_heatmap(fig, map_ax=ax_map, virus_ax=ax_virus, vector_ax=ax_mosquito, gdr_coords=gdr_coords, start_col=start_col, end_col=end_col, colour=GDR_COLOUR)

    # plot MVEV tree

    # save figure
    left_axes = [ax_mosquito, ax_virus]
    for ax in left_axes:
        x, y = ax.yaxis.label.get_position()
        ax.yaxis.set_label_coords(-0.34, y)
        ax.yaxis.label.set_verticalalignment("top")
    fig.align_ylabels(left_axes)
    fig.savefig(out_file, dpi=600)
    # save as pdf
    fig.savefig(out_file.replace(".png", ".pdf"), dpi=600)

def add_the_missing_MVEV_trap_as_dot(ax: plt.Axes, df: pl.DataFrame, trap_id: str, colour: str) -> None:
    # get lat and lon for the trap_id
    trap_row = df.filter(pl.col("sample") == trap_id).select(["latitude", "longitude"]).to_pandas()
    if trap_row.empty:
        print(f"Warning: trap_id {trap_id} not found in dataframe")
        return
    lat = trap_row["latitude"].values[0]
    lon = trap_row["longitude"].values[0]
    ax.scatter(lon, lat, color=colour, s=15, edgecolor=colour, marker="o", zorder=100)
