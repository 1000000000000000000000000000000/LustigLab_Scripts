#!/usr/bin/env python3

import argparse
import math
import sys
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_arguments():
	parser = argparse.ArgumentParser(
		description="Generate raw and normalized descriptor overlay plots from Nextflow descriptor and PyVkabat CSV outputs."
	)

	parser.add_argument(
		"--nextflow_csv",
		required=True,
		help="Path to <sequence_name>_nextflow.csv containing E6, isunstruct, Residue, and query_index columns."
	)

	parser.add_argument(
		"--vkabat_csv",
		required=True,
		help="Path to <sequence_name>_vkabat.csv containing vkabat and optional N columns."
	)

	parser.add_argument(
		"--outdir",
		default=".",
		help="Output directory for generated plot files. Default: current directory."
	)

	parser.add_argument(
		"--prefix",
		default=None,
		help="Output filename prefix and default plot-title prefix. Default: inferred from nextflow_csv filename."
	)

	parser.add_argument(
		"--title",
		default=None,
		help="Custom base plot title. If omitted, title is generated from prefix."
	)

	parser.add_argument(
		"--res_start",
		type=int,
		default=None,
		help="Optional first residue position to include in plots and output files. Default: include from first residue."
	)

	parser.add_argument(
		"--res_end",
		type=int,
		default=None,
		help="Optional last residue position to include in plots and output files. Default: include through last residue."
	)

	parser.add_argument(
		"--define_switch_regions",
		default=None,
		help="Optional comma-separated residue ranges to visualize as predefined/historical switch regions, e.g. '150-153,160-165,180-194'."
	)

	parser.add_argument(
		"--defined_switch_region_color",
		default="#7b3294",
		help="Color used to shade predefined/historical switch regions. Default: purple."
	)

	parser.add_argument(
		"--defined_switch_region_alpha",
		type=float,
		default=0.20,
		help="Transparency for predefined/historical switch-region shading. Default: 0.20."
	)

	parser.add_argument(
		"--defined_switch_region_label",
		default="Defined switch region",
		help="Legend label for predefined/historical switch regions. Default: Defined switch region."
	)

	parser.add_argument(
		"--e6_col",
		default="E6",
		help="Column name for E6 values in nextflow_csv. Default: E6."
	)

	parser.add_argument(
		"--isunstruct_col",
		default="isunstruct",
		help="Column name for IsUnstruct values in nextflow_csv. Default: isunstruct."
	)

	parser.add_argument(
		"--residue_col",
		default="Residue",
		help="Column name for one-letter residue abbreviations in nextflow_csv. Default: Residue."
	)

	parser.add_argument(
		"--position_col",
		default="query_index",
		help="Column name for residue position in nextflow_csv. Default: query_index."
	)

	parser.add_argument(
		"--vkabat_col",
		default="vkabat",
		help="Column name for Vkabat values in vkabat_csv. Default: vkabat."
	)

	parser.add_argument(
		"--n_col",
		default="N",
		help="Optional column name for number of predictions used in Vkabat calculation. Default: N."
	)

	parser.add_argument(
		"--N",
		dest="fallback_N",
		type=int,
		default=15,
		help="Fallback N value used when n_col is not present in vkabat_csv. Default: 15."
	)

	parser.add_argument(
		"--class_count",
		type=int,
		default=3,
		help="Number of mapped secondary-structure classes. Default: 3."
	)

	parser.add_argument(
		"--fig_width",
		type=float,
		default=10.0,
		help="Figure width in inches. Default: 10."
	)

	parser.add_argument(
		"--fig_height",
		type=float,
		default=5.5,
		help="Figure height in inches. Default: 5.5."
	)

	parser.add_argument(
		"--dpi",
		type=int,
		default=300,
		help="Output PNG resolution. Default: 300."
	)

	parser.add_argument(
		"--line_width",
		type=float,
		default=1.2,
		help="Line width for descriptor traces. Default: 1.2."
	)

	parser.add_argument(
		"--smoothing_window",
		type=int,
		default=1,
		help="Centered rolling-average window size in residues. Default: 1, no smoothing."
	)

	parser.add_argument(
		"--x_tick_interval",
		type=int,
		default=50,
		help="Spacing between x-axis tick marks in residue positions. Default: 50."
	)

	parser.add_argument(
		"--raw_y_tick_interval",
		type=float,
		default=1.0,
		help="Spacing between y-axis tick marks in the raw plot. Default: 1.0."
	)

	parser.add_argument(
		"--normalized_y_tick_interval",
		type=float,
		default=0.1,
		help="Spacing between y-axis tick marks in the normalized plot. Default: 0.1."
	)

	parser.add_argument(
		"--raw_ymax",
		type=float,
		default=None,
		help="Optional y-axis maximum for raw plot. Default: automatic."
	)

	parser.add_argument(
		"--e6_color",
		default="#1f77b4",
		help="Color for E6 trace. Default: blue."
	)

	parser.add_argument(
		"--isunstruct_color",
		default="#d62728",
		help="Color for IsUnstruct trace. Default: red."
	)

	parser.add_argument(
		"--vkabat_color",
		default="#ffb000",
		help="Color for Vkabat trace. Default: golden yellow."
	)

	parser.add_argument(
		"--legend_position",
		default="under_title",
		choices=["under_title", "right", "top", "best", "none"],
		help="Legend placement. Default: under_title."
	)

	parser.add_argument(
		"--legend_columns",
		type=int,
		default=4,
		help="Number of legend columns when legend is placed under the title or above the plot. Default: 4."
	)

	parser.add_argument(
		"--show_thresholds",
		action="store_true",
		help="Add horizontal dashed percentile-threshold lines and visualize threshold-derived switch-like regions."
	)

	parser.add_argument(
		"--show_threshold_percentiles_in_legend",
		action="store_true",
		help="Include threshold percentile values in legend labels. Default: omit percentile values from legend."
	)

	parser.add_argument(
		"--e6_threshold_percentile",
		type=float,
		default=50.0,
		help="Percentile used for the E6 threshold line when --show_thresholds is used. Residues must be <= this value. Default: 50."
	)

	parser.add_argument(
		"--isunstruct_threshold_percentile",
		type=float,
		default=50.0,
		help="Percentile used for the IsUnstruct threshold line when --show_thresholds is used. Residues must be <= this value. Default: 50."
	)

	parser.add_argument(
		"--vkabat_threshold_percentile",
		type=float,
		default=50.0,
		help="Percentile used for the Vkabat threshold line when --show_thresholds is used. Residues must be >= this value. Default: 50."
	)

	parser.add_argument(
		"--minimum_switch_region_length",
		type=int,
		default=3,
		help="Minimum number of consecutive residues required to define a threshold-derived switch-like region. Default: 3."
	)

	parser.add_argument(
		"--threshold_line_width",
		type=float,
		default=1.0,
		help="Line width for threshold lines. Default: 1.0."
	)

	parser.add_argument(
		"--threshold_alpha",
		type=float,
		default=0.8,
		help="Transparency for threshold lines. Default: 0.8."
	)

	parser.add_argument(
		"--threshold_linestyle",
		default="--",
		help="Matplotlib linestyle for threshold lines. Default: --."
	)

	parser.add_argument(
		"--switch_region_color",
		default="#ff0000",
		help="Color used to shade threshold-derived switch-like regions. Default: red."
	)

	parser.add_argument(
		"--switch_region_alpha",
		type=float,
		default=0.25,
		help="Transparency for threshold-derived switch-like regions. Default: 0.25."
	)

	parser.add_argument(
		"--no_grid",
		action="store_true",
		help="Disable plot grid."
	)

	return parser.parse_args()


def clean_dataframe(df):
	unnamed_columns = [column for column in df.columns if column.startswith("Unnamed:")]

	if unnamed_columns:
		df = df.drop(columns=unnamed_columns)

	return df


def infer_prefix(nextflow_csv_path):
	stem = Path(nextflow_csv_path).stem

	for suffix in ["_nextflow", "_final", ".nextflow", ".final"]:
		if stem.endswith(suffix):
			stem = stem[: -len(suffix)]

	return stem


def get_range_text(args):
	if args.res_start is None and args.res_end is None:
		return None

	if args.res_start is not None and args.res_end is not None:
		return f"Residues {args.res_start}-{args.res_end}"

	if args.res_start is not None:
		return f"Residues {args.res_start}-end"

	return f"Residues start-{args.res_end}"


def get_raw_title(args, prefix):
	if args.title is not None:
		title = args.title
	else:
		title = f"{prefix} Descriptor Overlay"

	range_text = get_range_text(args)

	if range_text is not None:
		title = f"{title} ({range_text})"

	if args.smoothing_window > 1:
		title = f"{title} ({args.smoothing_window}-Residue Rolling Average)"

	return title


def get_normalized_title(args, prefix):
	if args.title is not None:
		title = f"Normalized {args.title}"
	else:
		title = f"{prefix} Normalized Descriptor Overlay"

	range_text = get_range_text(args)

	if range_text is not None:
		title = f"{title} ({range_text})"

	if args.smoothing_window > 1:
		title = f"{title} ({args.smoothing_window}-Residue Rolling Average)"

	return title


def require_column(df, column_name, file_label):
	if column_name not in df.columns:
		available = ", ".join(df.columns)
		raise ValueError(
			f"Column '{column_name}' was not found in {file_label}. Available columns: {available}"
		)


def calculate_vkabat_max(n_value, class_count):
	n_value = int(n_value)

	if n_value < 1:
		raise ValueError(f"N must be at least 1. Encountered N={n_value}.")

	k_max = min(class_count, n_value)
	n1_min = math.ceil(n_value / k_max)

	return k_max * n_value / n1_min


def normalize_vkabat(vkabat_value, n_value, class_count):
	v_min = 1.0
	v_max = calculate_vkabat_max(n_value, class_count)

	if v_max == v_min:
		return 0.0

	return (vkabat_value - v_min) / (v_max - v_min)


def validate_percentile(value, argument_name):
	if value < 0 or value > 100:
		raise ValueError(f"{argument_name} must be between 0 and 100. Encountered {value}.")


def validate_arguments(args):
	if args.smoothing_window < 1:
		raise ValueError("--smoothing_window must be >= 1. Use 1 for no smoothing.")

	if args.x_tick_interval < 1:
		raise ValueError("--x_tick_interval must be >= 1.")

	if args.raw_y_tick_interval <= 0:
		raise ValueError("--raw_y_tick_interval must be > 0.")

	if args.normalized_y_tick_interval <= 0:
		raise ValueError("--normalized_y_tick_interval must be > 0.")

	if args.legend_columns < 1:
		raise ValueError("--legend_columns must be >= 1.")

	if args.threshold_line_width <= 0:
		raise ValueError("--threshold_line_width must be > 0.")

	if args.threshold_alpha < 0 or args.threshold_alpha > 1:
		raise ValueError("--threshold_alpha must be between 0 and 1.")

	if args.switch_region_alpha < 0 or args.switch_region_alpha > 1:
		raise ValueError("--switch_region_alpha must be between 0 and 1.")

	if args.defined_switch_region_alpha < 0 or args.defined_switch_region_alpha > 1:
		raise ValueError("--defined_switch_region_alpha must be between 0 and 1.")

	if args.minimum_switch_region_length < 1:
		raise ValueError("--minimum_switch_region_length must be >= 1.")

	if args.fallback_N < 1:
		raise ValueError("--N must be >= 1.")

	if args.class_count < 1:
		raise ValueError("--class_count must be >= 1.")

	if args.res_start is not None and args.res_start < 1:
		raise ValueError("--res_start must be >= 1.")

	if args.res_end is not None and args.res_end < 1:
		raise ValueError("--res_end must be >= 1.")

	if args.res_start is not None and args.res_end is not None and args.res_start > args.res_end:
		raise ValueError("--res_start must be less than or equal to --res_end.")

	validate_percentile(args.e6_threshold_percentile, "--e6_threshold_percentile")
	validate_percentile(args.isunstruct_threshold_percentile, "--isunstruct_threshold_percentile")
	validate_percentile(args.vkabat_threshold_percentile, "--vkabat_threshold_percentile")


def filter_residue_range(plot_df, args):
	if args.res_start is not None:
		plot_df = plot_df[plot_df["residue_position"] >= args.res_start]

	if args.res_end is not None:
		plot_df = plot_df[plot_df["residue_position"] <= args.res_end]

	plot_df = plot_df.reset_index(drop=True)

	if plot_df.empty:
		raise ValueError(
			"No residues remained after applying the selected residue range. "
			"Check --res_start, --res_end, and the residue positions in the input CSV file."
		)

	return plot_df


def add_smoothed_columns(plot_df, args):
	window = args.smoothing_window

	if window == 1:
		plot_df["E6_plot"] = plot_df["E6"]
		plot_df["IsUnstruct_plot"] = plot_df["IsUnstruct"]
		plot_df["Vkabat_plot"] = plot_df["Vkabat"]
		plot_df["E6_normalized_plot"] = plot_df["E6_normalized"]
		plot_df["IsUnstruct_normalized_plot"] = plot_df["IsUnstruct_normalized"]
		plot_df["Vkabat_normalized_plot"] = plot_df["Vkabat_normalized"]

		return plot_df

	plot_df["E6_plot"] = plot_df["E6"].rolling(window=window, center=True, min_periods=1).mean()
	plot_df["IsUnstruct_plot"] = plot_df["IsUnstruct"].rolling(window=window, center=True, min_periods=1).mean()
	plot_df["Vkabat_plot"] = plot_df["Vkabat"].rolling(window=window, center=True, min_periods=1).mean()

	plot_df["E6_normalized_plot"] = plot_df["E6_normalized"].rolling(window=window, center=True, min_periods=1).mean()
	plot_df["IsUnstruct_normalized_plot"] = plot_df["IsUnstruct_normalized"].rolling(window=window, center=True, min_periods=1).mean()
	plot_df["Vkabat_normalized_plot"] = plot_df["Vkabat_normalized"].rolling(window=window, center=True, min_periods=1).mean()

	return plot_df


def build_descriptor_dataframe(args):
	nextflow_df = clean_dataframe(pd.read_csv(args.nextflow_csv))
	vkabat_df = clean_dataframe(pd.read_csv(args.vkabat_csv))

	require_column(nextflow_df, args.e6_col, args.nextflow_csv)
	require_column(nextflow_df, args.isunstruct_col, args.nextflow_csv)
	require_column(nextflow_df, args.residue_col, args.nextflow_csv)
	require_column(nextflow_df, args.position_col, args.nextflow_csv)

	require_column(vkabat_df, args.vkabat_col, args.vkabat_csv)

	min_length = min(len(nextflow_df), len(vkabat_df))

	if len(nextflow_df) != len(vkabat_df):
		print(
			f"WARNING: Input files have different lengths. nextflow_csv has {len(nextflow_df)} rows; "
			f"vkabat_csv has {len(vkabat_df)} rows. Truncating to {min_length} rows.",
			file=sys.stderr
		)

	nextflow_df = nextflow_df.iloc[:min_length].reset_index(drop=True)
	vkabat_df = vkabat_df.iloc[:min_length].reset_index(drop=True)

	if args.n_col in vkabat_df.columns:
		n_values = pd.to_numeric(vkabat_df[args.n_col], errors="coerce").fillna(args.fallback_N)
		n_source = f"column '{args.n_col}'"
	else:
		n_values = pd.Series([args.fallback_N] * min_length)
		n_source = f"fallback --N value ({args.fallback_N})"
		print(
			f"WARNING: Column '{args.n_col}' was not found in {args.vkabat_csv}. "
			f"Using fallback --N value of {args.fallback_N} for all residues.",
			file=sys.stderr
		)

	plot_df = pd.DataFrame(
		{
			"residue_position": pd.to_numeric(nextflow_df[args.position_col], errors="coerce"),
			"residue": nextflow_df[args.residue_col].astype(str),
			"E6": pd.to_numeric(nextflow_df[args.e6_col], errors="coerce"),
			"IsUnstruct": pd.to_numeric(nextflow_df[args.isunstruct_col], errors="coerce"),
			"Vkabat": pd.to_numeric(vkabat_df[args.vkabat_col], errors="coerce"),
			"N": pd.to_numeric(n_values, errors="coerce"),
			"N_source": n_source
		}
	)

	plot_df["vkabat_row_index_position"] = plot_df.index + 1

	plot_df = plot_df.dropna(
		subset=[
			"residue_position",
			"E6",
			"IsUnstruct",
			"Vkabat",
			"N"
		]
	)

	plot_df["residue_position"] = plot_df["residue_position"].astype(int)
	plot_df["N"] = plot_df["N"].astype(int)
	plot_df = plot_df.sort_values("residue_position").reset_index(drop=True)

	plot_df["E6_max"] = math.log2(6)
	plot_df["E6_normalized"] = plot_df["E6"] / plot_df["E6_max"]

	plot_df["Vkabat_min"] = 1.0
	plot_df["Vkabat_max"] = plot_df["N"].apply(
		lambda n_value: calculate_vkabat_max(n_value, args.class_count)
	)

	plot_df["Vkabat_normalized"] = plot_df.apply(
		lambda row: normalize_vkabat(row["Vkabat"], row["N"], args.class_count),
		axis=1
	)

	plot_df["IsUnstruct_normalized"] = plot_df["IsUnstruct"]

	for column in ["E6_normalized", "Vkabat_normalized", "IsUnstruct_normalized"]:
		if (plot_df[column] < 0).any() or (plot_df[column] > 1).any():
			print(
				f"WARNING: Some values in {column} fall outside the expected 0-1 range.",
				file=sys.stderr
			)

	plot_df = filter_residue_range(plot_df, args)
	plot_df = add_smoothed_columns(plot_df, args)

	return plot_df


def get_empty_defined_region_dataframe():
	columns = [
		"defined_region_id",
		"input_region",
		"start_residue_position",
		"end_residue_position",
		"length",
		"visible_in_plot",
		"clipped_start_residue_position",
		"clipped_end_residue_position",
		"clipped_length"
	]

	return pd.DataFrame(columns=columns)


def parse_single_defined_region(region_text):
	region_text = region_text.strip()

	if not region_text:
		raise ValueError("Empty region definition encountered in --define_switch_regions.")

	if "-" in region_text:
		parts = region_text.split("-")

		if len(parts) != 2:
			raise ValueError(f"Invalid region definition '{region_text}'. Use formats like 150-153 or 150.")

		start = int(parts[0].strip())
		end = int(parts[1].strip())
	else:
		start = int(region_text)
		end = start

	if start < 1 or end < 1:
		raise ValueError(f"Defined switch-region positions must be >= 1. Encountered '{region_text}'.")

	if start > end:
		raise ValueError(f"Defined switch-region start must be <= end. Encountered '{region_text}'.")

	return start, end


def parse_defined_switch_regions(args, plot_df):
	if args.define_switch_regions is None or not args.define_switch_regions.strip():
		return get_empty_defined_region_dataframe()

	plot_min = int(plot_df["residue_position"].min())
	plot_max = int(plot_df["residue_position"].max())

	records = []
	region_tokens = [token.strip() for token in args.define_switch_regions.split(",") if token.strip()]

	for region_index, region_text in enumerate(region_tokens, start=1):
		start, end = parse_single_defined_region(region_text)
		clipped_start = max(start, plot_min)
		clipped_end = min(end, plot_max)
		visible_in_plot = clipped_start <= clipped_end
		clipped_length = clipped_end - clipped_start + 1 if visible_in_plot else 0

		records.append(
			{
				"defined_region_id": region_index,
				"input_region": region_text,
				"start_residue_position": start,
				"end_residue_position": end,
				"length": end - start + 1,
				"visible_in_plot": visible_in_plot,
				"clipped_start_residue_position": clipped_start if visible_in_plot else None,
				"clipped_end_residue_position": clipped_end if visible_in_plot else None,
				"clipped_length": clipped_length
			}
		)

	if not records:
		return get_empty_defined_region_dataframe()

	return pd.DataFrame(records)


def make_tick_list(start, stop, interval):
	ticks = []
	current = start

	while current <= stop:
		ticks.append(round(current, 6))
		current += interval

	if round(stop, 6) not in ticks:
		ticks.append(round(stop, 6))

	return ticks


def format_x_axis(ax, plot_df, args):
	max_position = int(plot_df["residue_position"].max())

	if args.res_start is not None:
		x_min = args.res_start
	else:
		x_min = 0

	if args.res_end is not None:
		x_max = args.res_end
	else:
		x_max = int(math.ceil(max_position / args.x_tick_interval) * args.x_tick_interval)

	if x_max < max_position:
		x_max = max_position

	first_interval_tick = int(math.ceil(x_min / args.x_tick_interval) * args.x_tick_interval)
	last_interval_tick = int(math.floor(x_max / args.x_tick_interval) * args.x_tick_interval)

	ticks = []

	if args.res_start is not None:
		ticks.append(x_min)
	else:
		ticks.append(0)

	if first_interval_tick <= last_interval_tick:
		ticks.extend(range(first_interval_tick, last_interval_tick + 1, args.x_tick_interval))

	if args.res_end is not None:
		ticks.append(x_max)

	ticks = sorted(set(ticks))

	ax.set_xlim(left=x_min, right=x_max)
	ax.set_xticks(ticks)
	ax.margins(x=0)


def format_y_axis(ax, plot_df, args, normalized=False):
	if normalized:
		y_min = 0.0
		y_max = 1.0
		ticks = make_tick_list(y_min, y_max, args.normalized_y_tick_interval)
	else:
		y_min = 0.0

		if args.raw_ymax is not None:
			y_max = float(args.raw_ymax)
		else:
			y_max = max(
				float(plot_df["E6"].max()),
				float(plot_df["IsUnstruct"].max()),
				float(plot_df["Vkabat"].max()),
				float(plot_df["Vkabat_max"].max())
			)
			y_max = math.ceil(y_max)

		ticks = make_tick_list(y_min, y_max, args.raw_y_tick_interval)

	ax.set_ylim(bottom=y_min, top=y_max)
	ax.set_yticks(ticks)
	ax.margins(y=0)


def get_unique_legend_items(ax):
	handles, labels = ax.get_legend_handles_labels()
	unique = {}

	for handle, label in zip(handles, labels):
		if label not in unique:
			unique[label] = handle

	return list(unique.values()), list(unique.keys())


def format_legend(fig, ax, args):
	if args.legend_position == "none":
		return None

	handles, labels = get_unique_legend_items(ax)

	if args.legend_position == "under_title":
		fig.legend(
			handles,
			labels,
			loc="upper center",
			bbox_to_anchor=(0.5, 0.925),
			ncol=args.legend_columns,
			frameon=False
		)

		return [0.0, 0.0, 1.0, 0.82]

	if args.legend_position == "right":
		ax.legend(
			handles,
			labels,
			loc="center left",
			bbox_to_anchor=(1.02, 0.5),
			ncol=1,
			frameon=False,
			borderaxespad=0
		)

		return None

	if args.legend_position == "top":
		fig.legend(
			handles,
			labels,
			loc="upper center",
			bbox_to_anchor=(0.5, 0.93),
			ncol=args.legend_columns,
			frameon=False
		)

		return [0.0, 0.0, 1.0, 0.84]

	ax.legend(
		handles,
		labels,
		loc="best",
		frameon=False
	)

	return None


def percentile_label(percentile):
	if float(percentile).is_integer():
		return str(int(percentile))

	return str(percentile)


def calculate_percentile_threshold(series, percentile):
	return float(series.quantile(percentile / 100.0))


def calculate_threshold_info(plot_df, args):
	if not args.show_thresholds:
		return None

	raw_thresholds = {
		"E6": calculate_percentile_threshold(plot_df["E6_plot"], args.e6_threshold_percentile),
		"IsUnstruct": calculate_percentile_threshold(plot_df["IsUnstruct_plot"], args.isunstruct_threshold_percentile),
		"Vkabat": calculate_percentile_threshold(plot_df["Vkabat_plot"], args.vkabat_threshold_percentile)
	}

	normalized_thresholds = {
		"E6": calculate_percentile_threshold(plot_df["E6_normalized_plot"], args.e6_threshold_percentile),
		"IsUnstruct": calculate_percentile_threshold(plot_df["IsUnstruct_normalized_plot"], args.isunstruct_threshold_percentile),
		"Vkabat": calculate_percentile_threshold(plot_df["Vkabat_normalized_plot"], args.vkabat_threshold_percentile)
	}

	threshold_info = {
		"raw": raw_thresholds,
		"normalized": normalized_thresholds
	}

	return threshold_info


def add_switch_like_flags(plot_df, threshold_info, args):
	if threshold_info is None:
		plot_df["E6_meets_threshold"] = False
		plot_df["IsUnstruct_meets_threshold"] = False
		plot_df["Vkabat_meets_threshold"] = False
		plot_df["meets_switch_like_criteria"] = False

		return plot_df

	raw_thresholds = threshold_info["raw"]

	plot_df["E6_threshold_value"] = raw_thresholds["E6"]
	plot_df["IsUnstruct_threshold_value"] = raw_thresholds["IsUnstruct"]
	plot_df["Vkabat_threshold_value"] = raw_thresholds["Vkabat"]

	plot_df["E6_meets_threshold"] = plot_df["E6_plot"] <= raw_thresholds["E6"]
	plot_df["IsUnstruct_meets_threshold"] = plot_df["IsUnstruct_plot"] <= raw_thresholds["IsUnstruct"]
	plot_df["Vkabat_meets_threshold"] = plot_df["Vkabat_plot"] >= raw_thresholds["Vkabat"]

	plot_df["meets_switch_like_criteria"] = (
		plot_df["E6_meets_threshold"]
		& plot_df["IsUnstruct_meets_threshold"]
		& plot_df["Vkabat_meets_threshold"]
	)

	return plot_df


def get_empty_region_dataframe():
	columns = [
		"region_id",
		"start_residue_position",
		"end_residue_position",
		"length",
		"residue_sequence",
		"mean_E6_plot",
		"mean_IsUnstruct_plot",
		"mean_Vkabat_plot",
		"mean_E6_normalized_plot",
		"mean_IsUnstruct_normalized_plot",
		"mean_Vkabat_normalized_plot",
		"max_Vkabat_plot",
		"min_E6_plot",
		"max_E6_plot",
		"min_IsUnstruct_plot",
		"max_IsUnstruct_plot",
		"E6_threshold_percentile",
		"IsUnstruct_threshold_percentile",
		"Vkabat_threshold_percentile",
		"E6_threshold_value",
		"IsUnstruct_threshold_value",
		"Vkabat_threshold_value",
		"minimum_switch_region_length",
		"smoothing_window"
	]

	return pd.DataFrame(columns=columns)


def add_region_record(records, region_df, threshold_info, args):
	if len(region_df) < args.minimum_switch_region_length:
		return

	raw_thresholds = threshold_info["raw"]

	record = {
		"region_id": len(records) + 1,
		"start_residue_position": int(region_df["residue_position"].iloc[0]),
		"end_residue_position": int(region_df["residue_position"].iloc[-1]),
		"length": int(len(region_df)),
		"residue_sequence": "".join(region_df["residue"].astype(str).tolist()),
		"mean_E6_plot": float(region_df["E6_plot"].mean()),
		"mean_IsUnstruct_plot": float(region_df["IsUnstruct_plot"].mean()),
		"mean_Vkabat_plot": float(region_df["Vkabat_plot"].mean()),
		"mean_E6_normalized_plot": float(region_df["E6_normalized_plot"].mean()),
		"mean_IsUnstruct_normalized_plot": float(region_df["IsUnstruct_normalized_plot"].mean()),
		"mean_Vkabat_normalized_plot": float(region_df["Vkabat_normalized_plot"].mean()),
		"max_Vkabat_plot": float(region_df["Vkabat_plot"].max()),
		"min_E6_plot": float(region_df["E6_plot"].min()),
		"max_E6_plot": float(region_df["E6_plot"].max()),
		"min_IsUnstruct_plot": float(region_df["IsUnstruct_plot"].min()),
		"max_IsUnstruct_plot": float(region_df["IsUnstruct_plot"].max()),
		"E6_threshold_percentile": float(args.e6_threshold_percentile),
		"IsUnstruct_threshold_percentile": float(args.isunstruct_threshold_percentile),
		"Vkabat_threshold_percentile": float(args.vkabat_threshold_percentile),
		"E6_threshold_value": float(raw_thresholds["E6"]),
		"IsUnstruct_threshold_value": float(raw_thresholds["IsUnstruct"]),
		"Vkabat_threshold_value": float(raw_thresholds["Vkabat"]),
		"minimum_switch_region_length": int(args.minimum_switch_region_length),
		"smoothing_window": int(args.smoothing_window)
	}

	records.append(record)


def find_switch_like_regions(plot_df, threshold_info, args):
	if threshold_info is None:
		return get_empty_region_dataframe()

	records = []
	current_indices = []
	previous_position = None

	for idx, row in plot_df.iterrows():
		current_position = int(row["residue_position"])
		meets_criteria = bool(row["meets_switch_like_criteria"])

		if meets_criteria:
			if not current_indices:
				current_indices = [idx]
			elif previous_position is not None and current_position == previous_position + 1:
				current_indices.append(idx)
			else:
				region_df = plot_df.loc[current_indices]
				add_region_record(records, region_df, threshold_info, args)
				current_indices = [idx]
		else:
			if current_indices:
				region_df = plot_df.loc[current_indices]
				add_region_record(records, region_df, threshold_info, args)
				current_indices = []

		previous_position = current_position

	if current_indices:
		region_df = plot_df.loc[current_indices]
		add_region_record(records, region_df, threshold_info, args)

	if not records:
		return get_empty_region_dataframe()

	return pd.DataFrame(records)


def add_defined_switch_region_spans(ax, defined_regions_df, args):
	if defined_regions_df.empty:
		return

	visible_regions = defined_regions_df[defined_regions_df["visible_in_plot"] == True]

	if visible_regions.empty:
		return

	label_added = False

	for _, row in visible_regions.iterrows():
		label = args.defined_switch_region_label if not label_added else None

		ax.axvspan(
			float(row["clipped_start_residue_position"]) - 0.5,
			float(row["clipped_end_residue_position"]) + 0.5,
			color=args.defined_switch_region_color,
			alpha=args.defined_switch_region_alpha,
			linewidth=0,
			label=label
		)

		label_added = True


def add_switch_like_spans(ax, regions_df, args):
	if regions_df.empty:
		return

	label_added = False

	for _, row in regions_df.iterrows():
		label = "Threshold-derived switch-like region" if not label_added else None

		ax.axvspan(
			float(row["start_residue_position"]) - 0.5,
			float(row["end_residue_position"]) + 0.5,
			color=args.switch_region_color,
			alpha=args.switch_region_alpha,
			linewidth=0,
			label=label
		)

		label_added = True


def build_threshold_legend_label(spec, args):
	if args.show_threshold_percentiles_in_legend:
		return f"{spec['descriptor']} {spec['operator']} P{percentile_label(spec['percentile'])}"

	return f"{spec['descriptor']} threshold"


def add_threshold_lines(ax, threshold_info, args, normalized=False):
	if threshold_info is None:
		return

	if normalized:
		thresholds = threshold_info["normalized"]
	else:
		thresholds = threshold_info["raw"]

	threshold_specs = [
		{
			"descriptor": "E6",
			"value": thresholds["E6"],
			"percentile": args.e6_threshold_percentile,
			"color": args.e6_color,
			"operator": "≤"
		},
		{
			"descriptor": "IsUnstruct",
			"value": thresholds["IsUnstruct"],
			"percentile": args.isunstruct_threshold_percentile,
			"color": args.isunstruct_color,
			"operator": "≤"
		},
		{
			"descriptor": "Vkabat",
			"value": thresholds["Vkabat"],
			"percentile": args.vkabat_threshold_percentile,
			"color": args.vkabat_color,
			"operator": "≥"
		}
	]

	for spec in threshold_specs:
		ax.axhline(
			y=spec["value"],
			color=spec["color"],
			linestyle=args.threshold_linestyle,
			linewidth=args.threshold_line_width,
			alpha=args.threshold_alpha,
			label=build_threshold_legend_label(spec, args)
		)


def plot_raw(plot_df, regions_df, defined_regions_df, threshold_info, args, output_path, prefix):
	fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))

	add_defined_switch_region_spans(ax, defined_regions_df, args)
	add_switch_like_spans(ax, regions_df, args)

	ax.plot(
		plot_df["residue_position"],
		plot_df["E6_plot"],
		label="E6",
		color=args.e6_color,
		linewidth=args.line_width
	)

	ax.plot(
		plot_df["residue_position"],
		plot_df["IsUnstruct_plot"],
		label="IsUnstruct",
		color=args.isunstruct_color,
		linewidth=args.line_width
	)

	ax.plot(
		plot_df["residue_position"],
		plot_df["Vkabat_plot"],
		label="Vkabat",
		color=args.vkabat_color,
		linewidth=args.line_width
	)

	add_threshold_lines(ax, threshold_info, args, normalized=False)

	fig.suptitle(get_raw_title(args, prefix), y=0.985)
	ax.set_xlabel("Residue Number")
	ax.set_ylabel("Raw Descriptor Value")

	if not args.no_grid:
		ax.grid(True, alpha=0.35)

	format_x_axis(ax, plot_df, args)
	format_y_axis(ax, plot_df, args, normalized=False)
	layout_rect = format_legend(fig, ax, args)

	if layout_rect is not None:
		fig.tight_layout(rect=layout_rect)
	else:
		fig.tight_layout()

	fig.savefig(output_path, dpi=args.dpi, bbox_inches="tight")
	plt.close(fig)


def plot_normalized(plot_df, regions_df, defined_regions_df, threshold_info, args, output_path, prefix):
	fig, ax = plt.subplots(figsize=(args.fig_width, args.fig_height))

	add_defined_switch_region_spans(ax, defined_regions_df, args)
	add_switch_like_spans(ax, regions_df, args)

	ax.plot(
		plot_df["residue_position"],
		plot_df["E6_normalized_plot"],
		label="E6 (Normalized)",
		color=args.e6_color,
		linewidth=args.line_width
	)

	ax.plot(
		plot_df["residue_position"],
		plot_df["IsUnstruct_normalized_plot"],
		label="IsUnstruct",
		color=args.isunstruct_color,
		linewidth=args.line_width
	)

	ax.plot(
		plot_df["residue_position"],
		plot_df["Vkabat_normalized_plot"],
		label="Vkabat (Normalized)",
		color=args.vkabat_color,
		linewidth=args.line_width
	)

	add_threshold_lines(ax, threshold_info, args, normalized=True)

	fig.suptitle(get_normalized_title(args, prefix), y=0.985)
	ax.set_xlabel("Residue Number")
	ax.set_ylabel("Normalized Descriptor Value")

	if not args.no_grid:
		ax.grid(True, alpha=0.35)

	format_x_axis(ax, plot_df, args)
	format_y_axis(ax, plot_df, args, normalized=True)
	layout_rect = format_legend(fig, ax, args)

	if layout_rect is not None:
		fig.tight_layout(rect=layout_rect)
	else:
		fig.tight_layout()

	fig.savefig(output_path, dpi=args.dpi, bbox_inches="tight")
	plt.close(fig)


def build_output_suffix(args):
	suffix = ""

	if args.res_start is not None or args.res_end is not None:
		start_label = args.res_start if args.res_start is not None else "start"
		end_label = args.res_end if args.res_end is not None else "end"
		suffix += f"_res{start_label}-{end_label}"

	if args.smoothing_window > 1:
		suffix += f"_smoothed_w{args.smoothing_window}"

	if args.show_thresholds:
		suffix += "_thresholds"

	if args.define_switch_regions is not None and args.define_switch_regions.strip():
		suffix += "_defined_regions"

	return suffix


def write_parameters_file(args, prefix, plot_df, regions_df, defined_regions_df, threshold_info, output_paths):
	lines = []

	lines.append("Descriptor Overlay Plotting Parameters")
	lines.append("=" * 40)
	lines.append(f"Generated: {datetime.now().isoformat(timespec='seconds')}")
	lines.append("")
	lines.append("Input files")
	lines.append("-" * 40)
	lines.append(f"Nextflow CSV: {args.nextflow_csv}")
	lines.append(f"Vkabat CSV: {args.vkabat_csv}")
	lines.append("")
	lines.append("Output files")
	lines.append("-" * 40)

	for label, path in output_paths.items():
		lines.append(f"{label}: {path}")

	lines.append("")
	lines.append("General settings")
	lines.append("-" * 40)
	lines.append(f"Prefix: {prefix}")
	lines.append(f"Title argument: {args.title}")
	lines.append(f"Output directory: {args.outdir}")
	lines.append(f"Figure width: {args.fig_width}")
	lines.append(f"Figure height: {args.fig_height}")
	lines.append(f"DPI: {args.dpi}")
	lines.append(f"Line width: {args.line_width}")
	lines.append(f"Smoothing window: {args.smoothing_window}")
	lines.append(f"Residue range start: {args.res_start}")
	lines.append(f"Residue range end: {args.res_end}")
	lines.append(f"Legend position: {args.legend_position}")
	lines.append(f"Legend columns: {args.legend_columns}")
	lines.append(f"Show threshold percentiles in legend: {args.show_threshold_percentiles_in_legend}")
	lines.append(f"X-axis tick interval: {args.x_tick_interval}")
	lines.append(f"Raw y-axis tick interval: {args.raw_y_tick_interval}")
	lines.append(f"Normalized y-axis tick interval: {args.normalized_y_tick_interval}")
	lines.append(f"Raw y-axis maximum: {args.raw_ymax}")
	lines.append("")
	lines.append("Input columns")
	lines.append("-" * 40)
	lines.append(f"E6 column: {args.e6_col}")
	lines.append(f"IsUnstruct column: {args.isunstruct_col}")
	lines.append(f"Residue column: {args.residue_col}")
	lines.append(f"Residue position column: {args.position_col}")
	lines.append(f"Vkabat column: {args.vkabat_col}")
	lines.append(f"N column: {args.n_col}")
	lines.append(f"Fallback N value: {args.fallback_N}")
	lines.append("")
	lines.append("Normalization")
	lines.append("-" * 40)
	lines.append("E6 normalized = E6 / log2(6)")
	lines.append("Vkabat normalized = (Vkabat - Vkabat_min) / (Vkabat_max - Vkabat_min)")
	lines.append("Vkabat_min = 1.0")
	lines.append("Vkabat_max is calculated dynamically from N and class_count.")
	lines.append(f"Class count: {args.class_count}")
	lines.append("")
	lines.append("Data summary")
	lines.append("-" * 40)
	lines.append(f"Number of residues plotted: {len(plot_df)}")
	lines.append(f"Minimum residue position plotted: {int(plot_df['residue_position'].min())}")
	lines.append(f"Maximum residue position plotted: {int(plot_df['residue_position'].max())}")
	lines.append(f"Unique N values: {sorted(plot_df['N'].unique().tolist())}")
	lines.append(f"N source: {plot_df['N_source'].iloc[0]}")
	lines.append(f"Minimum Vkabat_max: {float(plot_df['Vkabat_max'].min()):.6f}")
	lines.append(f"Maximum Vkabat_max: {float(plot_df['Vkabat_max'].max()):.6f}")
	lines.append("")
	lines.append("Predefined/historical switch-region settings")
	lines.append("-" * 40)
	lines.append(f"Defined switch regions argument: {args.define_switch_regions}")
	lines.append(f"Defined switch region color: {args.defined_switch_region_color}")
	lines.append(f"Defined switch region alpha: {args.defined_switch_region_alpha}")
	lines.append(f"Defined switch region label: {args.defined_switch_region_label}")
	lines.append(f"Number of defined switch regions: {len(defined_regions_df)}")

	if not defined_regions_df.empty:
		for _, row in defined_regions_df.iterrows():
			visible_text = "visible" if bool(row["visible_in_plot"]) else "not visible in plotted range"
			lines.append(
				f"Defined region {int(row['defined_region_id'])}: "
				f"{int(row['start_residue_position'])}-{int(row['end_residue_position'])} "
				f"({visible_text})"
			)

	lines.append("")
	lines.append("Threshold and switch-like region settings")
	lines.append("-" * 40)
	lines.append(f"Thresholds shown: {args.show_thresholds}")
	lines.append(f"E6 threshold percentile: {args.e6_threshold_percentile}")
	lines.append(f"IsUnstruct threshold percentile: {args.isunstruct_threshold_percentile}")
	lines.append(f"Vkabat threshold percentile: {args.vkabat_threshold_percentile}")
	lines.append(f"Minimum switch-like region length: {args.minimum_switch_region_length}")
	lines.append(f"Switch-like region alpha: {args.switch_region_alpha}")
	lines.append("Threshold-derived switch-like residue criterion:")
	lines.append("E6 <= E6 percentile threshold")
	lines.append("IsUnstruct <= IsUnstruct percentile threshold")
	lines.append("Vkabat >= Vkabat percentile threshold")
	lines.append("Consecutive switch-like residues are grouped into regions if their length is >= minimum_switch_region_length.")

	if threshold_info is not None:
		lines.append("")
		lines.append("Raw threshold values")
		lines.append("-" * 40)
		lines.append(f"E6 raw threshold: {threshold_info['raw']['E6']:.6f}")
		lines.append(f"IsUnstruct raw threshold: {threshold_info['raw']['IsUnstruct']:.6f}")
		lines.append(f"Vkabat raw threshold: {threshold_info['raw']['Vkabat']:.6f}")
		lines.append("")
		lines.append("Normalized threshold values")
		lines.append("-" * 40)
		lines.append(f"E6 normalized threshold: {threshold_info['normalized']['E6']:.6f}")
		lines.append(f"IsUnstruct normalized threshold: {threshold_info['normalized']['IsUnstruct']:.6f}")
		lines.append(f"Vkabat normalized threshold: {threshold_info['normalized']['Vkabat']:.6f}")

	lines.append("")
	lines.append("Threshold-derived switch-like region results")
	lines.append("-" * 40)
	lines.append(f"Number of threshold-derived switch-like regions: {len(regions_df)}")

	if not regions_df.empty:
		for _, row in regions_df.iterrows():
			lines.append(
				f"Region {int(row['region_id'])}: "
				f"{int(row['start_residue_position'])}-{int(row['end_residue_position'])} "
				f"(length {int(row['length'])})"
			)

	parameter_path = output_paths["parameters_txt"]
	Path(parameter_path).write_text("\n".join(lines) + "\n")


def main():
	args = parse_arguments()
	validate_arguments(args)

	outdir = Path(args.outdir)
	outdir.mkdir(parents=True, exist_ok=True)

	prefix = args.prefix if args.prefix else infer_prefix(args.nextflow_csv)
	suffix = build_output_suffix(args)

	plot_df = build_descriptor_dataframe(args)
	defined_regions_df = parse_defined_switch_regions(args, plot_df)
	threshold_info = calculate_threshold_info(plot_df, args)
	plot_df = add_switch_like_flags(plot_df, threshold_info, args)
	regions_df = find_switch_like_regions(plot_df, threshold_info, args)

	raw_output_path = outdir / f"{prefix}_descriptor_overlay_raw{suffix}.png"
	normalized_output_path = outdir / f"{prefix}_descriptor_overlay_normalized{suffix}.png"
	combined_csv_path = outdir / f"{prefix}_descriptor_overlay_values{suffix}.csv"
	switch_regions_csv_path = outdir / f"{prefix}_switch_like_regions{suffix}.csv"
	defined_switch_regions_csv_path = outdir / f"{prefix}_defined_switch_regions{suffix}.csv"
	parameters_txt_path = outdir / f"{prefix}_descriptor_overlay_parameters{suffix}.txt"

	output_paths = {
		"raw_plot_png": raw_output_path,
		"normalized_plot_png": normalized_output_path,
		"combined_values_csv": combined_csv_path,
		"threshold_derived_switch_like_regions_csv": switch_regions_csv_path,
		"defined_switch_regions_csv": defined_switch_regions_csv_path,
		"parameters_txt": parameters_txt_path
	}

	plot_raw(plot_df, regions_df, defined_regions_df, threshold_info, args, raw_output_path, prefix)
	plot_normalized(plot_df, regions_df, defined_regions_df, threshold_info, args, normalized_output_path, prefix)

	plot_df.to_csv(combined_csv_path, index=False)
	regions_df.to_csv(switch_regions_csv_path, index=False)
	defined_regions_df.to_csv(defined_switch_regions_csv_path, index=False)
	write_parameters_file(args, prefix, plot_df, regions_df, defined_regions_df, threshold_info, output_paths)

	print(f"Wrote raw overlay plot: {raw_output_path}")
	print(f"Wrote normalized overlay plot: {normalized_output_path}")
	print(f"Wrote combined descriptor CSV: {combined_csv_path}")
	print(f"Wrote threshold-derived switch-like regions CSV: {switch_regions_csv_path}")
	print(f"Wrote defined switch regions CSV: {defined_switch_regions_csv_path}")
	print(f"Wrote parameters TXT: {parameters_txt_path}")


if __name__ == "__main__":
	main()
