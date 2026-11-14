#!/usr/bin/env python3
"""Summarize transport event CSV files; plotting is optional."""

import argparse
from bisect import bisect_right
import csv
import math
from pathlib import Path
import statistics
import sys


BALANCE_LIMIT_MEV = 1e-6
TAIL_ANGLE_RAD = 0.005
COLUMNS = (
    "event_id", "initial_energy_MeV", "final_energy_MeV", "loss_MeV",
    "edep_MeV", "initial_charge", "final_charge", "exit_angle_rad",
    "exited", "energy_balance_MeV",
)
METRICS = ("loss_mean_MeV", "loss_variance_MeV2",
           "edep_mean_MeV", "edep_variance_MeV2")
ERRORS = ("loss_mean_se_MeV", "loss_variance_se_MeV2",
          "edep_mean_se_MeV", "edep_variance_se_MeV2")
ANGLE_QUANTILES = (("q68", 0.68), ("q95", 0.95))
SUMMARY_FIELDS = (
    "sample", "input_file", "events", "exited_events", "invalid_rows",
    "max_abs_energy_balance_MeV", "integrity", "energy_core_sigma",
    "loss_core_events", "edep_core_events", "charge_fractions_available",
    "loss_mean_MeV", "loss_mean_se_MeV", "loss_sigma_MeV",
    "loss_variance_MeV2", "loss_variance_se_MeV2",
    "edep_mean_MeV", "edep_mean_se_MeV", "edep_sigma_MeV",
    "edep_variance_MeV2", "edep_variance_se_MeV2",
    "exit_angle_q68_rad", "exit_angle_q95_rad",
    "exit_angle_above_5mrad_count", "exit_angle_above_5mrad_fraction",
)


def sample_names(paths):
    parts = [path.with_suffix("").parts for path in paths]
    names = [part[-1] for part in parts]
    for depth in range(2, max(map(len, parts)) + 1):
        duplicates = {name for name in names if names.count(name) > 1}
        if not duplicates:
            break
        names = ["/".join(part[-depth:]) if name in duplicates else name
                 for part, name in zip(parts, names)]
    return names


def read_events(path):
    values = {key: [] for key in ("loss", "edep", "angle", "charge")}
    invalid = 0
    largest_balance = 0.0
    event_ids = set()
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        missing = set(COLUMNS) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            try:
                if None in row:
                    raise ValueError("extra CSV fields")
                numbers = {key: float(row[key]) for key in COLUMNS}
                if not all(math.isfinite(value) for value in numbers.values()):
                    raise ValueError("non-finite value")
                event_id = int(numbers["event_id"])
                if event_id < 0 or event_id != numbers["event_id"]:
                    raise ValueError("event_id must be a nonnegative integer")
                if event_id in event_ids:
                    raise ValueError(f"duplicate event_id {event_id}")
                if numbers["exited"] not in (0, 1):
                    raise ValueError("exited must be 0 or 1")
            except (ValueError, TypeError) as error:
                invalid += 1
                if invalid <= 5:
                    print(f"{path}:{reader.line_num}: {error}", file=sys.stderr)
                continue
            event_ids.add(event_id)
            values["loss"].append(numbers["loss_MeV"])
            values["edep"].append(numbers["edep_MeV"])
            largest_balance = max(largest_balance, abs(numbers["energy_balance_MeV"]))
            if numbers["exited"]:
                values["angle"].append(numbers["exit_angle_rad"])
                values["charge"].append(numbers["final_charge"])
    return values, invalid, largest_balance


def moments(values):
    n = len(values)
    if not n:
        return dict.fromkeys(("mean", "mean_se", "sigma", "variance", "variance_se"))
    mean = statistics.fmean(values)
    if n < 2:
        return {"mean": mean, "mean_se": None, "sigma": None,
                "variance": None, "variance_se": None}
    differences = [value - mean for value in values]
    variance = math.fsum(value * value for value in differences) / (n - 1)
    fourth_moment = math.fsum(value ** 4 for value in differences) / n
    variance_se = math.sqrt(max(0.0, (
        fourth_moment - (n - 3) / (n - 1) * variance ** 2
    ) / n))
    return {"mean": mean, "mean_se": math.sqrt(variance / n),
            "sigma": math.sqrt(variance), "variance": variance,
            "variance_se": variance_se}


def quantile(sorted_values, probability):
    if not sorted_values:
        return None
    position = probability * (len(sorted_values) - 1)
    lower = int(position)
    upper = min(lower + 1, len(sorted_values) - 1)
    return sorted_values[lower] + (position - lower) * (
        sorted_values[upper] - sorted_values[lower])


def select_core(values, width):
    if not width or len(values) < 2:
        return values
    selected = values
    for _ in range(100):
        mean = statistics.fmean(selected)
        sigma = statistics.stdev(selected) if len(selected) > 1 else 0.0
        updated = [value for value in values if abs(value - mean) <= width * sigma]
        if updated == selected:
            return selected
        if not updated:
            raise ValueError("The core window contains no events")
        selected = updated
    raise ValueError("The core window did not converge")


def discrete_charges(values):
    if not values or any(abs(value - round(value)) > 1e-6 for value in values):
        return None
    counts = {}
    for value in values:
        charge = round(value)
        counts[charge] = counts.get(charge, 0) + 1
    return counts


def summarize(name, path, values, invalid, balance, core_sigma=0):
    n = len(values["loss"])
    row = {
        "sample": name, "input_file": str(path), "events": n,
        "exited_events": len(values["angle"]), "invalid_rows": invalid,
        "max_abs_energy_balance_MeV": balance if n else None,
        "integrity": "OK" if n and not invalid and balance <= BALANCE_LIMIT_MEV else "ERROR",
        "energy_core_sigma": core_sigma,
        "charge_fractions_available": int(discrete_charges(values["charge"]) is not None),
    }
    for observable in ("loss", "edep"):
        selected = select_core(values[observable], core_sigma)
        row[f"{observable}_core_events"] = len(selected)
        for key, value in moments(selected).items():
            unit = "MeV2" if key.startswith("variance") else "MeV"
            row[f"{observable}_{key}_{unit}"] = value
    angles = sorted(values["angle"])
    row["exit_angle_q68_rad"] = quantile(angles, 0.68)
    row["exit_angle_q95_rad"] = quantile(angles, 0.95)
    row["exit_angle_above_5mrad_count"] = sum(angle > TAIL_ANGLE_RAD for angle in angles)
    row["exit_angle_above_5mrad_fraction"] = row["exit_angle_above_5mrad_count"] / len(angles) if angles else None
    return row


def read_reference(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    reference = {}
    for row in rows:
        name = row.get("sample")
        if not name or name in reference:
            raise ValueError(f"{path}: missing or duplicate sample name")
        core = float(row.get("energy_core_sigma") or 0)
        if not math.isfinite(core) or core < 0:
            raise ValueError(f"{path}: invalid energy_core_sigma")
        reference[name] = row
    return reference


def read_reference_charges(path):
    if not path.exists():
        return None
    reference = {}
    with path.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            try:
                name = row["sample"]
                charge, count = float(row["final_charge"]), int(row["count"])
                counts = reference.setdefault(name, {})
                if not name or not math.isfinite(charge) or count < 0 or charge in counts:
                    raise ValueError("invalid or duplicate charge count")
                counts[charge] = count
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError(f"{path}: invalid reference charge row") from error
    return reference


def number(value):
    return "" if value is None else format(value, ".12g")


def standard_score(difference, uncertainty):
    if uncertainty:
        return difference / uncertainty
    return 0.0 if difference == 0 else math.copysign(math.inf, difference)


def fraction_score(count, n, previous_count, previous_n):
    pooled = (count + previous_count) / (n + previous_n)
    uncertainty = math.sqrt(pooled * (1 - pooled) * (1 / n + 1 / previous_n))
    return standard_score(count / n - previous_count / previous_n, uncertainty)


def compare(row, values, counts, reference, reference_charges):
    old = reference.get(row["sample"])
    status = "UNMATCHED" if old is None else "COMPARED"
    same_core = old is None or row["energy_core_sigma"] == float(old.get("energy_core_sigma") or 0)
    if not same_core:
        print("  Energy comparison unavailable: core windows differ")
        status = "UNAVAILABLE"
    for metric, error in zip(METRICS, ERRORS):
        z = None
        if old is not None and same_core:
            try:
                previous, previous_se = float(old[metric]), float(old[error])
                if not math.isfinite(previous) or not math.isfinite(previous_se) or previous_se < 0:
                    raise ValueError("invalid reference value or uncertainty")
            except (KeyError, TypeError, ValueError):
                status = "UNAVAILABLE" if status != "REVIEW" else status
                previous = previous_se = None
            current, current_se = row[metric], row[error]
            if current is None or current_se is None:
                status = "UNAVAILABLE" if status != "REVIEW" else status
            if previous is not None and current is not None and current_se is not None:
                difference = current - previous
                combined_se = math.hypot(current_se, previous_se)
                z = standard_score(difference, combined_se)
                flag = " REVIEW" if abs(z) > 4 else ""
                if flag:
                    status = "REVIEW"
                print(f"  {metric}: {number(current)} +/- {number(current_se)}; "
                      f"reference {number(previous)} +/- {number(previous_se)}; "
                      f"difference {number(difference)}, z={number(z)}{flag}")
        row[f"reference_{metric}_z"] = z
    n = row["exited_events"]
    try:
        previous_n = int(old["exited_events"]) if old is not None else 0
    except (KeyError, TypeError, ValueError):
        previous_n = 0
    for label, probability in ANGLE_QUANTILES:
        fraction = z = None
        if old is not None:
            try:
                cutoff = float(old[f"exit_angle_{label}_rad"])
                if not math.isfinite(cutoff) or not 0 <= cutoff <= math.pi or n <= 0 or previous_n <= 0:
                    raise ValueError("invalid angular reference")
                fraction = sum(angle <= cutoff for angle in values["angle"]) / n
                combined_se = math.sqrt(probability * (1 - probability) * (1 / n + 1 / previous_n))
                z = standard_score(fraction - probability, combined_se)
                flag = " REVIEW" if abs(z) > 4 else ""
                if flag:
                    status = "REVIEW"
                print(f"  exit_angle_{label}: fraction <= {number(cutoff)} rad "
                      f"is {number(fraction)}; reference {number(probability)}; "
                      f"z={number(z)}{flag}")
            except (KeyError, TypeError, ValueError):
                status = "UNAVAILABLE" if status != "REVIEW" else status
        row[f"reference_exit_angle_{label}_cdf"] = fraction
        row[f"reference_exit_angle_{label}_z"] = z
    tail_z = None
    if old is not None and old.get("exit_angle_above_5mrad_count") not in (None, ""):
        try:
            previous_count = int(old["exit_angle_above_5mrad_count"])
            if n <= 0 or previous_n <= 0 or not 0 <= previous_count <= previous_n:
                raise ValueError("invalid angular tail reference")
            tail_z = fraction_score(row["exit_angle_above_5mrad_count"], n, previous_count, previous_n)
            flag = " REVIEW" if abs(tail_z) > 4 else ""
            if flag:
                status = "REVIEW"
            print(f"  angle > 5 mrad: {row['exit_angle_above_5mrad_count']}/{n}; "
                  f"reference {previous_count}/{previous_n}; z={number(tail_z)}{flag}")
        except (TypeError, ValueError):
            status = "UNAVAILABLE" if status != "REVIEW" else status
    row["reference_exit_angle_above_5mrad_z"] = tail_z
    max_charge_z = None
    old_discrete = old is not None and old.get("charge_fractions_available", "1") == "1"
    if old is not None and old_discrete != bool(row["charge_fractions_available"]):
        status = "UNAVAILABLE" if status != "REVIEW" else status
    elif old_discrete and reference_charges is not None:
        previous_counts = reference_charges.get(row["sample"], {})
        if n > 0 and previous_n > 0 and sum(previous_counts.values()) == previous_n:
            scores = {}
            for charge in counts.keys() | previous_counts.keys():
                count, previous_count = counts.get(charge, 0), previous_counts.get(charge, 0)
                scores[charge] = abs(fraction_score(count, n, previous_count, previous_n))
            worst_charge = max(scores, key=scores.get)
            max_charge_z = scores[worst_charge]
            flag = " REVIEW" if max_charge_z > 4 else ""
            if flag:
                status = "REVIEW"
            print(f"  charge fractions: max |z|={number(max_charge_z)} "
                  f"at charge {number(worst_charge)}{flag}")
        else:
            status = "UNAVAILABLE" if status != "REVIEW" else status
    row["reference_charge_max_abs_z"] = max_charge_z
    row["reference_status"] = status
    if status in ("UNMATCHED", "UNAVAILABLE"):
        print(f"  Reference comparison: {status}")


def write_csv(path, fields, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_histograms(path):
    result = {}
    if path is not None and path.exists():
        with path.open(newline="", encoding="utf-8-sig") as stream:
            for row in csv.DictReader(stream):
                low, high, density = (float(row[key]) for key in ("bin_low", "bin_high", "density"))
                if not all(math.isfinite(x) for x in (low, high, density)) or high <= low or density < 0:
                    raise ValueError(f"{path}: invalid histogram bin")
                result.setdefault((row["sample"], row["observable"]), []).append((low, high, density))
    return {key: sorted(bins) for key, bins in result.items()}


def histogram(values, edges):
    counts = [0] * (len(edges) - 1)
    for value in values:
        index = min(bisect_right(edges, value) - 1, len(counts) - 1)
        if index >= 0 and value <= edges[-1]:
            counts[index] += 1
    return [count / (len(values) * (right - left))
            for count, left, right in zip(counts, edges, edges[1:])]


def plot_samples(samples, output, references):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise ValueError("--plot requires Matplotlib (python -m pip install matplotlib)") from error

    def bins_for(key):
        values = [value for _, data in samples for value in data[key]]
        if not values:
            return [0.0, 1.0]
        low, high = min(values), max(values)
        if key == "angle" and high > 0:
            positive_low = min(value for value in values if value > 0)
            if positive_low < high:
                log_low, log_span = math.log(positive_low), math.log(high / positive_low)
                return [0.0] + [math.exp(log_low + log_span * i / 119) for i in range(120)]
        if low == high:
            margin = max(abs(low) * 0.01, 1e-6)
            low, high = low - margin, high + margin
        return [low + (high - low) * i / 80 for i in range(81)]

    def plot_density(ax, name, data, key, bins):
        ref = references.get((name, key))
        edges = [ref[0][0]] + [item[1] for item in ref] if ref else bins
        line = ax.stairs(histogram(data[key], edges), edges, label=name)
        if ref:
            ax.stairs([item[2] for item in ref], edges, linestyle="--",
                      color=line.get_edgecolor(), label=name + " reference")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    for ax, key, label in zip(axes, ("loss", "edep"), ("Energy loss (MeV)", "Deposited energy (MeV)")):
        bins = bins_for(key)
        for name, data in samples:
            if data[key]:
                plot_density(ax, name, data, key, bins)
        ax.set(xlabel=label, ylabel="Probability density (1/MeV)")
        ax.legend(fontsize="small")
    energy_path = output.with_name(output.stem + "_energy.png")
    fig.savefig(energy_path, dpi=160)
    plt.close(fig)

    fig, (ax, tail) = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    bins = bins_for("angle")
    for name, data in samples:
        if data["angle"]:
            plot_density(ax, name, data, "angle", bins)
            angles = sorted(data["angle"])
            line, = tail.step(angles, [(len(angles)-i)/len(angles) for i in range(len(angles))],
                              where="post", label=name)
            ref = references.get((name, "angle"))
            if ref:
                cumulative = 0.0
                survival = []
                for low, high, density in ref:
                    survival.append(max(0.0, 1-cumulative))
                    cumulative += (high-low)*density
                tail.step([item[0] for item in ref], survival, where="post", linestyle="--",
                          color=line.get_color(), label=name + " reference")
    ax.set(xlabel="Exit angle (rad)", ylabel="Probability density (1/rad)")
    if any(data["angle"] for _, data in samples):
        limit = max(quantile(sorted(data["angle"]), 0.999)
                    for _, data in samples if data["angle"])
        if limit > 0:
            ax.set_xlim(0, limit)
        ax.legend(fontsize="small")
        tail.legend(fontsize="small")
    tail.set(xlabel="Exit angle (rad)", ylabel="Fraction above angle", xscale="log", yscale="log")
    angle_path = output.with_name(output.stem + "_angles.png")
    fig.savefig(angle_path, dpi=160)
    plt.close(fig)
    print(f"Plots: {energy_path}, {angle_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path, metavar="events.csv")
    parser.add_argument("--output", type=Path, default=Path("summary.csv"))
    parser.add_argument("--reference", type=Path, help="Compare energy moments, exit angles and adjacent charges.csv, matching sample names")
    parser.add_argument("--core-sigma", type=float, help="Iterative energy-core window; otherwise use the reference setting or all events")
    parser.add_argument("--plot", action="store_true", help="Save energy and exit-angle histograms; requires Matplotlib")
    args = parser.parse_args()
    if args.core_sigma is not None and (not math.isfinite(args.core_sigma) or args.core_sigma <= 0):
        parser.error("--core-sigma must be finite and positive")
    paths = [path.resolve() for path in args.inputs]
    if len(set(paths)) != len(paths):
        parser.error("each input file must be listed only once")
    output = args.output.resolve()
    charges_path = output.with_name("charges.csv")
    if output == charges_path or output in paths or charges_path in paths:
        parser.error("output filenames must differ from inputs and charges.csv")
    if args.reference and (output == args.reference.resolve()
                           or charges_path == args.reference.with_name("charges.csv").resolve()):
        parser.error("write analysis results outside the reference directory")
    reference = read_reference(args.reference) if args.reference else None
    reference_charges = read_reference_charges(args.reference.with_name("charges.csv")) if args.reference else None
    summaries, charges, samples = [], [], []
    for name, path in zip(sample_names(paths), paths):
        values, invalid, balance = read_events(path)
        old = reference.get(name, {}) if reference else {}
        core_sigma = args.core_sigma if args.core_sigma is not None else float(old.get("energy_core_sigma") or 0)
        row = summarize(name, path, values, invalid, balance, core_sigma)
        summaries.append(row)
        plotted = dict(values)
        for key in ("loss", "edep"):
            plotted[key] = select_core(values[key], core_sigma)
        samples.append((name, plotted))
        print(f"{name}: {row['events']} events, {row['exited_events']} exited; "
              f"max |energy balance| = {number(row['max_abs_energy_balance_MeV'])} MeV; "
              f"data integrity {row['integrity']}")
        for key in ("loss", "edep"):
            print(f"  {key}: mean {number(row[key + '_mean_MeV'])} +/- "
                  f"{number(row[key + '_mean_se_MeV'])} MeV, "
                  f"sigma {number(row[key + '_sigma_MeV'])} MeV")
        if core_sigma:
            print(f"  Energy core +/-{core_sigma:g} sigma: "
                  f"{row['loss_core_events']} loss and {row['edep_core_events']} deposit events")
        counts = discrete_charges(values["charge"])
        if counts is None:
            print("  Discrete charge fractions unavailable: no exits or noninteger effective charges")
        for charge, count in sorted((counts or {}).items()):
            fraction = count / row["exited_events"]
            charges.append({"sample": name, "final_charge": charge, "count": count,
                            "fraction": fraction,
                            "fraction_se": math.sqrt(fraction * (1 - fraction) / row["exited_events"])})
        if reference is not None:
            compare(row, values, counts or {}, reference, reference_charges)
    fields = list(SUMMARY_FIELDS)
    if reference is not None:
        fields += [f"reference_{metric}_z" for metric in METRICS]
        fields += [f"reference_exit_angle_{label}_{metric}"
                   for label, _ in ANGLE_QUANTILES for metric in ("cdf", "z")]
        fields += ["reference_exit_angle_above_5mrad_z", "reference_charge_max_abs_z", "reference_status"]
    output.parent.mkdir(parents=True, exist_ok=True)
    write_csv(output, fields, summaries)
    write_csv(charges_path, ("sample", "final_charge", "count", "fraction", "fraction_se"), charges)
    print(f"Written: {output}, {charges_path}")
    if args.plot:
        histograms = read_histograms(args.reference.with_name("histograms.csv") if args.reference else None)
        plot_samples(samples, output, histograms)
    return 1 if any(row["integrity"] != "OK" for row in summaries) else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, ArithmeticError, csv.Error) as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(2)
