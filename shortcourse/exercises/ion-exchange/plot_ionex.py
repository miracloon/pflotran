#!/usr/bin/env python3
"""Plot Voegelin Ca-Mg-Na breakthrough curves.

Experimental points come from exp1.dat. PFLOTRAN curves come from
voegelin-obs-0.pft, the column-outlet observation file produced by
running voegelin.in with the input prefix voegelin. Use --obs to select
another observation file, including older PFLOTRAN output formats.

Pore volume uses the Darcy velocity in voegelin.in:

    q = 0.61 * 5.6 cm/min = 5.69333e-4 m/s
    t_pv = 0.481 m * 0.61 / q

Run from this directory:

    pflotran -input_prefix voegelin
    python3 plot_ionex.py
    python3 plot_ionex.py --obs voegelin-obs-0.pft --output ionex.png
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
from matplotlib.ticker import FuncFormatter, MultipleLocator
import numpy as np

HERE = Path(__file__).resolve().parent

LENGTH_M = 0.481
POROSITY = 0.61
DARCY_M_PER_S = 5.69333e-4
SECONDS_PER_PORE_VOLUME = LENGTH_M * POROSITY / DARCY_M_PER_S

SPECIES = ('Ca++', 'Mg++', 'Na+')

# One Okabe–Ito hue per cation. Filled symbols are the measurements.
# Solid lines are PFLOTRAN.
STYLE = {
    'Ca++': {'color': '#0072B2', 'marker': 'o', 'label': r'Ca$^{2+}$'},
    'Mg++': {'color': '#D55E00', 'marker': 's', 'label': r'Mg$^{2+}$'},
    'Na+': {'color': '#009E73', 'marker': 'D', 'label': r'Na$^{+}$'},
}

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Avenir', 'Helvetica Neue', 'Arial', 'DejaVu Sans'],
    'axes.linewidth': 0.8,
    'axes.labelsize': 13,
    'axes.labelcolor': '#222222',
    'axes.edgecolor': '#333333',
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'xtick.color': '#333333',
    'ytick.color': '#333333',
    'text.color': '#222222',
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'savefig.facecolor': 'white',
})

# Legacy headerless observation columns (1-based: time, Na, Ca, Mg).
LEGACY_COLUMNS = {'time': 0, 'Na+': 2, 'Ca++': 3, 'Mg++': 4}

SECONDS_PER_UNIT = {
    's': 1.0,
    'sec': 1.0,
    'min': 60.0,
    'h': 3600.0,
    'hr': 3600.0,
    'd': 86400.0,
    'day': 86400.0,
    'y': 365.25 * 86400.0,
    'yr': 365.25 * 86400.0,
}


def read_experiment(path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Read exp1.dat. Columns are (PV, Ca), (PV, Mg), (PV, Na).

    Later rows omit Mg and Na. Those blank fields stay missing so the Ca
    series can continue on its own.
    """
    rows = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith('#'):
            continue
        fields = line.split('\t') if '\t' in line else stripped.split()
        if len(fields) < 6:
            fields = fields + [''] * (6 - len(fields))
        values = []
        for field in fields[:6]:
            text = field.strip()
            values.append(np.nan if text == '' else float(text))
        rows.append(values)
    if not rows:
        raise SystemExit(f'{path} has no experimental points')
    table = np.asarray(rows, dtype=float)
    if table.shape[1] < 6:
        raise SystemExit(f'{path} needs six columns: PV Ca PV Mg PV Na')
    return {
        'Ca++': (table[:, 0], table[:, 1]),
        'Mg++': (table[:, 2], table[:, 3]),
        'Na+': (table[:, 4], table[:, 5]),
    }


def _header_fields(line: str) -> list[str]:
    return [field.strip() for field in next(csv.reader([line.strip()]))]


def _time_unit(header: list[str]) -> str:
    name = header[0]
    if '[' in name and ']' in name:
        return name.split('[', 1)[1].split(']', 1)[0].strip().lower()
    return 's'


def _species_column(header: list[str], species: str) -> int | None:
    """Column of aqueous total molality, skipping bulk and sorbed totals."""
    matches = []
    for index, name in enumerate(header):
        lowered = name.lower()
        if species.lower() not in lowered:
            continue
        if 'bulk' in lowered or 'sorbed' in lowered or 'gamma' in lowered:
            continue
        if 'free ' in lowered:
            continue
        if not lowered.lstrip().split()[0].rstrip('-').isdigit() and \
                not name.lower().startswith('total'):
            # Quoted names look like '3-Total Ca++ [m] east (...)'.
            if 'total' not in lowered:
                continue
        if '[m]' not in lowered and 'total' not in lowered:
            continue
        matches.append(index)
    if not matches:
        return None
    # Prefer an explicit Total molality column over any other hit.
    for index in matches:
        if 'total' in header[index].lower() and '[m]' in header[index].lower():
            return index
    return matches[0]


def read_observation(path: Path) -> dict[str, np.ndarray]:
    """Return pore volume and concentration arrays keyed by species name."""
    lines = path.read_text().splitlines()
    header = None
    start = 0
    for i, line in enumerate(lines):
        stripped = line.strip()
        if not stripped or stripped.startswith('#'):
            continue
        if stripped.upper().startswith('ZONE') or stripped.upper().startswith('VARIABLES'):
            continue
        if '"' in stripped:
            header = _header_fields(stripped)
            start = i + 1
            break
        start = i
        break

    rows = []
    for line in lines[start:]:
        stripped = line.strip()
        if not stripped or stripped.startswith('#') or stripped.startswith('"'):
            continue
        if stripped.upper().startswith('ZONE') or stripped.upper().startswith('VARIABLES'):
            continue
        rows.append([float(value) for value in stripped.split()])
    if not rows:
        raise SystemExit(f'{path} has no observation data')
    data = np.asarray(rows, dtype=float)

    if header is None:
        needed = max(LEGACY_COLUMNS.values()) + 1
        if data.shape[1] < needed:
            raise SystemExit(
                f'{path} has no header and fewer than {needed} columns. '
                'Legacy observations require time, an unused column, Na+, Ca++, Mg++.')
        seconds = data[:, LEGACY_COLUMNS['time']]
        concentrations = {
            species: data[:, LEGACY_COLUMNS[species]] for species in SPECIES
        }
    else:
        if len(header) != data.shape[1]:
            raise SystemExit(
                f'{path} header has {len(header)} names but data rows have '
                f'{data.shape[1]} values')
        unit = _time_unit(header)
        if unit not in SECONDS_PER_UNIT:
            raise SystemExit(f'Unrecognized time unit "{unit}" in {path}')
        seconds = data[:, 0] * SECONDS_PER_UNIT[unit]
        concentrations = {}
        for species in SPECIES:
            column = _species_column(header, species)
            if column is None:
                raise SystemExit(
                    f'{path} has no Total {species} [m] column. Header names:\n'
                    + '\n'.join(header))
            concentrations[species] = data[:, column]

    pore_volume = seconds / SECONDS_PER_PORE_VOLUME
    concentrations['pore_volume'] = pore_volume
    return concentrations


def find_observation(explicit: Path | None) -> Path | None:
    if explicit is not None:
        if not explicit.is_file():
            raise SystemExit(f'Observation file not found: {explicit}')
        return explicit
    candidate = HERE / 'voegelin-obs-0.pft'
    return candidate if candidate.is_file() else None


def plot(experiment, observation, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 5.15), layout='constrained')
    handles = []
    labels = []
    for species in SPECIES:
        style = STYLE[species]
        pv, conc = experiment[species]
        measured, = ax.plot(
            pv, conc,
            linestyle='none',
            marker=style['marker'],
            color=style['color'],
            markerfacecolor=style['color'],
            markeredgecolor='white',
            markeredgewidth=0.7,
            markersize=6.5,
            zorder=3,
        )
        if observation is None:
            handles.append(measured)
        else:
            modeled, = ax.plot(
                observation['pore_volume'], observation[species],
                color=style['color'],
                linewidth=1.8,
                solid_capstyle='round',
                solid_joinstyle='round',
                zorder=2,
            )
            handles.append((measured, modeled))
        labels.append(style['label'])

    ax.set_xlim(0, 80)
    ax.set_ylim(0, 0.01)
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_major_locator(MultipleLocator(0.002))
    ax.yaxis.set_minor_locator(MultipleLocator(0.001))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _pos: f'{value:g}'))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _pos: f'{value:g}'))
    ax.set_xlabel('Pore volume')
    ax.set_ylabel(r'Concentration [mol/kg-H$_2$O]')
    ax.tick_params(width=0.8, length=4)
    ax.tick_params(which='minor', length=2.4, width=0.6)
    ax.yaxis.grid(True, color='#E4E4E4', linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    legend = ax.legend(
        handles, labels,
        handler_map={tuple: HandlerTuple(ndivide=None, pad=0.6)},
        title='Experiment  ·  PFLOTRAN',
        loc='upper left',
        frameon=False,
        fontsize=12,
        title_fontsize=11,
        handlelength=3.2,
        borderaxespad=0.2,
        labelspacing=0.45,
    )
    legend.get_title().set_color('#444444')
    fig.savefig(output, dpi=160)
    print(f'Wrote {output}')


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--exp', type=Path, default=HERE / 'exp1.dat',
        help='Experimental breakthrough file (default: exp1.dat beside this script)')
    parser.add_argument(
        '--obs', type=Path, default=None,
        help='PFLOTRAN observation file (default: voegelin-obs-0.pft beside this script)')
    parser.add_argument(
        '--output', type=Path, default=HERE / 'ionex.png',
        help='Figure path (default: ionex.png beside this script)')
    args = parser.parse_args(argv)

    if not args.exp.is_file():
        raise SystemExit(f'Experimental file not found: {args.exp}')
    observation_path = find_observation(args.obs)
    if observation_path is None:
        print(
            'No voegelin-obs-0.pft next to this script. '
            'Run pflotran -input_prefix voegelin or select a file with --obs. '
            'Plotting the experimental points only.',
            file=sys.stderr)
        observation = None
    else:
        print(f'Reading {observation_path}')
        observation = read_observation(observation_path)
    plot(read_experiment(args.exp), observation, args.output)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
